"""
    relative_algorithm_uncertainty(model, stations, events, differentials, result;
        event_ids, absolute_picks=[], relocate_options=(;),
        observation_covariance_s2=nothing, assume_independent=false,
        method=:auto, noise_samples=128, seed=0, step=.02, replay_tasks=1, progress=nothing,
        derivative_rtol=1e-3, derivative_position_atol_km=1e-7,
        derivative_origin_atol_s=1e-8, event_confidence_probability=.95)

Estimate how observation errors change the positions and origin times
returned by the complete relocation procedure. Perturb observations and rerun
all configured iterations, weights and constraints; these complete reruns are
called replays below. This is first-order error propagation through the
scaled, damped relocation, including finite schedules, rather than an estimate
of velocity-model error or location accuracy. `events` is the ORIGINAL starting
catalog, held exact; `relocate_options` must contain the actual relocation
keywords. A same-data replay must first reproduce `result`. No final-step
inverse or correction covariance is substituted for the returned catalog.

For observation covariance C=L*L', `:exact` streams each nonzero column of L.
`:randomized` instead streams `noise_samples` independent Rademacher vectors
L*z, E[z*z']=I. `:auto` uses exact propagation for at most 128 nonzero noise
directions and randomized propagation otherwise. Independent sigma uses a
Diagonal factor; neither method builds an observation-by-sample matrix.

Central replays at y +/- h*v and y +/- h/2*v give D_h and D_h/2. The factor is the Richardson
value (4D_h/2-D_h)/3; |D_h/2-D_h|/3 must meet the stated, separate km and s
budgets. `step` is dimensionless in standard-noise coordinates. Continuous
weights, column scales and directions are recomputed in every replay.
Changed selection, exits, retained rank, stages, stopping or ray branches,
non-smooth rays and active depth clipping are declined. This is a local
numerical check, not proof of differentiability for finite observational
errors. Full covariance follows ALL input differentials then absolute picks.

For randomized propagation the event covariance is mean(q*q'), without
subtracting the known zero directional mean. This is an unbiased covariance
estimator for an exactly linear response. `sampling_check` separately reports
the empirical Monte Carlo standard error of every covariance element, the
diagonal relative errors, and the two-half covariance difference. Those are
sampling diagnostics, not Richardson errors or a guaranteed accuracy bound.
The seed, actual sample count and RNG are saved. Only requested event factors
(four coordinates by the direction/sample count) are retained.

`replay_tasks` bounds concurrent +h/-h/+h/2/-h/2 solves within ONE noise
direction. Directions, random draws, failure checks and early stopping retain
their serial order. Choose either this parallelism or `relocate_options.tasks`
above one, to avoid nested solve workspaces. Perturbed runs retain only history
needed by the same discrete-path check; the baseline keeps full history.

The result describes Gaussian linearized fluctuations of the algorithm's
returned coordinates conditional on its initial catalog, model and settings.
It excludes initialization/reference/model error and finite-iteration or
truncation bias. Marginal normal half-widths are propagation intervals, not
verified coverage of true hypocenters. A zero response in an overall position
or common origin-time shift does not establish its physical accuracy.
Uncertainties of fixed reference-event values stay missing.
"""
function relative_algorithm_uncertainty(model,stations::AbstractDict,events,differentials,result;
        event_ids,absolute_picks=NamedTuple[],relocate_options::NamedTuple=(;),
        observation_covariance_s2=nothing,assume_independent::Bool=false,step::Real=.02,
        method::Symbol=:auto,noise_samples::Integer=128,seed::Integer=0,replay_tasks::Integer=1,progress=nothing,
        derivative_rtol::Real=1e-3,derivative_position_atol_km::Real=1e-7,
        derivative_origin_atol_s::Real=1e-8,event_confidence_probability::Real=.95)
    catalog,data,arrivals=collect(events),collect(differentials),collect(absolute_picks)
    requested=String.(collect(event_ids))
    final=Dict(String(event.event_id)=>event for event in result.events)
    all(haskey(final,id) for id in requested) || throw(ArgumentError("Unknown uncertainty event ID."))
    get(relocate_options,:solver,:damped)==:damped ||
        throw(ArgumentError("Algorithm replay describes the column-scaled damped relocation."))
    haskey(relocate_options,:absolute_picks) && throw(ArgumentError("Supply absolute_picks separately."))
    all(value->isfinite(value) && value>0,(step,derivative_rtol,derivative_position_atol_km,
        derivative_origin_atol_s)) || throw(ArgumentError("Replay steps and derivative budgets must be finite positive."))
    isfinite(event_confidence_probability) && 0<event_confidence_probability<1 ||
        throw(ArgumentError("Marginal propagation probability must lie strictly between zero and one."))
    method in (:auto,:exact,:randomized) || throw(ArgumentError("Replay method must be :auto, :exact or :randomized."))
    noise_samples>=2 || throw(ArgumentError("Randomized covariance requires at least two noise samples."))
    seed>=0 || throw(ArgumentError("Replay seed must be nonnegative."))
    replay_tasks>=1 || throw(ArgumentError("replay_tasks must be a positive integer."))
    replay_tasks>1 && get(relocate_options,:tasks,1)>1 &&
        throw(ArgumentError("Choose replay_tasks>1 or relocate_options.tasks>1, not both; each level needs its own solve workspaces."))
    perturbation_options=merge(relocate_options,(_history=:replay,))
    n=length(data)+length(arrivals)
    shared_covariance=observation_covariance_s2===nothing ? _catalog_ct_covariance_s2(data,arrivals) : nothing
    if shared_covariance!==nothing && observation_covariance_s2===nothing && assume_independent
        throw(ArgumentError("Catalog CT rows equal arrival_i-arrival_j and share those pick errors. assume_independent=true would understate the formal error. Omit it to use the catalog-difference covariance, or pass observation_covariance_s2. Relocation weights are unchanged."))
    end
    effective_covariance=observation_covariance_s2===nothing && shared_covariance!==nothing && !assume_independent ?
        shared_covariance : observation_covariance_s2
    supplied=effective_covariance!==nothing
    noise=supplied ? _observation_covariance_factor(effective_covariance,n) :
        Diagonal(assume_independent ? Float64[row.sigma_s for row in Iterators.flatten((data,arrivals))] : zeros(n))
    # The existing PSD factor preserves exact zero eigenvalues; omit only
    # identically zero directions, without imposing a statistical rank cutoff.
    active_noise_columns=noise isa Diagonal ? findall(!iszero,noise.diag) :
        [j for j in axes(noise,2) if any(!iszero,view(noise,:,j))]
    noise_rank=length(active_noise_columns)
    actual_method=method==:auto ? (noise_rank<=128 ? :exact : :randomized) : method
    directions=actual_method==:exact ? noise_rank : noise_rank==0 ? 0 : Int(noise_samples)
    rng=MersenneTwister(seed)
    probability=Float64(event_confidence_probability)
    scale=sqrt(_spatial_chisq_quantile(probability,1))
    assumptions=(quantity=:algorithm_returned_coordinates,estimator=:complete_scaled_relocation_algorithm,
        covariance_method=actual_method==:exact ? :directional_central_replay : :randomized_directional_central_replay,
        replay_method=actual_method,requested_replay_method=method,
        direction_distribution=actual_method==:randomized ? :rademacher : :coordinate_directions,
        random_seed=actual_method==:randomized ? seed : nothing,
        rng=actual_method==:randomized ? :MersenneTwister : nothing,
        rng_julia_version=actual_method==:randomized ? string(VERSION) : nothing,
        requested_noise_samples=actual_method==:randomized ? Int(noise_samples) : nothing,
        noise_factor_storage=noise isa Diagonal ? :diagonal : :dense_supplied_covariance,
        sampling_error_reported=actual_method==:randomized,monte_carlo_error_added_to_covariance=false,
        sampling_error_scope=:covariance_monte_carlo_standard_error,
        linearization=:first_order_delta_method,
        observation_covariance=observation_covariance_s2!==nothing ? :supplied_full :
            effective_covariance!==nothing ? :catalog_ct_difference :
            assume_independent ? :independent_sigma : :unspecified,
        covariance_order=:differentials_then_absolute_picks,initial_catalog_fixed=true,
        initial_catalog_errors_included=false,reference_errors_included=false,velocity_model_fixed=true,
        station_coordinates_fixed=true,phase_assignment_fixed=true,station_delays_fixed=true,
        continuous_weights_recomputed=true,discrete_path_conditioned=true,residual_rescaling=false,
        finite_iteration_bias_included=false,truncation_bias_included=false,model_bias_included=false,
        real_accuracy_claim=false,confidence_distribution=:gaussian,
        confidence_scope=:separate_marginal_algorithm_fluctuations,confidence_probability=probability,
        confidence_covers_true_coordinates=false,location_refitted=true)
    summaries=Dict(component.component_id=>component for component in result.components)
    components=unique([final[id].component_id for id in requested if
        get(final[id],:selection_status,:selected)==:selected && haskey(summaries,final[id].component_id)])
    reasons=Dict{Int,Union{Nothing,Symbol}}()
    signatures=Dict{Int,Any}()
    baseline_position_error,baseline_origin_error=0.,0.
    replay_count=0
    progress===nothing || progress((;stage=:relative_replay_baseline,processed=0,total=1,status=:started))
    baseline=relocate(model,stations,catalog,data;absolute_picks=arrivals,relocate_options...)
    replay_count+=1
    progress===nothing || progress((;stage=:relative_replay_baseline,processed=1,total=1,status=:completed))
    baseline_events=Dict(String(event.event_id)=>event for event in baseline.events)
    baseline_components=Dict(component.component_id=>component for component in baseline.components)
    same_observations=length(result.residuals)==length(data) && length(result.absolute_residuals)==length(arrivals) &&
        [(r.event_i,r.event_j,r.station_id,r.phase,r.dt_s,r.sigma_s,r.source) for r in result.residuals]==
        [(r.event_i,r.event_j,r.station_id,r.phase,r.dt_s,r.sigma_s,r.source) for r in baseline.residuals] &&
        [(r.event_id,r.station_id,r.phase,r.arrival_s,r.sigma_s) for r in result.absolute_residuals]==
        [(r.event_id,r.station_id,r.phase,r.arrival_s,r.sigma_s) for r in baseline.absolute_residuals] &&
        [[get(r,key,nothing) for key in (:original_quality,:quality,:quality_source)] for r in result.residuals]==
        [[get(r,key,nothing) for key in (:original_quality,:quality,:quality_source)] for r in baseline.residuals] &&
        isequal([r.equation_weight for r in result.residuals],[r.equation_weight for r in baseline.residuals]) &&
        isequal([r.equation_weight for r in result.absolute_residuals],[r.equation_weight for r in baseline.absolute_residuals])
    same_catalog=Set(keys(final))==Set(keys(baseline_events))
    same_initial=same_catalog && all(id->all(key->isequal(get(final[id],key,nothing),
        get(baseline_events[id],key,nothing)),(:catalog_source_km,:catalog_origin_s,:initial_source_km,:initial_origin_s)),keys(final))
    if same_catalog
        baseline_position_error=maximum((maximum(abs.(final[id].source_km.-baseline_events[id].source_km)) for id in keys(final));init=0.)
        baseline_origin_error=maximum((abs(final[id].origin_s-baseline_events[id].origin_s) for id in keys(final));init=0.)
    end
    for cid in components
        summary=summaries[cid]
        reason=if !(supplied || assume_independent)
            :observation_covariance_unspecified
        elseif get(summary,:solver,:damped)!=:damped
            :recorded_solver_is_not_damped
        elseif !same_initial || !same_observations || baseline_position_error>2e-10 || baseline_origin_error>2e-11 ||
               !haskey(baseline_components,cid) || result.depth_bounds!=baseline.depth_bounds
            :baseline_replay_mismatch
        elseif get(summary,:secant_gradient_count,0)>0
            :nondifferentiable_replay_path
        elseif !(summary.status in (:converged,:converged_regularized,:converged_minimum_norm,
                                    :converged_truncated,:completed_schedule,:max_iterations))
            :component_has_no_regular_returned_estimate
        elseif get(summary,:successful_iterations,0)==0 && !(summary.status in
                (:converged,:converged_regularized,:converged_minimum_norm,:converged_truncated))
            :no_successful_updates
        elseif any(row->row.accepted && row.step_fraction!=1.,summary.trace) ||
               any(event->event.component_id==cid && event.depth_at_bound && !(event.event_id in summary.reference_event_ids),result.events)
            :active_depth_bound
        else
            nothing
        end
        if reason===nothing
            signature=_dd_algorithm_signature(model,stations,result,cid)
            reference_signature=_dd_algorithm_signature(model,stations,baseline,cid)
            reason=signature===nothing ? :nondifferentiable_replay_path :
                signature!=reference_signature || !isequal(summary.final_equation_weights,baseline_components[cid].final_equation_weights) ||
                !isequal([r.row_weights for r in summary.trace],[r.row_weights for r in baseline_components[cid].trace]) ?
                :baseline_replay_mismatch : nothing
            signatures[cid]=signature
        end
        reasons[cid]=reason
    end
    sensitivity=Dict(id=>zeros(4,directions) for id in requested)
    numerical_error=Dict(id=>zeros(4) for id in requested)
    numerical_ratio=Dict(id=>0. for id in requested)
    completed_directions=0
    # Replays consume these working vectors before the next direction. Keep
    # inactive covariance columns zero without allocating per noise sample.
    needs_directions=directions>0 && any(isnothing,values(reasons))
    direction=needs_directions ? zeros(n) : Float64[]
    coefficients=needs_directions && actual_method==:randomized && !(noise isa Diagonal) ?
        zeros(size(noise,2)) : Float64[]
    progress===nothing || progress((;stage=:relative_replay_directions,processed=0,total=directions,status=:started))
    for j in 1:directions
        any(isnothing,values(reasons)) || break
        fill!(direction,0.)
        if actual_method==:exact
            column=active_noise_columns[j]
            if noise isa Diagonal
                direction[column]=noise.diag[column]
            else
                direction.=view(noise,:,column)
            end
        elseif noise isa Diagonal
            for i in active_noise_columns
                direction[i]=noise.diag[i]*(rand(rng,Bool) ? 1. : -1.)
            end
        else
            fill!(coefficients,0.)
            for i in active_noise_columns
                coefficients[i]=rand(rng,Bool) ? 1. : -1.
            end
            mul!(direction,noise,coefficients)
        end
        # Each complete rerun with perturbed observations includes the full
        # catalog and observation history.
        # Keep only requested coordinates and path comparisons between runs.
        values_by_event=Dict(id=>Vector{Float64}[] for id in requested)
        stable_path=Dict(cid=>true for cid in components if reasons[cid]===nothing)
        multipliers=(Float64(step),-Float64(step),Float64(step)/2,-Float64(step)/2)
        perturbations=_ordered_tasks(multipliers;tasks=replay_tasks) do multiplier
            shifted_data=[merge(row,(dt_s=row.dt_s+multiplier*direction[i],)) for (i,row) in enumerate(data)]
            shifted_arrivals=[merge(row,(arrival_s=row.arrival_s+multiplier*direction[length(data)+i],)) for (i,row) in enumerate(arrivals)]
            replay=relocate(model,stations,catalog,shifted_data;absolute_picks=shifted_arrivals,perturbation_options...)
            stable=Dict(cid=>_dd_algorithm_signature(model,stations,replay,cid)==signatures[cid] for cid in keys(stable_path))
            replay_events=Dict(event.event_id=>event for event in replay.events)
            values=Dict(id=>[replay_events[id].source_km...,replay_events[id].origin_s] for id in requested)
            (;stable,values)
        end
        replay_count+=length(perturbations)
        for perturbation in perturbations
            for cid in keys(stable_path)
                stable_path[cid] &= perturbation.stable[cid]
            end
            for id in requested
                push!(values_by_event[id],perturbation.values[id])
            end
        end
        completed_directions+=1
        for cid in components
            reasons[cid]===nothing || continue
            if !stable_path[cid]
                reasons[cid]=:unstable_discrete_replay_path
                continue
            end
            for id in requested
                final[id].component_id==cid || continue
                values_by_run=values_by_event[id]
                coarse=(values_by_run[1]-values_by_run[2])/(2step)
                fine=(values_by_run[3]-values_by_run[4])/step
                richardson=(4fine-coarse)/3
                error_estimate=abs.(fine-coarse)/3
                budget=[fill(Float64(derivative_position_atol_km),3);Float64(derivative_origin_atol_s)].+
                    Float64(derivative_rtol).*abs.(richardson)
                sensitivity[id][:,j]=richardson
                numerical_error[id].=max.(numerical_error[id],error_estimate)
                numerical_ratio[id]=max(numerical_ratio[id],maximum(error_estimate./budget))
                all(isfinite,richardson) && all(error_estimate.<=budget) ||
                    (reasons[cid]=:unstable_numerical_derivative)
            end
        end
        if progress!==nothing && completed_directions<directions && any(isnothing,values(reasons))
            progress((;stage=:relative_replay_directions,processed=completed_directions,total=directions,status=:running))
        end
    end
    progress===nothing || progress((;stage=:relative_replay_directions,processed=completed_directions,total=directions,
        status=completed_directions==directions ? :completed : :stopped))
    rows=NamedTuple[]
    for id in requested
        event=final[id]
        summary=get(summaries,event.component_id,nothing)
        references=summary===nothing ? String[] : summary.reference_event_ids
        reason=get(event,:selection_status,:selected)!=:selected ? :queried_event_not_selected :
            get(event,:exit_reason,nothing)!==nothing ? :queried_event_exited :
            id in references ? :reference_fixed : get(reasons,event.component_id,:no_component_uncertainty)
        factor=sensitivity[id]
        sampling_check=_dd_algorithm_sampling_check(factor,actual_method,completed_directions,reason===nothing)
        normalization=actual_method==:randomized && completed_directions>0 ? completed_directions : 1
        covariance=reason===nothing ? factor*factor'/normalization : nothing
        standard=reason===nothing ? sqrt.(max.(diag(covariance),0.)) : fill(nothing,4)
        se=(x_km=standard[1],y_km=standard[2],z_km=standard[3],origin_s=standard[4])
        halfwidths=map(value->value===nothing ? nothing : scale*value,se)
        numerical_check=(step=Float64(step),half_step=Float64(step)/2,
            derivative_rtol=Float64(derivative_rtol),position_atol_km=Float64(derivative_position_atol_km),
            origin_atol_s=Float64(derivative_origin_atol_s),maximum_error_budget_ratio=numerical_ratio[id],
            maximum_position_factor_error_km=maximum(numerical_error[id][1:3];init=0.),
            maximum_origin_factor_error_s=numerical_error[id][4],
            checked=reason===nothing,passed=reason===nothing)
        push!(rows,(;event_id=id,component_id=event.component_id,source_km=event.source_km,origin_s=event.origin_s,
            status=reason===nothing ? :estimated_algorithm_conditional : :not_estimated,reason,
            solver_status=summary===nothing ? nothing : summary.status,
            parameter_names=(:x_km,:y_km,:z_km,:origin_s),covariance,standard_deviations=se,
            marginal_confidence_half_widths=halfwidths,confidence_probability=probability,confidence_scale=scale,
            rank=summary===nothing ? nothing : summary.retained_rank,reference_event_ids=references,numerical_check,sampling_check,
            assumptions=merge(assumptions,(solver=summary===nothing ? nothing : summary.solver,
                datum=summary===nothing ? :unavailable : summary.datum))))
    end
    numerical_check=(step=Float64(step),half_step=Float64(step)/2,
        baseline_position_error_km=baseline_position_error,baseline_origin_error_s=baseline_origin_error,
        baseline_position_tolerance_km=2e-10,baseline_origin_tolerance_s=2e-11,
        derivative_rtol=Float64(derivative_rtol),derivative_position_atol_km=Float64(derivative_position_atol_km),
        derivative_origin_atol_s=Float64(derivative_origin_atol_s),method=actual_method,
        noise_factor_rank=noise_rank,completed_noise_directions=completed_directions,
        requested_noise_samples=actual_method==:randomized ? Int(noise_samples) : nothing,
        actual_noise_samples=actual_method==:randomized ? completed_directions : 0,
        seed=actual_method==:randomized ? seed : nothing)
    (;events=rows,observation_count=n,noise_direction_count=directions,replay_count,assumptions,numerical_check)
end

# The random directional mean is zero by construction. Estimate its second
# moment, and the Monte Carlo SE of that second moment, without centering q.
function _dd_algorithm_sampling_check(factor,method,samples,available)
    if method==:exact
        return (;method,actual_samples=0,checked=available,
            covariance_element_standard_error=nothing,diagonal_relative_standard_error=nothing,
            two_half_covariance_relative_difference=nothing,
            error_scope=:no_random_sampling,accuracy_bound=false)
    end
    if !available
        return (;method,actual_samples=samples,checked=false,
            covariance_element_standard_error=nothing,diagonal_relative_standard_error=nothing,
            two_half_covariance_relative_difference=nothing,
            error_scope=:covariance_monte_carlo_standard_error,accuracy_bound=false)
    end
    samples==0 && return (;method,actual_samples=0,checked=true,
        covariance_element_standard_error=zeros(4,4),diagonal_relative_standard_error=fill(nothing,4),
        two_half_covariance_relative_difference=nothing,
        error_scope=:known_zero_observation_covariance,accuracy_bound=false)
    covariance=factor*factor'/samples
    standard_error=[sqrt(sum((factor[a,j]*factor[b,j]-covariance[a,b])^2 for j in 1:samples)/
        (samples*(samples-1))) for a in 1:4,b in 1:4]
    relative=[covariance[i,i]>0 ? standard_error[i,i]/covariance[i,i] : nothing for i in 1:4]
    half=fld(samples,2)
    left=view(factor,:,1:half)
    right=view(factor,:,half+1:samples)
    difference=norm(left*left'/half-right*right'/(samples-half))
    denominator=norm(covariance)
    (;method,actual_samples=samples,checked=true,covariance_element_standard_error=standard_error,
        diagonal_relative_standard_error=relative,
        two_half_covariance_relative_difference=denominator>0 ? difference/denominator : 0.,
        error_scope=:covariance_monte_carlo_standard_error,accuracy_bound=false)
end

# Discrete computation choices are compared, while numerical weights/scales
# are allowed to vary. Ray identities are checked along the complete saved
# path; zero-gradient/secant endpoints are never assigned Gaussian errors.
function _dd_algorithm_signature(model,stations,result,cid)
    component=findfirst(c->c.component_id==cid,result.components)
    component===nothing && return nothing
    summary=result.components[component]
    summary.secant_gradient_count>0 && return nothing
    members=[event for event in result.events if event.component_id==cid]
    index=Dict(event.event_id=>i for (i,event) in enumerate(members))
    references=Set(summary.reference_event_ids)
    rays=sort!(unique([(String(id),String(row.station_id),String(row.phase))
        for row in [result.residuals;result.absolute_residuals] if row.component_id==cid &&
        get(row,:selection_status,:selected)==:selected
        for id in (hasproperty(row,:event_id) ? (row.event_id,) : (row.event_i,row.event_j)) if !(id in references)]))
    ray_identity(ray)=(get(ray,:branch,nothing),get(ray,:refractor_layer,nothing),get(ray,:turning_layer,nothing),
        get(ray,:region_status,nothing),hasproperty(ray,:components) ?
            [(c.model_id,ray_identity(c.prediction)) for c in ray.components] : nothing)
    initial=reduce(vcat,[reshape(Float64[event.initial_source_km...,event.initial_origin_s],1,4) for event in members])
    ray_path=Any[]
    for parameters in (initial,(row.parameters_after for row in summary.trace)...)
        identities=Any[]
        for (id,station,phase) in rays
            prediction=try
                travel_time(model,Tuple(parameters[index[id],1:3]),stations,station,phase)
            catch exception
                exception isa DomainError || rethrow()
                return nothing
            end
            prediction.dsource_s_km===nothing && return nothing
            push!(identities,ray_identity(prediction))
        end
        push!(ray_path,identities)
    end
    weights=[(row.attempt,row.input_row,row.row_kind,get(row,:source,nothing),row.removal_reason,
        row.empirical_base_weight===nothing ? nothing : row.empirical_base_weight*row.distance_weight>1e-6)
        for row in summary.weight_trace]
    # MAD is piecewise smooth: retain its actual ordering/absolute-value branch.
    mad_path=Any[]
    for attempt in unique([row.attempt for row in summary.weight_trace]),source in ("CC","CT")
        rows=[row for row in summary.weight_trace if row.attempt==attempt && row.source==source]
        isempty(rows) && continue
        scale=source=="CC" ? rows[1].residual_scale_cc : rows[1].residual_scale_ct
        scale.status==:zero_mad && return nothing
        if scale.status==:mad
            values=[row.residual_s for row in rows]
            centered=values.-_dd_weight_median(values)
            push!(mad_path,(attempt,source,sortperm(values),sign.(centered),sortperm(abs.(centered))))
        end
    end
    (;solver=summary.solver,row_weighting=summary.row_weighting,damp=summary.damp,
        cutoff=summary.svd_relative_cutoff,schedule=summary.weight_schedule,
        status=summary.status,stop_reason=summary.stop_reason,reference_event_ids=summary.reference_event_ids,
        parameter_keys=summary.parameter_keys,retained=summary.diagnostics===nothing ? nothing : summary.diagnostics.retained,
        events=[(e.event_id,e.status,e.exit_reason,get(e,:selection_status,:selected),e.depth_at_bound) for e in members],
        rows=[(get(r,:input_row,0),r.used,r.unused_reason,get(r,:selection_status,:selected))
            for r in [result.residuals;result.absolute_residuals] if r.component_id==cid],
        trace=[(r.accepted,r.reason,r.iteration,r.stage_index,r.applied_damp,r.retained_rank,r.solver_rank,
            r.step_fraction==1.,r.airquake_event_ids,r.removed_or_disconnected_ids,r.removed_by_weighting_ids,
            r.active_station_ids) for r in summary.trace],weights,mad_path,ray_path)
end
