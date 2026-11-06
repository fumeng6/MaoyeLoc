# Covariance of the full-rank, frozen-weight LOCAL least-squares problem.
# This is not the covariance of one spectrally damped iteration step.
function _absolute_frozen_weight_covariance(design,weights,variance_s2,noise_model)
    A=design.*weights
    decomposition,rank=_absolute_svd(A)
    rank==size(design,2) || return (covariance=nothing,rank)
    inverse=decomposition.V*Diagonal(inv.(decomposition.S))*decomposition.U'
    map=noise_model==:weighted_rows ? inverse : inverse.*permutedims(weights)
    covariance=variance_s2.*(map*map')
    (;covariance=(covariance+covariance')/2,rank)
end

"""
    absolute_formal_uncertainty(model,stations,picks,location;
        rderr_s=.15,ercof=1.,noise_model=:weighted_rows)

Hypoinverse's explicitly requested empirical timing variance
`q = RDERR^2 + ERCOF*RMS^2`, in s², and a formal covariance of the local
least-squares problem with the final fitted weights held fixed. This does not change the fitted location
or its sigma metadata. RMS² is sum(w²*r²)/sum(w²), without an n-p correction.
Reading identities and their event, station, phase, arrival and sigma must
match the fitted records; input order may change without merging readings.

`:weighted_rows` declares Cov(W*epsilon)=q*I, hence original-arrival variance
q/w² on selected positive-weight rows. `:arrival_rows` instead declares
Cov(epsilon)=q*I and uses the complete weighted-row sandwich. The weights are
the saved empirical equation factors.

The quantity is `:linearized_frozen_weight_least_squares`, conditional on q,
weights, selected rows, the model and branch. This Gauss-Newton approximation
omits residual-times-Hessian terms and is not the exact nonlinear refit
derivative when residual curvature is appreciable. It does not propagate stochastic
weight selection, estimation of q, nonlinear bias, or model uncertainty.
Spectral shift and iteration damping control parameter updates; they do not
replace the local least-squares map from data errors to parameter errors when
propagating covariance. An explicit legacy-rule normal stop is eligible for
this local frozen-weight calculation; its status remains recorded and does
not establish convergence under the stricter undamped-step criterion.
Incomplete/truncated solutions,
active inequality bounds and nonsmooth fitted rays do not receive a full
parameter covariance. Fixed parameters have missing standard deviations.
"""
function absolute_formal_uncertainty(model,stations::AbstractDict,picks,location;
        rderr_s=.15,ercof=1.,noise_model=:weighted_rows)
    all(x->x isa Real && isfinite(x) && x>=0,(rderr_s,ercof)) ||
        throw(ArgumentError("RDERR seconds and ERCOF must be finite and nonnegative."))
    noise_model in (:weighted_rows,:arrival_rows) || throw(ArgumentError("Noise model must be weighted_rows or arrival_rows."))
    observations=collect(picks)
    n=length(observations)
    n==location.input_count==length(location.residuals) || throw(ArgumentError("Picks must match the fitted absolute observation count."))
    fitted=Dict(_absolute_observation_key(row)=>row for row in location.residuals)
    length(fitted)==n || throw(ArgumentError("Fitted absolute residual identities must be unique."))
    event_ids=unique(pick.event_id for pick in observations if hasproperty(pick,:event_id))
    length(event_ids)<=1 && (isempty(event_ids) || location.event_id===nothing || only(event_ids)==location.event_id) ||
        throw(ArgumentError("Formal-error picks must belong to the fitted event."))
    constraint=get(location,:constraint,:free)
    fixed_position=constraint in (:fixed_position,:fixed_all)
    fixed_depth=fixed_position || constraint==:fixed_depth || location.depth_bounds[1]==location.depth_bounds[2]
    axes=constraint==:fixed_all ? () : fixed_position ? (4,) : fixed_depth ? (1,2,4) : (1,2,3,4)
    parameter_names=Tuple((:x_km,:y_km,:z_km,:origin_s)[i] for i in axes)
    get(location,:row_weighting,:empirical)==:empirical ||
        throw(ArgumentError("Formal uncertainty uses the saved empirical equation weights."))
    weights=Float64[]
    residual=Float64[]
    seen=Set{Tuple}()
    for pick in observations
        key=_absolute_observation_key(pick)
        key in seen && throw(ArgumentError("Duplicate absolute observation identity in formal-error input."))
        push!(seen,key)
        original=get(fitted,key,nothing)
        original!==nothing && String(pick.station_id)==String(original.station_id) &&
            String(pick.phase)==String(original.phase) &&
            pick.arrival_s==original.arrival_s && pick.sigma_s==original.sigma_s &&
            get(pick,:event_id,nothing)==get(original,:event_id,nothing) ||
            throw(ArgumentError("Observation identity, arrival and sigma must match the fitted absolute solution."))
        haskey(stations,pick.station_id) && pick.phase in ("P","S") && isfinite(pick.sigma_s) && pick.sigma_s>0 ||
            throw(ArgumentError("Formal-error phases need known stations and positive finite sigma metadata."))
        weight=original.equation_weight
        isfinite(weight) && weight>=0 && isfinite(original.residual_s) || throw(ArgumentError("Saved row weights and residuals must be finite."))
        push!(weights,weight)
        push!(residual,original.residual_s)
    end
    normalization=1.
    selected=findall(>(0),weights)
    denominator=sum(abs2,weights)
    rms=denominator==0 ? nothing : sqrt(sum(abs2,weights.*residual)/denominator)
    variance=rms===nothing ? (ercof==0 ? Float64(rderr_s)^2 : nothing) : Float64(rderr_s)^2+Float64(ercof)*rms^2
    assumptions=(quantity=:linearized_frozen_weight_least_squares,
        fitted_estimator=:empirical_absolute_rows,
        selected_search_minimum_conditioned=false,
        observation_noise_model=noise_model,variance_reference=noise_model==:weighted_rows ? :normalized_weighted_arrivals : :original_arrivals,
        scale_fixed=true,weights_fixed=true,selected_rows_fixed=true,branch_fixed=true,
        velocity_model_fixed=true,station_coordinates_fixed=true,station_delays_fixed=true,
        phase_assignment_fixed=true,first_order=true,residual_curvature_ignored=true,scale_estimation_uncertainty_included=false,
        weight_selection_uncertainty_included=false,model_uncertainty_included=false,
        sigma_metadata_unchanged=true,location_refitted=false,real_accuracy_claim=false,
        location_status=location.status,
        spectral_shift_role=:iteration_control,fixed_parameters_conditional=fixed_depth,
        covariance_order=:input_picks_order)
    base=(quantity=:linearized_frozen_weight_least_squares,parameter_names,
        variance_s2=variance,timing_error_s=variance===nothing ? nothing : sqrt(variance),
        weighted_rms_s=rms,rderr_s=Float64(rderr_s),ercof=Float64(ercof),noise_model,
        equation_weights=weights,row_normalization=normalization,included_mask=weights.>0,
        observation_variance_s2=[variance===nothing || w==0 ? nothing : noise_model==:weighted_rows ? variance/w^2 : variance for w in weights],
        observation_count=n,included_count=length(selected),assumptions)
    unavailable(reason;status=:not_estimated,rank=nothing)=merge(base,(;status,reason,covariance=nothing,rank,
        standard_deviations=(x_km=nothing,y_km=nothing,z_km=nothing,origin_s=nothing)))
    constraint==:fixed_all && return unavailable(:all_parameters_fixed;status=:scale_only)
    isempty(selected) && return unavailable(:no_positive_weight_observations;status=:scale_only)
    get(location,:weak_direction_truncation,false) && return unavailable(:weak_direction_truncation;rank=get(location,:retained_rank,nothing))
    location.status==:rank_deficient && return unavailable(:rank_deficient;rank=get(location,:free_parameter_rank,nothing))
    location.status in (:converged,:converged_at_depth_bound,:stopped_by_legacy_rule) || return unavailable(:location_not_converged)
    !fixed_depth && (get(location,:depth_at_bound,false) || get(location,:depth_constraint_active,false) ||
        location.source_km[3]==_source_top(model,location.source_km)) && return unavailable(:active_depth_bound)
    length(selected)>=length(axes) || return unavailable(:insufficient_observations)
    design=zeros(length(selected),length(axes))
    for (r,i) in enumerate(selected)
        pick=observations[i]
        if fixed_position
            design[r,1]=1.
        else
            ray=try
                travel_time(model,location.source_km,stations,pick.station_id,pick.phase)
            catch error
                error isa DomainError || rethrow()
                nothing
            end
            (ray===nothing || ray.dsource_s_km===nothing) && return unavailable(:nondifferentiable_travel_time)
            gradient=[ray.dsource_s_km...;1.]
            design[r,:] .= gradient[collect(axes)]
        end
    end
    propagated=_absolute_frozen_weight_covariance(design,weights[selected],variance,noise_model)
    propagated.covariance===nothing && return unavailable(:rank_deficient;rank=propagated.rank)
    standard=sqrt.(diag(propagated.covariance))
    values=Union{Nothing,Float64}[nothing,nothing,nothing,nothing]
    values[collect(axes)] .= standard
    merge(base,(status=:estimated_formal_conditional,reason=nothing,covariance=propagated.covariance,rank=propagated.rank,
        standard_deviations=(x_km=values[1],y_km=values[2],z_km=values[3],origin_s=values[4])))
end
