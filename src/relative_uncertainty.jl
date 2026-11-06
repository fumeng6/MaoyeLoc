# A catalog CT that is exactly one unique pair's arrival_i-arrival_j shares
# those pick errors. CC and any CT that does not match within 1e-9 s stay
# independent. Ambiguous multi-pick matches are not guessed. Detect these links
# separately: reporting an unavailable parameter covariance needs no dense C.
function _catalog_ct_links(differentials,absolute_picks)
    links = Tuple{Int,Int,Int}[]
    isempty(absolute_picks) && return links
    buckets = Dict{Tuple{String,String,String},Vector{Int}}()
    for (k,pick) in enumerate(absolute_picks)
        key = (String(pick.event_id),String(pick.station_id),String(pick.phase))
        push!(get!(buckets,key,Int[]),k)
    end
    for (r,row) in enumerate(differentials)
        String(row.source)=="CT" || continue
        left = get(buckets,(String(row.event_i),String(row.station_id),String(row.phase)),Int[])
        right = get(buckets,(String(row.event_j),String(row.station_id),String(row.phase)),Int[])
        (isempty(left) || isempty(right)) && continue
        matches = Tuple{Int,Int}[]
        for i in left, j in right
            abs(row.dt_s-(absolute_picks[i].arrival_s-absolute_picks[j].arrival_s))<=1e-9 &&
                push!(matches,(i,j))
        end
        length(matches)==1 || continue
        push!(links,(r,matches[1]...))
    end
    links
end

# The returned covariance is ordered differentials then absolute picks.
# Nothing means no CT shares an absolute pick's error.
function _catalog_ct_covariance_s2(differentials,absolute_picks)
    links = _catalog_ct_links(differentials,absolute_picks)
    isempty(links) && return nothing
    n,m = length(differentials),length(absolute_picks)
    linked = Dict(r=>(i,j) for (r,i,j) in links)
    nsource = m+n-length(linked)
    mapping = zeros(n+m,nsource)
    variance = zeros(nsource)
    for k in 1:m
        mapping[n+k,k] = 1
        variance[k] = absolute_picks[k].sigma_s^2
    end
    source = m
    for r in 1:n
        if haskey(linked,r)
            i,j = linked[r]
            mapping[r,i] = 1
            mapping[r,j] = -1
        else
            source += 1
            mapping[r,source] = 1
            variance[source] = differentials[r].sigma_s^2
        end
    end
    factor = mapping*Diagonal(sqrt.(variance))
    factor*factor'
end

"""
    relative_uncertainty(model, stations, events, differentials, result;
        event_pairs=[], event_ids=[], absolute_picks=[], assume_independent=false,
        observation_covariance_s2=nothing, event_confidence_probability=.95)

Report that an analytic covariance is not estimated for the column-scaled
damped relocation. `event_ids` keeps the requested order. A queried event that
was not selected, left the cluster, or is an explicit reference is labeled
with that reason. Cross-component pairs are not estimated. This is not a
parameter covariance and not a claim about absolute accuracy.

A CT row whose `dt_s` equals `arrival_i - arrival_j` for one unique pair of
absolute picks at the same station and phase (tolerance 1e-9 s) shares those
pick errors. `assume_independent=true` is rejected while any such row is
present. A supplied covariance is checked for shape, symmetry and positive
semidefiniteness, in the order of all differentials then all absolute picks.
Neither check is turned into a damped-step covariance. Observation-error
propagation of the finished relocation is `relative_algorithm_uncertainty`.
"""
function relative_uncertainty(model,stations::AbstractDict,events,differentials,result;
                              event_pairs=Tuple{String,String}[],event_ids=String[],absolute_picks=NamedTuple[],
                              assume_independent::Bool=false,observation_covariance_s2=nothing,
                              event_confidence_probability=.95)
    catalog, observations, arrivals = collect(events),collect(differentials),collect(absolute_picks)
    ids = String[event.event_id for event in catalog]
    length(unique(ids))==length(ids) || throw(ArgumentError("Duplicate event ID."))
    Set(ids)==Set(event.event_id for event in result.events) ||
        throw(ArgumentError("Events must match the relocated catalog."))
    final = Dict(event.event_id=>event for event in result.events)
    summaries = Dict(component.component_id=>component for component in result.components)
    length(observations)==length(result.residuals) && length(arrivals)==length(result.absolute_residuals) ||
        throw(ArgumentError("Observations must match the relocation counts."))
    differential_key(row) = (row.event_i,row.event_j,row.station_id,row.phase,row.dt_s,row.sigma_s,row.source)
    absolute_key(row) = (row.event_id,_absolute_observation_key(row),row.station_id,row.phase,row.arrival_s,row.sigma_s)
    fitted_differentials=Dict(differential_key(row)=>row for row in result.residuals)
    fitted_arrivals=Dict(absolute_key(row)=>row for row in result.absolute_residuals)
    selected_event(event)=get(event,:selection_status,:selected)==:selected
    selected_differential(row)=get(fitted_differentials[differential_key(row)],:selection_status,:selected)==:selected
    selected_arrival(row)=get(fitted_arrivals[absolute_key(row)],:selection_status,:selected)==:selected
    Set(differential_key(row) for row in observations)==Set(differential_key(row) for row in result.residuals) &&
        Set(absolute_key(row) for row in arrivals)==Set(absolute_key(row) for row in result.absolute_residuals) ||
        throw(ArgumentError("Observations and sigma weights must match the fitted relocation."))
    for event in result.events
        length(event.source_km)==3 && all(isfinite,event.source_km) && isfinite(event.origin_s) ||
            throw(ArgumentError("Final event coordinates and origin time must be finite."))
        !selected_event(event) || haskey(summaries,event.component_id) || throw(ArgumentError("Unknown result component."))
    end
    for row in observations
        haskey(final,row.event_i) && haskey(final,row.event_j) && row.event_i!=row.event_j &&
            haskey(stations,row.station_id) || throw(ArgumentError("Invalid differential identity."))
        !selected_differential(row) || final[row.event_i].component_id==final[row.event_j].component_id ||
            throw(ArgumentError("Differential spans reported components."))
        row.phase in ("P","S") && row.source in ("CT","CC") && isfinite(row.dt_s) &&
            isfinite(row.sigma_s) && row.sigma_s>0 || throw(ArgumentError("Invalid differential phase/value/sigma."))
    end
    for row in arrivals
        haskey(final,row.event_id) && haskey(stations,row.station_id) ||
            throw(ArgumentError("Invalid absolute-constraint identity."))
        row.phase in ("P","S") && isfinite(row.arrival_s) && isfinite(row.sigma_s) && row.sigma_s>0 ||
            throw(ArgumentError("Invalid absolute-constraint phase/value/sigma."))
    end
    pairs = Tuple{String,String}[]
    for pair in event_pairs
        length(pair)==2 || throw(ArgumentError("An event pair must contain two IDs."))
        i,j = String(pair[1]),String(pair[2])
        haskey(final,i) && haskey(final,j) && i!=j || throw(ArgumentError("Invalid query event pair."))
        push!(pairs,(i,j))
    end
    requested_events=String.(collect(event_ids))
    all(haskey(final,id) for id in requested_events) || throw(ArgumentError("Unknown event in individual uncertainty query."))
    event_confidence_probability isa Real && isfinite(event_confidence_probability) && 0<event_confidence_probability<1 ||
        throw(ArgumentError("Individual marginal confidence probability must lie strictly between zero and one."))
    event_confidence_scale=isempty(requested_events) ? nothing : sqrt(_spatial_chisq_quantile(Float64(event_confidence_probability),1))
    lower,upper = result.depth_bounds
    isfinite(lower) && lower>=model_top(model) && upper>=lower || throw(ArgumentError("Invalid result depth bounds."))
    all(lower<=event.source_km[3]<=upper for event in result.events) ||
        throw(ArgumentError("Final depth outside result bounds."))
    fixed_depth = lower==upper
    parameter_names = fixed_depth ? (:x_km,:y_km,:origin_s) : (:x_km,:y_km,:z_km,:origin_s)
    shared_picks = observation_covariance_s2===nothing && !isempty(_catalog_ct_links(observations,arrivals))
    if shared_picks && assume_independent
        throw(ArgumentError("Catalog CT rows equal arrival_i-arrival_j and share those pick errors. assume_independent=true would understate the formal error. Omit it to use the catalog-difference covariance, or pass observation_covariance_s2. Relocation weights are unchanged."))
    end
    observation_count = length(observations)+length(arrivals)
    # The catalog-difference C is known to be PSD by its construction L*L'.
    # Only a user-supplied C needs validation; no covariance is estimated here.
    observation_covariance_s2===nothing || _observation_covariance_factor(observation_covariance_s2,observation_count)
    scaled_estimators = !isempty(result.components) &&
        all(get(summary,:solver,:damped)==:damped for summary in result.components)
    assumptions = (quantity=isempty(requested_events) ? :event_i_minus_event_j : :query_specific,
        estimator=scaled_estimators ? :component_specific : :unavailable,
        observation_covariance=observation_covariance_s2!==nothing ? :supplied_full :
            shared_picks ? :catalog_ct_difference :
            assume_independent ? :independent_sigma : :unspecified,
        covariance_order=:differentials_then_absolute_picks,location_refitted=false,
        linearization=:first_order,residual_curvature_ignored=true,
        weight_selection_uncertainty_included=false,velocity_model_fixed=true,station_coordinates_fixed=true,
        phase_assignment_fixed=true,first_arrival_branch_fixed=true,
        station_delays_fixed=true,station_delay_errors_included=false,
        reference_errors_included=false,group_center_covariance=:not_computed,
        fixed_depth_conditional=fixed_depth,residual_rescaling=false,
        model_bias_included=false,real_accuracy_claim=false,
        covariance_method=scaled_estimators ? :component_specific : :not_computed)

    requested_components = unique(vcat(
        [final[i].component_id for (i,j) in pairs if selected_event(final[i]) && selected_event(final[j]) && final[i].component_id==final[j].component_id],
        [final[id].component_id for id in requested_events if selected_event(final[id])]))
    analyses = Dict{Int,NamedTuple}()
    for component_id in requested_components
        summary = summaries[component_id]
        members = [final[id] for id in ids if final[id].component_id==component_id]
        index = Dict(event.event_id=>i for (i,event) in enumerate(members))
        references = Set(summary.reference_event_ids)
        all(haskey(index,id) for id in references) || throw(ArgumentError("Reference outside its component."))
        solver = get(summary,:solver,:damped)
        solver==:damped || throw(ArgumentError("Relative uncertainty describes the column-scaled damped relocation."))
        analyses[component_id]=(;reason=:scaled_damping_covariance_not_implemented,rank=nothing,
            references=String.(summary.reference_event_ids),solver)
    end

    pair_results = NamedTuple[]
    for (i,j) in pairs
        component_i,component_j = final[i].component_id,final[j].component_id
        component_id = component_i==component_j ? component_i : nothing
        difference = (x_km=final[i].source_km[1]-final[j].source_km[1],
            y_km=final[i].source_km[2]-final[j].source_km[2],
            z_km=final[i].source_km[3]-final[j].source_km[3],
            origin_s=final[i].origin_s-final[j].origin_s)
        query_exited = any(final[id].status in (:airquake_removed,:disconnected_after_airquake) for id in (i,j))
        row = (;event_i=i,event_j=j,component_id,component_i,component_j,difference,
            status=:not_estimated,reason=query_exited ? :queried_event_exited : :disconnected_components,parameter_names,
            covariance=nothing,standard_deviations=(x_km=nothing,y_km=nothing,z_km=nothing,origin_s=nothing),
            rank=nothing,reference_event_ids=String[],assumptions=merge(assumptions,(quantity=:event_i_minus_event_j,)))
        if !selected_event(final[i]) || !selected_event(final[j])
            push!(pair_results,merge(row,(reason=:queried_event_not_selected,)))
            continue
        end
        if component_id===nothing
            push!(pair_results,row)
            continue
        end
        state = analyses[component_id]
        reason = query_exited ?
            :queried_event_exited : i in state.references && j in state.references ? :both_reference_fixed : state.reason
        pair_assumptions = merge(assumptions,(estimator=:scaled_tikhonov,
            quantity=:event_i_minus_event_j,covariance_method=:not_computed,solver=state.solver))
        push!(pair_results,merge(row,(;reason,rank=state.rank,reference_event_ids=state.references,assumptions=pair_assumptions)))
    end

    event_results=NamedTuple[]
    for id in requested_events
        event=final[id]
        component_id=event.component_id
        summary=get(summaries,component_id,nothing)
        state=get(analyses,component_id,nothing)
        solver=state===nothing ? nothing : state.solver
        references=state===nothing ? String[] : state.references
        event_assumptions=merge(assumptions,(quantity=:event_conditional_coordinates,
            datum=summary===nothing ? :unavailable : get(summary,:datum,:unknown),
            estimator=solver==:damped ? :scaled_tikhonov : :unavailable,
            covariance_method=:not_computed,
            confidence_distribution=:gaussian,confidence_scope=:separate_marginal_axes,
            confidence_probability=Float64(event_confidence_probability),
            reference_values_conditioned=!isempty(references),reference_true_errors_included=false))
        reason=if !selected_event(event)
            :queried_event_not_selected
        elseif event.status in (:airquake_removed,:disconnected_after_airquake,:low_weight_removed,:disconnected_after_weighting)
            :queried_event_exited
        elseif id in references
            :reference_fixed
        elseif state===nothing
            :no_component_uncertainty
        else
            state.reason
        end
        row=(;event_id=id,component_id,source_km=event.source_km,origin_s=event.origin_s,
            status=:not_estimated,reason,parameter_names,covariance=nothing,
            standard_deviations=(x_km=nothing,y_km=nothing,z_km=nothing,origin_s=nothing),
            marginal_confidence_half_widths=(x_km=nothing,y_km=nothing,z_km=nothing,origin_s=nothing),
            confidence_probability=Float64(event_confidence_probability),confidence_scale=event_confidence_scale,
            rank=state===nothing ? nothing : state.rank,reference_event_ids=references,assumptions=event_assumptions)
        push!(event_results,row)
    end
    return (;pairs=pair_results,events=event_results,observation_count,assumptions)
end
