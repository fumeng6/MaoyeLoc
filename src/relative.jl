function _validate_relative_controls(;solver=:damped,linear_solver=:auto,lsqr_rtol=1e-6,
        lsqr_max_iterations=0,lsqr_condition_limit=1e12,row_weighting=:empirical,
        empirical_weighting=nothing,airquake_policy=:bounds,airquake_depth_km=0.,
        damp=0.,svd_relative_cutoff=1e-7,row_weights=nothing)
    solver==:damped || throw(ArgumentError("Relative location uses one column-scaled damped path. qr and the sum-penalty svd path are not available."))
    linear_solver in (:auto,:dense,:lsqr) || throw(ArgumentError("linear_solver must be auto, dense or lsqr."))
    isfinite(lsqr_rtol) && 0<lsqr_rtol<1 && lsqr_max_iterations>=0 &&
        isfinite(lsqr_condition_limit) && lsqr_condition_limit>1 || throw(ArgumentError("Invalid LSQR controls."))
    row_weighting==:empirical || throw(ArgumentError("Relative location reweights rows by data class, distance and residual. sigma, unit and explicit row weights are not a separate path."))
    row_weights===nothing || throw(ArgumentError("Explicit row_weights are not a separate relative-location path."))
    airquake_policy in (:bounds,:exclude_repeat) && isfinite(airquake_depth_km) ||
        throw(ArgumentError("Use airquake_policy=:bounds/:exclude_repeat and a finite depth zero plane."))
    isfinite(damp) && damp>=0 && isfinite(svd_relative_cutoff) && 0<=svd_relative_cutoff<=1 ||
        throw(ArgumentError("Invalid DD DAMP or relative singular-value cutoff."))
    nothing
end

"""
    relocate(model, stations, events, differentials; reference_event_ids=[], ...)

Estimate event positions and origin times from differences between arrivals
at the same station and phase. CT denotes differences formed from catalog
picks; CC denotes differences supplied by an external waveform
cross-correlation measurement. This routine does not measure CC from waveforms.
Events
provide event_id, source_km and origin_s; differences provide event_i, event_j,
station_id, phase, dt_s, sigma_s and source (CT/CC). dt_s is the COMPLETE arrival
difference i-j, including origin times. The only path reweights those rows by
data class, inter-event distance and residual, column-scales the four
corrections per event, and damps that scaled step. A small system uses the
dense SVD of this same damped problem; a large system uses LSQR. No event is
fixed automatically. sigma_s remains the recorded differential scatter and is
not the equation weight.

Optional `absolute_picks` add arrival rows with event_id, station_id, phase,
arrival_s and sigma_s. A component with such rows needs no automatic reference:
absolute data determine its overall position and common origin-time shift
if the full design has sufficient rank.
This is a diagonal weighted objective; shared-pick or CT/CC error covariance
is not modeled, so row counts do not imply independent information.

Each connected component starts from its supplied catalog (the default) or,
with initialization=:centroid, one equal-weight spatial centroid. An explicit
initialization_frame selects geographic latitude/longitude/depth averaging;
otherwise coordinates are averaged in local km. Per-event origins and full
observed differences are unchanged. Original and actual initial values are
returned; a fixed reference event retains its position and origin time AFTER
this initialization.

Each connected component is solved separately. Isolated events stay in the
result with status no_observations. No event is fixed unless its ID is listed
in reference_event_ids.

An explicit `selection` table enables IDAT/IPHA/ID, DIST and OBSCC/OBSCT/CID.
DIST requires `selection_frame` for geographic center-to-station distances.
A positive OBSCC+OBSCT threshold uses strong links to name clusters; all
eligible intracluster rows are then reread at that cluster's center. A zero
threshold does not merge disconnected events: each connected piece of the
retained differential graph is its own component. CID
preserves the complete cluster numbering.
Unselected events and rows keep their original identities and exclusion
reasons. Selection applies to differential rows; explicitly supplied absolute
arrivals remain additional constraints for the selected cluster members.
The linear step updates all four parameters of each free event together.
Columns are scaled by their RMS, then damped. damp=0 is the minimum-norm
update in that scaled space. linear_solver=:auto uses dense SVD for a small
matrix and LSQR for a large one; :dense and :lsqr select either algorithm of
this same damped problem. LSQR reports its iteration count, normal residual
and estimated condition number.
The HypoDD airquake rule identifies a trial source above the declared depth
threshold (negative depth with the default zero threshold).
airquake_policy=:exclude_repeat rejects the whole full trial if any proposed
depth lies strictly below airquake_depth_km, removes that event's rows and
disconnected events, and retries the same successful iteration without
re-centering or re-clustering. The default :bounds is a separate constraint.

`tasks=1` solves components sequentially; a larger value bounds concurrent
component solves. Results retain their original order. An optional
`progress(done,total,index,result)` callback reports each completed component.
"""
function relocate(model, stations::AbstractDict, events, differentials;
                  absolute_picks=NamedTuple[],
                  initialization::Symbol=:catalog, initialization_frame=nothing,
                  reference_event_ids=String[], depth_bounds=(model_top(model), Inf),
                  max_iterations::Integer=50, location_tolerance_km::Real=1e-6,
                  time_tolerance_s::Real=1e-7, finite_difference_step_km::Real=1e-5,
                  solver::Symbol=:damped, row_weighting::Symbol=:empirical, row_weights=nothing,
                  damp::Real=0.,svd_relative_cutoff::Real=1e-7,
                  linear_solver::Symbol=:auto,lsqr_rtol::Real=1e-6,
                  lsqr_max_iterations::Integer=0,lsqr_condition_limit::Real=1e12,
                  airquake_policy::Symbol=:bounds,airquake_depth_km::Real=0.,
                  selection=nothing,selection_frame=nothing,empirical_weighting=nothing,
                  tasks::Integer=1,progress=nothing,_history::Symbol=:full)
    catalog = collect(events)
    _history in (:full,:replay) || throw(ArgumentError("Internal DD history must be full or replay."))
    initialization in (:catalog,:centroid) ||
        throw(ArgumentError("Relative initialization must be :catalog or :centroid."))
    (initialization_frame===nothing || initialization_frame isa LocalFrame) ||
        throw(ArgumentError("initialization_frame must be nothing or the catalog's LocalFrame."))
    observations = collect(differentials)
    _validate_relative_controls(;solver,linear_solver,lsqr_rtol,lsqr_max_iterations,lsqr_condition_limit,
        row_weighting,empirical_weighting,airquake_policy,airquake_depth_km,damp,svd_relative_cutoff,row_weights)
    ids = [String(event.event_id) for event in catalog]
    length(unique(ids)) == length(ids) || throw(ArgumentError("Duplicate event ID."))
    index = Dict(id => i for (i, id) in enumerate(ids))
    references = Set{String}(String(id) for id in reference_event_ids)
    all(haskey(index, id) for id in references) || throw(ArgumentError("Unknown reference event ID."))
    lower, upper = Float64.(depth_bounds)
    lower==upper && throw(ArgumentError("This relative path estimates four parameters per event. Equal depth bounds are not a fixed-depth solve."))
    isfinite(lower) && lower >= model_top(model) && upper >= lower ||
        throw(ArgumentError("Depth bounds must lie within the model, lower <= upper."))
    max_iterations >= 0 && location_tolerance_km > 0 && time_tolerance_s > 0 &&
        finite_difference_step_km > 0 || throw(ArgumentError("Invalid iteration settings."))
    for event in catalog
        length(event.source_km) == 3 && all(isfinite, event.source_km) && isfinite(event.origin_s) ||
            throw(ArgumentError("Event coordinates and origin time must be finite."))
        lower <= event.source_km[3] <= upper || throw(ArgumentError("Initial depth outside bounds."))
        event.source_km[3] >= _source_top(model,event.source_km) ||
            throw(ArgumentError("Initial source is above an active model top."))
    end
    neighbors = [Int[] for _ in catalog]
    arrivals = collect(absolute_picks)
    seen_arrivals = Set{Tuple}()
    for pick in arrivals
        haskey(index,pick.event_id) && haskey(stations,pick.station_id) ||
            throw(ArgumentError("Unknown event/station in absolute constraints."))
        pick.phase in ("P","S") && isfinite(pick.arrival_s) && isfinite(pick.sigma_s) && pick.sigma_s>0 ||
            throw(ArgumentError("Absolute constraints require P/S, finite arrival and positive sigma."))
        key=(String(pick.event_id),_absolute_observation_key(pick))
        key in seen_arrivals && throw(ArgumentError("Duplicate absolute constraint."))
        push!(seen_arrivals,key)
    end
    empirical_state=nothing
    weight_schedule=nothing
    empirical_default_quality=nothing
    min_weight=1e-5
    empirical_weighting===nothing && (empirical_weighting=Dict{String,Any}(
        "stages"=>_dd_default_stages(),"default_quality"=>1.))
    weights = if row_weighting==:empirical
        row_weights===nothing || throw(ArgumentError("Explicit row_weights cannot be combined with empirical quality weights."))
        settings=empirical_weighting
        stage_options=_dd_weight_option(settings,:stages,nothing)
        allowed=stage_options===nothing ? ("wtccp","wtccs","wtctp","wtcts","default_quality","absolute_row_weights") :
            ("stages","default_quality","absolute_row_weights","min_weight")
        all(String(key) in allowed for key in keys(settings)) || throw(ArgumentError("Use four static WT factors or explicit stages, with quality and absolute-row weights as needed."))
        empirical_default_quality=_dd_weight_option(settings,:default_quality,1.)
        if stage_options===nothing
            stage=only(_dd_weight_schedule([merge(Dict{String,Any}(String(k)=>v for (k,v) in pairs(settings)),Dict("iterations"=>1))]))
            empirical_state=_dd_empirical_weights(observations,Dict(row.event_id=>row.source_km for row in catalog),zeros(length(observations)),stage;
                default_quality=empirical_default_quality)
        else
            weight_schedule=_dd_weight_schedule(stage_options)
            foreach(row->_dd_row_quality(row,empirical_default_quality),observations)
            min_weight=_dd_weight_option(settings,:min_weight,1e-5)
            min_weight isa Real && isfinite(min_weight) && min_weight>0 || throw(ArgumentError("DD min_weight must be finite and positive."))
        end
        absolute_weights=_dd_weight_option(settings,:absolute_row_weights,ones(length(arrivals)))
        absolute_weights isa AbstractVector && length(absolute_weights)==length(arrivals) &&
            all(w->w isa Real && isfinite(w) && w>=0,absolute_weights) ||
            throw(ArgumentError("Joint absolute rows require explicit nonnegative absolute_row_weights in empirical mode."))
        [empirical_state===nothing ? ones(length(observations)) : empirical_state.effective_weights;Float64.(absolute_weights)]
    else
        throw(ArgumentError("Relative location reweights rows by data class, distance and residual."))
    end
    seen = Set{Tuple{String,String,String,String,String}}()
    for obs in observations
        haskey(index, obs.event_i) && haskey(index, obs.event_j) || throw(ArgumentError("Unknown event ID."))
        obs.event_i != obs.event_j || throw(ArgumentError("A differential requires two distinct events."))
        haskey(stations, obs.station_id) || throw(ArgumentError("Unknown station ID $(obs.station_id)."))
        obs.phase in ("P", "S") && obs.source in ("CT", "CC") ||
            throw(ArgumentError("Require phase P/S and source CT/CC."))
        isfinite(obs.dt_s) && isfinite(obs.sigma_s) && obs.sigma_s > 0 ||
            throw(ArgumentError("Differential and positive sigma must be finite."))
        a, b = minmax(String(obs.event_i), String(obs.event_j))
        key = (a, b, String(obs.station_id), String(obs.phase), String(obs.source))
        key in seen && throw(ArgumentError("Duplicate differential, including reversed pair: $key"))
        push!(seen, key)
        i, j = index[obs.event_i], index[obs.event_j]
        push!(neighbors[i], j)
        push!(neighbors[j], i)
    end
    initial_distance=weight_schedule===nothing ? nothing : _dd_initial_distance_selection(observations,
        Dict(row.event_id=>row.source_km for row in catalog),first(weight_schedule))
    effective_selection=selection===nothing && weight_schedule!==nothing ? (idat=3,ipha=3,obscc=0,obsct=0) : selection
    plan=effective_selection===nothing ? nothing : _dd_selection(catalog,stations,observations;
        selection=effective_selection,frame=selection_frame,
        initial_exclusion_reasons=initial_distance===nothing ? nothing : initial_distance.reasons)
    membership = zeros(Int, length(catalog))
    groups = Vector{Int}[]
    for i in eachindex(catalog)
        membership[i] == 0 || continue
        group = Int[]
        push!(groups, group)
        stack = [i]
        membership[i] = length(groups)
        while !isempty(stack)
            node = pop!(stack)
            push!(group, node)
            for other in neighbors[node]
                if membership[other] == 0
                    membership[other] = length(groups)
                    push!(stack, other)
                end
            end
        end
        sort!(group) # Preserve original catalog order within each component.
    end
    if plan!==nothing
        groups=plan.groups
        membership=plan.membership
    end
    selected_groups=plan===nothing ? collect(eachindex(groups)) : plan.selected_group_indices
    solved_components = _ordered_tasks(selected_groups;tasks,progress) do component_id
        group=groups[component_id]
        original_subset = catalog[group]
        center = initialization==:centroid ? _dd_centroid(original_subset,initialization_frame) : nothing
        subset = center===nothing ? original_subset : [merge(event,(source_km=center,)) for event in original_subset]
        for event in subset
            lower<=event.source_km[3]<=upper && event.source_km[3]>=_source_top(model,event.source_km) ||
                throw(ArgumentError("Initialized relative source is outside the depth/model domain; it is not clipped."))
        end
        differential_rows = plan===nothing ? findall(obs->membership[index[obs.event_i]]==component_id,observations) : plan.group_rows[component_id]
        absolute_rows = findall(pick->membership[index[pick.event_id]]==component_id,arrivals)
        # Infer a concrete record type when these rows share their fields.
        # map preserves every value and also accepts heterogeneous input rows.
        data = map(identity,observations[differential_rows])
        absolute_data = map(identity,arrivals[absolute_rows])
        fixed = String[event.event_id for event in subset if event.event_id in references]
        automatic = false
        result = _relocate_scaled_component(model,stations,subset,data,fixed,(lower,upper),
            max_iterations,location_tolerance_km,time_tolerance_s,finite_difference_step_km;
            absolute_observations=absolute_data,solver,row_weighting,
            weights=weights[[differential_rows;length(observations).+absolute_rows]],
            damp,svd_relative_cutoff,airquake_policy,airquake_depth_km,
            linear_solver,lsqr_rtol,lsqr_max_iterations,lsqr_condition_limit,
            weight_schedule,default_quality=empirical_default_quality,min_weight,_history)
        (;component_id,center,fixed,automatic,differential_rows,absolute_rows,result)
    end
    event_results = Dict{String,NamedTuple}()
    residuals, absolute_residuals, components = NamedTuple[], NamedTuple[], NamedTuple[]
    # Workers own the component state; merge only here, in original group order.
    for solved in solved_components
        (;component_id,center,fixed,automatic,differential_rows,absolute_rows,result)=solved
        for row in result.events
            original = catalog[index[row.event_id]]
            initial = center===nothing ? original.source_km : center
            event_results[row.event_id] = merge(row, (; component_id,initialization,selection_status=:selected,
                catalog_source_km=Tuple(original.source_km),catalog_origin_s=original.origin_s,
                initial_source_km=Tuple(initial),initial_origin_s=original.origin_s))
        end
        append!(residuals, [merge(row, (; component_id,input_row=input_row,selection_status=:selected))
            for (row,input_row) in zip(result.residuals,differential_rows)])
        append!(absolute_residuals,[merge(row,(;component_id,input_row=input_row,selection_status=:selected))
            for (row,input_row) in zip(result.absolute_residuals,absolute_rows)])
        push!(components, merge(result.summary, (; component_id,initialization,
               initialization_coordinates=initialization_frame===nothing ? :cartesian_km : :geographic_mean,
               initial_source_km=center,
               reference_event_ids=String.(fixed), automatic_reference=automatic,
               datum=isempty(fixed) ? :scaled_parameter_damping : :fixed_reference_and_scaled_parameter_damping)))
    end
    if plan!==nothing
        for (i,event) in enumerate(catalog)
            haskey(event_results,event.event_id) && continue
            reason=plan.event_statuses[i]
            event_results[event.event_id]=merge(event,(source_km=Tuple(event.source_km),status=reason,
                component_id=membership[i]==0 ? nothing : membership[i],iterations=0,input_count=0,used_count=0,
                depth_at_bound=event.source_km[3]==lower || event.source_km[3]==upper,
                initialization,selection_status=reason,catalog_source_km=Tuple(event.source_km),catalog_origin_s=event.origin_s,
                initial_source_km=Tuple(event.source_km),initial_origin_s=event.origin_s,solver,exit_reason=reason))
        end
        for row in plan.excluded
            push!(residuals,merge(row,(dt_pred_s=nothing,residual_s=nothing,normalized_residual=nothing,
                used=false,unused_reason=row.reason,component_id=row.cluster_id==0 ? nothing : row.cluster_id,
                equation_weight=0.,weighted_residual=nothing,selection_status=row.reason)))
        end
        selected_absolute=Set(row.input_row for row in absolute_residuals)
        for (input_row,row) in enumerate(arrivals)
            input_row in selected_absolute && continue
            reason=plan.event_statuses[index[row.event_id]]
            push!(absolute_residuals,merge(row,(arrival_pred_s=nothing,residual_s=nothing,normalized_residual=nothing,
                used=false,unused_reason=reason,component_id=nothing,equation_weight=0.,weighted_residual=nothing,
                input_row,selection_status=reason)))
        end
        sort!(residuals;by=row->row.input_row)
        sort!(absolute_residuals;by=row->row.input_row)
    end
    if empirical_state!==nothing
        residuals=[merge(row,(quality=empirical_state.quality[row.input_row],
            quality_source=empirical_state.quality_source[row.input_row],
            data_class=empirical_state.classes[row.input_row],class_weight=empirical_state.class_weight[row.input_row],
            empirical_base_weight=empirical_state.base[row.input_row],
            unused_reason=get(row,:selection_status,:selected)==:selected && empirical_state.reasons[row.input_row]!==nothing ?
                empirical_state.reasons[row.input_row] : row.unused_reason)) for row in residuals]
    end
    return (events=[event_results[id] for id in ids], residuals=residuals,
            absolute_residuals=absolute_residuals, components=components, uncertainty=nothing,
            observation_covariance=:diagonal_approximation,depth_bounds=(lower,upper),selection=plan,
            solver,row_weighting)
end

# ISTART=1 changes spatial starting values only. The original per-event
# origins remain the reference of the complete observed arrival differences.
function _dd_centroid(events,frame)
    depth = sum(event.source_km[3] for event in events)/length(events)
    if frame===nothing
        return (sum(event.source_km[1] for event in events)/length(events),
                sum(event.source_km[2] for event in events)/length(events),depth)
    end
    coordinates = [local_to_geographic(frame,event.source_km) for event in events]
    latitude = sum(row.latitude_deg for row in coordinates)/length(coordinates)
    longitude = sum(row.longitude_deg for row in coordinates)/length(coordinates)
    geographic_to_local(frame,latitude,longitude,-depth)
end

# Like the layered-ray candidates, a validity flag distinguishes an absent
# gradient. Its storage values are never used for a free-source equation
# before the required analytic or secant gradient has been obtained.
struct _DDRay
    time_s::Float64
    gradient_s_km::NTuple{3,Float64}
    has_gradient::Bool
end
_DDRay(time,::Nothing)=_DDRay(time,(0.,0.,0.),false)
_DDRay(time,gradient::NTuple{3,Float64})=_DDRay(time,gradient,true)
_dd_ray_cache() = Dict{Tuple{Int,String,String},_DDRay}()
_dd_ray_cache(count::Integer)=fill!(Vector{Union{Nothing,_DDRay}}(undef,count),nothing)

_dd_source_secant(model,source,stations,station,phase,step)=
    _absolute_secant(model,source,stations,station,phase,step)

function _dd_source_secant(model::LayeredModel,source,stations,station,phase,step)
    gradient=_absolute_secant(model,source,stations,station,phase,step)
    layer=searchsortedlast(model.layer_top_km,source[3])
    velocity=phase=="P" ? model.vp_km_s : model.vs_km_s
    if layer>1 && source[3]==model.layer_top_km[layer] && velocity[layer]!=velocity[layer-1]
        # At a velocity jump the vector of arrivals has two vertical limits.
        # Use its shallow-side secant consistently for all rows of this source.
        # The source itself stays at its supplied coordinate.
        before=copy(source)
        offset=min(step,(source[3]-model.layer_top_km[layer-1])/2)
        before[3]=min(source[3]-offset,prevfloat(source[3]))
        center_time=travel_time(model,Tuple(source),stations,station,phase).time_s
        before_time=travel_time(model,Tuple(before),stations,station,phase).time_s
        vertical=(center_time-before_time)/(source[3]-before[3])
        return (gradient[1],gradient[2],vertical)
    end
    gradient
end

function _dd_ray_indices(observations,arrivals,index)
    slots=Dict{Tuple{Int,String,String},Int}()
    station_names=String[]
    station_slots=Dict{String,Int}()
    function slot(event,station,phase)
        get!(slots,(event,String(station),String(phase))) do
            length(slots)+1
        end
    end
    function station_slot(station)
        get!(station_slots,String(station)) do
            push!(station_names,String(station))
            length(station_names)
        end
    end
    differential=Vector{NTuple{5,Int}}(undef,length(observations))
    for (row_index,row) in enumerate(observations)
        i,j=index[row.event_i],index[row.event_j]
        differential[row_index]=(i,j,slot(i,row.station_id,row.phase),slot(j,row.station_id,row.phase),station_slot(row.station_id))
    end
    absolute=Vector{NTuple{3,Int}}(undef,length(arrivals))
    for (row_index,row) in enumerate(arrivals)
        event=index[row.event_id]
        absolute[row_index]=(event,slot(event,row.station_id,row.phase),station_slot(row.station_id))
    end
    (;differential,absolute,station_names,ray_count=length(slots))
end

@inline function _dd_equation_ray!(rays,model,stations,event,station,phase,parameters,columns,step,build_design,
        key=(event,String(station),String(phase)))
    ray = get(rays,key,nothing)
    ray!==nothing && (!build_design || ray.has_gradient) && return ray,0
    needs_gradient = build_design && any(axis->columns[event,axis]>0,1:4)
    source = ntuple(axis->parameters[event,axis],3)
    if ray===nothing
        prediction = try
            travel_time(model,source,stations,station,phase)
        catch exception
            # Forward predictions at a fixed coincident endpoint remain valid;
            # a free source still requires the same gradient/domain checks.
            !needs_gradient && exception isa DomainError && Tuple(source)==Tuple(stations[station]) || rethrow()
            (time_s=_forward_time(model,source,stations,station,phase),dsource_s_km=nothing)
        end
        ray = _DDRay(prediction.time_s,prediction.dsource_s_km)
        rays[key] = ray
    end
    secants = 0
    if needs_gradient && !ray.has_gradient
        # A forward-only pass can cache a coincident endpoint. It must still
        # raise the original free-source domain error before trying secants.
        source==Tuple(stations[station]) && travel_time(model,source,stations,station,phase)
        derivative = _dd_source_secant(model,collect(source),stations,station,phase,step)
        ray = _DDRay(ray.time_s,derivative)
        rays[key] = ray
        secants = 1
    end
    ray,secants
end

function _dd_equations(model, stations, observations, index, parameters, columns, difference_step;
                       build_design=true,absolute_observations=NamedTuple[],equation_weights=nothing,ray_cache=nothing,
                       compact_design=false,row_indices=nothing,absolute_indices=nothing)
    n = length(observations)+length(absolute_observations)
    parameter_count=maximum(columns;init=0)
    use_compact=build_design && compact_design && parameter_count>128 &&
        parameter_count<=typemax(Int32) && _dd_compact_columns(columns)
    compact=use_compact ? _DDDesign(zeros(Int32,2,n),zeros(8,n),parameter_count) : nothing
    rows, cols, values = Int[], Int[], Float64[]
    if build_design && !use_compact
        capacity=8length(observations)+4length(absolute_observations)
        sizehint!(rows,capacity);sizehint!(cols,capacity);sizehint!(values,capacity)
    end
    residual = zeros(n)
    rays = ray_cache===nothing ? _dd_ray_cache() : ray_cache
    secants = 0
    @inbounds for (row, obs) in enumerate(observations)
        i,j = row_indices===nothing ? (index[obs.event_i],index[obs.event_j]) : row_indices[row][1:2]
        key_i=row_indices===nothing ? (i,String(obs.station_id),String(obs.phase)) : row_indices[row][3]
        key_j=row_indices===nothing ? (j,String(obs.station_id),String(obs.phase)) : row_indices[row][4]
        ray_i,si = _dd_equation_ray!(rays,model,stations,i,obs.station_id,obs.phase,parameters,columns,difference_step,build_design,key_i)
        ray_j,sj = _dd_equation_ray!(rays,model,stations,j,obs.station_id,obs.phase,parameters,columns,difference_step,build_design,key_j)
        secants += si+sj
        prediction = parameters[i, 4]-parameters[j, 4]+ray_i.time_s-ray_j.time_s
        multiplier=equation_weights===nothing ? inv(obs.sigma_s) : equation_weights[row]
        residual[row] = (obs.dt_s-prediction)*multiplier
        if use_compact
            _dd_set_design_row!(compact,row,columns[i,1],ray_i.gradient_s_km,
                columns[j,1],ray_j.gradient_s_km,multiplier)
        elseif build_design
            for (event, ray, sign) in ((i, ray_i, 1.0), (j, ray_j, -1.0))
                all(axis->columns[event,axis]==0,1:4) && continue
                gradient = (ray.gradient_s_km..., 1.0)
                for axis in 1:4
                    column = columns[event, axis]
                    column == 0 && continue
                    push!(rows, row); push!(cols, column)
                    push!(values, sign*gradient[axis]*multiplier)
                end
            end
        end
    end
    @inbounds for (offset,pick) in enumerate(absolute_observations)
        row=length(observations)+offset
        event=absolute_indices===nothing ? index[pick.event_id] : absolute_indices[offset][1]
        key=absolute_indices===nothing ? (event,String(pick.station_id),String(pick.phase)) : absolute_indices[offset][2]
        ray,s = _dd_equation_ray!(rays,model,stations,event,pick.station_id,pick.phase,parameters,columns,difference_step,build_design,key)
        secants+=s
        multiplier=equation_weights===nothing ? inv(pick.sigma_s) : equation_weights[row]
        residual[row]=(pick.arrival_s-parameters[event,4]-ray.time_s)*multiplier
        if use_compact
            _dd_set_design_row!(compact,row,columns[event,1],ray.gradient_s_km,0,nothing,multiplier)
        elseif build_design && any(axis->columns[event,axis]>0,1:4)
            gradient=(ray.gradient_s_km...,1.0)
            for axis in 1:4
                column=columns[event,axis]
                column==0 && continue
                push!(rows,row);push!(cols,column);push!(values,gradient[axis]*multiplier)
            end
        end
    end
    design = use_compact ? compact : sparse(rows,cols,values,n,parameter_count)
    return residual, design, secants
end

