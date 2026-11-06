# Four physical components are kept separate: XYZ in km and origin in seconds.
function _dd_history_values(prefix,values;step=false)
    axes = step ? (:dx_km,:dy_km,:dz_km,:dorigin_s) : (:x_km,:y_km,:z_km,:origin_s)
    names = Tuple(Symbol(prefix,axis) for axis in axes)
    NamedTuple{names}(values===nothing ? (nothing,nothing,nothing,nothing) : Tuple(values))
end

function _dd_history_corrections(step,indices)
    if step===nothing || isempty(indices)
        empty_values = _dd_history_values("mean_",nothing;step=true)
        return merge(empty_values,_dd_history_values("rms_",nothing;step=true),
            _dd_history_values("max_abs_",nothing;step=true),
            (mean_position_step_km=nothing,rms_position_step_km=nothing,max_position_step_km=nothing))
    end
    changes = step[indices,:]
    count = length(indices)
    means = [sum(changes[:,axis])/count for axis in 1:4]
    rms = [sqrt(sum(abs2,changes[:,axis])/count) for axis in 1:4]
    maxima = [maximum(abs,changes[:,axis]) for axis in 1:4]
    positions = [sqrt(sum(abs2,changes[i,1:3])) for i in 1:count]
    merge(_dd_history_values("mean_",means;step=true),
        _dd_history_values("rms_",rms;step=true),
        _dd_history_values("max_abs_",maxima;step=true),
        (mean_position_step_km=sum(positions)/count,
         rms_position_step_km=sqrt(sum(abs2,positions)/count),max_position_step_km=maximum(positions)))
end

"""
    relative_iteration_diagnostics(result)

Return `(sources, centroids, corrections)` for the actual saved DD history.
`sources` includes the true iteration-zero state of each selected component,
then each recorded attempt, retaining original event IDs. No iteration records
are invented for unselected or not-entered events. Raw/full steps unavailable
after a failed solve are `nothing`, whereas a rejected update has its actual
zero state change recorded explicitly.

Raw steps are the undamped column-scaled update. Full steps are the damped
update before a rejected trial is discarded. Actual steps are the accepted
XYZ and origin corrections. Origin is not profiled separately, so
`origin_profiled` stays false. Airquake rejection leaves all actual XYZ and
origin changes at zero. An event listed in `reference_event_ids`, and a fixed
depth, stay explicit conditions rather than inferred small changes.

Centroids are equal-weight means in local XYZ km and origin seconds, using the
same membership for initial, before and current values. `scope=:all` retains
the original selected component, `:active` uses survivors after this attempt,
and `:reference` retains the declared references even if they exited. Active
membership changes therefore cannot create a spurious displacement. Each
correction row gives componentwise mean, RMS and maximum absolute corrections,
plus spatial norm statistics, for `scope=:all/:free/:active` and
`step_kind=:raw/:full/:actual`. Zero-member statistics are missing.
"""
function relative_iteration_diagnostics(result)
    sources,centroids,corrections = NamedTuple[],NamedTuple[],NamedTuple[]
    fixed_depth = result.depth_bounds[1]==result.depth_bounds[2]
    for component in result.components
        events = [event for event in result.events if event.component_id==component.component_id &&
            get(event,:selection_status,:selected)==:selected]
        isempty(events) && continue
        ids = String[event.event_id for event in events]
        event_by_id = Dict(event.event_id=>event for event in events)
        references = Set(component.reference_event_ids)
        initial = reduce(vcat,[reshape(Float64[event.initial_source_km...,event.initial_origin_s],1,4)
            for event in events])
        previous = copy(initial)
        active = Set(ids)
        exit_reasons = Dict{String,Symbol}()
        solver = get(component,:solver,:damped)
        initial_row = (attempt=0,iteration=0,successful_iterations=0,solver,stage_index=nothing,
            accepted=nothing,reason=:initial,origin_profiled=false,
            parameters_after=initial,actual_step=zeros(length(ids),4),
            full_raw_step=nothing,full_step=nothing,event_ids=ids)
        for row in [initial_row;component.trace]
            # The scaled solver keeps original component membership and row
            # order. Use its explicit IDs rather than assuming an index map.
            row_indices=Dict{String,Int}()
            for (index,id) in enumerate(row.event_ids)
                get!(row_indices,id,index)
            end
            order = [row_indices[id] for id in ids]
            current = row.parameters_after[order,:]
            actual = row.actual_step[order,:]
            raw = row.full_raw_step===nothing ? nothing : row.full_raw_step[order,:]
            full = row.full_step===nothing ? nothing : row.full_step[order,:]
            for id in get(row,:removed_by_weighting_ids,String[])
                delete!(active,id)
                reason = get(event_by_id[id],:exit_reason,nothing)
                exit_reasons[id] = reason===nothing ? :disconnected_after_weighting : Symbol(reason)
            end
            airquakes = Set(get(row,:airquake_event_ids,String[]))
            for id in get(row,:removed_or_disconnected_ids,String[])
                delete!(active,id)
                exit_reasons[id] = id in airquakes ? :airquake_removed : :disconnected_after_airquake
            end
            base = (component_id=component.component_id,attempt=row.attempt,iteration=row.iteration,
                successful_iterations=row.successful_iterations,solver,stage_index=get(row,:stage_index,nothing),
                accepted=row.accepted,reason=row.reason,origin_profiled=get(row,:origin_profiled,false))
            for (i,id) in enumerate(ids)
                push!(sources,merge(base,(event_id=id,is_reference=id in references,
                    depth_fixed=fixed_depth || id in references,active=id in active,
                    exit_reason=get(exit_reasons,id,nothing)),
                    _dd_history_values("",current[i,:]),
                    _dd_history_values("raw_",raw===nothing ? nothing : raw[i,:];step=true),
                    _dd_history_values("full_",full===nothing ? nothing : full[i,:];step=true),
                    _dd_history_values("actual_",actual[i,:];step=true)))
            end
            for (scope,indices) in ((:all,collect(eachindex(ids))),
                    (:active,findall(id->id in active,ids)),
                    (:reference,findall(id->id in references,ids)))
                means = [isempty(indices) ? nothing :
                    [sum(values[indices,axis])/length(indices) for axis in 1:4]
                    for values in (initial,previous,current)]
                initial_mean,before_mean,current_mean = means
                drift = current_mean===nothing ? nothing : current_mean-initial_mean
                change = current_mean===nothing ? nothing : current_mean-before_mean
                push!(centroids,merge(base,(;scope,event_count=length(indices),event_ids=ids[indices]),
                    _dd_history_values("initial_",initial_mean),
                    _dd_history_values("before_",before_mean),_dd_history_values("current_",current_mean),
                    _dd_history_values("drift_",drift;step=true),_dd_history_values("actual_",change;step=true)))
            end
            if row.attempt>0
                for (scope,indices) in ((:all,collect(eachindex(ids))),
                        (:free,findall(id->!(id in references),ids)),
                        (:active,findall(id->id in active,ids)))
                    for (step_kind,step) in ((:raw,raw),(:full,full),(:actual,actual))
                        push!(corrections,merge(base,(;scope,step_kind,event_count=length(indices)),
                            _dd_history_corrections(step,indices)))
                    end
                end
            end
            previous = current
        end
    end
    (;sources,centroids,corrections)
end
