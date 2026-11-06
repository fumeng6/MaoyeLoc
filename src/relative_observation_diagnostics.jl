# Saved final rows are the observation authority. In particular, empirical
# weights are frozen solver multipliers, not recomputed from final residuals.
function _relative_diagnostic_weight(row,component)
    weighting=get(row,:selection_status,:selected)==:selected ? get(component,:row_weighting,:empirical) : :not_applied
    hasproperty(row,:equation_weight) ||
        throw(ArgumentError("Relative diagnostics require the saved final equation_weight."))
    weight=row.equation_weight
    residual=get(row,:residual_s,nothing)
    weighted=get(row,:weighted_residual,residual===nothing ? nothing : weight*residual)
    (;equation_weight=weight,weighted_residual=weighted,row_weighting=weighting,
      row_weight_definition=:saved_final_equation_weight,
      equation_weight_unit=weighting==:not_applied ? :not_applicable : :dimensionless,
      weighted_residual_unit=weighting==:not_applied ? :not_applicable : :seconds)
end

struct _RelativeEndpointGeometry
    distance_km::Float64
    distance_m::Float64
    azimuth_deg::Union{Nothing,Float64}
    source_angle_deg::Union{Nothing,Float64}
    angle_definition::Symbol
    angle_reason::Union{Nothing,Symbol}
    source_angle_components::Vector{NamedTuple}
end

function _relative_geometry_fields(geometry)
    (;geometry.distance_km,geometry.distance_m,geometry.azimuth_deg,geometry.source_angle_deg,
        geometry.angle_definition,geometry.angle_reason,geometry.source_angle_components)
end

function _relative_diagnostic_geometry(model,stations,event,row,evaluated;cache=nothing)
    if cache!==nothing
        key=(event.event_id,row.station_id,row.phase,evaluated)
        return get!(cache,key) do
            _relative_diagnostic_geometry(model,stations,event,row,evaluated)
        end
    end
    source=event.source_km
    station=stations[row.station_id]
    dx,dy=station[1]-source[1],station[2]-source[2]
    distance=hypot(dx,dy)
    azimuth=distance==0 ? nothing : mod(atand(dx,dy),360.)
    angle = if !evaluated
        (;source_angle_deg=nothing,angle_definition=:unavailable,
          angle_reason=:observation_not_evaluated,source_angle_components=NamedTuple[])
    else
        try
            _geometry_source_angle(travel_time(model,source,stations,row.station_id,row.phase))
        catch error
            error isa DomainError || rethrow()
            (;source_angle_deg=nothing,angle_definition=:unavailable,
              angle_reason=:undefined_ray,source_angle_components=NamedTuple[])
        end
    end
    _RelativeEndpointGeometry(distance,1000distance,azimuth,angle.source_angle_deg,
        angle.angle_definition,angle.angle_reason,angle.source_angle_components)
end

_relative_event_separation(a,b)=hypot((a.source_km[k]-b.source_km[k] for k in 1:3)...)

function _relative_differential_diagnostic(model,stations,row,result_row,events,components,geometry_cache)
    scale1000(value)=value===nothing ? nothing : 1000value
    input_row=get(row,:input_row,result_row)
    component=get(components,get(row,:component_id,nothing),(;))
    weight=_relative_diagnostic_weight(row,component)
    evaluated=get(row,:residual_s,nothing)!==nothing && get(row,:selection_status,:selected)==:selected
    a,b=events[row.event_i],events[row.event_j]
    ga=_relative_diagnostic_geometry(model,stations,a,row,evaluated;cache=geometry_cache)
    gb=_relative_diagnostic_geometry(model,stations,b,row,evaluated;cache=geometry_cache)
    separation=_relative_event_separation(a,b)
    class=_dd_observation_class(row)
    idx=class==:ccp ? 1 : class==:ccs ? 2 : class==:ctp ? 3 : 4
    quality_source=get(row,:quality_source,hasproperty(row,:original_quality) ? :original_quality :
        hasproperty(row,:quality) ? :quality : :not_supplied)
    original_quality=get(row,:original_quality,get(row,:quality,nothing))
    return merge(row,(;input_row,result_row,
        input_row_provenance=hasproperty(row,:input_row) ? :saved_input_row : :result_row_order,
        idx,original_quality,quality_source,dt_ms=scale1000(row.dt_s),sigma_ms=scale1000(row.sigma_s),
        dt_pred_ms=scale1000(get(row,:dt_pred_s,nothing)),residual_ms=scale1000(get(row,:residual_s,nothing)),
        weight...,separation_km=separation,separation_m=1000separation,
        distance_i_km=ga.distance_km,distance_i_m=ga.distance_m,azimuth_i_deg=ga.azimuth_deg,
        source_angle_i_deg=ga.source_angle_deg,angle_i_definition=ga.angle_definition,
        angle_i_reason=ga.angle_reason,source_angle_i_components=ga.source_angle_components,
        distance_j_km=gb.distance_km,distance_j_m=gb.distance_m,azimuth_j_deg=gb.azimuth_deg,
        source_angle_j_deg=gb.source_angle_deg,angle_j_definition=gb.angle_definition,
        angle_j_reason=gb.angle_reason,source_angle_j_components=gb.source_angle_components))
end

function _relative_absolute_diagnostic(model,stations,row,result_row,events,components,geometry_cache)
    scale1000(value)=value===nothing ? nothing : 1000value
    input_row=get(row,:input_row,result_row)
    component=get(components,get(row,:component_id,nothing),(;))
    weight=_relative_diagnostic_weight(row,component)
    evaluated=get(row,:residual_s,nothing)!==nothing && get(row,:selection_status,:selected)==:selected
    geometry=_relative_diagnostic_geometry(model,stations,events[row.event_id],row,evaluated;cache=geometry_cache)
    return merge(row,(;input_row,result_row,
        input_row_provenance=hasproperty(row,:input_row) ? :saved_input_row : :result_row_order,
        arrival_ms=scale1000(row.arrival_s),sigma_ms=scale1000(row.sigma_s),
        arrival_pred_ms=scale1000(get(row,:arrival_pred_s,nothing)),
        residual_ms=scale1000(get(row,:residual_s,nothing)),weight...,_relative_geometry_fields(geometry)...))
end

"""
    relative_observation_diagnostics(model, stations, result)

Return final double-difference `rows` and separate joint `absolute_rows`.
Rows retain saved input fields, identities, selection/exit reasons and final
residuals. `input_row` is the original index saved by `relocate`; results are
ordered by that index. A manually supplied result lacking it instead records
`:result_row_order` as `input_row_provenance`. No observations are constructed.

IDX (`idx`) is 1=CC-P, 2=CC-S, 3=CT-P, 4=CT-S. `original_quality` is the
supplied original_quality/quality; absent values stay missing,
and an explicitly assumed empirical default is identified by quality_source.
The actual saved final equation multiplier and weighted residual are reused.
`sigma_s` stays the recorded measurement scatter and is not that multiplier.

`separation_km` is the final three-dimensional event distance from result.events,
not a distance saved during reweighting. `distance_i_km`/`distance_j_km` are
epicentral DIST; AZM is clockwise from north to the station. Source angles
are measured from downward vertical (0 downward, 180 upward). Smooth isotropic
angles use the Fermat direction -grad_source(T); nonunique rays have no single
angle. Regional mixtures retain branch angles rather than invent a mixed ray;
a single actual branch can retain its unique angle. Interpolated table normals
are identified separately by angle_definition.

Rows without a saved evaluated residual do not trigger extra forward modeling.
Their coordinate distances and azimuths still describe the supplied final
catalog, with selection/exit status preserved. Seconds/milliseconds and
kilometers/meters have separate named fields; `units` declares their meanings.
"""
function relative_observation_diagnostics(model,stations::AbstractDict,result)
    events=Dict(event.event_id=>event for event in result.events)
    components=Dict(component.component_id=>component for component in get(result,:components,NamedTuple[]))
    # Final source positions are fixed for this call. Repeated DD endpoints
    # share their geometry, while unevaluated observations never request a ray.
    geometry_cache=Dict{Tuple{String,String,String,Bool},_RelativeEndpointGeometry}()
    rows=NamedTuple[]
    for (result_row,row) in enumerate(result.residuals)
        push!(rows,_relative_differential_diagnostic(model,stations,row,result_row,events,components,geometry_cache))
    end
    absolute_rows=NamedTuple[]
    for (result_row,row) in enumerate(get(result,:absolute_residuals,NamedTuple[]))
        push!(absolute_rows,_relative_absolute_diagnostic(model,stations,row,result_row,events,components,geometry_cache))
    end
    sort!(rows;by=row->row.input_row)
    sort!(absolute_rows;by=row->row.input_row)
    (;rows,absolute_rows,units=(;time=:seconds,time_milli=:milliseconds,distance=:kilometers,
        distance_milli=:meters,angle=:degrees,separation=:three_dimensional,DIST=:epicentral),
      source_angle_definition=:clockwise_azimuth_from_north_and_source_angle_from_downward_vertical)
end
