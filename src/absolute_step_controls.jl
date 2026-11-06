"""
Apply the Hypoinverse DAM step controls to `[dx, dy, dz, dt]`.

The controls select an iterative step; they do not change the arrival-time
objective. All spatial values are in km, time is in seconds, and z is positive
downward. The source is the position before the step. `damp` is halved (or
multiplied by the specified `late_damping_factor`) after floor(0.6 ITRLIM), as
in the reference source. Vertical limiting then changes only dz; the model-top
and full spatial-length controls scale all four components together.
"""
function _absolute_controlled_step(step, source;
        iteration::Integer, max_iterations::Integer, damp::Real=1.,
        late_damping_factor::Real=1., dzmax_km::Real=Inf, dzair::Real=.5,
        dxmax_km::Real=Inf, model_top_km::Real=0.)
    length(step)==4 && length(source)==3 ||
        throw(ArgumentError("Step must contain dx, dy, dz, dt and source must contain x, y, z."))
    all(isfinite,step) && all(isfinite,source) && isfinite(model_top_km) ||
        throw(ArgumentError("Step, source and model top must be finite."))
    iteration>=1 && max_iterations>=0 ||
        throw(ArgumentError("Iteration must be positive and max_iterations nonnegative."))
    isfinite(damp) && damp>0 && isfinite(late_damping_factor) && late_damping_factor>0 ||
        throw(ArgumentError("Damping factors must be finite and positive."))
    0<dzair<=1 || throw(ArgumentError("dzair must lie in (0, 1]."))
    dzmax_km>0 && dxmax_km>0 ||
        throw(ArgumentError("Step limits must be positive, with Inf meaning no limit."))
    source[3]>=model_top_km || throw(ArgumentError("Source lies above the model top."))

    damping_factor = Float64(damp) *
        (iteration>floor(Int,0.6 * max_iterations) ? Float64(late_damping_factor) : 1.)
    controlled = Float64[step...] .* damping_factor
    air_factor = 1.
    if source[3]+controlled[3]<model_top_km
        vertical = -Float64(dzair)*(source[3]-model_top_km)
        air_factor = vertical/controlled[3]
        controlled .*= air_factor
        controlled[3] = vertical
    end
    vertical_limited = abs(controlled[3])>dzmax_km
    if vertical_limited
        controlled[3] *= dzmax_km/(abs(controlled[3])+dzmax_km)
    end
    spatial_length = norm(view(controlled,1:3))
    spatial_factor = spatial_length>dxmax_km ? Float64(dxmax_km)/spatial_length : 1.
    controlled .*= spatial_factor
    return (;step=controlled,damping_factor,air_factor,vertical_limited,spatial_factor)
end

"""Second-nearest observed physical station's horizontal distance in km.

`station_site_ids` identifies components belonging to one physical station.
Unmapped IDs remain distinct sites, including sites with identical coordinates.
With component coordinates that differ, use that site's shortest distance.
"""
function _absolute_second_station_distance(stations, observations, source, station_site_ids)
    distances = Dict{String,Float64}()
    for observation in observations
        id = String(observation.station_id)
        haskey(stations,id) || throw(ArgumentError("Unknown station ID $id."))
        station = stations[id]
        distance = hypot(station[1]-source[1],station[2]-source[2])
        site = String(get(station_site_ids,id,id))
        distances[site] = min(get(distances,site,Inf),distance)
    end
    length(distances)<2 && return nothing
    return sort!(collect(values(distances)))[2]
end
