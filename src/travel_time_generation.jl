"""
    generate_table(model; distance_km, depth_km, phase="P", station_plane_km=0, title="")

Evaluate first-arrival travel times on explicit km axes using an existing
one-dimensional physical ray model. Receivers lie at z=station_plane_km;
source z is station_plane_km+depth_km. The selected P or S times are stored
directly as Float64 seconds, without rounding or reduced-time encoding.

Generating accurate nodes does not guarantee accurate interpolation between
nodes; choose the grid from the ray curvature and intended location accuracy.
No node-error bound is inferred from the solver or interpolation.
"""
function generate_table(model::Union{HomogeneousModel,LayeredModel,LinearGradientModel,
                                    FiniteGradientModel,PiecewiseGradientModel};
                        distance_km,depth_km,phase="P",station_plane_km=0.,title="")
    phase in ("P","S") || throw(ArgumentError("Phase must be P or S."))
    plane = Float64(station_plane_km)
    isfinite(plane) && plane >= model_top(model) ||
        throw(ArgumentError("Station plane must be finite and covered by the ray model."))
    distances,depths = _table_axis(distance_km),_table_axis(depth_km)
    times = [r == z == 0 ? 0.0 : travel_time(model,(r,0.,plane+z),(0.,0.,plane),phase).time_s
             for z in depths, r in distances]
    TravelTimeTable(distances,depths,times;title)
end
