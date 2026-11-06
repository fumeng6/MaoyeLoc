"""
    trace_ray(model::PiecewiseGradientModel, source_depth_km, launch_angle_deg; phase="P")

Follow one ray through the continuous, piecewise-linear speed profile to the
model's top plane. Depths are the model's common z coordinate (km, down).
The launch angle is measured from upward vertical: 0 is upward, 90 horizontal,
and 180 downward. Select the actual P or S speed profile with `phase`.

This is a prescribed launch ray, which need not be the first arrival. Return
surface distance/time only for `status=:surface`. `:halfspace` means the ray
continues downward without returning, `:horizontal` a ray confined to a
constant-speed layer, and `:waveguide` two buried turns. Turning quantities
describe actual turns, not the entry point into a non-returning half-space.
No reflection, P/S conversion, attenuation or waveform amplitude is computed.
"""
function trace_ray(model::PiecewiseGradientModel,source_depth_km::Real,
                   launch_angle_deg::Real;phase="P")
    depth,angle = Float64(source_depth_km),Float64(launch_angle_deg)
    isfinite(depth) && depth >= model_top(model) ||
        throw(ArgumentError("Ray source depth must be finite and covered by the model."))
    isfinite(angle) && 0 <= angle <= 180 ||
        throw(ArgumentError("Ray launch angle must be finite, from 0 to 180 degrees from upward vertical."))
    z=model.depth_km
    v=phase == "P" ? model.vp_km_s : phase == "S" ? model.vs_km_s :
        throw(ArgumentError("Phase must be P or S."))
    source_v=_piece_speed(z,v,depth)
    p=sind(angle)/source_v
    direction=angle <= 90 ? -1 : 1
    turns,turn_layer=0,nothing
    turn_depth,turn_speed=nothing,nothing
    current=depth
    distance,time=0.0,0.0
    function result(status)
        (;source_depth_km=depth,launch_angle_deg=angle,phase,
          ray_parameter_s_km=p,status,
          distance_km=status == :surface ? distance : nothing,
          time_s=status == :surface ? time : nothing,
          turning_layer=turn_layer,turning_depth_km=turn_depth,
          turning_speed_km_s=turn_speed,turn_count=turns)
    end
    if angle == 90
        current == first(z) && return result(:surface)
        above=searchsortedfirst(z,current)-1
        if above >= 1 && v[above] < source_v
            direction=-1
        else
            below=searchsortedlast(z,current)
            if below < length(z) && v[below+1] < source_v
                direction=1
                turns=1 # The source itself is an upper turn.
                turn_layer,turn_depth,turn_speed=below,current,source_v
            else
                return result(:horizontal)
            end
        end
    end
    while true
        current == first(z) && direction < 0 && return result(:surface)
        k=direction > 0 ? searchsortedlast(z,current) : searchsortedfirst(z,current)-1
        k == length(z) && direction > 0 && return result(:halfspace)
        boundary=direction > 0 ? z[k+1] : z[k]
        initial_v,final_v=_piece_speed(z,v,current),_piece_speed(z,v,boundary)
        invp=p == 0 ? Inf : inv(p)
        turning=final_v >= invp && final_v > initial_v
        destination=turning ? current+(invp-initial_v)*(boundary-current)/(final_v-initial_v) : boundary
        if destination != current
            legs=_piece_legs(z,v,current,destination)
            if p == 0
                time+=sum(leg[3]*_finite_logratio((leg[2]-leg[1])/leg[1])/leg[1] for leg in legs)
            else
                distance+=_piece_fixed_offset(invp,0.0,legs)
                time+=sum(_piece_leg_time(invp,0.0,leg) for leg in legs)
            end
        end
        current=destination
        current == first(z) && return result(:surface)
        if turning
            turns+=1
            if turn_depth === nothing || current > turn_depth
                turn_layer,turn_depth,turn_speed=k,current,invp
            end
            turns >= 2 && return result(:waveguide)
            direction=-direction
        end
    end
end

"""
    ray_geometric_density(model::PiecewiseGradientModel, source_depth_km,
                          launch_angle_deg; phase="P")

Return sin(phi)/(r*abs(dr/dphi)), in km^-2, for the surface ring reached by a
prescribed launch ray. phi is the launch angle in radians from upward
vertical; the public input is degrees. This geometric density is not a
measured seismic amplitude or a complete elastic spreading factor.

The vertical upward limit is analytic. Other angles use centered derivatives
of the same surface branch. Return `nothing` for non-returning/trapped rays,
the source on the top plane, branch changes or an unresolved caustic.
"""
function ray_geometric_density(model::PiecewiseGradientModel,source_depth_km::Real,
                               launch_angle_deg::Real;phase="P")
    ray=trace_ray(model,source_depth_km,launch_angle_deg;phase)
    ray.status == :surface && source_depth_km > model_top(model) || return nothing
    angle=Float64(launch_angle_deg)
    if angle == 0
        v=phase == "P" ? model.vp_km_s : model.vs_km_s
        source_v=_piece_speed(model.depth_km,v,source_depth_km)
        slope=sum(h*(a+b)/(2source_v) for (a,b,h) in
            _piece_legs(model.depth_km,v,model_top(model),source_depth_km))
        return inv(slope^2)
    end
    step=min(rad2deg(cbrt(eps(Float64))),angle/2,(180-angle)/2)
    step > 0 || return nothing
    left=trace_ray(model,source_depth_km,angle-step;phase)
    right=trace_ray(model,source_depth_km,angle+step;phase)
    left.status == right.status == :surface &&
        left.turning_layer == right.turning_layer == ray.turning_layer || return nothing
    slope=(right.distance_km-left.distance_km)/(2deg2rad(step))
    isfinite(slope) && abs(slope) > 32sqrt(eps(Float64))*max(ray.distance_km,1.0) &&
        ray.distance_km > 0 || return nothing
    sind(angle)/(ray.distance_km*abs(slope))
end
