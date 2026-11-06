"""
    LinearGradientModel(vp_at_top, vs_at_top, gradient_p, gradient_s; top_km)

A continuous linear-gradient half-space, with speeds in km/s and z down in km:
`v(z) = v_at_top + gradient * (z - top_km)`. Gradients have units s⁻¹.
Require `Vp_at_top > Vs_at_top > 0` and `gradient_p >= gradient_s >= 0`, so
P remains faster than S throughout the half-space. Zero gradients are allowed.
The model top, source and station use the same vertical reference; h gives z=-h.
"""
struct LinearGradientModel
    vp_at_top_km_s::Float64
    vs_at_top_km_s::Float64
    gradient_p_s_inv::Float64
    gradient_s_s_inv::Float64
    top_km::Float64

    function LinearGradientModel(vp_at_top::Real, vs_at_top::Real,
                                 gradient_p::Real, gradient_s::Real; top_km::Real)
        vp, vs, gp, gs, top = Float64.((vp_at_top, vs_at_top, gradient_p, gradient_s, top_km))
        all(isfinite, (vp, vs, gp, gs, top)) ||
            throw(ArgumentError("Linear-gradient model parameters must be finite."))
        vp > vs > 0 || throw(ArgumentError("Require Vp_at_top > Vs_at_top > 0 km/s."))
        gp >= gs >= 0 ||
            throw(ArgumentError("Require gradient_p >= gradient_s >= 0 s⁻¹ throughout the half-space."))
        new(vp, vs, gp, gs, top)
    end
end

model_top(model::LinearGradientModel) = model.top_km

"""
    travel_time(model::LinearGradientModel, source_km, station_km, phase)

First P/S arrival along the circular ray in a continuous positive linear gradient
(a vertical line for zero horizontal separation). For gradient g and endpoint
speeds v_s,v_r, `T = 2*asinh(g*r/(2*sqrt(v_s*v_r)))/g`.
The equivalent expression below retains the g→0 limit without cancellation.
The circular path lies on the deeper side of its straight chord and cannot rise
above the shallower endpoint. Its local curvature points toward the shallower,
lower-speed side of the medium.

Return `time_s` and analytic `dsource_s_km` in s and s/km, with the receiver
fixed. At the model top the normal derivative is from within the model.
Zero gradient gives the exact homogeneous formula. Coincident endpoints have
no unique source gradient and are rejected, as in the other travel-time models.
"""
function travel_time(model::LinearGradientModel, source_km, station_km,
                     phase::AbstractString)
    length(source_km) == 3 && length(station_km) == 3 ||
        throw(ArgumentError("Source and station must each contain (x, y, z) in km."))
    source = ntuple(i -> Float64(source_km[i]), 3)
    station = ntuple(i -> Float64(station_km[i]), 3)
    all(isfinite, source) && all(isfinite, station) ||
        throw(ArgumentError("Source and station coordinates must be finite."))
    min(source[3], station[3]) >= model.top_km ||
        throw(ArgumentError("Source and station must be at or below the model top."))
    speed_top, gradient = if phase == "P"
        model.vp_at_top_km_s, model.gradient_p_s_inv
    elseif phase == "S"
        model.vs_at_top_km_s, model.gradient_s_s_inv
    else
        throw(ArgumentError("Phase must be P or S."))
    end
    delta = ntuple(i -> source[i] - station[i], 3)
    distance = hypot(delta...)
    distance > 0 ||
        throw(DomainError(distance, "Coincident endpoints have no unique source gradient."))
    if gradient == 0
        return (time_s=distance/speed_top,
                dsource_s_km=ntuple(i -> (delta[i]/distance)/speed_top, 3))
    end
    speed_source = muladd(gradient, source[3]-model.top_km, speed_top)
    speed_station = muladd(gradient, station[3]-model.top_km, speed_top)
    geometric_speed = sqrt(speed_source)*sqrt(speed_station)
    half_gradient_distance = gradient*(distance/2)
    a = half_gradient_distance/geometric_speed
    time_s = (distance/geometric_speed)*(a == 0 ? 1.0 : asinh(a)/a)
    denominator = hypot(geometric_speed, half_gradient_distance)
    dsource_s_km = ((delta[1]/distance)/denominator,
                    (delta[2]/distance)/denominator,
                    (delta[3]/distance-half_gradient_distance/speed_source)/denominator)
    return (; time_s, dsource_s_km)
end
