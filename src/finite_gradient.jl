"""
    FiniteGradientModel(vp_top, vs_top, gp, gs, vp_half, vs_half; top_km, bottom_km)

One finite, nonnegative linear-gradient layer over a constant-speed half-space.
In the upper layer `v(z)=v_top+g*(z-top_km)`; at `bottom_km` the speed may
jump upward to `v_half`. Speeds are km/s, depths km (positive down), and g s⁻¹.
Require P>S>0 at both ends of the finite layer and in the half-space. Each
half-space speed must be at least its corresponding upper-layer bottom speed.
"""
struct FiniteGradientModel
    vp_at_top_km_s::Float64
    vs_at_top_km_s::Float64
    gradient_p_s_inv::Float64
    gradient_s_s_inv::Float64
    vp_halfspace_km_s::Float64
    vs_halfspace_km_s::Float64
    top_km::Float64
    bottom_km::Float64

    function FiniteGradientModel(vp_top::Real, vs_top::Real, gp::Real, gs::Real,
                                 vp_half::Real, vs_half::Real;
                                 top_km::Real, bottom_km::Real)
        vp, vs, g_p, g_s, vh_p, vh_s, top, bottom =
            Float64.((vp_top,vs_top,gp,gs,vp_half,vs_half,top_km,bottom_km))
        all(isfinite,(vp,vs,g_p,g_s,vh_p,vh_s,top,bottom)) ||
            throw(ArgumentError("Finite-gradient parameters must be finite."))
        bottom > top || throw(ArgumentError("Require bottom_km > top_km."))
        g_p >= 0 && g_s >= 0 || throw(ArgumentError("Require nonnegative P/S gradients."))
        vp > vs > 0 && vh_p > vh_s > 0 ||
            throw(ArgumentError("Require Vp > Vs > 0 at the top and in the half-space."))
        vp_bottom,vs_bottom = muladd(g_p,bottom-top,vp),muladd(g_s,bottom-top,vs)
        isfinite(vp_bottom) && isfinite(vs_bottom) && vp_bottom > vs_bottom > 0 ||
            throw(ArgumentError("Require finite Vp > Vs > 0 at the layer bottom."))
        vh_p >= vp_bottom && vh_s >= vs_bottom ||
            throw(ArgumentError("Each half-space speed must be >= its layer-bottom speed."))
        new(vp,vs,g_p,g_s,vh_p,vh_s,top,bottom)
    end
end

model_top(model::FiniteGradientModel) = model.top_km

_finite_logratio(x) = x == 0 ? 1.0 : log1p(x)/x

# One monotonic upper-layer leg. Its p is parameterized by sine/cosine in the
# fastest half-space, retaining small vertical slowness at grazing incidence.
function _finite_gradient_leg(va,vb,h,g,vhalf,sine_half,cosine_half)
    h == 0 && return (offset_km=0.0,time_s=0.0)
    ra,rb = va/vhalf,vb/vhalf
    ca = sqrt((1-ra)*(1+ra)+(ra*cosine_half)^2)
    cb = sqrt((1-rb)*(1+rb)+(rb*cosine_half)^2)
    ca+cb > 0 || return (offset_km=Inf,time_s=Inf)
    p = sine_half/vhalf
    offset_km = p*h*(va+vb)/(ca+cb)
    if g == 0
        return (;offset_km,time_s=h/(va*ca))
    end
    # Integrate dz/[v*sqrt(1-p²v²)]. Rationalization and log1p retain both
    # g→0 and critical incidence, without subtracting nearly equal cosines.
    factor = p^2*(va+vb)/((ca+cb)*(1+cb))
    k,ell = g*h/va,factor*g*h
    time_s = h*(_finite_logratio(k)/va+factor*_finite_logratio(ell))
    return (;offset_km,time_s)
end

# Circular ray in the untruncated gradient. The receiver's vertical derivative
# determines whether the actual arc, rather than its extended circle, turns.
function _finite_upper_ray(vtop,g,top,source,station,delta,distance,offset)
    va = muladd(g,source[3]-top,vtop)
    vb = muladd(g,station[3]-top,vtop)
    if g == 0
        grad = ntuple(k -> (delta[k]/distance)/vtop,3)
        return (time_s=distance/vtop,p=(offset/distance)/vtop,
                dsource_s_km=grad,receiver_dz=-grad[3])
    end
    geometric_speed = sqrt(va)*sqrt(vb)
    half_g_distance = g*(distance/2)
    a = half_g_distance/geometric_speed
    denominator = hypot(geometric_speed,half_g_distance)
    time_s = (distance/geometric_speed)*(a == 0 ? 1.0 : asinh(a)/a)
    grad = ((delta[1]/distance)/denominator,(delta[2]/distance)/denominator,
            (delta[3]/distance-half_g_distance/va)/denominator)
    return (;time_s,p=(offset/distance)/denominator,dsource_s_km=grad,
            receiver_dz=(-delta[3]/distance-half_g_distance/vb)/denominator)
end

"""
    travel_time(model::FiniteGradientModel, source_km, station_km, phase)

First P/S arrival in a finite gradient layer over a faster/equal half-space:
an upper-layer circular ray (`:gradient`), a critical head wave (`:head_down`),
a monotonic transmitted ray (`:transmitted`), or a half-space straight ray
(`:halfspace`). The circular ray is retained only if its actual arc stays
above the bottom. The head wave requires nonnegative along-interface length.

Return `time_s`, `dsource_s_km` (s/km), and `branch`. At a source velocity jump,
or a first-arrival tie with different gradients, `dsource_s_km=nothing`.
A continuous bottom speed retains its unique first derivative even if the
velocity gradient jumps. Endpoints above the model top and coincident points
are rejected; neither endpoint is moved to avoid the bottom interface.
"""
function travel_time(model::FiniteGradientModel, source_km, station_km,
                     phase::AbstractString)
    length(source_km) == 3 && length(station_km) == 3 ||
        throw(ArgumentError("Source and station must each contain (x, y, z) in km."))
    source = ntuple(i -> Float64(source_km[i]),3)
    station = ntuple(i -> Float64(station_km[i]),3)
    all(isfinite,source) && all(isfinite,station) ||
        throw(ArgumentError("Source and station coordinates must be finite."))
    min(source[3],station[3]) >= model.top_km ||
        throw(ArgumentError("Source and station must be at or below the model top."))
    vtop,g,vhalf = if phase == "P"
        model.vp_at_top_km_s,model.gradient_p_s_inv,model.vp_halfspace_km_s
    elseif phase == "S"
        model.vs_at_top_km_s,model.gradient_s_s_inv,model.vs_halfspace_km_s
    else
        throw(ArgumentError("Phase must be P or S."))
    end
    delta = ntuple(k -> source[k]-station[k],3)
    distance = hypot(delta...)
    distance > 0 || throw(DomainError(distance,"Coincident endpoints have no unique source gradient."))
    offset = hypot(delta[1],delta[2])
    bottom = model.bottom_km
    vbottom = muladd(g,bottom-model.top_km,vtop)
    source_interface = source[3] == bottom && vhalf != vbottom
    horizontal(p) = offset == 0 ? (0.0,0.0) : (p*delta[1]/offset,p*delta[2]/offset)
    candidates = NamedTuple[]

    if min(source[3],station[3]) >= bottom
        push!(candidates,(time_s=distance/vhalf,
              dsource_s_km=ntuple(k -> (delta[k]/distance)/vhalf,3),branch=:halfspace))
    elseif max(source[3],station[3]) <= bottom
        ray = _finite_upper_ray(vtop,g,model.top_km,source,station,delta,distance,offset)
        turns_inside = ray.dsource_s_km[3] < 0 && ray.receiver_dz < 0
        if !turns_inside || ray.p*vbottom >= 1
            push!(candidates,(time_s=ray.time_s,dsource_s_km=ray.dsource_s_km,branch=:gradient))
        end
        va = muladd(g,source[3]-model.top_km,vtop)
        vb = muladd(g,station[3]-model.top_km,vtop)
        leg_a = _finite_gradient_leg(va,vbottom,bottom-source[3],g,vhalf,1.0,0.0)
        leg_b = _finite_gradient_leg(vb,vbottom,bottom-station[3],g,vhalf,1.0,0.0)
        critical_offset = leg_a.offset_km+leg_b.offset_km
        if offset >= critical_offset
            q = sqrt((1-va/vhalf)*(1+va/vhalf))/va
            grad = (horizontal(1/vhalf)...,-q)
            time_s = leg_a.time_s+leg_b.time_s+(offset-critical_offset)/vhalf
            push!(candidates,(;time_s,dsource_s_km=grad,branch=:head_down))
        end
    else
        upper_z,lower_z = minmax(source[3],station[3])
        va = muladd(g,upper_z-model.top_km,vtop)
        hu,hl = bottom-upper_z,lower_z-bottom
        # Unique monotonic transmitted ray. Use its lower horizontal length,
        # rather than p alone, to retain a tiny lower-layer vertical slowness.
        lo,hi = 0.0,offset
        for _ in 1:80
            xh = (lo+hi)/2
            (xh == lo || xh == hi) && break
            length_half = hypot(xh,hl)
            leg = _finite_gradient_leg(va,vbottom,hu,g,vhalf,xh/length_half,hl/length_half)
            if xh+leg.offset_km < offset
                lo = xh
            else
                hi = xh
            end
        end
        xh = (lo+hi)/2
        length_half = hypot(xh,hl)
        sine_half,cosine_half = xh/length_half,hl/length_half
        p = sine_half/vhalf
        leg = _finite_gradient_leg(va,vbottom,hu,g,vhalf,sine_half,cosine_half)
        # Fermat form corrects only the rounding-sized final offset mismatch.
        time_s = muladd(p,offset-leg.offset_km-xh,leg.time_s+length_half/vhalf)
        q = if source[3] > bottom
            cosine_half/vhalf
        else
            ratio = va/vhalf
            -sqrt((1-ratio)*(1+ratio)+(ratio*cosine_half)^2)/va
        end
        push!(candidates,(;time_s,dsource_s_km=(horizontal(p)...,q),branch=:transmitted))
    end

    best = candidates[argmin(getproperty.(candidates,:time_s))]
    differentiable = !source_interface
    for ray in candidates
        if abs(ray.time_s-best.time_s) <= 1e-11*max(1.0,best.time_s) &&
           maximum(abs.(ray.dsource_s_km.-best.dsource_s_km)) > 1e-9
            differentiable = false
        end
    end
    return merge(best,(dsource_s_km=differentiable ? best.dsource_s_km : nothing,))
end
