"""
    LayeredModel(layer_top_km, vp_km_s, vs_km_s)

Horizontal, constant-velocity layers with strictly increasing tops (km, z down).
The last layer is a half-space. Low-velocity layers are allowed; gradients are
not represented by this model. Both source and station must lie within it.
"""
struct LayeredModel
    layer_top_km::Vector{Float64}
    vp_km_s::Vector{Float64}
    vs_km_s::Vector{Float64}

    function LayeredModel(layer_top_km, vp_km_s, vs_km_s)
        top, vp, vs = Float64.(layer_top_km), Float64.(vp_km_s), Float64.(vs_km_s)
        !isempty(top) && length(top) == length(vp) == length(vs) ||
            throw(ArgumentError("Layer tops, Vp and Vs must have equal nonzero length."))
        all(isfinite, top) && all(isfinite, vp) && all(isfinite, vs) ||
            throw(ArgumentError("Layer parameters must be finite."))
        all(diff(top) .> 0) || throw(ArgumentError("Layer tops must increase strictly."))
        all(vp .> vs) && all(vs .> 0) || throw(ArgumentError("Require Vp > Vs > 0 in every layer."))
        new(collect(top), collect(vp), collect(vs))
    end
end

model_top(model::HomogeneousModel) = model.top_km
model_top(model::LayeredModel) = first(model.layer_top_km)

# Vertical thickness of one layer between two endpoint depths. Zero-thickness
# layers are skipped by the caller; their absence does not change the sums.
@inline function _layer_span(top, za, zb, k)
    shallow, deep = minmax(za, zb)
    lower = k < length(top) ? top[k + 1] : deep
    max(0.0, min(deep, lower) - max(shallow, top[k]))
end

# Horizontal distance X(u) and dX/du in km, with u=tan(theta)
# in the fastest traversed layer. Newton steps stay inside a valid bracket.
@inline function _direct_horizontal_derivative(top,v,za,zb,vmax,u)
    total,slope=0.,0.
    for k in eachindex(top)
        span=_layer_span(top,za,zb,k)
        span==0 && continue
        ratio=v[k]/vmax
        denominator=hypot(1.,sqrt((1-ratio)*(1+ratio))*u)
        inverse=inv(denominator)
        total+=span*ratio*(u*inverse)
        slope+=span*ratio*inverse^3
    end
    total,slope
end
function _direct_arrival(top,v,za,zb,offset)
    vmax=-Inf
    thickness=0.
    for k in eachindex(top)
        span=_layer_span(top,za,zb,k)
        span==0 && continue
        thickness+=span
        vmax=max(vmax,v[k])
    end
    isfinite(vmax) || error("Direct ray has no traversed layer.")
    uniform=true
    for k in eachindex(top)
        if _layer_span(top,za,zb,k)>0 && v[k]!=vmax
            uniform=false
            break
        end
    end
    if uniform
        distance=hypot(offset,thickness)
        return (time_s=distance/vmax,p=(offset/distance)/vmax,cosine_fast=thickness/distance,vmax=vmax)
    end
    lo,hi=0.,max(1.,offset/thickness)
    while first(_direct_horizontal_derivative(top,v,za,zb,vmax,hi))<offset
        hi*=2
        isfinite(hi) || error("Direct-ray geometry exceeds Float64 range.")
    end
    u=offset/thickness
    for iteration=1:80
        distance,slope=_direct_horizontal_derivative(top,v,za,zb,vmax,u)
        mismatch=distance-offset
        abs(mismatch)<=4eps(offset) && break
        if mismatch<0
            lo=u
        else
            hi=u
        end
        proposal=u-mismatch/slope
        next=isfinite(proposal) && lo<proposal<hi ? proposal : lo+(hi-lo)/2
        next==u && break
        u=next
    end
    cosine_fast=inv(hypot(1.,u))
    p=(u*cosine_fast)/vmax
    time_s=p*offset
    for k in eachindex(top)
        span=_layer_span(top,za,zb,k)
        span==0 && continue
        ratio=v[k]/vmax
        time_s+=span*sqrt((1-ratio)*(1+ratio)+(ratio*cosine_fast)^2)/v[k]
    end
    (;time_s,p,cosine_fast,vmax)
end

function _head_arrival(top, v, z_source, z_station, boundary, refractor, offset)
    speed = v[refractor]
    xcritical = 0.0
    intercept = 0.0
    for k in eachindex(top)
        leg = _layer_span(top, z_source, boundary, k) + _layer_span(top, z_station, boundary, k)
        leg == 0.0 && continue
        v[k] < speed || return nothing
        sine = v[k] / speed
        cosine = sqrt((1 - sine) * (1 + sine))
        xcritical += leg * sine / cosine
        intercept += leg * cosine / v[k]
    end
    offset >= xcritical || return nothing
    (time_s=offset / speed + intercept, p=1 / speed)
end

function _gradient_components(offset, dx, dy, p, zsign, q)
    (offset == 0 ? 0.0 : p * dx / offset,
     offset == 0 ? 0.0 : p * dy / offset,
     zsign * q)
end

function _layered_source_gradient(source_interface, v_source, offset, dx, dy, p, zsign)
    source_interface && return nothing
    q = sqrt(max(0.0, (1 - p * v_source) * (1 + p * v_source))) / v_source
    _gradient_components(offset, dx, dy, p, zsign, q)
end

function _layered_source_gradient(source_interface, ::Any, offset, dx, dy, p, zsign, vertical_slowness::Float64)
    source_interface && return nothing
    _gradient_components(offset, dx, dy, p, zsign, vertical_slowness)
end

struct _LayeredCandidate
    time_s::Float64
    gx::Float64
    gy::Float64
    gz::Float64
    has_gradient::Bool
    branch::Symbol
    layer::Int
end

function _remember_ray!(rays, time_s, gradient, branch, layer)
    if gradient === nothing
        push!(rays, _LayeredCandidate(time_s, 0.0, 0.0, 0.0, false, branch, layer))
    else
        push!(rays, _LayeredCandidate(time_s, gradient[1], gradient[2], gradient[3], true, branch, layer))
    end
    rays
end

"""
    travel_time(model::LayeredModel, source_km, station_km, phase)

First P/S arrival among monotonic-depth transmitted rays and valid head waves
on faster layers above or below both endpoints. Both endpoint depths enter
each ray leg. Reflections and conversions cannot be first arrivals here.

Return `time_s`, `dsource_s_km`, `branch` (:direct, :head_down, :head_up), and
`refractor_layer` (0 for direct). The gradient is in s/km for the source with
the receiver fixed. At a source velocity discontinuity or a first-arrival
tie with differing gradients, it is `nothing`: time remains well defined.
No source coordinate is moved to avoid a layer boundary. At the model top,
the normal derivative is the derivative from within the model.
"""
function travel_time(model::LayeredModel, source_km, station_km, phase::AbstractString)
    length(source_km) == 3 && length(station_km) == 3 ||
        throw(ArgumentError("Source and station must each contain (x, y, z) in km."))
    source = ntuple(i -> Float64(source_km[i]), 3)
    station = ntuple(i -> Float64(station_km[i]), 3)
    all(isfinite, source) && all(isfinite, station) || throw(ArgumentError("Coordinates must be finite."))
    min(source[3], station[3]) >= model_top(model) || throw(ArgumentError("Endpoint above model top."))
    v = phase == "P" ? model.vp_km_s : phase == "S" ? model.vs_km_s :
        throw(ArgumentError("Phase must be P or S."))
    source != station || throw(DomainError(0.0, "Coincident endpoints have no unique source gradient."))
    top = model.layer_top_km
    source_layer = searchsortedlast(top, source[3])
    source_interface = source_layer > 1 && source[3] == top[source_layer] &&
                       v[source_layer] != v[source_layer-1]
    dx, dy = source[1]-station[1], source[2]-station[2]
    offset = hypot(dx, dy)
    v_source = v[source_layer]
    rays = _LayeredCandidate[]
    if source[3] == station[3]
        p = 1/v_source
        _remember_ray!(rays, offset*p,
                       _layered_source_gradient(source_interface, v_source, offset, dx, dy, p, 0),
                       :direct, 0)
    else
        ray = _direct_arrival(top, v, source[3], station[3], offset)
        ratio = v_source/ray.vmax
        q_source = source_interface ? 0.0 :
            sqrt((1-ratio)*(1+ratio)+(ratio*ray.cosine_fast)^2)/v_source
        _remember_ray!(rays, ray.time_s,
                       _layered_source_gradient(source_interface, v_source, offset, dx, dy, ray.p,
                                               sign(source[3]-station[3]), q_source),
                       :direct, 0)
    end
    for k in 2:length(top)
        boundary = top[k]
        # At zero-length legs, along-interface propagation on either side is allowed.
        if max(source[3], station[3]) <= boundary
            head = _head_arrival(top, v, source[3], station[3], boundary, k, offset)
            head === nothing || _remember_ray!(rays, head.time_s,
                _layered_source_gradient(source_interface, v_source, offset, dx, dy, head.p, -1),
                :head_down, k)
        end
        if min(source[3], station[3]) >= boundary
            head = _head_arrival(top, v, source[3], station[3], boundary, k-1, offset)
            head === nothing || _remember_ray!(rays, head.time_s,
                _layered_source_gradient(source_interface, v_source, offset, dx, dy, head.p, 1),
                :head_up, k-1)
        end
    end
    best_index = 1
    for i in 2:length(rays)
        rays[i].time_s < rays[best_index].time_s && (best_index = i)
    end
    best = rays[best_index]
    differentiable = best.has_gradient
    # A branch label can change without a physical kink when the gradients coincide.
    tolerance = 1e-11 * max(1.0, best.time_s)
    for ray in rays
        abs(ray.time_s - best.time_s) <= tolerance || continue
        if !ray.has_gradient || !best.has_gradient ||
           max(abs(ray.gx - best.gx), abs(ray.gy - best.gy), abs(ray.gz - best.gz)) > 1e-9
            differentiable = false
        end
    end
    return (time_s=best.time_s,
            dsource_s_km=differentiable ? (best.gx, best.gy, best.gz) : nothing,
            branch=best.branch, refractor_layer=best.layer)
end
