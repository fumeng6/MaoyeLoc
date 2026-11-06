"""
    PiecewiseGradientModel(depth_km, vp_km_s, vs_km_s)

Continuous piecewise-linear P/S speeds between strictly increasing depth nodes,
followed by a constant half-space. The bottom speed is a global maximum. As in
the TTGEN model, each phase may have one low-speed zone (one descending group,
possibly containing flat segments). No velocity jumps are implied.
All nodes and endpoints share the same vertical datum; depths are km, z down.
"""
struct PiecewiseGradientModel
    depth_km::Vector{Float64}
    vp_km_s::Vector{Float64}
    vs_km_s::Vector{Float64}

    function PiecewiseGradientModel(depth_km, vp_km_s, vs_km_s)
        z, vp, vs = collect(Float64.(depth_km)), collect(Float64.(vp_km_s)), collect(Float64.(vs_km_s))
        length(z) >= 2 && length(z) == length(vp) == length(vs) ||
            throw(ArgumentError("At least two equally sized depth/P/S node arrays are required."))
        all(isfinite,z) && all(isfinite,vp) && all(isfinite,vs) && all(diff(z) .> 0) ||
            throw(ArgumentError("Nodes must be finite with strictly increasing depths."))
        all(vp .> vs) && all(vs .> 0) || throw(ArgumentError("Require Vp > Vs > 0 at every node."))
        for v in (vp,vs)
            v[end] == maximum(v) || throw(ArgumentError("The bottom half-space speed must be a global maximum."))
            signs = filter(!iszero,sign.(diff(v)))
            groups = sum(signs[k] < 0 && (k == 1 || signs[k-1] > 0) for k in eachindex(signs);init=0)
            groups <= 1 || throw(ArgumentError("TTGEN permits only one descending low-speed group per phase."))
        end
        new(z,vp,vs)
    end
end

model_top(model::PiecewiseGradientModel) = first(model.depth_km)

function _piece_speed(z,v,depth)
    k = searchsortedlast(z,depth)
    k == length(z) && return v[end]
    muladd((v[k+1]-v[k])/(z[k+1]-z[k]),depth-z[k],v[k])
end

# Each fixed monotonic leg is represented by its endpoint speeds and thickness.
# Reordering the speeds leaves both integrals unchanged and keeps log1p positive.
function _piece_legs(z,v,za,zb)
    a,b = minmax(za,zb)
    legs = NTuple{3,Float64}[]
    k = searchsortedlast(z,a)
    while a < b
        right = k < length(z) ? min(b,z[k+1]) : b
        va,vb = minmax(_piece_speed(z,v,a),_piece_speed(z,v,right))
        push!(legs,(va,vb,right-a))
        a = right
        k += 1
    end
    legs
end

# V=1/p=hypot(m,u). Unlike p alone, u retains the small vertical slowness
# when V rounds to the maximum fixed-leg speed m at nearly grazing incidence.
_piece_q(m,u,v) = hypot(sqrt(max(0.0,(m-v)*(m+v))),u)

function _piece_offset(m,u,leg)
    a,b,h = leg
    h*(a+b)/(_piece_q(m,u,a)+_piece_q(m,u,b))
end

function _piece_leg_time(m,u,leg)
    a,b,h = leg
    qa,qb = _piece_q(m,u,a),_piece_q(m,u,b)
    V = hypot(m,u)
    a == b && return h*V/(a*qa)
    factor = (a+b)/((qa+qb)*(V+qb))
    h*(_finite_logratio((b-a)/a)/a+factor*_finite_logratio(factor*(b-a)))
end

_piece_fixed_offset(m,u,legs) = sum((_piece_offset(m,u,leg) for leg in legs);init=0.0)

# Bounds on the magnitude of dX/du for a fixed leg. Do not assume that this
# magnitude is monotone: bound the numerator and denominator separately.
function _piece_slope_bounds(m,lo,hi,leg)
    a,b,h = leg
    if a == b == m
        return (h*m/hi^2,h*m/lo^2)
    elseif b == m
        qa_lo,qa_hi = _piece_q(m,lo,a),_piece_q(m,hi,a)
        c = h*(a+b)
        return (c/(qa_hi*(hi+qa_hi)),c/(qa_lo*(lo+qa_lo)))
    end
    qa_lo,qb_lo = _piece_q(m,lo,a),_piece_q(m,lo,b)
    qa_hi,qb_hi = _piece_q(m,hi,a),_piece_q(m,hi,b)
    c = h*(a+b)
    (lo*c/(qa_hi*qb_hi*(qa_hi+qb_hi)),hi*c/(qa_lo*qb_lo*(qa_lo+qb_lo)))
end

# Isolate all roots on a turning segment. A(u) grows and B(u) decreases, so
# [A(lo)+B(hi), A(hi)+B(lo)] encloses X. Derivative bounds decide when a root
# is unique; ambiguous intervals are divided, rather than sampled on a grid.
function _piece_turn_roots(m,base_speed,g,legs,upper,offset)
    a(u) = 2*_piece_q(m,u,base_speed)/g
    b(u) = _piece_fixed_offset(m,u,legs)
    x(u) = a(u)+b(u)
    da(u) = base_speed == m ? 2/g : 2u/(g*_piece_q(m,u,base_speed))
    roots = Float64[]
    stack = [(0.0,upper)]
    while !isempty(stack)
        lo,hi = pop!(stack)
        alo,ahi,blo,bhi = a(lo),a(hi),b(lo),b(hi)
        lower,upper_x = alo+bhi,ahi+blo
        margin = 32eps(max(1.0,offset,isfinite(lower) ? lower : 0.0))
        lower > offset+margin && continue
        upper_x < offset-margin && continue
        f_lo,f_hi = alo+blo-offset,ahi+bhi-offset
        f_lo == 0 && push!(roots,lo)
        f_hi == 0 && push!(roots,hi)
        d_lo,d_hi = da(lo),da(hi)
        slope_scale = max(1.0,abs(d_lo),abs(d_hi))
        for leg in legs
            small,large = _piece_slope_bounds(m,lo,hi,leg)
            d_lo -= large
            d_hi -= small
            isfinite(large) && (slope_scale += large)
        end
        slope_margin = 32eps(slope_scale)
        increasing,decreasing = d_lo > slope_margin,d_hi < -slope_margin
        if increasing || decreasing
            (f_lo < 0 < f_hi || f_hi < 0 < f_lo) || continue
            for _ in 1:160
                mid = (lo+hi)/2
                (mid == lo || mid == hi) && break
                if (x(mid) < offset) == increasing
                    lo = mid
                else
                    hi = mid
                end
            end
            push!(roots,(lo+hi)/2)
            continue
        end
        mid = (lo+hi)/2
        if mid == lo || mid == hi || upper_x-lower <= margin
            # At a fold or interval endpoint only rounding-sized distance
            # uncertainty remains. The Fermat correction below preserves T.
            choices = (lo,mid,hi)
            nearest = choices[argmin(abs.((x.(choices)).-offset))]
            abs(x(nearest)-offset) <= margin && push!(roots,nearest)
        else
            push!(stack,(mid,hi),(lo,mid))
        end
    end
    sort!(unique!(roots))
end

"""
    travel_time(model::PiecewiseGradientModel, source_km, station_km, phase)

First arrival among monotonic-depth rays, all accessible upward/downward
gradient turns, and critical propagation on constant layers or velocity peaks.
Critical propagation is also considered when a peak lies between the endpoints.
Return time_s, dsource_s_km, branch and turning_layer. At first-arrival crossings
with different endpoint slownesses the gradient is nothing. Continuous gradient
changes alone do not remove a unique derivative. Neither endpoint is moved.
"""
function travel_time(model::PiecewiseGradientModel,source_km,station_km,phase::AbstractString)
    length(source_km) == 3 && length(station_km) == 3 ||
        throw(ArgumentError("Source and station must contain (x,y,z) in km."))
    source,station = Tuple(Float64.(source_km)),Tuple(Float64.(station_km))
    all(isfinite,source) && all(isfinite,station) || throw(ArgumentError("Coordinates must be finite."))
    min(source[3],station[3]) >= model_top(model) || throw(ArgumentError("Endpoint above model top."))
    source != station || throw(DomainError(0.0,"Coincident endpoints have no unique source gradient."))
    v = phase == "P" ? model.vp_km_s : phase == "S" ? model.vs_km_s :
        throw(ArgumentError("Phase must be P or S."))
    z = model.depth_km
    vs = _piece_speed(z,v,source[3])
    dx,dy = source[1]-station[1],source[2]-station[2]
    offset = hypot(dx,dy)
    function candidate(m,u,time,vertical_sign,branch,layer)
        V = hypot(m,u)
        grad = (offset == 0 ? 0.0 : dx/(offset*V),
                offset == 0 ? 0.0 : dy/(offset*V),
                vertical_sign*_piece_q(m,u,vs)/(V*vs))
        (;time_s=time,dsource_s_km=grad,branch,turning_layer=layer)
    end
    candidates = NamedTuple[]
    legs = _piece_legs(z,v,source[3],station[3])
    if !isempty(legs)
        m = maximum(leg[2] for leg in legs)
        if offset == 0
            time = sum(leg[3]*_finite_logratio((leg[2]-leg[1])/leg[1])/leg[1] for leg in legs)
            push!(candidates,(time_s=time,dsource_s_km=(0.0,0.0,sign(source[3]-station[3])/vs),
                              branch=:direct,turning_layer=0))
        elseif _piece_fixed_offset(m,0.0,legs) >= offset
            lo,hi = 0.0,m
            while _piece_fixed_offset(m,hi,legs) > offset
                hi *= 2
                isfinite(hi) || error("Ray geometry exceeds Float64 range.")
            end
            for _ in 1:160
                mid = (lo+hi)/2
                (mid == lo || mid == hi) && break
                _piece_fixed_offset(m,mid,legs) > offset ? (lo=mid) : (hi=mid)
            end
            u = (lo+hi)/2
            time = sum(_piece_leg_time(m,u,leg) for leg in legs)+
                   (offset-_piece_fixed_offset(m,u,legs))/hypot(m,u)
            push!(candidates,candidate(m,u,time,sign(source[3]-station[3]),:direct,0))
        end
    else
        k = searchsortedlast(z,source[3])
        if k == length(z) || v[k] == v[k+1] || (k > 1 && source[3] == z[k] && v[k-1] == v[k])
            push!(candidates,candidate(vs,0.0,offset/vs,0,:horizontal,0))
        end
    end

    for k in 1:length(z)-1
        g = (v[k+1]-v[k])/(z[k+1]-z[k])
        g == 0 && continue
        down = g > 0
        base = down ? max(source[3],station[3],z[k]) : min(source[3],station[3],z[k+1])
        (down ? base < z[k+1] : base > z[k]) || continue
        fixed = vcat(_piece_legs(z,v,source[3],base),_piece_legs(z,v,station[3],base))
        vb = _piece_speed(z,v,base)
        m = maximum((leg[2] for leg in fixed);init=vb)
        end_speed = down ? v[k+1] : v[k]
        end_speed >= m || continue
        upper = sqrt((end_speed-m)*(end_speed+m))
        for u in _piece_turn_roots(m,vb,abs(g),fixed,upper,offset)
            qbase = _piece_q(m,u,vb)
            x = _piece_fixed_offset(m,u,fixed)+2qbase/abs(g)
            time = sum((_piece_leg_time(m,u,leg) for leg in fixed);init=0.0)+
                   2asinh(qbase/vb)/abs(g)+(offset-x)/hypot(m,u)
            push!(candidates,candidate(m,u,time,down ? -1 : 1,down ? :turn_down : :turn_up,k))
        end
    end

    # Peaks include both ends of a flat maximum. No same-side endpoint rule:
    # a buried source and a surface receiver can meet at the low-zone top.
    for k in eachindex(z)
        left = k == 1 ? 0.0 : v[k]-v[k-1]
        right = k == length(z) ? 0.0 : v[k+1]-v[k]
        left >= 0 && right <= 0 || continue
        fixed = vcat(_piece_legs(z,v,source[3],z[k]),_piece_legs(z,v,station[3],z[k]))
        all(leg[2] <= v[k] for leg in fixed) || continue
        critical = _piece_fixed_offset(v[k],0.0,fixed)
        offset >= critical || continue
        time = sum((_piece_leg_time(v[k],0.0,leg) for leg in fixed);init=0.0)+(offset-critical)/v[k]
        push!(candidates,candidate(v[k],0.0,time,sign(source[3]-z[k]),:head,k))
    end
    isempty(candidates) && error("No accessible first-arrival ray found.")
    best = candidates[argmin(getproperty.(candidates,:time_s))]
    differentiable = all(ray -> abs(ray.time_s-best.time_s) > 1e-11*max(1.0,best.time_s) ||
        maximum(abs.(ray.dsource_s_km.-best.dsource_s_km)) <= 1e-9,candidates)
    merge(best,(dsource_s_km=differentiable ? best.dsource_s_km : nothing,))
end
