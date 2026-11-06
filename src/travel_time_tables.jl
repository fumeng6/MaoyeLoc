"""
    TravelTimeTable(distance_km, depth_km, time_s; title="", node_error_s=0)

Single-phase travel time T(r,z), in seconds, to a receiver on a fixed plane.
Both axes start at zero and contain at least three increasing nodes in km.
`time_s` has one row per source depth and one column per horizontal distance.
Values are stored as Float64, without time reduction or integer quantization.

`node_error_s` is an optional absolute bound supplied by the caller for each
node. Zero means no node-error assumption was supplied; it is not a claim of
zero interpolation, numerical, or physical-model error.
"""
struct TravelTimeTable
    title::String
    distance_km::Vector{Float64}
    depth_km::Vector{Float64}
    time_s::Matrix{Float64}
    node_error_s::Float64

    function TravelTimeTable(distance_km,depth_km,time_s;title="",node_error_s=0.)
        distances,depths = _table_axis(distance_km),_table_axis(depth_km)
        values = Matrix{Float64}(time_s)
        size(values) == (length(depths),length(distances)) ||
            throw(ArgumentError("Travel-time rows must match depth_km and columns must match distance_km."))
        all(isfinite,values) || throw(ArgumentError("Travel-time values must be finite seconds."))
        bound = Float64(node_error_s)
        isfinite(bound) && bound >= 0 || throw(ArgumentError("Node error bound must be finite and nonnegative seconds."))
        new(String(title),distances,depths,values,bound)
    end
end

function _table_axis(values)
    axis = Float64.(collect(values))
    length(axis) >= 3 && all(isfinite,axis) && first(axis) == 0 && all(diff(axis) .> 0) ||
        throw(ArgumentError("A travel-time axis must start at zero and contain at least three increasing finite nodes."))
    axis
end

"""Read a native TOML travel-time table. `time_s` is an array of depth rows."""
function read_travel_time_table(path::AbstractString)
    data = TOML.parsefile(path)
    rows = data["time_s"]
    distances,depths = data["distance_km"],data["depth_km"]
    length(rows) == length(depths) && all(length(row) == length(distances) for row in rows) ||
        throw(ArgumentError("Travel-time rows must match depth_km and columns must match distance_km."))
    times = reduce(vcat,(permutedims(Float64.(row)) for row in rows))
    TravelTimeTable(distances,depths,times;title=get(data,"title",""),node_error_s=get(data,"node_error_s",0.))
end

"""Write a native TOML table, retaining Float64 times and explicit axes."""
function write_travel_time_table(path::AbstractString,table::TravelTimeTable)
    data = Dict("title"=>table.title,"distance_km"=>table.distance_km,
        "depth_km"=>table.depth_km,"time_s"=>[collect(row) for row in eachrow(table.time_s)],
        "node_error_s"=>table.node_error_s)
    open(path,"w") do io
        TOML.print(io,data)
    end
    path
end

# Use the three nodes around the nearest interior node, and differentiate the
# same Lagrange polynomial. Beyond the final node use the final two nodes.
function _table_stencil(nodes,x;side=0)
    if x > last(nodes) || (x == last(nodes) && side > 0)
        ids = (length(nodes)-1):length(nodes)
    else
        lower = clamp(searchsortedlast(nodes,x),1,length(nodes)-1)
        midpoint = (nodes[lower]+nodes[lower+1])/2
        nearest = x < midpoint || (x == midpoint && side < 0) ? lower : lower+1
        center = clamp(nearest,2,length(nodes)-1)
        ids = (center-1):(center+1)
    end
    w,dw = zeros(length(ids)),zeros(length(ids))
    for (i,k) in enumerate(ids)
        others = filter(!=(k),collect(ids))
        denominator = prod(nodes[k]-nodes[j] for j in others)
        w[i] = prod(x-nodes[j] for j in others)/denominator
        dw[i] = sum(prod((x-nodes[j] for j in others if j != l);init=1.0) for l in others)/denominator
    end
    (;ids,w,dw)
end

function _table_value(table,r,z;radial_side=0,depth_side=0,node_error_s=table.node_error_s)
    rw = _table_stencil(table.distance_km,r;side=radial_side)
    zw = _table_stencil(table.depth_km,z;side=depth_side)
    values = table.time_s[zw.ids,rw.ids]
    time = dot(zw.w,values*rw.w)
    dr = dot(zw.w,values*rw.dw)
    dz = dot(zw.dw,values*rw.w)
    bounds = (time_s=node_error_s*sum(abs,zw.w)*sum(abs,rw.w),
              dtdr_s_km=node_error_s*sum(abs,zw.w)*sum(abs,rw.dw),
              dtdz_s_km=node_error_s*sum(abs,zw.dw)*sum(abs,rw.w))
    (;time_s=time,dtdr_s_km=dr,dtdz_s_km=dz,node_error_bounds=bounds)
end

"""
    evaluate_table(table, distance_km, depth_km; node_error_s=table.node_error_s)

Interpolate T and the exact derivatives of the same local polynomial. Beyond
either maximum axis use linear extrapolation in that direction. Negative
distance/depth is rejected. `differentiable=false` marks a stencil boundary
where the time or either derivative has unequal one-sided limits.

`node_error_bounds` propagates only the caller's bound on table-node errors.
It does not bound interpolation error or model error. The default zero means
no node-error assumption has been supplied. `source_angle_deg` is the angle
of the inferred ray direction from downward vertical, defined from the source
gradient as atan2(abs(dT/dr),-dT/dz); it is not an independently traced angle.
"""
function evaluate_table(table::TravelTimeTable,distance_km::Real,depth_km::Real;node_error_s=table.node_error_s)
    r,z = Float64(distance_km),Float64(depth_km)
    all(isfinite,(r,z)) && min(r,z) >= 0 ||
        throw(ArgumentError("Table distance and source depth must be finite and nonnegative km."))
    isfinite(node_error_s) && node_error_s >= 0 || throw(ArgumentError("Node error bound must be finite and nonnegative seconds."))
    result = _table_value(table,r,z;node_error_s)
    differentiable = true
    for (axis,coordinate,nodes) in ((:r,r,table.distance_km),(:z,z,table.depth_km))
        at_switch = coordinate == last(nodes) ||
            any(i->coordinate == (nodes[i]+nodes[i+1])/2,1:length(nodes)-1)
        if coordinate > 0 && at_switch
            left = axis == :r ? _table_value(table,r,z;radial_side=-1) : _table_value(table,r,z;depth_side=-1)
            right = axis == :r ? _table_value(table,r,z;radial_side=1) : _table_value(table,r,z;depth_side=1)
            for key in (:time_s,:dtdr_s_km,:dtdz_s_km)
                abs(getproperty(left,key)-getproperty(right,key)) > 64eps(max(1.0,abs(getproperty(result,key)))) && (differentiable=false)
            end
        end
    end
    merge(result,(;differentiable,source_angle_deg=atand(abs(result.dtdr_s_km),-result.dtdz_s_km),
        extrapolated=r>last(table.distance_km) || z>last(table.depth_km)))
end

"""P/S tables at one explicit receiver plane. The table cannot infer station
elevation corrections: off-plane stations require a ray model. Source depth
in the table is source_z-station_plane_km. Independent S tables are used
directly; `vp_vs_ratio` instead specifies T_S=vp_vs_ratio*T_P.
"""
struct TravelTimeTableModel
    p::TravelTimeTable
    s::TravelTimeTable
    station_plane_km::Float64
    s_time_scale::Float64
    function TravelTimeTableModel(p::TravelTimeTable,s::TravelTimeTable;station_plane_km=0.0)
        isfinite(station_plane_km) || throw(ArgumentError("Station plane must be finite."))
        new(p,s,Float64(station_plane_km),1.0)
    end
    function TravelTimeTableModel(p::TravelTimeTable;vp_vs_ratio,station_plane_km=0.0)
        factor = Float64(vp_vs_ratio)
        factor > 1 && isfinite(factor) || throw(ArgumentError("An explicit Vp/Vs ratio > 1 is required."))
        isfinite(station_plane_km) || throw(ArgumentError("Station plane must be finite."))
        new(p,p,Float64(station_plane_km),factor)
    end
end
model_top(model::TravelTimeTableModel) = model.station_plane_km

function travel_time(model::TravelTimeTableModel,source_km,station_km,phase::AbstractString)
    length(source_km) == 3 && length(station_km) == 3 || throw(ArgumentError("Source and station must contain (x,y,z)."))
    source,station = Tuple(Float64.(source_km)),Tuple(Float64.(station_km))
    all(isfinite,source) && all(isfinite,station) || throw(ArgumentError("Coordinates must be finite."))
    station[3] == model.station_plane_km || throw(ArgumentError("Table receiver is fixed on station_plane_km; off-plane stations require a ray model."))
    table = phase == "P" ? model.p : phase == "S" ? model.s : throw(ArgumentError("Phase must be P or S."))
    depth = source[3]-model.station_plane_km
    depth >= 0 || throw(ArgumentError("Source lies above the table receiver plane."))
    source != station || throw(DomainError(0.,"Coincident endpoints have no unique source gradient."))
    dx,dy = source[1]-station[1],source[2]-station[2]
    r = hypot(dx,dy)
    result = evaluate_table(table,r,depth)
    if phase == "S" && model.s_time_scale != 1
        factor = model.s_time_scale
        bounds = result.node_error_bounds
        result = merge(result,(time_s=factor*result.time_s,dtdr_s_km=factor*result.dtdr_s_km,
            dtdz_s_km=factor*result.dtdz_s_km,
            node_error_bounds=(time_s=factor*bounds.time_s,dtdr_s_km=factor*bounds.dtdr_s_km,
                               dtdz_s_km=factor*bounds.dtdz_s_km)))
    end
    radial_tolerance = result.node_error_bounds.dtdr_s_km + 64eps(max(1.,abs(result.dtdr_s_km)))
    result.dtdr_s_km < -radial_tolerance && throw(DomainError(result.dtdr_s_km,
        "Nonphysical first-arrival radial slowness: negative dT/dr exceeds the supplied table-node error bound; inspect or refine the table."))
    gradient = if !result.differentiable || (r == 0 && abs(result.dtdr_s_km) > 1e-12)
        nothing
    else
        (r == 0 ? 0.0 : result.dtdr_s_km*dx/r,r == 0 ? 0.0 : result.dtdr_s_km*dy/r,result.dtdz_s_km)
    end
    merge(result,(dsource_s_km=gradient,branch=:table))
end
