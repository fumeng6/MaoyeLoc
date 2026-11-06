# Resolve bundled data at run time so compiled applications can be moved.
function _default_data_directory()
    bundled=normpath(joinpath(Sys.BINDIR,"..","share","maoyeloc","data"))
    isfile(joinpath(bundled,"geographic_regions.toml")) ? bundled : joinpath(@__DIR__,"data")
end

# Geographic names are reference data, not velocity-model regions.
# NET 1 is a piecewise set of ordered geometric predicates. NET 2/3 use the
# original polygon tables and ordered search, with depth bands for NET 3.
const _REGION_ORDER_2 = [8,7,12,43,11,4,22,23,24,21,31,90,96,98,79,25,9,5,10,6,
    13,14,26,19,16,18,27,17,20,29,28,30,15,2,1,3,48,42,46,45,
    39,37,38,36,35,34,49,44,47,41,57,58,56,62,93,92,50,52,51,33,
    32,59,60,61,65,53,55,54,40,72,73,69,70,71,74,75,76,77,78,80,
    81,82,83,84,85,86,87,88,89,91,94,95,97,99,100,101,102,103]
const _REGION_ORDER_3 = [27,28,29,15,5,2,3,26,16,9,10,30,1,6,7,25,35,11,31,8,
    14,32,4,36,44,12,13,24,33,34,23,17,18,19,51,21,22,20,38,39,
    37,42,40,43,41,50,45,46,47,48,49]

function _old_hawaii_region(y,x,z)
    closed(bounds)=bounds[1]<=y<=bounds[2] && bounds[3]<=x<=bounds[4]
    open(bounds)=bounds[1]<y<bounds[2] && bounds[3]<x<bounds[4]
    # Each row is an ordered geographic domain and its depth/subregion rule.
    domains=(
        ((19.367,19.467,-155.317,-155.233),true,
            ()->z>13 ? 21 : z>5 ? 19 : y>19.409 ? 1 : x< -155.275 ? 2 : 3),
        ((19.167,19.433,-155.283,-155.),true,
            ()->z>13 ? 21 : z>5 ? (x< -155.242 ? 9 : x< -155.205 ? 10 : x< -155.152 ? 11 : x< -155.088 ? 12 : 13) :
                y<19.342 || y<70.1+.3271x ? 7 : x< -155.233 ? 6 : x< -155.12 ? 4 : 5),
        ((19.167,19.367,-155.417,-155.283),true,()->z<13 ? 18 : 21),
        ((19.317,19.5,-155.533,-155.317),true,()->z<13 ? 20 : 23),
        ((19.267,19.533,-155.,-154.667),true,()->z<5 ? 8 : 14),
        ((18.667,19.317,-155.717,-155.417),true,()->z<13 ? 16 : 22),
        ((19.317,19.583,-155.717,-155.317),true,()->z<13 ? 15 : 23),
        ((19.,19.650,-156.333,-155.717),false,()->25),
        ((19.650,19.917,-156.333,-155.717),false,()->26),
        ((19.917,20.417,-156.333,-155.567),false,()->27),
        ((18.667,19.167,-155.417,-155.),false,()->24),
        ((19.533,19.783,-155.15,-154.667),false,()->29),
        ((19.43,19.583,-155.317,-155.),true,()->z<13 ? 17 : 21),
        ((19.583,20.417,-155.717,-154.667),false,()->28))
    for (bounds,inclusive,region) in domains
        (inclusive ? closed(bounds) : open(bounds)) && return region()
    end
    30
end

# Winding number in the longitude/latitude plane; a point ON an edge belongs
# to the polygon. This direct geometric formulation handles concave polygons
# without copying the old signed half-crossing implementation.
function _inside_region_polygon(x,y,vertices)
    winding=0
    for i in eachindex(vertices)
        a,b=vertices[i],vertices[mod1(i+1,length(vertices))]
        cross=(b[1]-a[1])*(y-a[2])-(x-a[1])*(b[2]-a[2])
        cross==0 && min(a[1],b[1])<=x<=max(a[1],b[1]) && min(a[2],b[2])<=y<=max(a[2],b[2]) && return true
        a[2]<=y<b[2] && cross>0 && (winding+=1)
        b[2]<=y<a[2] && cross<0 && (winding-=1)
    end
    winding!=0
end

function _new_hawaii_depth_allowed(id,z)
    z<4.5 && return !(25<=id<=50 || id>=61)
    z<14 && return !(id<=14 || id==23 || 36<=id<=50 || id>=61)
    z<20 && return !(id<=23 || 25<=id<=43 || id in (51,52) || id>=61)
    !(id<=35 || 44<=id<=56 || id in (58,59))
end

# Prepare the selected reference once for a batch. The returned data is only
# read by queries; keeping it local to the run also rereads external tables on
# the next run instead of retaining a process-wide cached copy.
function _prepare_geographic_regions(net::Integer;reference_directory=_default_data_directory())
    net in 0:3 || throw(ArgumentError("Geographic NET must be 0, 1, 2 or 3."))
    net==0 && return nothing
    path=joinpath(reference_directory,"geographic_regions.toml")
    native=TOML.parsefile(path)["net$(net)"]
    source=abspath(path)
    if net==1
        codes,names=native["codes"],native["names"]
        length(codes)==length(names)==30 || throw(ArgumentError("Incomplete old Hawaii name reference."))
        return (;net=Int(net),reference_source=source,codes,names)
    end
    polygons=[(;region_id=Int(r["region_id"]),code=String(r["code"]),name=String(r["name"]),
            vertices=[(Float64(v[1]),Float64(v[2])) for v in r["vertices"]]) for r in native["regions"]]
    (;net=Int(net),reference_source=source,polygons)
end

"""
    geographic_region(latitude_deg,longitude_deg,calculation_depth_km;
        net=0,reference_directory=.../src/data)

Geographic NET names from the USGS HYPOINVERSE 1.40 reference tables:
0 disabled, 1 old Hawaii, 2 northern California, 3 new Hawaii.
Longitude is positive east. Query uses the original
calculation depth, before any report-only GEO/datum conversion; it never
selects a velocity model. Names and polygons are bundled as native TOML data.
An explicit reference_directory supplies the same native TOML format.
Search order, polygon-edge inclusion and
depth bands match the original KLAS definition. NET 3 intentionally searches
the original 51-region list, even though its table has 65 polygons. The old
outside-region IDs (NET 2:99, NET 3:52) are retained with code DIS, so identity
is the full (net,region_id,code) tuple, not the number alone.

Geometry is evaluated in Float64; the original default REAL input/table
rounding can differ at boundaries by its coordinate quantization.
"""
function geographic_region(latitude_deg,longitude_deg,calculation_depth_km;
        net::Integer=0,reference_directory=_default_data_directory(),_prepared=nothing)
    all(x->x isa Real && isfinite(x),(latitude_deg,longitude_deg,calculation_depth_km)) &&
        -90<=latitude_deg<=90 && -180<=longitude_deg<=180 || throw(ArgumentError("Region query needs finite geographic degrees and calculation depth km."))
    net in 0:3 || throw(ArgumentError("Geographic NET must be 0, 1, 2 or 3."))
    y,x,z=Float64.((latitude_deg,longitude_deg,calculation_depth_km))
    base=(;net=Int(net),latitude_deg=y,longitude_deg=x,calculation_depth_km=z)
    net==0 && return merge(base,(;status=:disabled,region_id=0,code="",name="",reference_source=nothing))
    data=_prepared===nothing ? _prepare_geographic_regions(net;reference_directory) : _prepared
    data.net==net || throw(ArgumentError("Prepared geographic reference NET does not match the query."))
    if net==1
        id=_old_hawaii_region(y,x,z)
        return merge(base,(;status=id==30 ? :outside_regions : :named_region,region_id=id,
            code=data.codes[id],name=data.names[id],reference_source=data.reference_source))
    end
    order=net==2 ? _REGION_ORDER_2 : _REGION_ORDER_3
    for id in order
        net==3 && !_new_hawaii_depth_allowed(id,z) && continue
        region=data.polygons[id]
        _inside_region_polygon(x,y,region.vertices) && return merge(base,
            (;status=:named_region,region_id=id,code=region.code,name=region.name,reference_source=data.reference_source))
    end
    merge(base,(;status=:outside_regions,region_id=length(order)+1,code="DIS",name="Distant",reference_source=data.reference_source))
end
