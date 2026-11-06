"""
    LocalFrame(latitude_deg, longitude_deg; vertical_datum)

Local east/north orthographic projection of the WGS84 ellipsoid at zero
ellipsoidal height. z=-elevation_km is supplied independently in the explicitly
named vertical reference plane. A geoid is the declared reference surface for
physical elevations; this projection does not calculate geoid undulation or
convert ellipsoidal heights to elevations above it.
Horizontal distances are planar approximations to geographic surface distances,
appropriate to a documented local study area, not global/geodesic distances.
"""
struct LocalFrame
    latitude_deg::Float64
    longitude_deg::Float64
    vertical_datum::String
    origin::NTuple{3,Float64}
    east::NTuple{3,Float64}
    north::NTuple{3,Float64}
    up::NTuple{3,Float64}

    function LocalFrame(latitude_deg::Real, longitude_deg::Real; vertical_datum::AbstractString)
        lat, lon = Float64(latitude_deg), Float64(longitude_deg)
        isfinite(lat) && isfinite(lon) && abs(lat) <= 90 && abs(lon) <= 180 ||
            throw(ArgumentError("Invalid WGS84 latitude/longitude in degrees."))
        !isempty(strip(vertical_datum)) || throw(ArgumentError("An explicit vertical datum is required."))
        a, b = 6378.137, 6356.752314245179
        s, c = sind(lat), cosd(lat)
        radius = a/sqrt(1-(1-(b/a)^2)*s^2)
        origin = (radius*c*cosd(lon), radius*c*sind(lon), radius*(b/a)^2*s)
        east = (-sind(lon), cosd(lon), 0.0)
        north = (-s*cosd(lon), -s*sind(lon), c)
        up = (c*cosd(lon), c*sind(lon), s)
        new(lat, lon, String(vertical_datum), origin, east, north, up)
    end
end

"Project WGS84 lat/lon (degrees), elevation (km, named datum) to east/north/down km."
function geographic_to_local(frame::LocalFrame, latitude_deg::Real, longitude_deg::Real, elevation_km::Real)
    lat, lon, height = Float64(latitude_deg), Float64(longitude_deg), Float64(elevation_km)
    all(isfinite, (lat, lon, height)) && abs(lat) <= 90 && abs(lon) <= 180 ||
        throw(ArgumentError("Geographic coordinates must be finite and in degrees/km."))
    a, b = 6378.137, 6356.752314245179
    s, c = sind(lat), cosd(lat)
    radius = a/sqrt(1-(1-(b/a)^2)*s^2)
    surface = (radius*c*cosd(lon), radius*c*sind(lon), radius*(b/a)^2*s)
    sum(surface[i]*frame.up[i]/(i == 3 ? b^2 : a^2) for i in 1:3) > 0 ||
        throw(ArgumentError("Point lies outside the local projection hemisphere."))
    delta = surface .- frame.origin
    return (dot(delta, frame.east), dot(delta, frame.north), -height)
end

"Inverse of the local WGS84 surface projection; returns lat/lon and datum elevation in km."
function local_to_geographic(frame::LocalFrame, source_km)
    length(source_km) == 3 && all(isfinite, source_km) || throw(ArgumentError("Require finite x,y,z in km."))
    a, b = 6378.137, 6356.752314245179
    tangent = frame.origin .+ source_km[1] .* frame.east .+ source_km[2] .* frame.north
    scale = (a^2, a^2, b^2)
    aa = sum(frame.up[i]^2/scale[i] for i in 1:3)
    bb = 2*sum(tangent[i]*frame.up[i]/scale[i] for i in 1:3)
    cc = sum(tangent[i]^2/scale[i] for i in 1:3)-1
    discriminant = bb^2-4*aa*cc
    discriminant >= 0 || throw(ArgumentError("Point lies outside the local projection ellipse."))
    normal = -2*cc/(bb+sqrt(discriminant))
    surface = tangent .+ normal .* frame.up
    longitude_deg = atand(surface[2], surface[1])
    latitude_deg = atand(surface[3]/b^2, hypot(surface[1], surface[2])/a^2)
    return (; latitude_deg, longitude_deg, elevation_km=-Float64(source_km[3]))
end
