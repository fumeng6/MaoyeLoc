function _origin_utc(reference::DateTime, relative_s)
    isfinite(relative_s) || throw(ArgumentError("Origin time must be finite."))
    whole = floor(Int, relative_s)
    nanoseconds = 1_000_000*millisecond(reference) + round(Int, (relative_s-whole)*1e9)
    carry, fraction = divrem(nanoseconds, 1_000_000_000)
    calendar = reference - Millisecond(millisecond(reference)) + Second(whole+carry)
    return Dates.format(calendar, dateformat"yyyy-mm-ddTHH:MM:SS") * "." * lpad(string(fraction),9,'0') * "Z"
end

"""
    write_geographic_results(directory, result, frame; time_reference_utc)

Write inverse-projected WGS84 lat/lon, common-datum depth and explicit UTC
origin labels alongside the authoritative local km/relative-second catalog.
UTC labels retain nine fractional digits; this formatting is not an estimate
of time accuracy. Native origin_s retains its Float64 value.
"""
function write_geographic_results(directory, result, frame::LocalFrame; time_reference_utc)
    endswith(time_reference_utc,"Z") || throw(ArgumentError("UTC reference must end in Z."))
    reference = DateTime(time_reference_utc[1:end-1])
    rows = Tuple[]
    for event in result.events
        geographic = local_to_geographic(frame,event.source_km)
        push!(rows,(event.event_id,geographic.latitude_deg,geographic.longitude_deg,event.source_km[3],
                    event.origin_s,_origin_utc(reference,event.origin_s),event.status,
                    hasproperty(event,:component_id) ? event.component_id : nothing,
                    _depth_result_values(event)...))
    end
    mkpath(directory)
    _write_result_csv(joinpath(directory,"geographic.csv"),
        ("event_id","latitude_deg","longitude_deg","depth_km","origin_s","origin_utc","status","component_id",
         _depth_result_header(result.events)...),rows)
    return joinpath(directory,"geographic.csv")
end
