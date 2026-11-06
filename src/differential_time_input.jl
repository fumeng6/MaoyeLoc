"""
    catalog_differences(picks, event_pairs)

Form CT differences only for supplied `(event_i, event_j)` string-ID pairs and
common station/phase picks. Arrivals already use one shared time reference;
the complete difference is `arrival_i - arrival_j`, without subtracting an
origin-time term. Pair direction is preserved and repeated/reversed pairs are
rejected rather than counted twice. Unmatched pairs are reported explicitly.

The propagated sigma is `hypot(sigma_i, sigma_j)` under the stated assumption
of independent pick errors. Differences that share a pick are correlated;
this function does not estimate or remove that covariance. No distance or
quality gate, automatic pair selection, or relaxed matching is performed.
"""
function catalog_differences(picks, event_pairs)
    by_event = Dict{String,Dict{Tuple{String,String},NamedTuple}}()
    for pick in picks
        pick.event_id isa AbstractString && pick.station_id isa AbstractString ||
            throw(ArgumentError("Pick event and station IDs must be strings."))
        pick.phase in ("P", "S") || throw(ArgumentError("Phase must be P or S."))
        isfinite(pick.arrival_s) && isfinite(pick.sigma_s) && pick.sigma_s > 0 ||
            throw(ArgumentError("Pick arrivals must be finite and sigma_s finite and positive."))
        event = get!(by_event, String(pick.event_id), Dict{Tuple{String,String},NamedTuple}())
        key = (String(pick.station_id), String(pick.phase))
        haskey(event, key) && throw(ArgumentError("CT differences require one reading per event/station/phase; select the intended observation_id before pairing."))
        event[key] = pick
    end
    observations, unused_pairs = NamedTuple[], NamedTuple[]
    seen = Set{Tuple{String,String}}()
    for pair in event_pairs
        length(pair) == 2 || throw(ArgumentError("Each requested event pair must have two IDs."))
        event_i, event_j = pair
        event_i isa AbstractString && event_j isa AbstractString && event_i != event_j ||
            throw(ArgumentError("Requested pair IDs must be distinct strings."))
        event_i, event_j = String(event_i), String(event_j)
        key = minmax(event_i, event_j)
        key in seen && throw(ArgumentError("Repeated or reversed requested pair $key."))
        push!(seen, key)
        if !haskey(by_event, event_i) || !haskey(by_event, event_j)
            push!(unused_pairs, (; event_i, event_j, reason=:no_picks_for_event))
            continue
        end
        left, right = by_event[event_i], by_event[event_j]
        common = sort!(collect(intersect(keys(left), keys(right))))
        if isempty(common)
            push!(unused_pairs, (; event_i, event_j, reason=:no_common_station_phase))
            continue
        end
        for (station_id, phase) in common
            i, j = left[(station_id, phase)], right[(station_id, phase)]
            push!(observations, (; event_i, event_j, station_id, phase,
                dt_s=Float64(i.arrival_s-j.arrival_s),
                sigma_s=hypot(Float64(i.sigma_s), Float64(j.sigma_s)), source="CT",
                arrival_i_s=Float64(i.arrival_s), arrival_j_s=Float64(j.arrival_s),
                pick_sigma_i_s=Float64(i.sigma_s), pick_sigma_j_s=Float64(j.sigma_s),
                sigma_assumption=:independent_pick_errors))
        end
    end
    return (; observations, unused_pairs)
end
