"""
    pair_catalog(events, stations, picks; max_separation_km,
                 max_separation_ratio, min_common_stations)

Make CT observations from catalog picks using explicit geometry criteria.
Candidate events lie within max_separation_km (3-D). A common station/phase
is retained only when separation / min(source-station distances) is no larger
than max_separation_ratio. Require min_common_stations distinct retained station
coordinates; colocated components do not increase this geometry count.
Thus one fixed separation in km is not the sole common-path criterion.

A spatial grid limits candidate searches; criteria are never relaxed to form
more pairs. Returned observations follow complete-arrival i-j and independent
pick sigma propagation. Shared-pick differences remain correlated. Metadata
reports candidate/retained pairs and every event's pairing status. The geometry
criteria alone do not establish cancellation of unknown 3-D path errors.
"""
function pair_catalog(events, stations, picks; max_separation_km::Real,
                      max_separation_ratio::Real, min_common_stations::Integer)
    isfinite(max_separation_km) && max_separation_km > 0 &&
        isfinite(max_separation_ratio) && 0 < max_separation_ratio < 1 && min_common_stations > 0 ||
        throw(ArgumentError("Require positive separation, 0 < separation ratio < 1, and positive common-station count."))
    catalog, observations = collect(events), collect(picks)
    ids = [String(event.event_id) for event in catalog]
    length(unique(ids)) == length(ids) || throw(ArgumentError("Duplicate catalog event ID."))
    index = Dict(id=>i for (i,id) in enumerate(ids))
    available = [Set{Tuple{String,String}}() for _ in catalog]
    for event in catalog
        length(event.source_km)==3 && all(isfinite,event.source_km) || throw(ArgumentError("Invalid event coordinate."))
    end
    for pick in observations
        haskey(index,pick.event_id) && haskey(stations,pick.station_id) || throw(ArgumentError("Unknown event/station in picks."))
        key = (String(pick.station_id),String(pick.phase))
        key in available[index[pick.event_id]] && throw(ArgumentError("CT pairing requires one reading per event/station/phase; select the intended observation_id before pairing."))
        push!(available[index[pick.event_id]],key)
    end
    grid = Dict{NTuple{3,Int},Vector{Int}}()
    candidate_counts, retained_counts = zeros(Int,length(catalog)),zeros(Int,length(catalog))
    pair_records, selected_pairs = NamedTuple[],Tuple{String,String}[]
    eligible = Set{NTuple{4,String}}()
    for (i,event) in enumerate(catalog)
        cell = ntuple(axis->floor(Int,event.source_km[axis]/max_separation_km),3)
        for ax in -1:1, ay in -1:1, az in -1:1
            neighbor = cell .+ (ax,ay,az)
            for j in get(grid,neighbor,Int[])
                other = catalog[j]
                separation = norm(collect(event.source_km) .- collect(other.source_km))
                separation <= max_separation_km || continue
                candidate_counts[i]+=1; candidate_counts[j]+=1
                common = intersect(available[j],available[i])
                paths = Tuple{String,String}[]
                ratios = Float64[]
                for (station_id,phase) in sort!(collect(common))
                    station = stations[station_id]
                    distance = min(norm(collect(event.source_km) .- collect(station)),
                                   norm(collect(other.source_km) .- collect(station)))
                    if distance > 0 && separation/distance <= max_separation_ratio
                        push!(paths,(station_id,phase))
                        push!(ratios,separation/distance)
                    end
                end
                station_count = length(unique(Tuple(stations[sid]) for (sid,_) in paths))
                retained = station_count >= min_common_stations
                reason = retained ? :retained : isempty(common) ? :no_common_station_phase : :insufficient_common_path_geometry
                push!(pair_records,(event_i=ids[j],event_j=ids[i],separation_km=separation,
                    common_phase_count=length(common),retained_phase_count=retained ? length(paths) : 0,
                    common_station_count=station_count,max_retained_separation_ratio=isempty(ratios) ? nothing : maximum(ratios),
                    retained,reason))
                if retained
                    push!(selected_pairs,(ids[j],ids[i]))
                    retained_counts[i]+=1; retained_counts[j]+=1
                    for (station_id,phase) in paths
                        push!(eligible,(ids[j],ids[i],station_id,phase))
                    end
                end
            end
        end
        push!(get!(grid,cell,Int[]),i)
    end
    generated = catalog_differences(observations,selected_pairs)
    selected = [row for row in generated.observations if (row.event_i,row.event_j,row.station_id,row.phase) in eligible]
    event_status = [(event_id=ids[i],candidate_pair_count=candidate_counts[i],retained_pair_count=retained_counts[i],
                    status=retained_counts[i]>0 ? :paired : isempty(available[i]) ? :no_picks :
                           candidate_counts[i]==0 ? :no_nearby_event : :no_pair_passed_geometry)
                    for i in eachindex(catalog)]
    return (observations=selected,pairs=pair_records,events=event_status,
            criteria=(;max_separation_km,max_separation_ratio,min_common_stations))
end

"""
    differential_geometry(events, stations, differentials; max_separation_km,
        max_separation_ratio, min_common_stations)

Reassess the ORIGINAL differential paths at the supplied result coordinates.
Report each pair's separation and distinct station coordinates still satisfying
the declared common-path geometry. Do not form new pairs, remove observations
or alter numerical solution states. Initial pairing can cease to apply after
a large relocation, even when the numerical iteration reports convergence.
"""
function differential_geometry(events, stations, differentials; max_separation_km,
                               max_separation_ratio,min_common_stations)
    positions = Dict(row.event_id=>row.source_km for row in events)
    paths = Dict{Tuple{String,String},Set{String}}()
    for row in differentials
        key = minmax(String(row.event_i),String(row.event_j))
        push!(get!(paths,key,Set{String}()),String(row.station_id))
    end
    rows = NamedTuple[]
    for (left,right) in sort!(collect(keys(paths)))
        a,b = positions[left],positions[right]
        separation = hypot((a[k]-b[k] for k in 1:3)...)
        kept = Set{NTuple{3,Float64}}()
        max_ratio = 0.0
        for id in paths[(left,right)]
            station = stations[id]
            distance = min(hypot((a[k]-station[k] for k in 1:3)...),
                           hypot((b[k]-station[k] for k in 1:3)...))
            ratio = distance>0 ? separation/distance : Inf
            max_ratio = max(max_ratio,ratio)
            ratio<=max_separation_ratio && push!(kept,Tuple(station))
        end
        satisfied = separation<=max_separation_km && length(kept)>=min_common_stations
        push!(rows,(event_i=left,event_j=right,separation_km=separation,
            max_original_path_ratio=max_ratio,common_path_station_count=length(kept),
            meets_declared_geometry=satisfied))
    end
    return rows
end
