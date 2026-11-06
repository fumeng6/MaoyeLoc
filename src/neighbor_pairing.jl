"""
    pair_neighbors(events, stations, picks; frame::LocalFrame, min_quality=0.,
        max_station_distance_km, max_separation_km, max_neighbors,
        min_links, min_observations, max_observations,
        coarse_velocity_km_s=(P=4., S=2.3), coarse_padding_s=.5)

Select neighbor event pairs and bounded common-station observations with the
strong-link search of Waldhauser's ph2dt (2001). This is an alternative to
`pair_catalog`'s common-path-ratio selection, not a different location solver.
Inputs are native records: finite arrivals and positive sigmas, nonnegative
`pair_quality` (dimensionless) and optional Boolean `important`. Importance
bypasses only min_quality, gives priority in the observation budget and adds
to that budget; it does not alter the measurement or its sigma. No negative
weight encoding is used. Mean pair quality is retained as original_quality
for an explicitly requested empirical DD estimator, never inferred as sigma.

Event separation is 3-D local km; station distance is from the geographic
pair center using the documented spheroidal approximation. frame is explicit.
Coarse exclusion uses abs((arrival_i-origin_i)-(arrival_j-origin_j)) >
separation/coarse_velocity + coarse_padding, an empirical filter whose P/S
speed and padding are stated, not a guarantee of physical common paths.
Strong links count common eligible readings before MAXOBS; directed searches
also count already visited reverse strong links. CT is the complete arrival
difference i-j, with independent-pick sigma propagation. Shared-pick CT rows
remain correlated. All pair, reading and event dispositions are returned.
"""
function pair_neighbors(events,stations,picks;frame::LocalFrame,min_quality=0.,
        max_station_distance_km,max_separation_km,max_neighbors,min_links,min_observations,max_observations,
        coarse_velocity_km_s=(P=4.,S=2.3),coarse_padding_s=.5)
    min_quality isa Real && isfinite(min_quality) && min_quality>=0 || throw(ArgumentError("min_quality must be nonnegative."))
    max_station_distance_km isa Real && !isnan(max_station_distance_km) && max_station_distance_km>=0 &&
        max_separation_km isa Real && isfinite(max_separation_km) && max_separation_km>=0 ||
        throw(ArgumentError("Pairing distances must be nonnegative km; source separation must be finite."))
    all(value->value isa Integer && value>=0,(max_neighbors,min_links,min_observations,max_observations)) ||
        throw(ArgumentError("MAXNGH/MINLNK/MINOBS/MAXOBS must be nonnegative integers."))
    all(v isa Real && isfinite(v) && v>0 for v in (coarse_velocity_km_s.P,coarse_velocity_km_s.S)) &&
        coarse_padding_s isa Real && isfinite(coarse_padding_s) && coarse_padding_s>=0 ||
        throw(ArgumentError("Coarse P/S speeds must be positive km/s and padding nonnegative seconds."))
    catalog,phase_rows=collect(events),collect(picks)
    ids=String[event.event_id for event in catalog]
    length(unique(ids))==length(ids) || throw(ArgumentError("Duplicate pairing event ID."))
    index=Dict(id=>i for (i,id) in enumerate(ids))
    eligible=[Int[] for _ in catalog]
    input_count=zeros(Int,length(catalog))
    quality=Float64[]
    relative_time=Float64[]
    pick_records=NamedTuple[]
    known=Set{Tuple{String,String,String}}()
    for (p,pick) in enumerate(phase_rows)
        haskey(index,pick.event_id) || throw(ArgumentError("Unknown event in pairing picks."))
        pick.phase in ("P","S") && isfinite(pick.arrival_s) && isfinite(pick.sigma_s) && pick.sigma_s>0 ||
            throw(ArgumentError("Pairing picks need P/S, finite complete arrivals and positive independent sigmas."))
        key=(String(pick.event_id),String(pick.station_id),String(pick.phase))
        key in known && throw(ArgumentError("Neighbor pairing requires one selected reading per event/station/phase."))
        push!(known,key)
        q=get(pick,:pair_quality,nothing)
        q isa Real && !(q isa Bool) && isfinite(q) && q>=0 ||
            throw(ArgumentError("Neighbor pairing requires explicit nonnegative dimensionless pair_quality, independent of sigma and quality_code."))
        important=get(pick,:important,false)
        important isa Bool || throw(ArgumentError("important must be a Boolean."))
        push!(quality,Float64(q))
        e=index[pick.event_id]
        input_count[e]+=1
        tt=pick.arrival_s-catalog[e].origin_s
        isfinite(tt) || throw(ArgumentError("Relative phase travel time must be finite seconds."))
        push!(relative_time,Float64(tt))
        retained=important || q>=min_quality
        retained && push!(eligible[e],p)
        push!(pick_records,(input_row=p,event_id=ids[e],station_id=pick.station_id,phase=pick.phase,
            pair_quality=Float64(q),important,retained_by_quality=retained,
            reason=retained ? nothing : :below_min_quality))
    end
    selected=[length(rows)>=min_observations for rows in eligible]
    active=findall(selected)
    geography=[local_to_geographic(frame,event.source_km) for event in catalog]
    coordinates=[collect(event.source_km) for event in catalog]
    phase_index=[Dict((phase_rows[p].station_id,phase_rows[p].phase)=>p for p in rows) for rows in eligible]
    station_geography=Dict(id=>local_to_geographic(frame,xyz) for (id,xyz) in stations)
    # 0=not examined, 1=strong, 2=weak. Symmetric storage avoids a second
    # evaluation/write, while directed visits still count reverse strong links.
    examined=Dict{Tuple{Int,Int},Bool}()
    pairs,readings,observations,visits=NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    neighbor_count,candidate_count,retained_count,strong_pair_count=zeros(Int,length(catalog)),zeros(Int,length(catalog)),zeros(Int,length(catalog)),zeros(Int,length(catalog))
    stop_reasons=fill(:not_selected,length(catalog))
    for i in active
        # Compute each distance once and keep the catalog index as its tie-break.
        # MAXNGH still counts strong links, so every neighbor remains available.
        neighbors=sort!([(norm(coordinates[i]-coordinates[j]),j) for j in active if j!=i])
        stop_reasons[i]=:all_neighbors_examined
        for (separation,j) in neighbors
            if neighbor_count[i]>=max_neighbors
                stop_reasons[i]=:max_strong_neighbors
                break
            end
            key=minmax(i,j)
            if haskey(examined,key)
                strong=examined[key]
                strong && (neighbor_count[i]+=1)
                push!(visits,(event_i=ids[i],event_j=ids[j],separation_km=separation,
                    action=strong ? :reverse_strong_counted : :reverse_weak_skipped,neighbor_count=neighbor_count[i]))
                continue
            end
            if separation>max_separation_km
                stop_reasons[i]=:max_pair_separation
                push!(visits,(event_i=ids[i],event_j=ids[j],separation_km=separation,
                    action=:pair_separation_stop,neighbor_count=neighbor_count[i]))
                break
            end
            candidate_count[i]+=1
            candidate_count[j]+=1
            center=(latitude_deg=(geography[i].latitude_deg+geography[j].latitude_deg)/2,
                longitude_deg=(geography[i].longitude_deg+geography[j].longitude_deg)/2)
            right=phase_index[j]
            row_indices,valid=Int[],Int[]
            for a in eligible[i]
                pick=phase_rows[a]
                b=get(right,(pick.station_id,pick.phase),nothing)
                b===nothing && continue
                other=phase_rows[b]
                distance=haskey(stations,pick.station_id) ? _dd_legacy_station_distance(center,station_geography[pick.station_id]) : nothing
                cutoff=separation/(pick.phase=="P" ? coarse_velocity_km_s.P : coarse_velocity_km_s.S)+coarse_padding_s
                time_difference=relative_time[a]-relative_time[b]
                important=get(pick,:important,false) || get(other,:important,false)
                reason=distance===nothing ? :unknown_station : distance>max_station_distance_km ? :station_distance_cutoff :
                    abs(time_difference)>cutoff ? :differential_time_outlier : nothing
                push!(readings,(event_i=ids[i],event_j=ids[j],station_id=pick.station_id,phase=pick.phase,
                    input_pick_i=a,input_pick_j=b,travel_time_i_s=relative_time[a],travel_time_j_s=relative_time[b],
                    original_quality_i=quality[a],original_quality_j=quality[b],quality=(quality[a]+quality[b])/2,
                    important,station_distance_km=distance,separation_km=separation,coarse_cutoff_s=cutoff,
                    relative_time_difference_s=time_difference,eligible=reason===nothing,used=false,reason))
                r=length(readings)
                push!(row_indices,r)
                reason===nothing && push!(valid,r)
            end
            count_valid=length(valid)
            important_count=count(r->readings[r].important,valid)
            strong=count_valid>=min_links
            retained=count_valid>=min_observations
            budget=count_valid>max_observations ? min(max_observations+important_count,count_valid) : count_valid
            # Importance gives the original zero sorting distance. Resolve an
            # exact zero-distance tie in favor of the important row, so a
            # station exactly at the pair center cannot displace it.
            order=sort(valid;by=r->(readings[r].important ? 0. : readings[r].station_distance_km,
                readings[r].important ? 0 : 1,r))
            kept=retained ? Set(order[1:budget]) : Set{Int}()
            for r in valid
                row=readings[r]
                used=r in kept
                reason=used ? nothing : retained ? :observation_budget : :insufficient_pair_observations
                readings[r]=merge(row,(;used,reason))
            end
            for r in order[1:(retained ? budget : 0)]
                row=readings[r]
                a,b=row.input_pick_i,row.input_pick_j
                push!(observations,(event_i=ids[i],event_j=ids[j],station_id=row.station_id,phase=row.phase,
                    dt_s=(catalog[i].origin_s-catalog[j].origin_s)+row.relative_time_difference_s,
                    sigma_s=hypot(phase_rows[a].sigma_s,phase_rows[b].sigma_s),source="CT",
                    original_quality=row.quality,quality=row.quality,important=row.important,
                    travel_time_i_s=row.travel_time_i_s,travel_time_j_s=row.travel_time_j_s,
                    original_quality_i=row.original_quality_i,original_quality_j=row.original_quality_j,
                    input_pick_i=a,input_pick_j=b))
            end
            examined[key]=strong
            strong && (neighbor_count[i]+=1;strong_pair_count[i]+=1;strong_pair_count[j]+=1)
            retained && (retained_count[i]+=1;retained_count[j]+=1)
            push!(pairs,(event_i=ids[i],event_j=ids[j],separation_km=separation,
                center_latitude_deg=center.latitude_deg,center_longitude_deg=center.longitude_deg,
                common_phase_count=length(row_indices),eligible_phase_count=count_valid,important_count,
                retained_phase_count=retained ? budget : 0,strong,retained,
                reason=retained ? nothing : :insufficient_pair_observations))
            push!(visits,(event_i=ids[i],event_j=ids[j],separation_km=separation,
                action=strong ? :new_strong_pair : :new_weak_pair,neighbor_count=neighbor_count[i]))
        end
    end
    event_records=[(event_id=ids[i],input_phase_count=input_count[i],eligible_phase_count=length(eligible[i]),
        selected=selected[i],candidate_pair_count=candidate_count[i],retained_pair_count=retained_count[i],
        strong_pair_count=strong_pair_count[i],strong_neighbor_count=neighbor_count[i],search_stop_reason=stop_reasons[i],
        status=!selected[i] ? :below_event_min_observations : retained_count[i]>0 ? :paired : :no_retained_pair)
        for i in eachindex(catalog)]
    (;observations,pairs,readings,events=event_records,pick_records,neighbor_visits=visits,
      criteria=(;min_quality,max_station_distance_km,max_separation_km,max_neighbors,min_links,min_observations,max_observations,coarse_velocity_km_s,coarse_padding_s),
      source_distance_definition=:symmetric_local_cartesian_3d_km,
      station_distance_definition=:geographic_pair_center_spheroidal_approximation,
      observation_covariance=:shared_pick_correlations_not_removed)
end
