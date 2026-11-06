# HypoDD 1.3 getdata/cluster1 selection. Strong edges establish cluster identity;
# all eligible observations within that cluster are then read again, including
# weak links and stations excluded by the first, whole-catalog DIST selection.
function _dd_legacy_station_distance(center, station)
    rad, halfpi = 1.745329e-2, 1.570796
    alat, blat = center.latitude_deg*rad, station.latitude_deg*rad
    acol = halfpi-atan(.993231*tan(alat))
    bcol = halfpi-atan(.993231*tan(blat))
    dlon = (station.longitude_deg-center.longitude_deg)*rad
    cosine = sin(acol)*sin(bcol)*cos(dlon)+cos(acol)*cos(bcol)
    # atan2 of cross/dot is the same spherical angle, without acos's
    # cancellation at the short station distances used by ph2dt.
    sine=hypot(sin(bcol)*sin(dlon),sin(acol)*cos(bcol)-cos(acol)*sin(bcol)*cos(dlon))
    angle = atan(sine,cosine)
    colat = halfpi-(alat+blat)/2
    # The original integer expression 1/3 is zero in Fortran. Use the intended
    # geographic radius correction rather than reproduce that truncation.
    radius = 6371.227*(1+3.37853e-3*(1/3-cos(colat)^2))
    angle*radius
end

function _dd_selection(catalog, stations, observations; selection, frame=nothing,initial_exclusion_reasons=nothing)
    initial_exclusion_reasons===nothing || length(initial_exclusion_reasons)==length(observations) ||
        throw(ArgumentError("Initial DD distance selection must match all observations."))
    values = Dict(String(key)=>value for (key,value) in pairs(selection))
    allowed = ("idat","ipha","max_station_distance_km","event_ids","obscc","obsct","cluster_id")
    all(key in allowed for key in keys(values)) || throw(ArgumentError("Unknown DD selection option."))
    idat,ipha = get(values,"idat",3),get(values,"ipha",3)
    idat isa Integer && idat in (1,2,3) || throw(ArgumentError("DD selection idat must be 1, 2 or 3; synthetic IDAT=0 is separate."))
    ipha isa Integer && ipha in (1,2,3) || throw(ArgumentError("DD selection ipha must be 1, 2 or 3."))
    obscc,obsct,cid = get(values,"obscc",0),get(values,"obsct",0),get(values,"cluster_id",0)
    all(value isa Integer && value>=0 for value in (obscc,obsct,cid)) ||
        throw(ArgumentError("DD OBSCC/OBSCT and cluster_id must be nonnegative integers."))
    maxdist = get(values,"max_station_distance_km",Inf)
    maxdist isa Real && !isnan(maxdist) && maxdist>=0 || throw(ArgumentError("DD station distance must be nonnegative km."))
    (frame===nothing || frame isa LocalFrame) && (!isfinite(maxdist) || frame isa LocalFrame) ||
        throw(ArgumentError("Finite DD station distance requires the catalog's LocalFrame."))
    requested = get(values,"event_ids",String[])
    requested isa AbstractVector && all(id isa AbstractString for id in requested) ||
        throw(ArgumentError("DD event_ids must be an explicit vector of strings."))
    ids = String[event.event_id for event in catalog]
    length(unique(ids))==length(ids) || throw(ArgumentError("Duplicate DD catalog event ID."))
    index = Dict(id=>i for (i,id) in enumerate(ids))
    length(unique(requested))==length(requested) && all(haskey(index,id) for id in requested) ||
        throw(ArgumentError("DD selected event IDs must be unique and present in the catalog."))
    selected_ids = isempty(requested) ? Set(ids) : Set(String.(requested))
    selected_events = findall(id->id in selected_ids,ids)
    geographic_events = frame===nothing ? nothing : [local_to_geographic(frame,event.source_km) for event in catalog]
    geographic_stations = frame===nothing ? nothing : Dict(id=>local_to_geographic(frame,xyz) for (id,xyz) in stations)
    function center(indices)
        (geographic_events===nothing || isempty(indices)) && return nothing
        (latitude_deg=sum(geographic_events[i].latitude_deg for i in indices)/length(indices),
         longitude_deg=sum(geographic_events[i].longitude_deg for i in indices)/length(indices))
    end
    function station_distance(center_value,row)
        !isfinite(maxdist) && return nothing
        _dd_legacy_station_distance(center_value,geographic_stations[row.station_id])
    end
    global_center = center(selected_events)
    base_reasons = Union{Nothing,Symbol}[nothing for _ in observations]
    initial_reasons = Union{Nothing,Symbol}[nothing for _ in observations]
    initial_rows = Int[]
    pair_counts = Dict{Tuple{Int,Int},NamedTuple{(:cc,:ct),Tuple{Int,Int}}}()
    matched = falses(length(catalog))
    for (r,row) in enumerate(observations)
        haskey(index,row.event_i) && haskey(index,row.event_j) || throw(ArgumentError("Unknown DD observation event."))
        row.event_i!=row.event_j || throw(ArgumentError("A DD observation needs distinct events."))
        row.source in ("CC","CT") && row.phase in ("P","S") || throw(ArgumentError("DD rows require CT/CC and P/S."))
        reason = if !(row.event_i in selected_ids && row.event_j in selected_ids)
            :event_not_selected
        elseif (idat==1 && row.source!="CC") || (idat==2 && row.source!="CT")
            :data_type_not_selected
        elseif (ipha==1 && row.phase!="P") || (ipha==2 && row.phase!="S")
            :phase_not_selected
        elseif !haskey(stations,row.station_id)
            :unknown_station
        elseif initial_exclusion_reasons!==nothing
            initial_exclusion_reasons[r]
        else
            nothing
        end
        base_reasons[r] = initial_reasons[r] = reason
        reason===nothing || continue
        distance = station_distance(global_center,row)
        if distance!==nothing && distance>maxdist
            initial_reasons[r] = :station_distance_initial
            continue
        end
        push!(initial_rows,r)
        i,j = minmax(index[row.event_i],index[row.event_j])
        matched[i] = matched[j] = true
        counts = get(pair_counts,(i,j),(cc=0,ct=0))
        pair_counts[(i,j)] = row.source=="CC" ? (cc=counts.cc+1,ct=counts.ct) : (cc=counts.cc,ct=counts.ct+1)
    end

    threshold = (idat==2 ? 0 : obscc)+(idat==1 ? 0 : obsct)
    clustering = threshold>0
    groups = Vector{Int}[]
    active = findall(matched)
    if !clustering
        # Zero OBSCC/OBSCT means "no minimum link count", not "one datum for
        # every event that has any observation". Union-find with threshold 0
        # would also join pairs that have no observation. Link only pairs that
        # retain at least one row, in catalog order, so each component can
        # carry its own reference.
        parent = collect(eachindex(catalog))
        function component_root(i)
            while parent[i]!=i
                parent[i] = parent[parent[i]]
                i = parent[i]
            end
            i
        end
        for ((i,j),counts) in pair_counts
            counts.cc+counts.ct>0 || continue
            parent[component_root(j)] = component_root(i)
        end
        roots = unique(component_root(i) for i in active)
        sort!(roots;by=r->minimum(i for i in active if component_root(i)==r))
        for r in roots
            group = [i for i in active if component_root(i)==r]
            sort!(group)
            push!(groups,group)
        end
    else
        # Original pair scanning uses numerically sorted CUSP IDs. General
        # string IDs use lexical order; numeric spellings/leading zeros survive.
        numbers = tryparse.(BigInt,ids)
        ordered = all(!isnothing,numbers) ? sort(active;by=i->numbers[i]) : sort(active;by=i->ids[i])
        parent = collect(eachindex(catalog))
        seed = fill(typemax(Int),length(catalog))
        function root(i)
            while parent[i]!=i
                parent[i] = parent[parent[i]]
                i = parent[i]
            end
            i
        end
        edge_order = 0
        for a in 2:length(ordered), b in 1:a-1
            i,j = ordered[a],ordered[b]
            counts = get(pair_counts,minmax(i,j),(cc=0,ct=0))
            counts.cc+counts.ct>=threshold || continue
            edge_order += 1
            ri,rj = root(i),root(j)
            parent[rj] = ri
            seed[ri] = min(seed[ri],seed[rj],edge_order)
        end
        roots = unique(root(i) for i in active if seed[root(i)]<typemax(Int))
        sort!(roots;by=i->seed[i])
        for r in roots
            push!(groups,[i for i in active if root(i)==r])
        end
        # Preserve the original largest-first ordering, including its <= tie
        # swaps, so CID remains reproducible for equal-size numeric-ID groups.
        for i in 1:length(groups)-1, j in i+1:length(groups)
            if length(groups[i])<=length(groups[j])
                groups[i],groups[j] = groups[j],groups[i]
            end
        end
    end
    cid<=length(groups) || throw(ArgumentError("DD cluster_id exceeds the number of selected clusters."))
    selected_group_indices = cid==0 ? collect(eachindex(groups)) : [Int(cid)]
    selected_cluster_set = Set(selected_group_indices)
    membership = zeros(Int,length(catalog))
    for (g,group) in enumerate(groups),i in group
        membership[i] = g
    end
    cluster_centers = center.(groups)
    # hypoDD skips getdata again when clustering is disabled, or when one
    # cluster already contains every matched event. The initial center can
    # include unmatched catalog entries, so skipping is physically observable.
    reread=[clustering && !(length(groups)==1 && length(group)==length(active)) for group in groups]
    distance_centers=[reread[g] ? cluster_centers[g] : global_center for g in eachindex(groups)]
    group_rows = [Int[] for _ in groups]
    excluded = NamedTuple[]
    for (r,row) in enumerate(observations)
        i,j = index[row.event_i],index[row.event_j]
        g = membership[i]
        reason = base_reasons[r]
        if reason===nothing
            if g==0 || membership[j]==0
                reason = initial_reasons[r]===nothing ? :not_clustered : initial_reasons[r]
            elseif g!=membership[j]
                reason = :different_clusters
            else
                distance = station_distance(distance_centers[g],row)
                if !reread[g] && initial_reasons[r]!==nothing
                    reason=initial_reasons[r]
                elseif distance!==nothing && distance>maxdist
                    reason = :station_distance_cluster
                else
                    push!(group_rows[g],r)
                    !(g in selected_cluster_set) && (reason=:cluster_not_selected)
                end
            end
        end
        reason===nothing || push!(excluded,merge(row,(input_row=r,reason,cluster_id=g>0 && g==membership[j] ? g : 0)))
    end
    event_statuses = [!(ids[i] in selected_ids) ? :event_not_selected : !matched[i] ? :no_selected_observations :
        membership[i]==0 ? :not_clustered : !(membership[i] in selected_cluster_set) ? :cluster_not_selected : :selected
        for i in eachindex(catalog)]
    settings = (idat=Int(idat),ipha=Int(ipha),max_station_distance_km=Float64(maxdist),
        event_ids=String.(requested),obscc=Int(obscc),obsct=Int(obsct),cluster_id=Int(cid),
        effective_threshold=threshold,clustering,
        distance_formula="geocentric_latitude_radius_with_float_one_third")
    counts = [(event_i=ids[i],event_j=ids[j],cc=n.cc,ct=n.ct,total=n.cc+n.ct,
        strong=clustering && n.cc+n.ct>=threshold) for ((i,j),n) in sort!(collect(pair_counts);by=first)]
    (;groups,group_rows,selected_group_indices,membership,event_statuses,excluded,catalog_event_ids=ids,
      initial_rows,initial_reasons,pair_counts=counts,reread,
      centers=(global_center=global_center,clusters=cluster_centers,distance_centers),settings)
end
