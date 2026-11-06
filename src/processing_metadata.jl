_processing_value(table,key,default) = table isa AbstractDict ? get(table,String(key),default) : get(table,key,default)
_processing_text(table,key,default="") = string(_processing_value(table,key,default))

"""
    processing_metadata(events,picks;metadata=nothing,known_station_ids=nothing)

Preserve event and reading labels and count supplied P first-motion labels.
Metadata may be a TOML table or file. An observation_id identifies each
reading; an override without one requires a unique event/station/phase match.
Metadata-only records never create an arrival or a solver row. A nonblank P
first-motion label counts when the station is known and quality_code < 4,
independently of sigma or final fit weight. Labels are preserved without
inferring compression/dilatation or performing focal-mechanism inversion.
Input IDs, arrival times and uncertainties are not modified.
"""
function processing_metadata(events,picks;metadata=nothing,known_station_ids=nothing)
    table=metadata isa AbstractString ? TOML.parsefile(metadata) : metadata===nothing ? (;) : metadata
    table=_processing_value(table,:processing_metadata,table)
    event_overrides=_processing_value(table,:events,[])
    pick_overrides=_processing_value(table,:picks,[])
    event_overrides isa AbstractVector && pick_overrides isa AbstractVector ||
        throw(ArgumentError("processing_metadata.events and .picks must be arrays of tables."))
    catalog,input_rows=collect(events),NamedTuple[row for row in picks]
    ids=String[event.event_id for event in catalog]
    length(unique(ids))==length(ids) || throw(ArgumentError("Duplicate event ID in processing metadata."))
    idset=Set(ids)
    station_set=known_station_ids===nothing ? nothing : Set(String.(collect(known_station_ids)))
    event_tables=Dict{String,Any}()
    for item in event_overrides
        id=String(_processing_value(item,:event_id,""))
        id in idset || throw(ArgumentError("Processing metadata event ID $id is absent from the catalog."))
        haskey(event_tables,id) && throw(ArgumentError("Duplicate event processing metadata for $id."))
        event_tables[id]=item
    end
    station_phase(row)=(String(_processing_value(row,:event_id,"")),String(_processing_value(row,:station_id,"")),String(_processing_value(row,:phase,"")))
    key(row)=(String(_processing_value(row,:event_id,"")),_absolute_observation_key(row))
    original_count=length(input_rows)
    positions=Dict{Tuple,Int}()
    station_phases=Dict{Tuple{String,String,String},Vector{Int}}()
    for (i,row) in enumerate(input_rows)
        identity=key(row)
        haskey(positions,identity) && throw(ArgumentError("Duplicate input processing reading $identity."))
        positions[identity]=i
        push!(get!(station_phases,station_phase(row),Int[]),i)
    end
    pick_tables=Dict{Tuple,Any}()
    unkeyed_tables=Dict{Tuple{String,String,String},Any}()
    unkeyed_order=Tuple{String,String,String}[]
    # Resolve explicit IDs first, independently of override list order.
    for item in pick_overrides
        identity=station_phase(item)
        identity[1] in idset && !isempty(identity[2]) && identity[3] in ("P","S") ||
            throw(ArgumentError("Reading metadata needs existing event_id, explicit station_id and P/S phase."))
        observation_id=_processing_value(item,:observation_id,nothing)
        reading_key=key(item)
        if observation_id===nothing
            haskey(unkeyed_tables,identity) && throw(ArgumentError("Duplicate unkeyed reading processing metadata $identity."))
            unkeyed_tables[identity]=item
            push!(unkeyed_order,identity)
            continue
        end
        haskey(pick_tables,reading_key) && throw(ArgumentError("Duplicate reading processing metadata $reading_key."))
        pick_tables[reading_key]=item
        if haskey(positions,reading_key)
            station_phase(input_rows[positions[reading_key]])==identity ||
                throw(ArgumentError("Processing observation_id $observation_id disagrees with its input station/phase."))
        else
            push!(input_rows,(event_id=identity[1],station_id=identity[2],phase=identity[3],observation_id=String(observation_id),metadata_only=true))
            positions[reading_key]=length(input_rows)
            push!(get!(station_phases,identity,Int[]),length(input_rows))
        end
    end
    for identity in unkeyed_order
        matches=get(station_phases,identity,Int[])
        length(matches)>1 && throw(ArgumentError("Ambiguous unkeyed reading processing metadata $identity; specify observation_id."))
        if isempty(matches)
            push!(input_rows,(event_id=identity[1],station_id=identity[2],phase=identity[3],metadata_only=true))
            positions[key(last(input_rows))]=length(input_rows)
            station_phases[identity]=[length(input_rows)]
        end
    end
    text_fields=(:remark,:analyst_remark,:program_remark,:data_source,:duration_data_source,:amplitude_data_source,:source,:analyst,:authority,:version,
        :human_review,:processing_domain,:processing_version,:review_status)
    function attributes(row,override;reading=false)
        fields=(;(name=>_processing_text(override,name,get(row,name,"")) for name in text_fields)...)
        tags=_processing_value(override,:review_tags,get(row,:review_tags,String[]))
        tags isa AbstractVector && all(tag->tag isa AbstractString,tags) || throw(ArgumentError("review_tags must be an array of text labels."))
        (;fields...,review_tags=String.(tags))
    end
    reading_records=NamedTuple[]
    for (index,row) in enumerate(input_rows)
        identity=station_phase(row)
        identity[1] in idset && identity[3] in ("P","S") || throw(ArgumentError("Invalid event/phase identity in processing metadata."))
        override=get(pick_tables,key(row),get(unkeyed_tables,identity,(;)))
        phase_remark=_processing_text(override,:phase_remark,get(row,:phase_remark,identity[3]))
        explicit_motion=_processing_value(override,:first_motion,get(row,:first_motion,nothing))
        motion=explicit_motion===nothing ? "" : string(explicit_motion)
        quality=_processing_value(override,:quality_code,get(row,:quality_code,0))
        quality isa Integer && 0<=quality<=9 || throw(ArgumentError("Processing quality_code must be an integer from 0 through 9."))
        station_known=station_set===nothing ? _processing_value(override,:station_known,get(row,:station_known,true)) : identity[2] in station_set
        station_known isa Bool || throw(ArgumentError("station_known must be an explicit boolean when supplied."))
        metadata_accepted=_processing_value(override,:metadata_accepted,get(row,:metadata_accepted,true))
        metadata_accepted isa Bool || throw(ArgumentError("metadata_accepted must be boolean."))
        metadata_exclusion_reason=_processing_value(override,:metadata_exclusion_reason,get(row,:metadata_exclusion_reason,nothing))
        reason=!metadata_accepted ? :excluded_metadata : identity[3]!="P" ? :not_p_phase : !station_known ? :unknown_station : quality>=4 ? :quality_code_disabled :
            isempty(strip(phase_remark)) ? :p_remark_blank : isempty(strip(motion)) ? :first_motion_blank : nothing
        metadata_only=get(row,:metadata_only,hasproperty(row,:arrival_s) && row.arrival_s===nothing)
        push!(reading_records,(event_id=identity[1],observation_id=get(row,:observation_id,nothing),station_id=identity[2],phase=identity[3],
            input_row=index<=original_count ? index : nothing,metadata_only,metadata_accepted,metadata_exclusion_reason,
            first_motion=motion,phase_remark,quality_code=Int(quality),station_known,
            nfm_counted=reason===nothing,nfm_exclusion_reason=reason,attributes(row,override;reading=true)...,
            equation_weight=metadata_only ? nothing : get(row,:equation_weight,nothing),raw_line=_processing_text(override,:raw_line,get(row,:raw_line,"")),
            source_file=_processing_text(override,:source_file,get(row,:source_file,"")),
            source_line=_processing_value(override,:source_line,get(row,:source_line,nothing))))
    end
    event_records=NamedTuple[]
    for event in catalog
        id=String(event.event_id)
        rows=filter(row->row.event_id==id,reading_records)
        counts=Dict{String,Int}()
        source_order=String[]
        for row in rows
            haskey(counts,row.data_source) || push!(source_order,row.data_source)
            counts[row.data_source]=get(counts,row.data_source,0)+1
        end
        weighted_rows=filter(row->row.equation_weight!==nothing && !row.metadata_only,rows)
        weighted=isempty(weighted_rows) ? nothing : Dict{String,Int}()
        primary=nothing
        order=String[]
        if weighted!==nothing
            for row in weighted_rows
                row.station_known && row.equation_weight>=.1 || continue
                haskey(weighted,row.data_source) || push!(order,row.data_source)
                weighted[row.data_source]=get(weighted,row.data_source,0)+1
            end
            !isempty(order) && (primary=order[argmax([weighted[source] for source in order])])
        end
        push!(event_records,(event_id=id,nfm=count(row->row.nfm_counted,rows),reported_reading_count=length(rows),
            metadata_only_reading_count=count(row->row.metadata_only,rows),attributes(event,get(event_tables,id,(;)))...,
            data_source_counts=counts,data_source_order=source_order,weighted_data_source_counts=weighted,
            weighted_data_source_order=order,weighted_source_unavailable_reading_count=count(row->!row.metadata_only && row.equation_weight===nothing,rows),
            primary_weighted_data_source=primary,
            raw_header=_processing_text(get(event_tables,id,(;)),:raw_header,get(event,:raw_header,""))))
    end
    (;events=event_records,picks=reading_records,
      assumptions=(nfm_rule=:matched_p_nonblank_first_motion_quality_below_four,
        nfm_uses_sigma=false,nfm_uses_final_weight=false,nfm_requires_arrival_time=false,
        source_count_basis=:reported_phase_records,weighted_source_threshold=.1,
        first_motion_interpreted=false,focal_mechanism_inverted=false,event_ids_preserved=true,
        reading_identity=:event_observation_id_with_unique_station_phase_match))
end
