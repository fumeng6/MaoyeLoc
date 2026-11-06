"Select and aggregate station measurements into the four original event magnitude slots."
function event_magnitudes(events,rows,settings;frame=nothing)
    selection=get(settings,"selection",Dict())
    aggregation=get(settings,"aggregation",Dict())
    method=Symbol(get(aggregation,"method","median"))
    quantization=Symbol(get(aggregation,"quantization","native"))
    channels=get(settings,"channels",Dict())
    event_ids=Set(event.event_id for event in events)
    all(row.event_id in event_ids for row in rows) ||
        throw(ArgumentError("Magnitude measurements must retain an event ID from the supplied catalog."))
    haskey(selection,"fmag2_region") && frame===nothing &&
        throw(ArgumentError("FMAG2 regional selection requires an explicit geographic frame."))
    selected_rows=NamedTuple[]
    summaries=NamedTuple[]
    preferred=NamedTuple[]
    display=NamedTuple[]
    external_records=NamedTuple[]
    slots=(:coda1,:coda2,:amplitude1,:amplitude2)
    rows_by_event=Dict{String,Vector{Int}}()
    for (index,row) in enumerate(rows)
        push!(get!(rows_by_event,row.event_id,Int[]),index)
    end
    for event in events
        geographic=frame===nothing ? nothing : local_to_geographic(frame,event.source_km)
        event_selected=NamedTuple[]
        event_summaries=NamedTuple[]
        for index in get(rows_by_event,event.event_id,Int[])
            row=rows[index]
            choices=magnitude_selection(row,get(channels[row.station_id],"component",""),selection;
                latitude_deg=geographic===nothing ? nothing : geographic.latitude_deg,
                longitude_deg=geographic===nothing ? nothing : geographic.longitude_deg)
            for slot in slots
                choice=getproperty(choices,slot)
                weight=magnitude_observation_weight(choice.quality_code,choice.station_weight)
                used=choice.selected && weight.status==:valid && weight.weight>0
                reason=!choice.selected ? choice.reason : weight.status!=:valid || weight.weight==0 ? weight.reason : :used
                push!(event_selected,merge(choice,(;event_id=row.event_id,station_id=row.station_id,input_row=row.input_row,
                    region_status=choices.region_status,observation_weight=weight.weight,
                    weight_status=weight.status,weight_reason=weight.reason,used,unused_reason=used ? nothing : reason)))
            end
        end
        append!(selected_rows,event_selected)
        for slot in slots
            candidates=[r for r in event_selected if r.slot==slot]
            used=[r for r in candidates if r.used]
            result=aggregate_magnitudes([r.magnitude for r in used],[r.observation_weight for r in used];method,quantization)
            configured=get(selection,String(slot),Dict())
            types=unique([r.magnitude_type for r in used])
            push!(event_summaries,merge(result,(;event_id=event.event_id,slot,label=get(configured,"label",""),
                configured_magnitude_type=get(configured,"magnitude_type",slot in (:coda1,:coda2) ? "md1" : "all"),
                magnitude_type=isempty(types) ? nothing : length(types)==1 ? only(types) : :mixed,
                n_selected=count(r->r.selected,candidates),n_unique_station_ids=length(unique([r.station_id for r in used])))))
        end
        append!(summaries,event_summaries)
        external=get(get(settings,"external",Dict()),event.event_id,nothing)
        legacy_external_amplitude=get(get(settings,"legacy_external_amplitude",Dict()),event.event_id,nothing)
        push!(external_records,(;event_id=event.event_id,external,legacy_external_amplitude))
        choice=preferred_magnitude(event_summaries;external,rules=get(get(settings,"preferred",Dict()),"rules",[]))
        push!(preferred,merge((;event_id=event.event_id),choice))
        in_region=_magsel_region(selection,geographic===nothing ? nothing : geographic.latitude_deg,
            geographic===nothing ? nothing : geographic.longitude_deg)==:inside
        displayed=magnitude_display(event_summaries;external,settings=get(settings,"display",Dict()),
            in_fmag2_region=in_region,legacy_external_amplitude)
        append!(display,[merge((;event_id=event.event_id),r) for r in displayed])
    end
    (;events=summaries,measurements=selected_rows,preferred,display,external_records,assumptions=(;
        method,quantization,weight_definition=:rounded_phase_quality_times_independent_station_magnitude_weight,
        count_definition=:positive_weight_original_measurement_cards,
        duplicate_station_measurements=:retained,location_weights_used=false,
        scatter_is_statistical_uncertainty=false,
        preferred_computed_before_display=true,external_scatter_is_not_assumed_zero=true))
end

function _write_event_magnitudes(directory,result)
    _write_result_csv(joinpath(directory,"event_magnitudes.csv"),
        ("event_id","slot","label","configured_magnitude_type","magnitude_type","magnitude","scatter_mag","scatter_kind",
         "method","quantization","n_measurements","n_unique_station_ids","n_selected","total_weight","status","reason"),
        ((r.event_id,r.slot,r.label,r.configured_magnitude_type,r.magnitude_type,r.magnitude,r.scatter_mag,r.scatter_kind,
          r.method,r.quantization,r.n_measurements,r.n_unique_station_ids,r.n_selected,r.total_weight,r.status,r.reason) for r in result.events))
    _write_result_csv(joinpath(directory,"magnitude_selection.csv"),
        ("input_row","event_id","station_id","slot","label","magnitude_type","magnitude","selected","status","reason",
         "quality_code","station_weight","observation_weight","weight_status","weight_reason","used","unused_reason","region_status"),
        ((r.input_row,r.event_id,r.station_id,r.slot,r.label,r.magnitude_type,r.magnitude,r.selected,r.status,r.reason,
          r.quality_code,r.station_weight,r.observation_weight,r.weight_status,r.weight_reason,r.used,r.unused_reason,r.region_status)
         for r in result.measurements))
    _write_result_csv(joinpath(directory,"preferred_magnitudes.csv"),
        ("event_id","slot","rule_index","magnitude","label","n_readings","scatter_mag","scatter_status","status","reason"),
        ((r.event_id,r.slot,r.rule_index,r.magnitude,r.label,r.n_readings,r.scatter_mag,r.scatter_status,r.status,r.reason)
         for r in result.preferred))
    _write_result_csv(joinpath(directory,"magnitude_display.csv"),
        ("event_id","position","slot","source_slot","magnitude","label","n_readings","scatter_mag","scatter_status",
         "status","reason","display_magnitude","display_label","calculated_magnitude","calculated_label"),
        ((r.event_id,r.position,r.slot,r.source_slot,r.magnitude,r.label,r.n_readings,r.scatter_mag,r.scatter_status,
          r.status,r.reason,r.display_magnitude,r.display_label,r.calculated_magnitude,r.calculated_label) for r in result.display))
    open(joinpath(directory,"event_magnitudes.toml"),"w") do io
        TOML.print(io,_result_metadata(result))
    end
end
