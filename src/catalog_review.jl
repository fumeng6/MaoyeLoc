# Native event/reading edits, actual recalculation and explicit catalog decisions.
function _review_subset(config,ids)
    selected=Set(ids)
    case=deepcopy(config)
    case["task"]="absolute"
    for key in ("absolute_constraints",)
        haskey(case,key) && (case[key]=Dict(id=>row for (id,row) in case[key] if id in selected))
    end
    for table_key in ("processing_metadata",)
        table=get(case,table_key,nothing)
        table===nothing && continue
        for key in ("events","picks")
            haskey(table,key) && (table[key]=[row for row in table[key] if row["event_id"] in selected])
        end
    end
    depth=get(case,"depth_reporting",Dict())
    haskey(depth,"event_stations") && (depth["event_stations"]=Dict(id=>row for (id,row) in depth["event_stations"] if id in selected))
    magnitude=get(case,"magnitudes",Dict())
    haskey(magnitude,"measurements") && (magnitude["measurements"]=[row for row in magnitude["measurements"] if row["event_id"] in selected])
    for key in ("external",)
        haskey(magnitude,key) && (magnitude[key]=Dict(id=>row for (id,row) in magnitude[key] if id in selected))
    end
    case
end

function _review_write_case(directory,config,events,picks)
    mkpath(directory)
    ids=[event.event_id for event in events]
    case=_review_subset(config,ids)
    _write_result_csv(joinpath(directory,"events.csv"),("event_id","x_km","y_km","z_km","origin_s"),
        [(event.event_id,event.source_km...,event.origin_s) for event in events])
    present=[hasproperty(row,:quality_code) for row in picks]
    if any(present) && !all(present)
        throw(ArgumentError("Review cannot save picks when only some arrivals have quality_code. A missing code is not quality 0. Give every arrival an explicit code, or omit the column for the whole file."))
    end
    header=any(present) ? ("event_id","station_id","phase","arrival_s","sigma_s","quality_code","observation_id") :
        ("event_id","station_id","phase","arrival_s","sigma_s","observation_id")
    rows=any(present) ?
        [(row.event_id,row.station_id,row.phase,row.arrival_s,row.sigma_s,row.quality_code,row.observation_id) for row in picks] :
        [(row.event_id,row.station_id,row.phase,row.arrival_s,row.sigma_s,row.observation_id) for row in picks]
    _write_result_csv(joinpath(directory,"picks.csv"),header,rows)
    case["files"]["events"]="events.csv"
    case["files"]["picks"]="picks.csv"
    case["output_directory"]="results"
    path=joinpath(directory,"case.toml")
    open(io->TOML.print(io,case),path,"w")
    path
end

"""
    review_catalog(config_path, steps; output_directory)

Review a native absolute catalog. Each step names the unchanged
event_id and action recalculate/accept/skip/delete. Reading edits require an
explicit observation_id; arrival_s, sigma_s, quality_code, phase_remark,
first_motion and remark may be changed. Quality codes must be explicit
integers from 0 to 9 and require empirical absolute weighting mode; editing
sigma_s applies directly to inverse-sigma weighting.

measurement_edits name a GLOBAL one-based measurement_index and edit coda or
amplitude quality, duration_s/amplitude/period_s; current event must agree.
Optional trial gives x_km/y_km/z_km/origin_s, or constraint gives a supported
event constraint. These explicit manual changes are the only data changes.

Recalculate writes and runs just this event, retaining every actual attempt.
Accept requires a current calculated attempt; pending edits cannot be accepted
as an old calculation. Skip retains input and previous attempts, while delete
removes the event from the saved working catalog without erasing source files.
Save all working edits, decisions, per-event attempts and accepted_catalog.csv.
A picks file with no quality_code column is written without that column; a
missing code is not stored as 0. Mixing arrivals with and without a code is
rejected. This is a plain scripted/terminal workflow, with no waveform processing.
"""
function review_catalog(config_path,steps;output_directory)
    base=dirname(abspath(config_path))
    config=_batch_paths!(TOML.parsefile(config_path),base)
    config["task"] in ("absolute","absolute_then_relative") ||
        throw(ArgumentError("Review requires a native absolute or absolute_then_relative case."))
    metadata=get(config,"processing_metadata",nothing)
    if metadata isa AbstractString
        table=TOML.parsefile(normpath(joinpath(base,metadata)))
        config["processing_metadata"]=get(table,"processing_metadata",table)
    end
    events=NamedTuple[read_events(config["files"]["events"])...]
    picks=NamedTuple[read_picks(config["files"]["picks"])...]
    event_ids=Set(event.event_id for event in events)
    last_results,dirty,decisions=Dict{String,Any}(),Dict(id=>true for id in event_ids),Dict(id=>:unreviewed for id in event_ids)
    deleted=Set{String}()
    history,attempts=NamedTuple[],NamedTuple[]
    output=abspath(output_directory)
    mkpath(output)
    empirical=get(get(config,"absolute_weighting",Dict()),"mode","sigma")=="empirical"
    function phase_quality!(index,requested)
        empirical || throw(ArgumentError("P/S quality edits require explicit absolute_weighting.mode=empirical; sigma-mode edits may change sigma_s."))
        row=picks[index]
        requested isa Integer && !(requested isa Bool) && 0<=requested<=9 ||
            throw(ArgumentError("Quality code must be an explicit integer from 0 to 9."))
        q=Int(requested)
        picks[index]=merge(row,(quality_code=q,))
        table=get!(config,"processing_metadata",Dict("events"=>Any[],"picks"=>Any[]))
        rows=get!(table,"picks",Any[])
        found=findfirst(item->get(item,"event_id",nothing)==row.event_id && get(item,"observation_id",nothing)==row.observation_id,rows)
        if found===nothing
            push!(rows,Dict("event_id"=>row.event_id,"station_id"=>row.station_id,"phase"=>row.phase,
                "observation_id"=>row.observation_id,"quality_code"=>q))
        else
            rows[found]["quality_code"]=q
        end
    end
    for (step_index,step) in enumerate(steps)
        all(Symbol(key) in (:event_id,:action,:edits,:measurement_edits,:trial,:constraint) for key in keys(step)) ||
            throw(ArgumentError("Review step accepts event_id, action, edits, measurement_edits, trial and constraint."))
        id=_report_value(step,:event_id)
        id in event_ids && !(id in deleted) || throw(ArgumentError("Review event $id is unavailable or deleted."))
        action=Symbol(_report_value(step,:action,"recalculate"))
        action in (:recalculate,:accept,:skip,:delete) || throw(ArgumentError("Review action must be recalculate, accept, skip or delete."))
        changes=NamedTuple[]
        for edit in _report_value(step,:edits,[])
            observation_id=_report_value(edit,:observation_id)
            found=findall(row->row.event_id==id && row.observation_id==observation_id,picks)
            length(found)==1 || throw(ArgumentError("Review needs one explicit event/observation_id: $id / $observation_id."))
            index=only(found)
            before=picks[index]
            for (name,value) in pairs(edit)
                key=Symbol(name)
                key==:observation_id && continue
                if key==:quality_code
                    phase_quality!(index,value)
                elseif key in (:arrival_s,:sigma_s)
                    value isa Real && isfinite(value) && (key!=:sigma_s || value>0) || throw(ArgumentError("Review times must be finite seconds; sigma positive."))
                    picks[index]=merge(picks[index],(;key=>Float64(value)))
                elseif key in (:phase_remark,:first_motion,:remark)
                    value isa AbstractString || throw(ArgumentError("Review reading labels must be text."))
                    table=get!(config,"processing_metadata",Dict("events"=>Any[],"picks"=>Any[]))
                    rows=get!(table,"picks",Any[])
                    match=findfirst(row->get(row,"event_id",nothing)==id && get(row,"observation_id",nothing)==observation_id,rows)
                    if match===nothing
                        push!(rows,Dict("event_id"=>id,"station_id"=>before.station_id,"phase"=>before.phase,"observation_id"=>observation_id))
                        match=length(rows)
                    end
                    rows[match][String(key)]=value
                else
                    throw(ArgumentError("Unsupported review reading edit $key."))
                end
            end
            push!(changes,(;kind=:reading,observation_id,before,after=picks[index],requested=edit))
        end
        magnitude=get(config,"magnitudes",Dict())
        measurements=get(magnitude,"measurements",[])
        for edit in _report_value(step,:measurement_edits,[])
            index=_report_value(edit,:measurement_index)
            index isa Integer && 1<=index<=length(measurements) && measurements[index]["event_id"]==id ||
                throw(ArgumentError("Measurement review index must identify this event's original config row."))
            before=deepcopy(measurements[index])
            for (name,value) in pairs(edit)
                key=String(name)
                key=="measurement_index" && continue
                if key in ("coda_quality_code","amplitude_quality_code")
                    value isa Integer && !(value isa Bool) && 0<=value<=9 ||
                        throw(ArgumentError("Magnitude quality code must be an explicit integer from 0 to 9."))
                    measurements[index][key]=Int(value)
                elseif key in ("duration_s","amplitude","period_s")
                    value isa Real && isfinite(value) && value>=0 || throw(ArgumentError("Manual magnitude measurement must be finite nonnegative."))
                    measurements[index][key]=Float64(value)
                else
                    throw(ArgumentError("Unsupported review measurement edit $key."))
                end
            end
            push!(changes,(;kind=:measurement,index,before,after=deepcopy(measurements[index])))
        end
        trial=_report_value(step,:trial)
        if trial!==nothing
            index=findfirst(event->event.event_id==id,events)
            before=events[index]
            allowed=(:x_km,:y_km,:z_km,:origin_s)
            all(Symbol(key) in allowed && value isa Real && isfinite(value) for (key,value) in pairs(trial)) || throw(ArgumentError("Review trial fields are finite x/y/z km and origin_s."))
            xyz=ntuple(i->_report_value(trial,allowed[i],before.source_km[i]),3)
            events[index]=merge(before,(source_km=xyz,origin_s=_report_value(trial,:origin_s,before.origin_s)))
            push!(changes,(;kind=:trial,before,after=events[index]))
        end
        constraint=_report_value(step,:constraint)
        if constraint!==nothing
            constraint in ("free","fixed_depth","fixed_position","fixed_all") || throw(ArgumentError("Unknown review constraint."))
            get!(config,"absolute_constraints",Dict())[id]=Dict("mode"=>constraint,"source"=>"explicit manual review step $step_index")
            push!(changes,(;kind=:constraint,constraint))
        end
        !isempty(changes) && (dirty[id]=true;decisions[id]=:unreviewed)
        if action==:recalculate
            selected_events=[event for event in events if event.event_id==id]
            selected_picks=[row for row in picks if row.event_id==id]
            try_index=length(attempts)+1
            path=_review_write_case(joinpath(output,"try_"*lpad(string(try_index),4,'0')),config,selected_events,selected_picks)
            result=run_case(path)
            last_results[id]=result
            dirty[id]=false
            event=only(result.absolute.events)
            push!(attempts,(;try_index,event_id=id,step_index,config_path=path,output_directory=result.output_directory,
                status=event.status,source_km=event.source_km,origin_s=event.origin_s))
        elseif action==:accept
            haskey(last_results,id) && !dirty[id] || throw(ArgumentError("Recalculate current edits before accepting event $id."))
            decisions[id]=:accepted
        elseif action==:skip
            decisions[id]=:skipped
        else
            decisions[id]=:deleted
            push!(deleted,id)
        end
        push!(history,(;step_index,event_id=id,action,changes,decision=decisions[id],pending_changes=dirty[id]))
    end
    working=_review_write_case(joinpath(output,"working"),config,[event for event in events if !(event.event_id in deleted)],
        [row for row in picks if !(row.event_id in deleted)])
    accepted=[only(last_results[event.event_id].absolute.events) for event in events if decisions[event.event_id]==:accepted]
    _write_result_csv(joinpath(output,"accepted_catalog.csv"),("event_id","x_km","y_km","z_km","origin_s","status"),
        [(event.event_id,event.source_km...,event.origin_s,event.status) for event in accepted])
    report=(;source_config=abspath(config_path),working_config=working,decisions,history,attempts,
        accepted_results=[last_results[event.event_id].absolute for event in accepted],
        assumptions=(;source_files_overwritten=false,reading_identity=:explicit_event_observation_id,
            measurement_index=:original_global_native_configuration_row,deleted_events=:omitted_from_working_catalog_source_retained,
            skipped_events=:working_input_retained_not_accepted,quality_definition=:explicit_zero_to_nine))
    open(io->TOML.print(io,_report_metadata(report)),joinpath(output,"review.toml"),"w")
    (;working_config=working,output_directory=output,decisions,attempts,accepted_catalog=joinpath(output,"accepted_catalog.csv"))
end

function run_catalog_review(path)
    config=TOML.parsefile(path)
    base=dirname(abspath(path))
    review_catalog(normpath(joinpath(base,config["case"])),config["steps"];
        output_directory=normpath(joinpath(base,config["output_directory"])))
end
