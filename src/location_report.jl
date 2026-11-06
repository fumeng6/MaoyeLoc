# Native scientific results in readable text, full-precision TOML and long-form CSV.
_report_value(row,key,default=nothing)=row isa AbstractDict ? get(row,String(key),get(row,key,default)) : get(row,key,default)
function _report_metadata(value)
    value===nothing || value===missing || value isa Symbol || value isa DateTime || value isa Date ||
        value isa AbstractString || value isa Number || value isa Bool || value isa NamedTuple ||
        value isa AbstractDict || value isa AbstractArray || value isa Tuple ||
        return Dict("type"=>string(typeof(value)),"fields"=>Dict(string(key)=>_report_metadata(getfield(value,key)) for key in fieldnames(typeof(value))))
    (value===nothing || value===missing) && return "not_available"
    (value isa Symbol || value isa DateTime || value isa Date) && return string(value)
    value isa NamedTuple || value isa AbstractDict ?
        Dict(string(key)=>_report_metadata(item) for (key,item) in pairs(value)) :
        value isa AbstractMatrix ? Dict("nrow"=>size(value,1),"ncol"=>size(value,2),"rows"=>[_report_metadata(Tuple(row)) for row in eachrow(value)]) :
        value isa AbstractArray || value isa Tuple ? [_report_metadata(item) for item in value] : value
end
function _report_lines(io,label,value)
    println(io,label)
    serialized=_report_metadata(value)
    if serialized isa AbstractDict
        TOML.print(io,serialized;sorted=true)
    elseif serialized isa AbstractVector && all(item->item isa AbstractDict,serialized)
        for (index,item) in enumerate(serialized)
            println(io,"record ",index)
            TOML.print(io,item;sorted=true)
            println(io)
        end
    else
        println(io,repr(serialized))
    end
    println(io)
end
function _report_companion(path,run;append=false)
    previous=append && isfile(path) ? TOML.parsefile(path) : Dict()
    runs=get(previous,"runs",Any[])
    push!(runs,_report_metadata(run))
    open(io->TOML.print(io,Dict("runs"=>runs)),path,"w")
    length(runs)
end
function _report_forward_readings(event,model,stations;source=event.source_km,origin=event.origin_s,iteration=nothing)
    rows=NamedTuple[]
    for input in event.residuals
        ray=try
            travel_time(model,source,stations,input.station_id,input.phase)
        catch error
            error isa DomainError || rethrow()
            nothing
        end
        # Coincidence removes a unique source gradient, not the forward time
        # of a fixed event (including an affine phase offset/known delay).
        travel=if ray===nothing
            try
                _forward_time(model,source,stations,input.station_id,input.phase)
            catch error
                error isa DomainError || rethrow()
                nothing
            end
        else
            ray.time_s
        end
        prediction=travel===nothing ? nothing : origin+travel
        get(event,:row_weighting,:empirical)==:empirical ||
            throw(ArgumentError("Location reports describe the empirical row-factor location."))
        saved=get(event,:weight_trace,NamedTuple[])
        history=iteration===nothing ? NamedTuple[] : [r for r in saved if r.iteration==iteration &&
            get(r,:observation_id,nothing)==get(input,:observation_id,nothing) && r.station_id==input.station_id && r.phase==input.phase &&
            get(r,:selected_for_equations,true)]
        weight=iteration===nothing ? get(input,:equation_weight,nothing) :
            isempty(history) ? nothing : last(history).equation_weight
        derivative=ray===nothing ? nothing : get(ray,:dsource_s_km,nothing)
        gradient=derivative===nothing ? nothing : (;x_km=derivative[1],y_km=derivative[2],z_km=derivative[3],origin_s=1.)
        weighted_gradient=gradient===nothing || weight===nothing ? nothing : (;(key=>value*weight for (key,value) in pairs(gradient))...)
        push!(rows,merge(input,(;event_id=event.event_id,iteration,source_km=source,origin_s=origin,
            travel_time_s=travel,ray_time_s=ray===nothing ? nothing : get(ray,:ray_time_s,ray.time_s),
            phase_offset_s=ray===nothing ? nothing : get(ray,:phase_offset_s,0.),
            station_delay_s=ray===nothing ? nothing : get(ray,:station_delay_s,0.),
            arrival_pred_s=prediction,residual_s=prediction===nothing ? nothing : input.arrival_s-prediction,
            equation_weight=weight,weight_source=iteration===nothing ? :saved_final : :saved_iteration_selection,
            prediction_gradient=gradient,weighted_prediction_gradient=weighted_gradient,
            ray)))
    end
    rows
end

"""
    write_magnitude_details(path,result; config=Dict(),source_metadata=nothing,append=false)

Long-form CSV: each actual station magnitude's additive term
occupies a row, retaining event/input-row/station identity, formula type,
value, status and reason. Their sum reconstructs the station magnitude.
The TOML companion stores full measurements, selected historical calibration,
station/component corrections, distance and response inputs, component/event
selection, aggregation, preferred and external magnitudes. Missing values
remain unavailable rather than fabricated zero corrections.
"""
function write_magnitude_details(path,result;config=Dict(),source_metadata=nothing,append=false)
    preparation=get(result,:magnitude_preparation,nothing)
    rows=preparation===nothing ? NamedTuple[] : preparation.rows
    companion=String(path)*".toml"
    previous=append && isfile(companion) ? TOML.parsefile(companion) : Dict()
    run_index=length(get(previous,"runs",[]))+1
    entries=Tuple[]
    for row in rows
        magnitudes=Pair{Any,Any}[]
        if row.coda!==nothing
            append!(magnitudes,[kind=>getproperty(row.coda,kind) for kind in (:md1,:md2,:mt)])
        end
        row.amplitude===nothing || push!(magnitudes,row.amplitude.magnitude_type=>row.amplitude)
        for (kind,magnitude) in magnitudes
            terms=get(magnitude,:terms,nothing)
            if terms===nothing
                push!(entries,(run_index,row.event_id,row.input_row,row.station_id,kind,magnitude.magnitude,magnitude.status,magnitude.reason,"",nothing))
            else
                for (name,value) in pairs(terms)
                    push!(entries,(run_index,row.event_id,row.input_row,row.station_id,kind,magnitude.magnitude,magnitude.status,magnitude.reason,name,value))
                end
            end
        end
    end
    header=("run_index","event_id","input_row","station_id","magnitude_type","magnitude","status","reason","term","term_mag")
    if append && isfile(path) && filesize(path)>0
        first(readlines(path))==join(header,',') || throw(ArgumentError("Existing magnitude-detail CSV has a different header."))
        open(path,"a") do io
            for entry in entries
                fields=[value===nothing ? "" : string(value) for value in entry]
                println(io,join((occursin(r"[,\"\r\n]",field) ? "\""*replace(field,"\""=>"\"\"")*"\"" : field for field in fields),','))
            end
        end
    else
        _write_result_csv(path,header,entries)
    end
    source=source_metadata isa AbstractString ? TOML.parsefile(source_metadata) : source_metadata
    _report_companion(companion,(;run_index,preparation,config,source_metadata=source,term_definition=:additive_magnitude_units,
        identity=:original_input_row_and_event_station,negative_magnitudes=:retained);append)
    (;path=abspath(path),metadata_path=abspath(companion),run_index,station_measurements=length(rows),term_rows=length(entries))
end

"""
    write_location_report(path,result,model,stations; config=Dict(),frame=nothing,
        time_reference_utc,source_metadata=nothing,detail=1,append=false)

Readable report plus a full-precision TOML companion and magnitude term CSV. The real model fields, station geometry/metadata, complete supplied
configuration, input metadata, all final rows including zero
weight, constraints, uncertainty and magnitude calibration remain available.
The writer evaluates the saved locations; it never fits or reweights them.

detail 0 prints final event results; 1 adds every final reading and magnitude;
2 adds saved iteration locations/adjustments and weight history; 3 adds the
saved singular spectrum, covariance, error axes, INFO and calibration history.
4 adds station predictions at each
saved iteration with its saved weights. Missing iteration weights stay blank.
No intermediate iterations are invented for an analytic fixed-position solve.
All detail levels retain the complete same native data in the companion.
append=true adds a separate run to text, TOML and magnitude details.
"""
function write_location_report(path,result,model,stations;config=Dict(),frame=nothing,
        time_reference_utc,source_metadata=nothing,detail=1,append=false,magnitude_path=nothing)
    detail isa Integer && !(detail isa Bool) && 0<=detail<=4 || throw(ArgumentError("Report detail must be an integer from 0 to 4."))
    endswith(time_reference_utc,"Z") || throw(ArgumentError("Report reference time must explicitly end in Z."))
    source=source_metadata isa AbstractString ? TOML.parsefile(source_metadata) : source_metadata
    companion=String(path)*".toml"
    previous=append && isfile(companion) ? TOML.parsefile(companion) : Dict()
    run_index=length(get(previous,"runs",[]))+1
    final_rows=NamedTuple[]
    iteration_rows=NamedTuple[]
    for event in result.events
        append!(final_rows,_report_forward_readings(event,model,stations))
        if detail>=4
            for state in get(event,:trace,NamedTuple[])
                append!(iteration_rows,_report_forward_readings(event,model,stations;source=state.source_km,
                    origin=state.origin_s,iteration=state.iteration))
            end
        end
    end
    preparation=get(result,:magnitude_preparation,nothing)
    magnitude=write_magnitude_details(magnitude_path===nothing ? String(path)*".magnitudes.csv" : magnitude_path,
        result;config,source_metadata=source,append)
    run=(;run_index,detail,time_reference_utc,frame,config,model,stations,source_metadata=source,
        result,final_readings=final_rows,iteration_readings=iteration_rows,
        assumptions=(;units="km, seconds, km/s and additive magnitude units",residual=:observed_minus_predicted,
            weights=:actual_saved_equation_rows_or_inverse_sigma,fit_changed=false,
            missing_values=:not_available,uncertainty=:same_estimator_conditions_as_saved_result,
            jacobian_columns=:full_prediction_before_fixed_column_selection,
            iteration_weights=:saved_selection_for_update_not_reweighted_at_displayed_point,
            iteration_locations=:actual_saved_trace_no_reconstructed_solver_steps))
    open(path,append ? "a" : "w") do io
        println(io,"MAOYELOC DETAILED LOCATION REPORT / RUN ",run_index)
        println(io,"UTC reference: ",time_reference_utc,"; detail: ",detail)
        println(io,"Units: km, seconds, km/s, magnitude. Residual = observed - predicted.")
        println(io,"Equation weights and statistical sigma are separate; unavailable values are not zero.\n")
        _report_lines(io,"ACTUAL MODEL",model)
        _report_lines(io,"STATIONS (geometry; original station parameters in source metadata)",stations)
        _report_lines(io,"SUPPLIED COMPUTATION PARAMETERS AND FILES",config)
        for event in result.events
            fields=(:event_id,:source_km,:origin_s,:status,:stop_reason,:constraint,:constraint_source,
                :iterations,:input_count,:used_count,:rms_residual_s,:weighted_sum_squares,
                :catalog_eligible,:catalog_exclusion_reason)
            _report_lines(io,"FINAL EVENT",(; (key=>get(event,key,nothing) for key in fields)...))
            frame===nothing || _report_lines(io,"GEOGRAPHIC LOCATION",local_to_geographic(frame,event.source_km))
            if detail>=2
                _report_lines(io,"SAVED ITERATIONS / ADJUSTMENTS",get(event,:trace,NamedTuple[]))
                _report_lines(io,"SAVED ITERATION WEIGHTS",get(event,:weight_trace,NamedTuple[]))
            end
            if detail>=3
                for key in (:spectral,:uncertainty,:formal_uncertainty,:error_ellipsoids,:diagnostics,:geometry,:depth_report,:geographic_region)
                    _report_lines(io,string(key),get(event,key,nothing))
                end
            end
        end
        if detail>=1
            _report_lines(io,"ALL FINAL READINGS (INCLUDING ZERO WEIGHT)",final_rows)
            _report_lines(io,"SUPPLIED INPUT METADATA / UNUSED OR EXCLUDED READINGS",source)
            _report_lines(io,"STATION MAGNITUDE TERMS AND CALIBRATION HISTORY",preparation)
            _report_lines(io,"PROCESSING / FIRST MOTIONS / SOURCE LABELS",get(result,:processing_metadata,nothing))
        end
        detail>=4 && _report_lines(io,"STATION PREDICTIONS AT EACH SAVED ITERATION",iteration_rows)
        println(io,"END RUN ",run_index,"\n")
    end
    _report_companion(companion,run;append)
    (;path=abspath(path),metadata_path=abspath(companion),run_index,detail,magnitude_details=magnitude,
        events=length(result.events),final_readings=length(final_rows),iteration_readings=length(iteration_rows))
end
