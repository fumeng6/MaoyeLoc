# Group once without sorting or merging repeated readings. Keep concrete native
# row types in each group; heterogeneous user records retain their extra fields.
function _group_by_event(rows; event_ids=nothing)
    row_type = eltype(rows)
    grouped = Dict{String,Vector{row_type}}()
    if event_ids !== nothing
        for id in event_ids
            grouped[id] = row_type[]
        end
    end
    counts = Dict{String,Int}()
    for row in rows
        if event_ids === nothing
            get!(() -> row_type[], grouped, row.event_id)
        else
            haskey(grouped, row.event_id) || throw(ArgumentError("Pick refers to unknown event."))
        end
        counts[row.event_id] = get(counts, row.event_id, 0) + 1
    end
    for (id, count) in counts
        sizehint!(grouped[id], count)
    end
    for row in rows
        push!(grouped[row.event_id], row)
    end
    grouped
end

"""
    initialize_catalog(model, stations, picks; depth_km)

Starting guesses from the horizontal centroid of each event's contributing
stations at an explicitly chosen depth. Initial origin time is the weighted
least-squares optimum with that position fixed. These are starting guesses,
not locations; no truth, old catalog, or preferred-looking branch is used.
"""
function initialize_catalog(model, stations, picks; depth_km::Real, picks_by_event=nothing)
    isfinite(depth_km) && depth_km >= model_top(model) || throw(ArgumentError("Initial depth outside model."))
    observations = collect(picks)
    ids = unique(pick.event_id for pick in observations)
    by_event = picks_by_event === nothing ? _group_by_event(observations; event_ids=ids) : picks_by_event
    events = NamedTuple{(:event_id,:source_km,:origin_s),Tuple{String,NTuple{3,Float64},Float64}}[]
    for id in ids
        selected = by_event[id]
        station_ids = unique(pick.station_id for pick in selected)
        all(haskey(stations, sid) for sid in station_ids) || throw(ArgumentError("Unknown station in initial catalog."))
        x = sum(stations[sid][1] for sid in station_ids)/length(station_ids)
        y = sum(stations[sid][2] for sid in station_ids)/length(station_ids)
        source = (x, y, Float64(depth_km))
        origin_sum, weight_sum = 0.0, 0.0
        for pick in selected
            isfinite(pick.arrival_s) && isfinite(pick.sigma_s) && pick.sigma_s > 0 ||
                throw(ArgumentError("Picks require finite arrivals and positive sigma."))
            weight = 1/pick.sigma_s^2
            time = _forward_time(model,source,stations,pick.station_id,pick.phase)
            origin_sum += weight*(pick.arrival_s-time)
            weight_sum += weight
        end
        push!(events, (event_id=String(id), source_km=source, origin_s=origin_sum/weight_sum))
    end
    return events
end

"Locate each catalog event without changing IDs or dropping events with no picks."
function locate_catalog(model, stations, events, picks;tasks::Integer=1,progress=nothing,picks_by_event=nothing,kwargs...)
    catalog, observations = collect(events), collect(picks)
    ids = [event.event_id for event in catalog]
    length(unique(ids)) == length(ids) || throw(ArgumentError("Duplicate event ID."))
    by_event = picks_by_event === nothing ? _group_by_event(observations; event_ids=ids) : picks_by_event
    if picks_by_event !== nothing
        valid_ids = Set(ids)
        all(pick.event_id in valid_ids for pick in observations) || throw(ArgumentError("Pick refers to unknown event."))
    end
    results = _ordered_tasks(catalog;tasks,progress) do event
        selected = by_event[event.event_id]
        event_options = merge((;kwargs...),
            (constraint=get(event,:constraint,get(kwargs,:constraint,:free)),
             constraint_source=get(event,:constraint_source,get(kwargs,:constraint_source,"unconstrained"))))
        location = try
            locate(model, stations, selected; initial_source_km=event.source_km,
                   initial_origin_s=event.origin_s, event_options...)
        catch exception
            exception isa DomainError && any(Tuple(event.source_km)==Tuple(stations[pick.station_id]) for pick in selected) || rethrow()
            # This is a failure of this event's initial geometry, not missing
            # data and not a reason to stop locating the remaining catalog.
            empty_result = locate(model,stations,selected[1:0];initial_source_km=event.source_km,
                                  initial_origin_s=event.origin_s,event_options...)
            rows = NamedTuple[]
            for pick in selected
                prediction = event.origin_s+_forward_time(model,event.source_km,stations,pick.station_id,pick.phase)
                residual = pick.arrival_s-prediction
                push!(rows,merge(pick,(arrival_pred_s=prediction,residual_s=residual,
                    normalized_residual=residual/pick.sigma_s,used=false,unused_reason=:coincident_source_station)))
            end
            merge(empty_result,(status=:coincident_source_station,residuals=rows,
                input_count=length(rows),used_count=0,unused_count=length(rows),
                forward_evaluated_count=length(rows),unused_reason=:coincident_source_station,
                weighted_sum_squares=sum(row.normalized_residual^2 for row in rows),
                rms_residual_s=sqrt(sum(row.residual_s^2 for row in rows)/length(rows)),
                jacobian_rank=nothing,free_parameter_rank=nothing,condition_number=nothing,
                unconstrained_condition_number=nothing))
        end
        location = merge(location, (event_id=String(event.event_id),))
        location
    end
    residuals=NamedTuple[]
    for location in results
        append!(residuals,location.residuals)
    end
    return (; events=results, residuals)
end

# Values held fixed are the event's supplied calculation coordinates and
# origin, not its optional GEO report. Constraints here apply to ABS only.
function _absolute_catalog_constraints(events,definitions)
    definitions isa AbstractDict || throw(ArgumentError("absolute_constraints must be event tables."))
    ids = Set(event.event_id for event in events)
    all(id -> id in ids,keys(definitions)) || throw(ArgumentError("Absolute constraint refers to an unknown event."))
    [if haskey(definitions,event.event_id)
        entry = definitions[event.event_id]
        entry isa AbstractDict && haskey(entry,"mode") && all(key -> key in ("mode","source"),keys(entry)) ||
            throw(ArgumentError("Each absolute constraint needs mode and optional source; fixed values come from the event coordinates/origin."))
        mode = Symbol(entry["mode"])
        mode in (:free,:fixed_depth,:fixed_position,:fixed_all) || throw(ArgumentError("Unknown absolute constraint mode."))
        source = get(entry,"source","absolute_constraints")
        source isa AbstractString && !isempty(strip(source)) || throw(ArgumentError("Constraint source must be nonempty text."))
        merge(event,(constraint=mode,constraint_source=String(source)))
    else
        event
    end for event in events]
end

function _write_csv_field(io,value)
    value===nothing && return
    field=string(value)
    if any(character->character in (',','"','\r','\n'),field)
        print(io,'"',replace(field,"\""=>"\"\""),'"')
    else
        print(io,field)
    end
end
_write_csv_field(io,value::Union{Integer,AbstractFloat})=print(io,value)
_write_csv_field(io,value,scratch)=_write_csv_field(io,value)
function _write_csv_field(io,value::Symbol,scratch)
    field=get!(scratch.symbols,value) do
        String(value)
    end
    _write_csv_field(io,field)
end
function _write_csv_field(io,value::Union{Float16,Float32,Float64},scratch)
    # Use Julia's usual shortest round-trip representation, reusing its workspace.
    last=Base.Ryu.writeshortest(scratch.bytes,1,value)-1
    write(io,view(scratch.bytes,1:last))
end
function _write_csv_row(io,row,scratch)
    for (index,value) in enumerate(row)
        index>1 && write(io,',')
        _write_csv_field(io,value,scratch)
    end
    write(io,'\n')
end
@generated function _write_csv_row(io,row::T,scratch) where {T<:Tuple}
    # Long heterogeneous result rows need constant field indices too.
    statements=Expr[]
    for index in 1:fieldcount(T)
        index>1 && push!(statements,:(write(io,',')))
        push!(statements,:(_write_csv_field(io,getfield(row,$index),scratch)))
    end
    push!(statements,:(write(io,'\n')))
    Expr(:block,statements...)
end
_write_csv_row(io,row,scratch,::Nothing)=_write_csv_row(io,row,scratch)
@generated function _write_csv_row(io,row::NamedTuple{names},scratch,::Val{fields}) where {names,fields}
    statements=Expr[]
    for (index,field) in enumerate(fields)
        index>1 && push!(statements,:(write(io,',')))
        value=field in names ? :(getfield(row,$(QuoteNode(field)))) : :(nothing)
        push!(statements,:(_write_csv_field(io,$value,scratch)))
    end
    push!(statements,:(write(io,'\n')))
    Expr(:block,statements...)
end
_result_csv_fields(row,::Val{fields}) where {fields} =
    ntuple(@inline(i->get(row,fields[i],nothing)),Val(length(fields)))
@generated function _result_csv_fields(row::NamedTuple{names},::Val{fields}) where {names,fields}
    values=[field in names ? :(getfield(row,$(QuoteNode(field)))) : :(nothing) for field in fields]
    Expr(:tuple,values...)
end
function _write_result_csv(path,header,rows;fields=nothing)
    open(path,"w") do io
        println(io,join(header,','))
        buffer=IOBuffer(sizehint=262144)
        scratch=(bytes=Vector{UInt8}(undef,Base.Ryu.neededdigits(Float64)),symbols=Dict{Symbol,String}())
        for row in rows
            _write_csv_row(buffer,row,scratch,fields)
            position(buffer)>=262144 && write(io,take!(buffer))
        end
        write(io,take!(buffer))
    end
    nothing
end

function _with_absolute_errors(result,model,stations,picks,settings; picks_by_event=nothing)
    formal=get(settings,"formal",nothing)
    formal===nothing || (formal isa AbstractDict && all(key in ("rderr_s","ercof","noise_model") for key in keys(formal))) ||
        throw(ArgumentError("uncertainty.formal accepts rderr_s, ercof and noise_model."))
    probability=get(settings,"ellipsoid_probability",.95)
    by_event = picks_by_event === nothing ? _group_by_event(picks) : picks_by_event
    events=NamedTuple[]
    for event in result.events
        selected=get(() -> eltype(picks)[], by_event, event.event_id)
        known=get(settings,"assume_independent_picks",false) ?
            absolute_uncertainty(model,stations,selected,event;assume_independent=true) : event.uncertainty
        working=formal===nothing ? nothing : absolute_formal_uncertainty(model,stations,selected,event;
            rderr_s=get(formal,"rderr_s",.15),ercof=get(formal,"ercof",1.),
            noise_model=Symbol(get(formal,"noise_model","weighted_rows")))
        ellipsoids=(sigma=error_ellipsoid(known;probability),formal=error_ellipsoid(working;probability))
        push!(events,merge(event,(uncertainty=known,formal_uncertainty=working,error_ellipsoids=ellipsoids)))
    end
    merge(result,(;events))
end

function _write_absolute_error_results(directory,events)
    working=[(event,estimate) for event in events for estimate in (get(event,:formal_uncertainty,nothing),) if estimate!==nothing]
    _write_result_csv(joinpath(directory,"formal_errors.csv"),
        ("event_id","status","reason","quantity","noise_model","rderr_s","ercof","weighted_rms_s",
         "variance_s2","timing_error_s","included_count","row_normalization","rank",
         "formal_se_x_km","formal_se_y_km","formal_se_z_km","formal_se_origin_s"),
        ((event.event_id,u.status,u.reason,u.quantity,u.noise_model,u.rderr_s,u.ercof,u.weighted_rms_s,
          u.variance_s2,u.timing_error_s,u.included_count,u.row_normalization,u.rank,Tuple(u.standard_deviations)...)
          for (event,u) in working))
    _write_result_csv(joinpath(directory,"formal_covariance.csv"),
        ("event_id","parameter_i","parameter_j","covariance"),
        ((event.event_id,u.parameter_names[i],u.parameter_names[j],u.covariance[i,j])
         for (event,u) in working if u.covariance!==nothing
         for i in eachindex(u.parameter_names) for j in eachindex(u.parameter_names)))
    ellipsoids=[(event.event_id,source,ellipse) for event in events
        for (source,ellipse) in pairs(get(event,:error_ellipsoids,NamedTuple()))]
    _write_result_csv(joinpath(directory,"error_ellipsoids.csv"),
        ("event_id","covariance_source","status","reason","probability","spatial_dimension","stochastic_rank",
         "confidence_scale","enclosed_probability","horizontal_radius_km","vertical_extent_km","legacy_erh_km","legacy_erz_km"),
        ((id,source,e.status,e.reason,e.probability,e.spatial_dimension,e.stochastic_rank,e.confidence_scale,
          get(e,:enclosed_probability,nothing),e.horizontal_radius_km,e.vertical_extent_km,e.legacy_erh_km,e.legacy_erz_km)
          for (id,source,e) in ellipsoids))
    _write_result_csv(joinpath(directory,"error_axes.csv"),
        ("event_id","covariance_source","axis","standard_error_km","semiaxis_km","east_component","north_component",
         "down_component","azimuth_deg","dip_deg","orientation_identifiable"),
        ((id,source,Tuple(axis)...) for (id,source,e) in ellipsoids for axis in e.axes))
end

function _with_absolute_diagnostics(result,model,stations,picks,settings,depth_settings; picks_by_event=nothing)
    get(settings,"absolute",false) || return result
    source=Symbol(get(settings,"quality_covariance_source","formal"))
    source in (:sigma,:formal) || throw(ArgumentError("Quality covariance source must be sigma or formal."))
    sites=get(settings,"station_site_ids",get(depth_settings,"station_site_ids",Dict()))
    by_event = picks_by_event === nothing ? _group_by_event(picks) : picks_by_event
    events=[merge(event,(
        diagnostics=absolute_diagnostics(model,stations,get(() -> eltype(picks)[], by_event, event.event_id),event),
        geometry=merge(absolute_geometry(model,stations,event;
            reported_depth_km=hasproperty(event,:depth_report) ? event.depth_report.reported_depth_km : event.source_km[3],
            ellipsoid=get(event.error_ellipsoids,source,nothing),station_site_ids=sites),
            (quality_covariance_source=source,)))) for event in result.events]
    merge(result,(;events))
end

function _write_absolute_diagnostics(directory,events)
    selected=[event for event in events if hasproperty(event,:diagnostics)]
    _write_result_csv(joinpath(directory,"information.csv"),
        ("event_id","status","reason","solution_status","retained_rank","numerical_rank","information_sum",
         "information_sum_minus_rank","conditioned_depth","row_normalization"),
        ((e.event_id,d.status,d.reason,d.solution_status,d.retained_rank,d.numerical_rank,d.information_sum,
          d.information_sum_minus_rank,d.conditioned_depth,d.row_normalization) for e in selected for d in (e.diagnostics,)))
    _write_result_csv(joinpath(directory,"reading_information.csv"),
        ("event_id","station_id","phase","equation_weight","included","used","info","info_milli","retained_projections","observation_id"),
        ((e.event_id,r.station_id,r.phase,r.equation_weight,r.included,r.used,r.info,r.info_milli,
          r.retained_projections===nothing ? nothing : join(r.retained_projections,'|'),get(r,:observation_id,nothing)) for e in selected for r in e.diagnostics.rows))
    _write_result_csv(joinpath(directory,"station_information.csv"),
        ("event_id","station_id","phase_count","included_count","p_info","s_info","info"),
        ((e.event_id,Tuple(r)...) for e in selected for r in e.diagnostics.station_contributions))
    _write_result_csv(joinpath(directory,"information_directions.csv"),
        ("event_id","direction","singular_value","retained","x_component","y_component","z_component","origin_component"),
        ((e.event_id,r.direction,r.singular_value,r.retained,Tuple(r.vector4)...) for e in selected for r in e.diagnostics.directions))
    _write_result_csv(joinpath(directory,"event_geometry.csv"),
        ("event_id","gap_deg","nearest_station_distance_km","weighted_reading_count","weighted_p_count","weighted_s_count",
         "valid_reading_count","effective_station_count","azimuth_station_count","weighted_rms_s","legacy_erh_km","legacy_erz_km",
         "reported_depth_km","row_weight_definition","row_normalization","geometry_reason","quality_code","error_quality_code",
         "geometry_quality_code","quality_reason","quality_covariance_source"),
        ((e.event_id,(getproperty(e.geometry,k) for k in (:gap_deg,:nearest_station_distance_km,:weighted_reading_count,:weighted_p_count,
          :weighted_s_count,:valid_reading_count,:effective_station_count,:azimuth_station_count,:weighted_rms_s,:legacy_erh_km,:legacy_erz_km,
          :reported_depth_km,:row_weight_definition,:row_normalization,:geometry_reason,:quality_code,:error_quality_code,
          :geometry_quality_code,:quality_reason,:quality_covariance_source))...) for e in selected))
    _write_result_csv(joinpath(directory,"reading_geometry.csv"),
        ("event_id","station_id","phase","site_id","distance_km","azimuth_deg","source_angle_deg","angle_definition","angle_reason",
         "equation_weight","used_for_fit","weighted_reading","geometry_eligible","valid_assigned_reading","quality_code","quality_source","observation_id"),
        ((e.event_id,(getproperty(r,k) for k in (:station_id,:phase,:site_id,:distance_km,:azimuth_deg,:source_angle_deg,:angle_definition,:angle_reason,
          :equation_weight,:used_for_fit,:weighted_reading,:geometry_eligible,:valid_assigned_reading,:quality_code,:quality_source))...,get(r,:observation_id,nothing))
          for e in selected for r in e.geometry.rows))
    _write_result_csv(joinpath(directory,"ray_components.csv"),
        ("event_id","station_id","phase","model_id","weight","source_angle_deg","angle_definition","angle_reason","observation_id"),
        ((e.event_id,r.station_id,r.phase,c.model_id,c.weight,c.source_angle_deg,c.angle_definition,c.angle_reason,get(r,:observation_id,nothing))
          for e in selected for r in e.geometry.rows for c in r.source_angle_components))
end

function _with_geographic_regions(result,frame,settings,base)
    isempty(settings) && return result
    get(settings,"net",0)==0 && return result
    frame===nothing && throw(ArgumentError("Geographic NET queries need an explicit geographic frame."))
    net=get(settings,"net",0)
    directory=get(settings,"reference_directory",_default_data_directory())
    directory=normpath(joinpath(base,directory))
    isempty(result.events) && return result
    prepared=_prepare_geographic_regions(net;reference_directory=directory)
    events=map(result.events) do event
        geo=local_to_geographic(frame,event.source_km)
        region=geographic_region(geo.latitude_deg,geo.longitude_deg,event.source_km[3];net,reference_directory=directory,_prepared=prepared)
        merge(event,(geographic_region=region,))
    end
    merge(result,(;events))
end

function _write_geographic_regions(directory,events)
    selected=[merge((event_id=event.event_id,),event.geographic_region) for event in events if hasproperty(event,:geographic_region)]
    isempty(selected) && return
    open(joinpath(directory,"regions.toml"),"w") do io
        TOML.print(io,Dict("events"=>_result_metadata(selected)))
    end
end

function _with_processing_metadata(result,stations,settings)
    settings===nothing && return result
    weighted=NamedTuple[]
    for event in result.events
        rows=event.residuals
        saved_weights=get(event,:row_weighting,:empirical)==:empirical
        normalization=isempty(rows) ? 1. : length(rows)/sum(inv(r.sigma_s) for r in rows)
        for row in rows
            weight=saved_weights ? row.equation_weight : normalization/row.sigma_s
            push!(weighted,merge(row,(equation_weight=weight,)))
        end
    end
    metadata=processing_metadata(result.events,weighted;metadata=settings,known_station_ids=keys(stations))
    merge(result,(processing_metadata=metadata,))
end

# Time actual phases without repeating work or consuming random numbers.
function _case_phase(f,stage,timings,progress)
    progress===nothing || progress((;kind=:phase,stage,status=:started))
    started=time_ns()
    status=:failed
    try
        value=f()
        status=:completed
        value
    finally
        elapsed=(time_ns()-started)/1e9
        timings[String(stage)]=get(timings,String(stage),0.)+elapsed
        progress===nothing || progress((;kind=:phase,stage,status,elapsed_seconds=elapsed))
    end
end

# Execute already prepared inputs; public path/dictionary entry points live in case.jl.
function _execute_case(prepared;progress=nothing,write_results=true,calculation=effective_config(prepared),
        timings=Dict{String,Float64}(),completed_stages=String[],stage_summaries=Dict{String,Any}())
    (;config_path,config,base,files,task,reference,frame,model,joint,absolute_failure_policy,
      stations,events,picks,output,
      bounds,options,absolute_controls,absolute_options,relative_controls,relative_initialization,
      relative_selection,initialization_coordinates,tasks,replay_tasks,output_detail,output_diagnostics)=prepared
    input_path(key)=normpath(joinpath(base,files[key]))
    picks_by_event=picks===nothing ? nothing : _group_by_event(picks;event_ids=(e.event_id for e in events))
    # Remove only our named products for the requested stages. A failed rerun
    # must not leave an old relative catalog beside new absolute locations.
    stages = task == "absolute_then_relative" ? ("absolute", "relative") : (task,)
    input_files = _case_input_paths(config,base,config_path)
    if write_results
        for stage in stages, name in _CASE_STAGE_FILES
            path = joinpath(output, stage, name)
            _case_path_key(path) in input_files && throw(ArgumentError("An output file is also a configured input: $path"))
            isfile(path) && rm(path)
        end
        if task != "absolute"
            for name in ("generated_ct.csv","pairing_events.csv","pairing_pairs.csv","final_pair_geometry.csv")
                path = joinpath(output,name)
                _case_path_key(path) in input_files && throw(ArgumentError("An output file is also a configured input: $path"))
                isfile(path) && rm(path)
            end
        end
        mkpath(output)
        open(joinpath(output,"calculation.toml"),"w") do io
            TOML.print(io,_result_metadata(calculation))
        end
    end
    absolute, relative = nothing, nothing
    if task in ("absolute", "absolute_then_relative")
        absolute_progress=progress===nothing ? nothing : (done,total,index,row)->progress(
            (;stage=:absolute,processed=done,total,event_id=row.event_id,status=row.status))
        absolute = _case_phase(:absolute_solve,timings,progress) do
            locate_catalog(model,stations,events,picks;tasks,picks_by_event,progress=absolute_progress,absolute_options...)
        end
        absolute = _case_phase(:absolute_uncertainty,timings,progress) do
            _with_absolute_errors(absolute,model,stations,picks,get(config,"uncertainty",Dict());picks_by_event)
        end
        absolute = _case_phase(:absolute_analysis,timings,progress) do
            report_stations = Dict(event.event_id=>[pick.station_id for pick in picks_by_event[event.event_id]]
                                   for event in absolute.events)

            absolute = _with_depth_reports(absolute,model,stations,report_stations,get(config,"depth_reporting",Dict()))
            processing_settings=get(config,"processing_metadata",nothing)
            processing_settings isa AbstractString && (processing_settings=normpath(joinpath(base,processing_settings)))
            absolute = _with_processing_metadata(absolute,stations,processing_settings)
            if output_diagnostics
                diagnostic_settings=merge(Dict("absolute"=>true),get(config,"diagnostics",Dict()))
                absolute = _with_absolute_diagnostics(absolute,model,stations,picks,diagnostic_settings,get(config,"depth_reporting",Dict());picks_by_event)
            end
            absolute = _with_geographic_regions(absolute,frame,get(config,"geographic_regions",Dict()),base)
            if haskey(config,"magnitudes")
                absolute=merge(absolute,(magnitude_preparation=magnitude_preparation(absolute.events,stations,config["magnitudes"];time_reference_utc=reference,model,frame),))
            end
            absolute
        end
        if write_results
            _case_phase(:absolute_export,timings,progress) do
                write_absolute_results(joinpath(output, "absolute"), absolute;detail=output_detail,diagnostics=output_diagnostics)
                frame !== nothing && write_geographic_results(joinpath(output,"absolute"),absolute,frame;time_reference_utc=reference)

            end
        end
        stage_summaries["absolute"]=_case_stage_summary(absolute)
        push!(completed_stages,"absolute")
    end
    if task in ("relative", "absolute_then_relative")
        failed=NamedTuple[]
        if absolute !== nothing
            failed=[event for event in absolute.events if !get(event,:catalog_eligible,true) ||
                !(event.status in (:converged,:converged_at_depth_bound,
                :converged_truncated,:stopped_by_legacy_rule,:fully_fixed))]
            isempty(failed) || absolute_failure_policy=="exclude" ||
                error("Some absolute locations failed; inspect outputs before starting relative relocation, or explicitly set relative.absolute_failure_policy=\"exclude\".")
            failed_ids=Set(event.event_id for event in failed)
            any(id in failed_ids for id in get(config,"reference_event_ids",String[])) &&
                error("An explicitly requested reference failed absolute location; its constraint is not silently replaced.")
            events = [(event_id=row.event_id, source_km=row.source_km, origin_s=row.origin_s)
                      for row in absolute.events if !(row.event_id in failed_ids)]
        end
        differences = if prepared.differences!==nothing
            prepared.differences
        else
            haskey(config,"pairing") || throw(ArgumentError("Supply files.differentials or explicit [pairing] criteria and files.picks."))
            criteria = config["pairing"]
            original_ids=Set(row.event_id for row in (absolute===nothing ? events : absolute.events))
            all(pick.event_id in original_ids for pick in picks) || throw(ArgumentError("Pick refers to unknown event."))
            retained_ids=Set(row.event_id for row in events)
            pairing_picks=[pick for pick in picks if pick.event_id in retained_ids]
            paired = pair_catalog(events,stations,pairing_picks;
                max_separation_km=criteria["max_separation_km"],
                max_separation_ratio=criteria["max_separation_ratio"],
                min_common_stations=criteria["min_common_stations"])
            for event in failed
                push!(paired.events,(event_id=event.event_id,candidate_pair_count=0,retained_pair_count=0,
                    status=:absolute_location_failed))
            end
            if write_results && output_diagnostics
                _write_result_csv(joinpath(output,"generated_ct.csv"),
                    ("event_i","event_j","station_id","phase","dt_s","sigma_s","source"),
                    ((row.event_i,row.event_j,row.station_id,row.phase,row.dt_s,row.sigma_s,row.source) for row in paired.observations))
                _write_result_csv(joinpath(output,"pairing_events.csv"),
                    ("event_id","candidate_pair_count","retained_pair_count","status"),
                    ((row.event_id,row.candidate_pair_count,row.retained_pair_count,row.status) for row in paired.events))
                _write_result_csv(joinpath(output,"pairing_pairs.csv"),
                    ("event_i","event_j","separation_km","common_phase_count","retained_phase_count",
                     "common_station_count","max_retained_separation_ratio","retained","reason"),
                    ((row.event_i,row.event_j,row.separation_km,row.common_phase_count,row.retained_phase_count,
                      row.common_station_count,row.max_retained_separation_ratio,row.retained,row.reason) for row in paired.pairs))
            end
            paired.observations
        end
        excluded=NamedTuple[]
        relative_picks=joint ? picks : NamedTuple[]
        differential_count,absolute_count=length(differences),length(relative_picks)
        differential_rows=eachindex(differences)
        absolute_rows=eachindex(relative_picks)
        if !isempty(failed)
            all_ids=Set(row.event_id for row in absolute.events)
            all(row.event_i in all_ids && row.event_j in all_ids for row in differences) ||
                throw(ArgumentError("Difference refers to an unknown catalog event."))
            failed_ids=Set(event.event_id for event in failed)
            excluded=[merge(row,(input_row=i,)) for (i,row) in enumerate(differences)
                if row.event_i in failed_ids || row.event_j in failed_ids]
            differential_rows=findall(row->!(row.event_i in failed_ids || row.event_j in failed_ids),differences)
            absolute_rows=findall(row->!(row.event_id in failed_ids),relative_picks)
            differences=differences[differential_rows]
            relative_picks=relative_picks[absolute_rows]
            # User arrays refer to the original DD rows, then original ABS rows.
            # Use these same retained indices for the solve and uncertainty replays.
            if get(relative_controls,:row_weighting,:sigma)==:explicit
                weights=relative_controls.row_weights
                length(weights)==differential_count+absolute_count ||
                    throw(ArgumentError("Explicit DD row_weights must match all original differential then absolute-pick rows before ABS exclusions."))
                relative_controls=merge(relative_controls,
                    (row_weights=weights[[differential_rows;differential_count.+absolute_rows]],))
            elseif get(relative_controls,:row_weighting,:sigma)==:empirical
                weighting=relative_controls.empirical_weighting
                weights=get(weighting,"absolute_row_weights",nothing)
                if weights!==nothing
                    weights isa AbstractVector && length(weights)==absolute_count ||
                        throw(ArgumentError("Empirical absolute_row_weights must match all original absolute-pick rows before ABS exclusions."))
                    relative_controls=merge(relative_controls,(empirical_weighting=
                        merge(weighting,Dict("absolute_row_weights"=>weights[absolute_rows])),))
                end
            end
        end
        relative_solve_options=merge(relative_controls,options,(tasks,
            selection=relative_selection,selection_frame=frame,
            initialization=relative_initialization,initialization_frame=initialization_coordinates=="geographic" ? frame : nothing,
            reference_event_ids=get(config,"reference_event_ids",String[])))
        relative_progress=progress===nothing ? nothing : (done,total,index,row)->progress(
            (;stage=:relative,processed=done,total,component_id=row.component_id,status=row.result.summary.status))
        relative = _case_phase(:relative_solve,timings,progress) do
            relocate(model,stations,events,differences;absolute_picks=relative_picks,
                progress=relative_progress,relative_solve_options...)
        end
        if output_diagnostics && haskey(config,"pairing")
            criteria = config["pairing"]
            geometry = differential_geometry(relative.events,stations,differences;
                max_separation_km=criteria["max_separation_km"],
                max_separation_ratio=criteria["max_separation_ratio"],
                min_common_stations=criteria["min_common_stations"])
            write_results && _write_result_csv(joinpath(output,"final_pair_geometry.csv"),
                ("event_i","event_j","separation_km","max_original_path_ratio",
                 "common_path_station_count","meets_declared_geometry"),(Tuple(row) for row in geometry))
        end
        uncertainty_settings=get(config,"uncertainty",Dict())
        if haskey(uncertainty_settings,"relative_pairs") || haskey(uncertainty_settings,"relative_events")
            requested=get(uncertainty_settings,"relative_pairs",[])
            all(length(pair)==2 for pair in requested) || throw(ArgumentError("Each uncertainty pair needs two IDs."))
            failed_ids=Set(row.event_id for row in failed)
            retained_pairs=[pair for pair in requested if !(pair[1] in failed_ids || pair[2] in failed_ids)]
            requested_events=get(uncertainty_settings,"relative_events",String[])
            requested_events=="all" && (requested_events=vcat([e.event_id for e in relative.events],[e.event_id for e in failed]))
            requested_events isa AbstractVector || throw(ArgumentError("relative_events must be an ID array or all."))
            retained_events=[id for id in requested_events if !(id in failed_ids)]
            estimates = _case_phase(:relative_uncertainty,timings,progress) do
                relative_uncertainty(model,stations,events,differences,relative;event_pairs=retained_pairs,
                    event_ids=retained_events,event_confidence_probability=get(uncertainty_settings,"relative_confidence_probability",.95),
                    absolute_picks=relative_picks,
                    assume_independent=get(uncertainty_settings,"assume_independent_relative_observations",false))
            end
            method=get(uncertainty_settings,"relative_covariance_method","analytic")
            method in ("analytic","algorithm_replay") || throw(ArgumentError("Relative covariance method must be analytic or algorithm_replay."))
            if method=="algorithm_replay"
                get(relative_solve_options,:solver,:damped)==:damped ||
                    throw(ArgumentError("Algorithm replay describes the column-scaled damped relocation."))
                isempty(retained_events) || begin
                    algorithm = _case_phase(:relative_algorithm_uncertainty,timings,progress) do
                        relative_algorithm_uncertainty(model,stations,events,differences,relative;event_ids=retained_events,
                            absolute_picks=relative_picks,relocate_options=relative_solve_options,
                            assume_independent=get(uncertainty_settings,"assume_independent_relative_observations",false),
                            event_confidence_probability=get(uncertainty_settings,"relative_confidence_probability",.95),
                            step=get(uncertainty_settings,"relative_replay_step",.02),
                            method=Symbol(get(uncertainty_settings,"relative_algorithm_method","auto")),
                            noise_samples=get(uncertainty_settings,"relative_algorithm_noise_samples",128),
                            seed=get(uncertainty_settings,"relative_algorithm_seed",0),replay_tasks,progress)
                    end
                    estimates=merge(estimates,(events=algorithm.events,algorithm_uncertainty=algorithm,
                        assumptions=merge(estimates.assumptions,(individual_events=algorithm.assumptions,))))
                end
            end
            reported=NamedTuple[]
            retained_index=0
            all_ids=Set{String}(row.event_id for row in events)
            union!(all_ids,failed_ids)
            for pair in requested
                all(id in all_ids for id in pair) && pair[1]!=pair[2] || throw(ArgumentError("Invalid uncertainty event pair."))
                if pair[1] in failed_ids || pair[2] in failed_ids
                    push!(reported,(event_i=String(pair[1]),event_j=String(pair[2]),component_id=nothing,
                        status=:not_estimated,reason=:event_not_entered_relative,parameter_names=(),covariance=nothing,
                        standard_deviations=(x_km=nothing,y_km=nothing,z_km=nothing,origin_s=nothing),rank=nothing))
                else
                    retained_index+=1
                    push!(reported,estimates.pairs[retained_index])
                end
            end
            event_errors=NamedTuple[]
            retained_index=0
            probability=get(uncertainty_settings,"relative_confidence_probability",.95)
            scale=sqrt(_spatial_chisq_quantile(probability,1))
            for id in requested_events
                id in all_ids || throw(ArgumentError("Invalid individual uncertainty event ID."))
                if id in failed_ids
                    event=only(e for e in failed if e.event_id==id)
                    push!(event_errors,(event_id=id,component_id=nothing,source_km=event.source_km,origin_s=event.origin_s,
                        status=:not_estimated,reason=:event_not_entered_relative,parameter_names=(),covariance=nothing,
                        standard_deviations=(x_km=nothing,y_km=nothing,z_km=nothing,origin_s=nothing),
                        marginal_confidence_half_widths=(x_km=nothing,y_km=nothing,z_km=nothing,origin_s=nothing),
                        confidence_probability=probability,confidence_scale=scale,rank=nothing,reference_event_ids=String[],
                        assumptions=(quantity=:event_conditional_coordinates,relative_stage_entered=false)))
                else
                    retained_index+=1
                    push!(event_errors,estimates.events[retained_index])
                end
            end
            relative=merge(relative,(uncertainty=merge(estimates,(pairs=reported,events=event_errors)),))
        end
        if !isempty(failed)
            # Numerical routines/replays use their retained input order. Map final
            # reading indices back only after those calculations are complete.
            relative=merge(relative,(events=NamedTuple[relative.events...],
                residuals=NamedTuple[merge(row,(input_row=differential_rows[row.input_row],)) for row in relative.residuals],
                absolute_residuals=NamedTuple[merge(row,(input_row=absolute_rows[row.input_row],)) for row in relative.absolute_residuals]))
            for event in failed
                count_diff=count(row->row.event_i==event.event_id || row.event_j==event.event_id,excluded)
                count_abs=joint ? count(row->row.event_id==event.event_id,picks) : 0
                push!(relative.events,(event_id=event.event_id,source_km=event.source_km,origin_s=event.origin_s,
                    status=:not_entered_relative,component_id=nothing,iterations=0,input_count=count_diff+count_abs,
                    used_count=0,depth_at_bound=event.depth_at_bound,previous_absolute_status=event.status))
            end
            for row in excluded
                push!(relative.residuals,merge(row,(dt_pred_s=nothing,residual_s=nothing,normalized_residual=nothing,
                    used=false,unused_reason="absolute_location_failed",component_id=nothing,
                    equation_weight=0.,weighted_residual=nothing,selection_status=:not_entered_relative)))
            end
            if joint
                failed_ids=Set(row.event_id for row in failed)
                for (input_row,row) in enumerate(picks)
                    row.event_id in failed_ids || continue
                    push!(relative.absolute_residuals,merge(row,(arrival_pred_s=nothing,residual_s=nothing,normalized_residual=nothing,
                        used=false,unused_reason="absolute_location_failed",component_id=nothing,input_row,
                        equation_weight=0.,weighted_residual=nothing,selection_status=:not_entered_relative)))
                end
            end
            sort!(relative.residuals;by=row->row.input_row)
            sort!(relative.absolute_residuals;by=row->row.input_row)
            # Keep the source catalog order even when failed events are appended.
            order=Dict(row.event_id=>i for (i,row) in enumerate(absolute.events))
            sort!(relative.events;by=row->order[row.event_id])
        end
        relative = _case_phase(:relative_analysis,timings,progress) do
            depth_settings=get(config,"depth_reporting",Dict())
            if !isempty(depth_settings)
                report_stations = Dict(event.event_id=>String[] for event in relative.events)
                for row in relative.residuals, id in (row.event_i,row.event_j)
                    push!(report_stations[id],row.station_id)
                end
                for row in relative.absolute_residuals
                    push!(report_stations[row.event_id],row.station_id)
                end
                if absolute !== nothing
                    for row in absolute.residuals
                        push!(get!(report_stations,row.event_id,String[]),row.station_id)
                    end
                end
                relative = _with_depth_reports(relative,model,stations,report_stations,depth_settings)
            end
            relative = _with_geographic_regions(relative,frame,get(config,"geographic_regions",Dict()),base)
            if output_diagnostics
                relative=merge(relative,(observation_diagnostics=relative_observation_diagnostics(model,stations,relative),))
            end
            relative
        end
        if write_results
            _case_phase(:relative_export,timings,progress) do
                write_relative_results(joinpath(output, "relative"), relative;detail=output_detail,diagnostics=output_diagnostics)
                output_diagnostics && _write_geographic_regions(joinpath(output,"relative"),relative.events)
                frame !== nothing && write_geographic_results(joinpath(output,"relative"),relative,frame;time_reference_utc=reference)

            end
        end
    end
    if task in ("relative","absolute_then_relative")
        stage_summaries["relative"]=_case_stage_summary(relative)
        push!(completed_stages,"relative")
    end
    return (;absolute,relative,output_directory=write_results ? output : nothing)
end
