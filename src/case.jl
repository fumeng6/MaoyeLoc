const _CASE_STAGE_FILES = ("search.toml", "search_candidates.csv", "output_detail.toml", "events.csv", "residuals.csv", "iterations.csv", "spectral_directions.csv", "components.csv", "geographic.csv", "datum_stations.csv", "absolute_residuals.csv", "conditional_covariance.csv", "pair_uncertainty.csv", "conditional_pair_covariance.csv", "source_iterations.csv", "column_scales.csv", "solver_directions.csv", "weight_iterations.csv","selection_events.csv","selection_clusters.csv","selection_pairs.csv","selection_settings.toml","formal_errors.csv","formal_covariance.csv","error_ellipsoids.csv","error_axes.csv","information.csv","reading_information.csv","station_information.csv","information_directions.csv","event_geometry.csv","reading_geometry.csv","ray_components.csv","regions.toml","processing_metadata.toml","event_uncertainty.csv","conditional_event_covariance.csv","uncertainty_assumptions.toml","event_statistics.csv","station_statistics.csv","component_statistics.csv","unused_statistics.csv","statistics.toml","observation_diagnostics.csv","absolute_observation_diagnostics.csv","observation_ray_components.csv","centroid_iterations.csv","correction_iterations.csv","run_summary.toml","algorithm_uncertainty_numerics.toml","magnitude_calibration.csv","magnitude_calibration.toml","coda_station_magnitudes.csv","coda_station_magnitudes.toml","amplitude_station_magnitudes.csv","amplitude_station_magnitudes.toml","event_magnitudes.csv","event_magnitudes.toml","magnitude_selection.csv","preferred_magnitudes.csv","magnitude_display.csv")

"""
    prepare_case(path)
    prepare_case(config; base_directory=pwd(), inputs=(;))

Read a native case and resolve its actual inputs and options without solving,
creating directories or removing outputs. Treat the returned data as read-only;
prepare again after changing scientific inputs or configuration.
`inputs` may replace native CSV tables with `stations` (ID => XYZ dictionary)
and `events`, `picks`, `differentials` (NamedTuple record vectors). Explicit
reading IDs survive; missing pick IDs use `memory_row_N` in original input order.
"""
function prepare_case(path::AbstractString)
    source=abspath(path)
    _prepare_case(TOML.parsefile(source);base_directory=dirname(source),config_path=source)
end

function prepare_case(config::AbstractDict;base_directory=pwd(),inputs=(;))
    _prepare_case(deepcopy(config);base_directory,config_path=nothing,inputs)
end

function _case_keys(table,allowed,name)
    table isa AbstractDict || throw(ArgumentError("$name must be a TOML table."))
    unknown=sort!([String(key) for key in keys(table) if !(key in allowed)])
    isempty(unknown) || throw(ArgumentError("Unknown $name key(s): $(join(unknown, ", ")). Valid keys: $(join(allowed, ", "))."))
end

# Windows treats these spellings as the same file. Use the same comparison
# for preparation and output cleanup.
function _case_path_key(path)
    absolute=normpath(abspath(path))
    Sys.iswindows() ? lowercase(absolute) : absolute
end

function _case_input_paths(config,base,source_config)
    paths=[joinpath(base,path) for path in values(config["files"])]
    source_config===nothing || push!(paths,source_config)
    processing=get(config,"processing_metadata",nothing)
    processing isa AbstractString && push!(paths,joinpath(base,processing))
    function model_paths(table)
        for (key,value) in table
            if key in ("p_table","s_table") && value isa AbstractString
                push!(paths,joinpath(base,value))
            elseif value isa AbstractDict
                model_paths(value)
            end
        end
    end
    model_paths(config["model"])
    regions=get(config,"geographic_regions",Dict())
    net=get(regions,"net",0)
    if net!=0
        directory=joinpath(base,get(regions,"reference_directory",_default_data_directory()))
        push!(paths,joinpath(directory,"geographic_regions.toml"))
    end
    Set(_case_path_key(path) for path in paths)
end

# One absolute path: the HypoInverse-style empirical weights and SVD step.
# Explicit keys override the same controls; they do not select another estimator.
function _absolute_method_config(settings,weighting,search)
    search===nothing || throw(ArgumentError("absolute_search was removed. Absolute location has one empirical SVD path."))
    haskey(settings,"loss") && throw(ArgumentError("solver.loss was removed. Absolute location does not select a separate least-squares or Huber loss."))
    haskey(settings,"huber_c") && throw(ArgumentError("solver.huber_c was removed with the Huber path."))
    method=get(settings,"method","hypoinverse")
    method=="hypoinverse" || throw(ArgumentError("Absolute location has one path. solver.method must be omitted or hypoinverse; least_squares, sigma_spectral and huber are not available."))
    weighting!==nothing && get(weighting,"mode","empirical")!="empirical" &&
        throw(ArgumentError("Absolute location uses empirical weights. absolute_weighting.mode=sigma is not available."))
    resolved=merge(Dict{String,Any}("eigtol"=>.012,"singular_value_shift"=>.006,"dxfix_km"=>7.,
        "dzmax_km"=>30.,"dzair"=>.5,"damp"=>.9,"late_damping_factor"=>.5,
        "dxmax_km"=>50.,"rback_s"=>.02,"bacfac"=>.6,"d2far_km"=>250.,
        "stopping_rule"=>"legacy","dquit_km"=>.04,"drqt_s"=>.001,
        "stop_start_iteration"=>1,"legacy_iteration_limit"=>(get(settings,"stopping_rule","legacy")=="step" ? nothing : 20)),settings)
    weighting=merge(Dict{String,Any}("mode"=>"empirical","quality_weights"=>[1.,.75,.5,.25],"s_weight"=>.5,
        "distance"=>Dict("start_iteration"=>4,"minimum_radius_km"=>50.,"inner_factor"=>1.,"outer_factor"=>3.),
        "initial_distance"=>Dict("start_iteration"=>100,"minimum_radius_km"=>50.,"inner_factor"=>3.,"outer_factor"=>6.),
        "residual"=>Dict("start_iteration"=>4,"floor_s"=>.16,"inner_factor"=>1.5,"outer_factor"=>3.,"passes"=>2)),
        something(weighting,Dict()))
    resolved,weighting,method
end

function _prepare_case(config; base_directory,config_path,inputs=(;))
    _case_keys(config,("task","time_reference_utc","coordinates","vertical_reference",
        "output_directory","files","model","frame","initial_depth_km","solver",
        "absolute_constraints","absolute_weighting","relative","reference_event_ids",
        "pairing","uncertainty","diagnostics","geographic_regions","magnitudes",
        "processing_metadata","depth_reporting","output","execution","absolute_search"),"case")
    base=abspath(base_directory)
    inputs isa NamedTuple && all(key in (:stations,:events,:picks,:differentials) for key in keys(inputs)) ||
        throw(ArgumentError("inputs must be a NamedTuple with stations, events, picks and/or differentials."))
    _case_keys(get(config,"solver",Dict()),(
        "depth_min_km","depth_max_km","max_iterations","location_tolerance_km","time_tolerance_s",
        "eigtol","singular_value_shift","dxfix_km","dzmax_km","dzair","damp","late_damping_factor",
        "dxmax_km","rback_s","bacfac","d2far_km","station_site_ids","stopping_rule","dquit_km",
        "drqt_s","stop_start_iteration","legacy_iteration_limit","loss","huber_c","method"),"solver")
    _case_keys(get(config,"relative",Dict()),(
        "use_absolute_picks","absolute_failure_policy","initialization","initialization_coordinates",
        "solver","row_weighting","row_weights","damp","svd_relative_cutoff","airquake_policy",
        "airquake_depth_km","linear_solver","lsqr_rtol","lsqr_max_iterations","lsqr_condition_limit",
        "weighting","selection"),"relative")
    _case_keys(get(config,"execution",Dict()),("tasks","replay_tasks"),"execution")
    _case_keys(get(config,"output",Dict()),("detail","diagnostics"),"output")
    tasks=get(get(config,"execution",Dict()),"tasks",1)
    tasks isa Integer && !(tasks isa Bool) && tasks>=1 ||
        throw(ArgumentError("execution.tasks must be a positive integer. Start Julia with --threads=N to enable concurrent tasks."))
    replay_tasks=get(get(config,"execution",Dict()),"replay_tasks",1)
    replay_tasks isa Integer && !(replay_tasks isa Bool) && replay_tasks>=1 ||
        throw(ArgumentError("execution.replay_tasks must be a positive integer."))
    tasks>1 && replay_tasks>1 &&
        throw(ArgumentError("Choose execution.tasks>1 or execution.replay_tasks>1, not both, to avoid nested solve workspaces."))
    output_detail=Symbol(get(get(config,"output",Dict()),"detail","full"))
    output_detail in (:core,:full) || throw(ArgumentError("output.detail must be core or full."))
    output_diagnostics=get(get(config,"output",Dict()),"diagnostics",false)
    output_diagnostics isa Bool || throw(ArgumentError("output.diagnostics must be true or false."))
    task = config["task"]
    task in ("absolute", "relative", "absolute_then_relative") || throw(ArgumentError("Unknown task."))
    haskey(config,"absolute_search") &&
        throw(ArgumentError("absolute_search was removed. Absolute location has one empirical SVD path."))
    task=="relative" && any(haskey(get(config,"solver",Dict()),key) for key in ("loss","huber_c","method")) &&
        throw(ArgumentError("solver.method, loss and huber_c apply to absolute location only."))
    task=="relative" && haskey(config,"absolute_constraints") &&
        throw(ArgumentError("absolute_constraints applies to an absolute stage; relative constraints use reference_event_ids and depth bounds."))
    reference = config["time_reference_utc"]
    endswith(reference, "Z") || throw(ArgumentError("time_reference_utc must use UTC with a trailing Z."))
    DateTime(reference[1:end-1]) # Validate the calendar reference, not just its label.
    for key in ("coordinates", "vertical_reference")
        isempty(strip(config[key])) && throw(ArgumentError("Missing coordinate metadata: $key"))
    end
    frame = if haskey(config,"frame")
        info = config["frame"]
        info["geographic_crs"] == "EPSG:4326" && info["horizontal_projection"] == "WGS84_local_orthographic" ||
            throw(ArgumentError("Geographic output requires the declared WGS84 local orthographic frame."))
        info["vertical_datum"] == config["vertical_reference"] || throw(ArgumentError("Inconsistent vertical datum declarations."))
        LocalFrame(info["latitude_deg"],info["longitude_deg"];vertical_datum=info["vertical_datum"])
    else
        nothing
    end
    configured_files=get(config,"files",Dict{String,String}())
    configured_files isa AbstractDict || throw(ArgumentError("files must be a TOML table."))
    # An explicit memory table replaces its file, including in path metadata
    # and collision checks; the caller's dictionary itself is never changed.
    files=Dict{String,Any}(key=>value for (key,value) in configured_files if !(Symbol(key) in keys(inputs)))
    config["files"]=files
    input_path(key) = normpath(joinpath(base, files[key]))
    model = _read_model_table(config["model"];base)
    joint = get(get(config,"relative",Dict()),"use_absolute_picks",false)
    absolute_failure_policy=get(get(config,"relative",Dict()),"absolute_failure_policy","stop")
    absolute_failure_policy in ("stop","exclude") || throw(ArgumentError("absolute_failure_policy must be stop or exclude."))
    stations = haskey(inputs,:stations) ? _memory_stations(inputs.stations) : read_stations(input_path("stations"))
    supplied_differences=haskey(inputs,:differentials) || haskey(files,"differentials")
    picks = task == "relative" && supplied_differences && !joint ? nothing :
        haskey(inputs,:picks) ? _memory_picks(inputs.picks) : read_picks(input_path("picks"))
    events = if haskey(inputs,:events)
        _memory_events(inputs.events)
    elseif haskey(files, "events")
        read_events(input_path("events"))
    elseif task != "relative" && haskey(config, "initial_depth_km")
        initialize_catalog(model, stations, picks; depth_km=config["initial_depth_km"])
    else
        throw(ArgumentError("Provide files.events or inputs.events, or initial_depth_km for an absolute task."))
    end
    task!="relative" && (events=_absolute_catalog_constraints(events,get(config,"absolute_constraints",Dict())))
    station_sets = model isa RegionalModels ?
        [branch for branch in values(model.models) if branch isa StationModels] :
        model isa StationModels ? [model] : StationModels[]
    for station_set in station_sets
        all(haskey(stations,id) for id in keys(station_set.station_model_ids)) ||
            throw(ArgumentError("Station-model assignment refers to an unknown station ID."))
        all(haskey(stations,key[2]) for key in keys(station_set.station_delays_s)) ||
            throw(ArgumentError("Station delay refers to an unknown station ID."))
    end
    output = normpath(joinpath(base, config["output_directory"]))
    settings = get(config, "solver", Dict())
    provided_solver_keys=sort!(collect(keys(settings)))
    weighting=get(config,"absolute_weighting",nothing)
    absolute_method="not_applied"
    if task!="relative"
        settings,weighting,absolute_method=_absolute_method_config(settings,weighting,get(config,"absolute_search",nothing))
        config["solver"]=settings
        weighting!==nothing && (config["absolute_weighting"]=weighting)
    end
    bounds = (get(settings, "depth_min_km", model_top(model)), get(settings, "depth_max_km", Inf))
    options = (depth_bounds=bounds, max_iterations=get(settings, "max_iterations", 50),
               location_tolerance_km=get(settings, "location_tolerance_km", 1e-6),
               time_tolerance_s=get(settings, "time_tolerance_s", 1e-7))
    options.max_iterations isa Integer && options.max_iterations>=0 &&
        options.location_tolerance_km>0 && options.time_tolerance_s>0 ||
        throw(ArgumentError("solver requires nonnegative integer max_iterations and positive step tolerances."))
    isfinite(bounds[1]) && bounds[1]>=model_top(model) && bounds[2]>=bounds[1] ||
        throw(ArgumentError("Depth bounds must lie within the model, lower <= upper."))
    # A fixed keyword layout avoids compiling a new orchestration method for
    # every subset of optional TOML keys. Values are the numerical API defaults.
    absolute_defaults=(eigtol=.012,singular_value_shift=.006,dxfix_km=7.,dzmax_km=30.,dzair=.5,
        damp=.9,late_damping_factor=.5,dxmax_km=50.,rback_s=.02,bacfac=.6,d2far_km=250.,
        station_site_ids=Dict{String,String}(),stopping_rule=:legacy,dquit_km=.04,drqt_s=.001,
        stop_start_iteration=1,legacy_iteration_limit=20)
    absolute_controls=(;(key=>(key==:stopping_rule ? Symbol(get(settings,String(key),value)) :
        get(settings,String(key),value)) for (key,value) in pairs(absolute_defaults))...)
    absolute_options = merge(options,absolute_controls)
    if weighting!==nothing
        weighting isa AbstractDict && get(weighting,"mode",nothing)=="empirical" ||
            throw(ArgumentError("absolute_weighting.mode must be empirical."))
    end
    w=something(weighting,Dict{String,Any}())
    absolute_options=merge(absolute_options,(row_weighting=Symbol(get(w,"mode","empirical")),
        quality_weights=get(w,"quality_weights",[1.,.75,.5,.25]),
        station_time_weights=get(w,"station_time_weights",Dict{String,Any}()),s_weight=get(w,"s_weight",.5),
        distance_weighting=get(w,"distance",nothing),initial_distance_weighting=get(w,"initial_distance",nothing),
        residual_weighting=get(w,"residual",nothing),event_rules=get(w,"event_rules",nothing),
        default_quality=get(w,"default_quality",nothing)))
    relative_settings = get(config,"relative",Dict())
    relative_initialization = Symbol(get(relative_settings,"initialization","catalog"))
    relative_initialization in (:catalog,:centroid) ||
        throw(ArgumentError("relative.initialization must be catalog or centroid."))
    dd_weighting=get(relative_settings,"weighting",nothing)
    relative_parameter_defaults=(solver=:damped,row_weighting=:empirical,
        row_weights=nothing,damp=0.,svd_relative_cutoff=1e-7,airquake_policy=:bounds,airquake_depth_km=0.,
        linear_solver=:auto,lsqr_rtol=1e-6,lsqr_max_iterations=0,lsqr_condition_limit=1e12,
        empirical_weighting=nothing)
    relative_controls=(;(key=>(value isa Symbol ? Symbol(get(relative_settings,String(key),value)) :
        get(relative_settings,String(key),value)) for (key,value) in pairs(relative_parameter_defaults))...)
    if dd_weighting===nothing
        stages=_dd_default_stages()
        haskey(relative_settings,"damp") && foreach(stage -> stage["damp"]=relative_settings["damp"], stages)
        relative_controls=merge(relative_controls,(;empirical_weighting=Dict{String,Any}(
            "stages"=>stages,"default_quality"=>1.)))
    else
        dd_weighting isa AbstractDict && get(dd_weighting,"mode",nothing)=="empirical" ||
            throw(ArgumentError("relative.weighting needs explicit mode=empirical."))
        get(relative_controls,:row_weighting,:empirical)==:empirical || throw(ArgumentError("Empirical weighting conflicts with relative.row_weighting."))
        relative_controls=merge(relative_controls,(row_weighting=:empirical,
            empirical_weighting=Dict{String,Any}(k=>v for (k,v) in dd_weighting if k!="mode")))
    end
    relative_selection=get(relative_settings,"selection",nothing)
    if relative_selection===nothing && haskey(something(relative_controls.empirical_weighting,Dict()),"stages")
        relative_selection=Dict{String,Any}("idat"=>3,"ipha"=>3,"obscc"=>0,"obsct"=>0)
    end
    initialization_coordinates = get(relative_settings,"initialization_coordinates","cartesian")
    initialization_coordinates in ("cartesian","geographic") ||
        throw(ArgumentError("relative.initialization_coordinates must be cartesian or geographic."))
    initialization_coordinates=="geographic" && frame===nothing &&
        throw(ArgumentError("Geographic centroid initialization requires the catalog's explicit frame."))
    _case_keys(files,("stations","events","picks","differentials"),"files")
    differences=task=="absolute" ? nothing : haskey(inputs,:differentials) ?
        _memory_differentials(inputs.differentials) : haskey(files,"differentials") ?
        read_differentials(input_path("differentials")) : nothing
    task!="absolute" && differences===nothing && !haskey(config,"pairing") &&
        throw(ArgumentError("Supply files.differentials or inputs.differentials, or explicit [pairing] criteria with picks."))
    ids=Set(String(event.event_id) for event in events)
    for pick in something(picks,NamedTuple[])
        pick.event_id in ids || throw(ArgumentError("Pick refers to unknown event $(pick.event_id)."))
        haskey(stations,pick.station_id) || throw(ArgumentError("Pick refers to unknown station $(pick.station_id)."))
    end
    for row in something(differences,NamedTuple[])
        row.event_i in ids && row.event_j in ids ||
            throw(ArgumentError("Difference refers to an unknown catalog event."))
        haskey(stations,row.station_id) || throw(ArgumentError("Difference refers to unknown station $(row.station_id)."))
    end
    all(id in ids for id in get(config,"reference_event_ids",String[])) ||
        throw(ArgumentError("A reference_event_id is not in the input catalog."))
    _case_validate_options(config,task,absolute_options,relative_controls,events,frame,stations)
    input_files=_case_input_paths(config,base,config_path)
    stages=task=="absolute_then_relative" ? ("absolute","relative") : (task,)
    native_outputs=[joinpath(output,stage,name) for stage in stages for name in _CASE_STAGE_FILES]
    append!(native_outputs,[joinpath(output,name) for name in ("calculation.toml","run_summary.toml",
        "generated_ct.csv","pairing_events.csv",
        "pairing_pairs.csv","final_pair_geometry.csv")])
    for path in native_outputs
        _case_path_key(path) in input_files &&
            throw(ArgumentError("An output file is also a configured input: $path"))
    end
    solver_metadata = merge(Dict{String,Any}("method"=>absolute_method,"depth_min_km"=>bounds[1],"depth_max_km"=>bounds[2],
        "max_iterations"=>options.max_iterations,"location_tolerance_km"=>options.location_tolerance_km,
        "time_tolerance_s"=>options.time_tolerance_s,
        "absolute_origin_handling"=>"full physical XYZ/origin singular subspace; no separate origin profiling",
        "absolute_spatial_step"=>"row-normalized full SVD; absolute EIGTOL cutoff; gain 1/(s+shift); explicit step controls"),
        Dict(String(key)=>value isa Symbol ? String(value) : value
            for (key,value) in pairs(absolute_controls) if value!==nothing))
    metadata = Dict("task" => task, "time_reference_utc" => reference,
            "coordinates" => config["coordinates"], "vertical_reference" => config["vertical_reference"],
            "model" => config["model"], "solver" => solver_metadata,
            "depth_reporting"=>get(config,"depth_reporting",Dict()),
            "absolute_constraints"=>get(config,"absolute_constraints",Dict()),
            "absolute_weighting"=>weighting===nothing ? Dict("mode"=>"empirical") : weighting,
            "absolute_constraint_scope"=>"absolute stage only; relative stage uses its own reference_event_ids and depth bounds",
            "input_config" => config_path===nothing ? "in-memory" : config_path,
            "input_format"=>"native_csv_or_memory", "input_files"=>Dict(key=>input_path(key) for key in keys(files)),
            "event_initialization" => haskey(files, "events") || haskey(inputs,:events) ? "supplied catalog" : "station centroid at explicit depth",
            "initial_depth_km" => get(config, "initial_depth_km",
                haskey(inputs,:events) ? "supplied in event records" : "supplied in event file"),
            "frame" => get(config,"frame",Dict()),
            "differential_source" => haskey(inputs,:differentials) ? "supplied in-memory records" :
                haskey(files,"differentials") ? "supplied file" : "catalog pairing if relative requested",
            "pairing" => haskey(files,"differentials") || haskey(inputs,:differentials) ?
                Dict("applied"=>false) : get(config,"pairing",Dict()),
            "relative"=>merge(Dict{String,Any}("use_absolute_picks"=>joint,"absolute_failure_policy"=>absolute_failure_policy,
                "selection"=>relative_selection===nothing ? Dict("enabled"=>false) : relative_selection,
                "initialization"=>String(relative_initialization),"initialization_coordinates"=>initialization_coordinates,
                "observation_covariance"=>"relocation keeps diagonal row weights; an exact catalog CT difference uses the shared pick covariance only when formal errors are requested",
                "solver"=>"damped","row_weighting"=>"empirical",
                "origin_handling"=>"joint full XYZ/origin update in scaled parameter space; damping holds the cluster near the start; no automatic hard reference"),
                Dict(String(key)=>value isa Symbol ? String(value) : value for (key,value) in pairs(relative_controls) if value!==nothing)),
            "uncertainty"=>get(config,"uncertainty",Dict("assume_independent_picks"=>false)),
            "diagnostics"=>get(config,"diagnostics",Dict("absolute"=>false)),
            "geographic_regions"=>get(config,"geographic_regions",Dict("net"=>0)),
            "magnitudes"=>get(config,"magnitudes",Dict()),
            "notes" => "Physical input metadata only; stage CSV statuses report outcomes. Residual=observed-predicted; km, km/s, s. Conditional uncertainty, when requested and estimable, is reported per event; it excludes velocity/station/model bias. DD datum is specified per component.")
    relative_defaults=Dict{String,Any}("solver"=>"damped","row_weighting"=>"empirical","damp"=>0.,
        "svd_relative_cutoff"=>1e-7,"linear_solver"=>"auto","lsqr_rtol"=>1e-6,
        "lsqr_max_iterations"=>0,"lsqr_condition_limit"=>1e12,"airquake_policy"=>"bounds",
        "airquake_depth_km"=>0.)
    metadata["relative"]=merge(relative_defaults,metadata["relative"])
    metadata["relative"]["damp_applied"]=true
    metadata["reference_event_ids"]=get(config,"reference_event_ids",String[])
    metadata["output_directory"]=output
    metadata["output"]=Dict("detail"=>String(output_detail),"diagnostics"=>output_diagnostics)
    metadata["execution"]=Dict("tasks"=>tasks,"replay_tasks"=>replay_tasks,
        "available_julia_threads"=>Threads.nthreads(),
        "maximum_concurrent_tasks"=>min(tasks,Threads.nthreads()),
        "maximum_concurrent_replays"=>min(replay_tasks,4,Threads.nthreads()),"blas_threads"=>BLAS.get_num_threads())
    metadata["provided_solver_keys"]=provided_solver_keys
    metadata["provided_relative_keys"]=sort!(collect(keys(relative_settings)))
    if !isempty(inputs)
        used_file_keys=[key for key in keys(files) if key in ("stations","events") ||
            key=="picks" && picks!==nothing || key=="differentials" && differences!==nothing]
        metadata["input_files"]=Dict(key=>input_path(key) for key in used_file_keys)
        metadata["input_sources"]=merge(Dict{String,Any}(key=>Dict("kind"=>"file","path"=>input_path(key)) for key in used_file_keys),
            Dict{String,Any}(String(key)=>Dict("kind"=>"memory","order"=>"original input order",
                "used"=>key==:picks ? picks!==nothing : key==:differentials ? differences!==nothing : true,
                "reading_identity"=>key==:picks ? "explicit observation_id; otherwise memory_row_N (one-based input index)" : "supplied string IDs")
                for key in keys(inputs)))
    end
    (;config_path,config,base,files,task,reference,frame,model,joint,absolute_failure_policy,
      stations,events,picks,differences,output,
      bounds,options,absolute_controls,absolute_options,relative_controls,relative_initialization,
      relative_selection,initialization_coordinates,metadata,tasks=Int(tasks),replay_tasks=Int(replay_tasks),
      output_detail,output_diagnostics)
end

"""Return prepared scientific settings and the current execution resources."""
function effective_config(prepared::NamedTuple)
    metadata=deepcopy(prepared.metadata)
    metadata["execution_at_preparation"]=deepcopy(metadata["execution"])
    merge!(metadata["execution"],Dict("available_julia_threads"=>Threads.nthreads(),
        "maximum_concurrent_tasks"=>min(prepared.tasks,Threads.nthreads()),
        "maximum_concurrent_replays"=>min(prepared.replay_tasks,4,Threads.nthreads()),
        "blas_threads"=>BLAS.get_num_threads()))
    metadata
end
effective_config(path::AbstractString)=effective_config(prepare_case(path))
effective_config(config::AbstractDict;kwargs...)=effective_config(prepare_case(config;kwargs...))

"""Check native inputs and effective settings without solving or writing outputs.
This does not establish ray reachability, numerical convergence or field accuracy.
"""
function check_case(prepared::NamedTuple)
    p=prepared
    requested=p.task=="absolute_then_relative" ? ["absolute","relative"] : [p.task]
    stale=[stage for stage in ("absolute","relative") if !(stage in requested) &&
        isdir(joinpath(p.output,stage))]
    Dict("status"=>"inputs_checked","task"=>p.task,"event_count"=>length(p.events),
        "station_count"=>length(p.stations),"pick_count"=>length(something(p.picks,NamedTuple[])),
        "differential_count"=>length(something(p.differences,NamedTuple[])),
        "output_directory"=>p.output,"requested_stages"=>requested,
        "existing_output_directory"=>isdir(p.output),"previous_stage_directories"=>stale,
        "tasks"=>p.tasks,"maximum_concurrent_tasks"=>min(p.tasks,Threads.nthreads()),
        "replay_tasks"=>p.replay_tasks,"maximum_concurrent_replays"=>min(p.replay_tasks,4,Threads.nthreads()),
        "scope"=>"Parsed inputs and settings only; rays and convergence are not evaluated.")
end
check_case(path::AbstractString)=check_case(prepare_case(path))
check_case(config::AbstractDict;kwargs...)=check_case(prepare_case(config;kwargs...))

"""
    run_case(path; progress=nothing, write_results=true)
    run_case(config; base_directory=pwd(), inputs=(;), progress=nothing, write_results=true)
    run_case(prepared; progress=nothing, write_results=true)

Run an explicit absolute, relative or absolute_then_relative case. File paths
are relative to the TOML file or the dictionary's explicit base directory.
Preparation is shared with `check_case` and `effective_config`. A progress
callback receives completed event/component records and phase notifications
(`kind=:phase`), including directions used to perturb observations and rerun
the complete relocation for error propagation. With
`write_results=false`, the same scientific path runs without creating, cleaning
or writing outputs, including scientific reports. Results include `calculation`
settings and `run_summary`; `output_directory` is then `nothing`.
"""
run_case(path::AbstractString;progress=nothing,write_results::Bool=true)=run_case(prepare_case(path);progress,write_results)
run_case(config::AbstractDict;base_directory=pwd(),inputs=(;),progress=nothing,write_results::Bool=true)=
    run_case(prepare_case(config;base_directory,inputs);progress,write_results)

function _case_stage_summary(result)
    result===nothing && return Dict{String,Any}("executed"=>false)
    counts=Dict{String,Int}()
    for event in result.events
        key=String(event.status)
        counts[key]=get(counts,key,0)+1
    end
    absolute_rows=get(result,:absolute_residuals,NamedTuple[])
    used=count(row->get(row,:used,false),result.residuals)+count(row->get(row,:used,false),absolute_rows)
    input_count=length(result.residuals)+length(absolute_rows)
    Dict{String,Any}("executed"=>true,"event_count"=>length(result.events),"status_counts"=>counts,
        "catalog_excluded"=>count(event->!get(event,:catalog_eligible,true),result.events),
        "input_observations"=>input_count,"used_observations"=>used,"unused_observations"=>input_count-used)
end

function run_case(prepared::NamedTuple;progress=nothing,write_results::Bool=true)
    completed=String[]
    stage_summaries=Dict{String,Any}()
    timings=Dict{String,Float64}()
    calculation=effective_config(prepared)
    calculation["output"]=merge(calculation["output"],Dict("write_results"=>write_results))
    summary=merge(check_case(prepared),Dict{String,Any}("status"=>"running",
        "completed_stages"=>completed,"output_detail"=>String(prepared.output_detail),
        "write_results"=>write_results,"output_diagnostics"=>prepared.output_diagnostics,"timings_seconds"=>timings,
        "execution"=>deepcopy(calculation["execution"]),
        "scope"=>"Execution and stage outcomes; stage completion does not imply every event converged.",
        "absolute"=>_case_stage_summary(nothing),"relative"=>_case_stage_summary(nothing)))
    # Preparation has checked input/output collisions before any file is changed.
    write_results && mkpath(prepared.output)
    function save_summary()
        write_results || return nothing
        open(joinpath(prepared.output,"run_summary.toml"),"w") do io
            TOML.print(io,_result_metadata(summary))
        end
    end
    save_summary()
    started=time_ns()
    try
        result=_execute_case(prepared;progress,write_results,calculation,timings,completed_stages=completed,stage_summaries)
        summary["status"]="completed"
        summary["absolute"]=_case_stage_summary(result.absolute)
        summary["relative"]=_case_stage_summary(result.relative)
        timings["total"]=(time_ns()-started)/1e9
        save_summary()
        merge(result,(;calculation,run_summary=summary))
    catch exception
        summary["status"]="failed"
        summary["error"]=sprint(showerror,exception)
        merge!(summary,stage_summaries)
        timings["total"]=(time_ns()-started)/1e9
        save_summary()
        rethrow()
    end
end
