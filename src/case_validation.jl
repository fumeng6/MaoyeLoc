# Static case controls only: do not trace rays, select observations or solve.
# Reuse the numerical entry points' constraints rather than inventing defaults.
function _case_validate_options(config,task,absolute_options,relative_controls,events,frame,stations)
    uncertainty=get(config,"uncertainty",Dict())
    _case_keys(uncertainty,("assume_independent_picks","ellipsoid_probability","formal",
        "relative_pairs","relative_events","relative_confidence_probability",
        "assume_independent_relative_observations","relative_covariance_method","relative_replay_step",
        "relative_algorithm_method","relative_algorithm_noise_samples","relative_algorithm_seed"),"uncertainty")
    diagnostics=get(config,"diagnostics",Dict())
    _case_keys(diagnostics,("absolute","quality_covariance_source","station_site_ids"),"diagnostics")
    for (table,key) in ((uncertainty,"assume_independent_picks"),
            (uncertainty,"assume_independent_relative_observations"),(diagnostics,"absolute"),
            (get(config,"relative",Dict()),"use_absolute_picks"))
        !haskey(table,key) || table[key] isa Bool || throw(ArgumentError("$key must be true or false."))
    end
    get(diagnostics,"quality_covariance_source","formal") in ("sigma","formal") ||
        throw(ArgumentError("Quality covariance source must be sigma or formal."))
    for key in ("ellipsoid_probability","relative_confidence_probability")
        probability=get(uncertainty,key,.95)
        probability isa Real && isfinite(probability) && 0<probability<1 ||
            throw(ArgumentError("uncertainty.$key must lie strictly between zero and one."))
    end
    formal=get(uncertainty,"formal",nothing)
    if formal!==nothing
        _case_keys(formal,("rderr_s","ercof","noise_model"),"uncertainty.formal")
        all(value->value isa Real && isfinite(value) && value>=0,
            (get(formal,"rderr_s",.15),get(formal,"ercof",1.))) ||
            throw(ArgumentError("RDERR seconds and ERCOF must be finite and nonnegative."))
        get(formal,"noise_model","weighted_rows") in ("weighted_rows","arrival_rows") ||
            throw(ArgumentError("Noise model must be weighted_rows or arrival_rows."))
    end
    get(uncertainty,"relative_covariance_method","analytic") in ("analytic","algorithm_replay") ||
        throw(ArgumentError("Relative covariance method must be analytic or algorithm_replay."))
    get(uncertainty,"relative_algorithm_method","auto") in ("auto","exact","randomized") ||
        throw(ArgumentError("Replay method must be auto, exact or randomized."))
    step=get(uncertainty,"relative_replay_step",.02)
    step isa Real && isfinite(step) && step>0 || throw(ArgumentError("Replay step must be finite positive."))
    for (key,minimum) in (("relative_algorithm_noise_samples",2),("relative_algorithm_seed",0))
        value=get(uncertainty,key,key=="relative_algorithm_seed" ? 0 : 128)
        value isa Integer && value>=minimum || throw(ArgumentError("uncertainty.$key must be an integer >= $minimum."))
    end
    ids=Set(row.event_id for row in events)
    for pair in get(uncertainty,"relative_pairs",[])
        pair isa AbstractVector && length(pair)==2 && pair[1]!=pair[2] && all(id in ids for id in pair) ||
            throw(ArgumentError("Each uncertainty pair needs two distinct catalog event IDs."))
    end
    requested=get(uncertainty,"relative_events",String[])
    requested=="all" || requested isa AbstractVector && all(id in ids for id in requested) ||
        throw(ArgumentError("relative_events must be a catalog ID array or all."))
    weighting=get(config,"absolute_weighting",nothing)
    if weighting!==nothing
        _case_keys(weighting,("mode","quality_weights","station_time_weights","s_weight","distance",
            "initial_distance","residual","event_rules","default_quality"),"absolute_weighting")
        for key in ("distance","initial_distance","residual","event_rules")
            haskey(weighting,key) || continue
            allowed=key=="event_rules" ? ("minimum_readings","insufficient_readings") :
                key=="residual" ? ("start_iteration","floor_s","inner_factor","outer_factor","passes") :
                ("start_iteration","minimum_radius_km","inner_factor","outer_factor")
            _case_keys(weighting[key],allowed,"absolute_weighting.$key")
        end
        quality=get(weighting,"default_quality",nothing)
        quality===nothing || (quality isa Integer && !(quality isa Bool) && 0<=quality<=9 &&
            get(weighting,"mode",nothing)=="empirical") ||
            throw(ArgumentError("absolute_weighting.default_quality must be an integer from 0 through 9 and requires mode=empirical."))
    end
    if task!="relative"
        o=absolute_options
        defaults=(eigtol=.012,singular_value_shift=.006,dxfix_km=7.,dzmax_km=30.,dzair=.5,
            damp=.9,late_damping_factor=.5,dxmax_km=50.,rback_s=.02,bacfac=.6,d2far_km=250.,
            stopping_rule=:legacy,dquit_km=.04,drqt_s=.001,stop_start_iteration=1)
        _validate_absolute_controls(;(key=>get(o,key,value) for (key,value) in pairs(defaults))...)
        get(o,:stop_start_iteration,1) isa Integer || throw(ArgumentError("stop_start_iteration must be an integer."))
        limit=get(o,:legacy_iteration_limit,nothing)
        limit===nothing || limit isa Integer && limit>=0 && get(o,:stopping_rule,:step)==:legacy ||
            throw(ArgumentError("legacy_iteration_limit is a nonnegative ITRLIM and requires stopping_rule=legacy."))
        if weighting!==nothing
            _validate_absolute_weight_factors(o.quality_weights,o.station_time_weights,o.s_weight)
            all(haskey(stations,id) for id in keys(o.station_time_weights)) ||
                throw(ArgumentError("Unknown station in station_time_weights."))
            _absolute_distance_options(o.distance_weighting)
            _absolute_distance_options(o.initial_distance_weighting;initial=true)
            _absolute_residual_options(o.residual_weighting)
            _absolute_event_rules(o.event_rules)
            o.row_weighting==:empirical ||
                throw(ArgumentError("Absolute location uses empirical quality, station, phase, distance and residual weights."))
        end
    end
    if task!="absolute"
        _validate_relative_controls(;relative_controls...)
        row_weights=get(relative_controls,:row_weights,nothing)
        row_weighting=get(relative_controls,:row_weighting,:empirical)
        row_weighting==:empirical && row_weights===nothing ||
            throw(ArgumentError("Relative location uses empirical class, distance and residual weights."))
        get(relative_controls,:lsqr_max_iterations,0) isa Integer || throw(ArgumentError("lsqr_max_iterations must be an integer."))
        absolute_options.depth_bounds[1]==absolute_options.depth_bounds[2] &&
            throw(ArgumentError("This relative path estimates four parameters per event and does not fix depth by equal bounds."))
        dd=get(get(config,"relative",Dict()),"weighting",nothing)
        if dd!==nothing
            staged=haskey(dd,"stages")
            _case_keys(dd,staged ? ("mode","stages","default_quality","absolute_row_weights","min_weight") :
                ("mode","wtccp","wtccs","wtctp","wtcts","default_quality","absolute_row_weights"),"relative.weighting")
            stages=staged ? dd["stages"] : [merge(dd,Dict("iterations"=>1))]
            if staged
                for stage in stages
                    _case_keys(stage,("iterations","wtccp","wtccs","wtctp","wtcts","wrcc","wrct","wdcc_km","wdct_km","damp"),"relative.weighting.stages")
                end
            end
            _dd_weight_schedule(stages)
            min_weight=get(dd,"min_weight",1e-5)
            min_weight isa Real && isfinite(min_weight) && min_weight>0 || throw(ArgumentError("DD min_weight must be finite and positive."))
            quality=get(dd,"default_quality",nothing)
            quality===nothing || quality isa Real && isfinite(quality) && quality>=0 || throw(ArgumentError("Default DD quality must be finite nonnegative."))
            # Row counts and required per-reading quality depend on any ABS
            # failures excluded later. Validate their final alignment in relocate.
            absolute_weights=get(dd,"absolute_row_weights",nothing)
            absolute_weights===nothing || absolute_weights isa AbstractVector &&
                all(w->w isa Real && isfinite(w) && w>=0,absolute_weights) ||
                throw(ArgumentError("Joint absolute_row_weights must be finite nonnegative."))
        end
        selection=get(get(config,"relative",Dict()),"selection",nothing)
        if selection!==nothing
            _case_keys(selection,("idat","ipha","max_station_distance_km","event_ids","obscc","obsct","cluster_id"),"relative.selection")
            for key in ("idat","ipha")
                value=get(selection,key,3)
                value isa Integer && value in (1,2,3) || throw(ArgumentError("DD selection $key must be 1, 2 or 3."))
            end
            for key in ("obscc","obsct","cluster_id")
                value=get(selection,key,0)
                value isa Integer && value>=0 || throw(ArgumentError("DD selection $key must be a nonnegative integer."))
            end
            distance=get(selection,"max_station_distance_km",Inf)
            distance isa Real && !isnan(distance) && distance>=0 && (!isfinite(distance) || frame!==nothing) ||
                throw(ArgumentError("DD station distance must be nonnegative km; a finite limit requires frame."))
            selected=get(selection,"event_ids",String[])
            selected isa AbstractVector && length(unique(selected))==length(selected) && all(id in ids for id in selected) ||
                throw(ArgumentError("DD selected event IDs must be unique catalog IDs."))
        end
    end
    if haskey(config,"pairing")
        criteria=config["pairing"]
        _case_keys(criteria,("max_separation_km","max_separation_ratio","min_common_stations"),"pairing")
        separation,ratio,count=criteria["max_separation_km"],criteria["max_separation_ratio"],criteria["min_common_stations"]
        isfinite(separation) && separation>0 && isfinite(ratio) && 0<ratio<1 && count isa Integer && count>0 ||
            throw(ArgumentError("Require positive separation, 0 < separation ratio < 1, and positive common-station count."))
    end
    regions=get(config,"geographic_regions",Dict())
    _case_keys(regions,("net","reference_directory"),"geographic_regions")
    net=get(regions,"net",0)
    net isa Integer && net in (0,1,2,3) && (net==0 || frame!==nothing) ||
        throw(ArgumentError("NET must be 0/1/2/3; enabled geographic queries require frame."))
    nothing
end
