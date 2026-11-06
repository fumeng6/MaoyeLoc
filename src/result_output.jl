function _relative_residual_csv_row(row)
    (row.event_i, row.event_j, row.station_id, row.phase, row.source, row.dt_s, row.sigma_s,
          row.dt_pred_s, row.residual_s, row.normalized_residual, row.used,
          row.unused_reason, row.component_id,get(row,:equation_weight,nothing),
          get(row,:weighted_residual,nothing),get(row,:input_row,nothing),
          get(row,:selection_status,:selected),
          _result_csv_fields(row,Val((:quality,:quality_source,:data_class,:class_weight,:empirical_base_weight,
              :distance_weight,:residual_weight,:raw_weight,:stage_index,:last_used_equation_weight)))...)
end

function _write_output_detail(directory, detail, omitted)
    open(joinpath(directory,"output_detail.toml"),"w") do io
        TOML.print(io,Dict("detail"=>String(detail),
            "scope"=>"file_export_only", "in_memory_history_retained"=>true,
            "not_exported_files"=>collect(omitted),
            "preexisting_not_exported_files"=>[name for name in omitted if isfile(joinpath(directory,name))]))
    end
end

function _write_location_catalog(directory,events;event_errors=NamedTuple[])
    errors=Dict(row.event_id=>row for row in event_errors)
    has_eligibility=any(event->hasproperty(event,:catalog_eligible),events)
    _write_result_csv(joinpath(directory,"events.csv"),
        ("event_id","x_km","y_km","z_km","origin_s","status","input_count","used_count",
         (has_eligibility ? ("catalog_eligible","catalog_exclusion_reason") : ())...,"uncertainty_status",
         "conditional_se_x_km","conditional_se_y_km","conditional_se_z_km","conditional_se_origin_s",
         _depth_result_header(events)...),
        (begin
            estimate=get(errors,event.event_id,get(event,:uncertainty,nothing))
            deviations=estimate===nothing ? (nothing,nothing,nothing,nothing) : Tuple(estimate.standard_deviations)
            (event.event_id,event.source_km...,event.origin_s,event.status,event.input_count,event.used_count,
                (has_eligibility ? (get(event,:catalog_eligible,true),get(event,:catalog_exclusion_reason,nothing)) : ())...,
                estimate===nothing ? :not_estimated : estimate.status,deviations...,_depth_result_values(event)...)
        end for event in events))
end

"""Write locations and estimated errors. `diagnostics=true` adds detailed records."""
function write_absolute_results(directory, result; detail::Symbol=:full,diagnostics::Bool=false)
    detail in (:core,:full) || throw(ArgumentError("Output detail must be :core or :full."))
    not_exported_files=detail==:full ? String[] :
        ["weight_iterations.csv","iterations.csv","spectral_directions.csv"]
    mkpath(directory)
    if !diagnostics
        _write_location_catalog(directory,result.events)
        if any(event->get(event,:uncertainty,nothing)!==nothing && event.uncertainty.covariance!==nothing,result.events)
            _write_result_csv(joinpath(directory,"conditional_covariance.csv"),
                ("event_id","parameter_i","parameter_j","covariance"),
                ((event.event_id,event.uncertainty.parameter_names[i],event.uncertainty.parameter_names[j],event.uncertainty.covariance[i,j])
                 for event in result.events if event.uncertainty!==nothing && event.uncertainty.covariance!==nothing
                 for i in eachindex(event.uncertainty.parameter_names) for j in eachindex(event.uncertainty.parameter_names)))
        end
        if any(event->get(event,:formal_uncertainty,nothing)!==nothing,result.events)
            _write_absolute_error_results(directory,result.events)
        end
        hasproperty(result,:magnitude_preparation) && _write_magnitude_preparation(directory,result.magnitude_preparation)
        return directory
    end
    uncertainty_status(event)=event.uncertainty===nothing ? :not_estimated : event.uncertainty.status
    uncertainty_reason(event)=event.uncertainty===nothing ? :not_requested : event.uncertainty.reason
    standard_errors(event)=event.uncertainty===nothing ? (nothing,nothing,nothing,nothing) :
        Tuple(event.uncertainty.standard_deviations)
    _write_result_csv(joinpath(directory, "events.csv"),
        ("event_id", "x_km", "y_km", "z_km", "origin_s", "status", "iterations", "input_count", "used_count",
         "rms_residual_s", "jacobian_rank", "free_parameter_rank", "condition_number", "depth_at_bound",
         "secant_gradient_count", "uncertainty_status", "uncertainty_reason",
         "conditional_se_x_km","conditional_se_y_km","conditional_se_z_km","conditional_se_origin_s",
         "constraint","constraint_source","free_parameter_count","forward_evaluated_count",
         "solver_method","retained_rank","weak_direction_truncation","stop_reason","weighted_rms_s",
         "temporarily_locked","second_station_distance_km",
         "row_weighting","positive_weight_count","weighted_reading_count","weighted_s_count","weight_normalization","weights_iteration","diagnostic_reason",
         "reweighting_cancelled","cancellation_iteration","minimum_readings","catalog_eligible","catalog_exclusion_reason","legacy_iteration_limit",
         "objective_value",
         _depth_result_header(result.events)...),
        ((event.event_id, event.source_km..., event.origin_s, event.status, event.iterations,
          event.input_count, event.used_count, event.rms_residual_s, event.jacobian_rank,
          event.free_parameter_rank, event.condition_number, event.depth_at_bound,
          event.secant_gradient_count,uncertainty_status(event),uncertainty_reason(event),
          standard_errors(event)...,event.constraint,event.constraint_source,event.free_parameter_count,
          event.forward_evaluated_count,get(event,:solver_method,:empirical_full_spectral),
          get(event,:retained_rank,nothing),get(event,:weak_direction_truncation,false),
          get(event,:stop_reason,nothing),get(event,:weighted_rms_s,nothing),
          get(event,:temporarily_locked,false),get(event,:second_station_distance_km,nothing),
          (get(event,key,nothing) for key in (:row_weighting,:positive_weight_count,:weighted_reading_count,:weighted_s_count,
              :weight_normalization,:weights_iteration,:diagnostic_reason,:reweighting_cancelled,:cancellation_iteration,
              :minimum_readings,:catalog_eligible,:catalog_exclusion_reason,:legacy_iteration_limit))...,
          get(event,:objective_value,.5*get(event,:weighted_sum_squares,NaN)),
          _depth_result_values(event)...) for event in result.events))
    _write_result_csv(joinpath(directory,"conditional_covariance.csv"),
        ("event_id","parameter_i","parameter_j","covariance"),
        ((event.event_id,event.uncertainty.parameter_names[i],event.uncertainty.parameter_names[j],event.uncertainty.covariance[i,j])
         for event in result.events if event.uncertainty!==nothing && event.uncertainty.covariance!==nothing
         for i in eachindex(event.uncertainty.parameter_names) for j in eachindex(event.uncertainty.parameter_names)))
    _write_absolute_error_results(directory,result.events)
    _write_depth_stations(directory,result.events)
    _write_absolute_diagnostics(directory,result.events)
    _write_geographic_regions(directory,result.events)
    hasproperty(result,:magnitude_preparation) && _write_magnitude_preparation(directory,result.magnitude_preparation)
    if hasproperty(result,:processing_metadata)
        open(joinpath(directory,"processing_metadata.toml"),"w") do io
            TOML.print(io,_result_metadata(result.processing_metadata))
        end
    end
    _write_result_csv(joinpath(directory, "residuals.csv"),
        ("event_id", "station_id", "phase", "arrival_s", "sigma_s", "arrival_pred_s", "residual_s",
         "normalized_residual", "used", "unused_reason","quality_code","quality_code_source",
         "quality_weight","station_weight","phase_weight","base_weight","distance_weight","residual_weight","raw_weight","equation_weight","weighted_residual","observation_id","robust_weight","robust_downweighted"),
        ((row.event_id, row.station_id, row.phase, row.arrival_s, row.sigma_s, row.arrival_pred_s,
          row.residual_s, row.normalized_residual, row.used, row.unused_reason,
          (get(row,key,nothing) for key in (:quality_code,:quality_code_source,:quality_weight,:station_weight,:phase_weight,
              :base_weight,:distance_weight,:residual_weight,:raw_weight,:equation_weight,:weighted_residual))...,get(row,:observation_id,nothing),get(row,:robust_weight,nothing),get(row,:robust_downweighted,nothing)) for row in result.residuals))
    if detail==:full
        _write_result_csv(joinpath(directory,"weight_iterations.csv"),
            ("event_id","iteration","station_id","phase","quality_code","quality_weight","station_weight","phase_weight",
             "base_weight","distance_weight","residual_weight","raw_weight","normalization","equation_weight",
             "rms_references_s","positive_count","weighted_reading_count","unused_reason",
             "weight_pass","selected_for_equations","reweighting_enabled","reweighting_cancelled","observation_id"),
            ((event.event_id,row.iteration,row.station_id,row.phase,row.quality_code,row.quality_weight,row.station_weight,
              row.phase_weight,row.base_weight,row.distance_weight,row.residual_weight,row.raw_weight,row.normalization,
              row.equation_weight,join(row.rms_references_s,'|'),row.positive_count,row.weighted_reading_count,row.unused_reason,
              row.weight_pass,row.selected_for_equations,row.reweighting_enabled,row.reweighting_cancelled,get(row,:observation_id,nothing))
              for event in result.events for row in get(event,:weight_trace,NamedTuple[])))
        _write_result_csv(joinpath(directory, "iterations.csv"),
            ("event_id", "iteration", "x_km", "y_km", "z_km", "origin_s", "weighted_sum_squares", "objective_value",
             "position_step_km", "origin_step_s", "active_depth", "numerical_damping",
             "weighted_rms_s","retained_rank","row_normalization","temporarily_locked","depth_released",
             "step_fraction","backtracks","damping_factor","air_factor","vertical_limited","spatial_factor",
             "second_station_distance_km","raw_dx_km","raw_dy_km","raw_dz_km","raw_dt_s",
             "actual_dx_km","actual_dy_km","actual_dz_km","actual_dt_s","iteration_action"),
            ((event.event_id, row.iteration, row.source_km..., row.origin_s, row.weighted_sum_squares, get(row,:objective_value,nothing),
              row.position_step_km, row.origin_step_s, row.active_depth, row.damping,
              (get(row,key,nothing) for key in (:weighted_rms_s,:retained_rank,:row_normalization,
                  :temporarily_locked,:depth_released,:step_fraction,:backtracks,:damping_factor,
                  :air_factor,:vertical_limited,:spatial_factor,:second_station_distance_km))...,
              get(row,:spectral_raw_step,(nothing,nothing,nothing,nothing))...,
              get(row,:actual_step,(nothing,nothing,nothing,nothing))...,get(row,:iteration_action,:solve))
              for event in result.events for row in event.trace))
        _write_result_csv(joinpath(directory,"spectral_directions.csv"),
            ("event_id","direction","normalized_singular_value","retained","x_component","y_component",
             "z_component","origin_component","row_normalization","eigtol","singular_value_shift"),
            ((event.event_id,j,event.spectral.singular_values[j],event.spectral.retained[j],
              event.spectral.vectors4[:,j]...,event.spectral.row_normalization,event.spectral.eigtol,
              event.spectral.singular_value_shift) for event in result.events if hasproperty(event,:spectral)
              for j in eachindex(event.spectral.singular_values)))
    end
    _write_output_detail(directory,detail,not_exported_files)
    return directory
end

# Candidate objectives use one declared loss and the original complete readings.
function _write_relative_statistics(directory,result)
    s=relative_statistics(result)
    fields=(:input_count,:used_count,:unused_count,:cc_p_count,:cc_s_count,:ct_p_count,:ct_s_count,:cc_count,:ct_count,
        :cc_rms_s,:ct_rms_s,:absolute_input_count,:absolute_used_count,:absolute_unused_count,:absolute_p_count,:absolute_s_count,:absolute_rms_s)
    maxima=(:max_station_cc_rms_s,:max_station_cc_station_id,:max_station_cc_component_id,
        :max_station_ct_rms_s,:max_station_ct_station_id,:max_station_ct_component_id)
    _write_result_csv(joinpath(directory,"event_statistics.csv"),
        ("event_id","component_id","status","selection_status","exit_reason",String.(fields)...),
        ((r.event_id,r.component_id,r.status,r.selection_status,r.exit_reason,(getproperty(r,k) for k in fields)...) for r in s.events))
    _write_result_csv(joinpath(directory,"station_statistics.csv"),("component_id","station_id",String.(fields)...),
        ((r.component_id,r.station_id,(getproperty(r,k) for k in fields)...) for r in s.stations))
    _write_result_csv(joinpath(directory,"component_statistics.csv"),("component_id","status","event_count",String.(fields)...,String.(maxima)...),
        ((r.component_id,r.status,r.event_count,(getproperty(r,k) for k in (fields...,maxima...))...) for r in s.components))
    _write_result_csv(joinpath(directory,"unused_statistics.csv"),
        ("scope","component_id","event_id","station_id","source","phase","reason","count"),(Tuple(r) for r in s.unused_reasons))
    open(joinpath(directory,"statistics.toml"),"w") do io
        TOML.print(io,_result_metadata((summary=s.summary,assumptions=s.assumptions)))
    end
    s
end

function _write_relative_observation_diagnostics(directory,diagnostics)
    fields=(:input_row,:input_row_provenance,:event_i,:event_j,:station_id,:phase,:source,:idx,:original_quality,:quality_source,
        :dt_s,:dt_ms,:sigma_s,:sigma_ms,:dt_pred_s,:dt_pred_ms,:residual_s,:residual_ms,:normalized_residual,
        :equation_weight,:weighted_residual,:equation_weight_unit,:weighted_residual_unit,:row_weight_definition,
        :used,:unused_reason,:component_id,:selection_status,:last_used_equation_weight,:stage_index,
        :separation_km,:separation_m,:distance_i_km,:distance_i_m,:azimuth_i_deg,:source_angle_i_deg,:angle_i_definition,:angle_i_reason,
        :distance_j_km,:distance_j_m,:azimuth_j_deg,:source_angle_j_deg,:angle_j_definition,:angle_j_reason)
    differential_fields=Val(fields)
    _write_result_csv(joinpath(directory,"observation_diagnostics.csv"),String.(fields),
        diagnostics.rows;fields=differential_fields)
    fields=(:input_row,:input_row_provenance,:event_id,:station_id,:phase,:arrival_s,:arrival_ms,:arrival_pred_s,:arrival_pred_ms,
        :sigma_s,:sigma_ms,:residual_s,:residual_ms,:normalized_residual,:equation_weight,:weighted_residual,:equation_weight_unit,
        :weighted_residual_unit,:row_weight_definition,:used,:unused_reason,:component_id,:selection_status,
        :distance_km,:distance_m,:azimuth_deg,:source_angle_deg,:angle_definition,:angle_reason,:observation_id)
    absolute_fields=Val(fields)
    _write_result_csv(joinpath(directory,"absolute_observation_diagnostics.csv"),String.(fields),
        diagnostics.absolute_rows;fields=absolute_fields)
    components=Tuple[]
    for r in diagnostics.rows,(endpoint,id) in ((:i,r.event_i),(:j,r.event_j))
        for c in getproperty(r,endpoint==:i ? :source_angle_i_components : :source_angle_j_components)
            push!(components,(r.input_row,endpoint,id,r.station_id,r.phase,c.model_id,c.weight,c.source_angle_deg,c.angle_definition,c.angle_reason,nothing))
        end
    end
    for r in diagnostics.absolute_rows,c in r.source_angle_components
        push!(components,(r.input_row,:absolute,r.event_id,r.station_id,r.phase,c.model_id,c.weight,c.source_angle_deg,c.angle_definition,c.angle_reason,get(r,:observation_id,nothing)))
    end
    _write_result_csv(joinpath(directory,"observation_ray_components.csv"),
        ("input_row","endpoint","event_id","station_id","phase","model_id","weight","source_angle_deg","angle_definition","angle_reason","observation_id"),components)
end

function _write_relative_error_results(directory,result;diagnostics=false)
    result.uncertainty===nothing && !diagnostics && return
    event_errors=result.uncertainty===nothing ? NamedTuple[] : get(result.uncertainty,:events,NamedTuple[])
    pairs=result.uncertainty===nothing ? NamedTuple[] : result.uncertainty.pairs
    if diagnostics || !isempty(pairs)
        _write_result_csv(joinpath(directory,"pair_uncertainty.csv"),
            ("event_i","event_j","component_id","status","reason","conditional_se_dx_km","conditional_se_dy_km",
             "conditional_se_dz_km","conditional_se_dorigin_s","jacobian_rank"),
            ((row.event_i,row.event_j,row.component_id,row.status,row.reason,row.standard_deviations...,
              row.rank) for row in pairs))
        _write_result_csv(joinpath(directory,"conditional_pair_covariance.csv"),
            ("event_i","event_j","parameter_i","parameter_j","covariance"),
            ((row.event_i,row.event_j,row.parameter_names[i],row.parameter_names[j],row.covariance[i,j])
             for row in pairs if row.covariance!==nothing
             for i in eachindex(row.parameter_names) for j in eachindex(row.parameter_names)))
    end
    component_status=Dict(c.component_id=>c.status for c in result.components)
    if diagnostics || !isempty(event_errors)
        _write_result_csv(joinpath(directory,"event_uncertainty.csv"),
            ("event_id","component_id","status","reason","conditional_se_x_km","conditional_se_y_km","conditional_se_z_km","conditional_se_origin_s",
             "confidence_probability","marginal_confidence_scale","marginal_half_width_x_km","marginal_half_width_y_km","marginal_half_width_z_km",
             "marginal_half_width_origin_s","rank","reference_event_ids","quantity","solver_status","covariance_method"),
            ((r.event_id,r.component_id,r.status,r.reason,Tuple(r.standard_deviations)...,r.confidence_probability,r.confidence_scale,
              Tuple(r.marginal_confidence_half_widths)...,r.rank,join(r.reference_event_ids,'|'),
              get(r.assumptions,:quantity,:event_conditional_coordinates),get(r,:solver_status,get(component_status,r.component_id,nothing)),
              get(r.assumptions,:covariance_method,:not_computed)) for r in event_errors))
        _write_result_csv(joinpath(directory,"conditional_event_covariance.csv"),
            ("event_id","parameter_i","parameter_j","covariance"),
            ((r.event_id,r.parameter_names[i],r.parameter_names[j],r.covariance[i,j]) for r in event_errors if r.covariance!==nothing
              for i in eachindex(r.parameter_names) for j in eachindex(r.parameter_names)))
    end
    if result.uncertainty!==nothing
        open(joinpath(directory,"uncertainty_assumptions.toml"),"w") do io
            TOML.print(io,_result_metadata((assumptions=result.uncertainty.assumptions,
                events=[(event_id=r.event_id,assumptions=r.assumptions) for r in event_errors])))
        end
        algorithm=get(result.uncertainty,:algorithm_uncertainty,nothing)
        if diagnostics && algorithm!==nothing
            # TOML uses arrays of rows for the 4-by-4 sampling-error matrix.
            sampling_record(check) = check===nothing ? nothing : merge(check,
                (covariance_element_standard_error=check.covariance_element_standard_error===nothing ? nothing :
                    [collect(row) for row in eachrow(check.covariance_element_standard_error)],))
            open(joinpath(directory,"algorithm_uncertainty_numerics.toml"),"w") do io
                TOML.print(io,_result_metadata((noise_direction_count=algorithm.noise_direction_count,
                    replay_count=algorithm.replay_count,observation_count=algorithm.observation_count,
                    numerical_check=algorithm.numerical_check,assumptions=algorithm.assumptions,
                    events=[(event_id=r.event_id,status=r.status,reason=r.reason,numerical_check=r.numerical_check,
                        sampling_check=sampling_record(get(r,:sampling_check,nothing))) for r in algorithm.events])))
            end
        end
    end
    nothing
end

"""Write DD locations and requested errors. `diagnostics=true` adds detailed records."""
function write_relative_results(directory, result; detail::Symbol=:full,diagnostics::Bool=false)
    detail in (:core,:full) || throw(ArgumentError("Output detail must be :core or :full."))
    not_exported_files=detail==:full ? String[] :
        ["weight_iterations.csv","iterations.csv","source_iterations.csv",
         "centroid_iterations.csv","correction_iterations.csv","solver_directions.csv"]
    mkpath(directory)
    if !diagnostics
        errors=result.uncertainty===nothing ? NamedTuple[] : get(result.uncertainty,:events,NamedTuple[])
        _write_location_catalog(directory,result.events;event_errors=errors)
        _write_relative_error_results(directory,result)
        return directory
    end
    statistics=_write_relative_statistics(directory,result)
    _write_depth_stations(directory,result.events)
    hasproperty(result,:observation_diagnostics) && _write_relative_observation_diagnostics(directory,result.observation_diagnostics)
    event_errors=result.uncertainty===nothing ? NamedTuple[] : get(result.uncertainty,:events,NamedTuple[])
    uncertainty_status=Dict(row.event_id=>row.status for row in event_errors)
    _write_result_csv(joinpath(directory, "events.csv"),
        ("event_id", "x_km", "y_km", "z_km", "origin_s", "status", "component_id", "iterations",
         "input_count", "used_count", "depth_at_bound", "uncertainty_status", "previous_absolute_status",
         "initialization","catalog_x_km","catalog_y_km","catalog_z_km","catalog_origin_s",
         "initial_x_km","initial_y_km","initial_z_km","initial_origin_s",
         "solver","exit_reason","selection_status",
         _depth_result_header(result.events)...),
        ((event.event_id, event.source_km..., event.origin_s, event.status, event.component_id,
          event.iterations, event.input_count, event.used_count, event.depth_at_bound,
          get(uncertainty_status,event.event_id,:not_estimated),get(event,:previous_absolute_status,nothing),get(event,:initialization,nothing),
          get(event,:catalog_source_km,(nothing,nothing,nothing))...,get(event,:catalog_origin_s,nothing),
          get(event,:initial_source_km,(nothing,nothing,nothing))...,get(event,:initial_origin_s,nothing),
          get(event,:solver,:damped),get(event,:exit_reason,nothing),get(event,:selection_status,:selected),
          _depth_result_values(event)...) for event in result.events))
    _write_result_csv(joinpath(directory, "residuals.csv"),
        ("event_i", "event_j", "station_id", "phase", "source", "dt_s", "sigma_s", "dt_pred_s",
         "residual_s", "normalized_residual", "used", "unused_reason", "component_id",
         "equation_weight","weighted_residual","input_row","selection_status",
         "quality","quality_source","data_class","class_weight","empirical_base_weight",
         "distance_weight","residual_weight","raw_weight","stage_index","last_used_equation_weight"),
        (_relative_residual_csv_row(row) for row in result.residuals))
    if detail==:full
        weight_keys=(:attempt,:iteration,:successful_updates,:input_row,:row_kind,:event_i,:event_j,:station_id,:phase,:source,
            :quality,:quality_source,:data_class,:class_weight,:empirical_base_weight,:distance_weight,:residual_weight,
            :raw_weight,:equation_weight,:stage_index,:residual_s,:separation_km,:removal_reason)
        stage_keys=(:wtccp,:wtccs,:wtctp,:wtcts,:wrcc,:wrct,:wdcc_km,:wdct_km,:damp)
        scale_keys=(:median_s,:mad_s,:scale_s,:cutoff_s,:status)
        _write_result_csv(joinpath(directory,"weight_iterations.csv"),
            ("component_id",String.(weight_keys)...,String.(stage_keys)...,
             ("cc_".*String.(scale_keys))...,("ct_".*String.(scale_keys))...),
            ((component.component_id,(getproperty(row,key) for key in weight_keys)...,
                (getproperty(row.stage,key) for key in stage_keys)...,
                (getproperty(row.residual_scale_cc,key) for key in scale_keys)...,
                (getproperty(row.residual_scale_ct,key) for key in scale_keys)...)
                for component in result.components for row in get(component,:weight_trace,NamedTuple[])))
    end
    plan=get(result,:selection,nothing)
    if plan!==nothing
        selection_metadata=Dict{String,Any}(String(k)=>v for (k,v) in Base.pairs(plan.settings))
        if plan.centers.global_center!==nothing
            selection_metadata["initial_center"]=Dict(String(k)=>v for (k,v) in Base.pairs(plan.centers.global_center))
        end
        open(joinpath(directory,"selection_settings.toml"),"w") do io
            TOML.print(io,selection_metadata)
        end
        _write_result_csv(joinpath(directory,"selection_events.csv"),
            ("event_id","cluster_id","selection_status"),
            ((id,plan.membership[i]==0 ? nothing : plan.membership[i],plan.event_statuses[i])
                for (i,id) in enumerate(plan.catalog_event_ids)))
        _write_result_csv(joinpath(directory,"selection_clusters.csv"),
            ("cluster_id","selected","event_count","reread_observation_count","center_latitude_deg","center_longitude_deg",
             "reread","distance_center_latitude_deg","distance_center_longitude_deg"),
            ((i,i in plan.selected_group_indices,length(group),length(plan.group_rows[i]),
                plan.centers.clusters[i]===nothing ? nothing : plan.centers.clusters[i].latitude_deg,
                plan.centers.clusters[i]===nothing ? nothing : plan.centers.clusters[i].longitude_deg,
                plan.reread[i],plan.centers.distance_centers[i]===nothing ? nothing : plan.centers.distance_centers[i].latitude_deg,
                plan.centers.distance_centers[i]===nothing ? nothing : plan.centers.distance_centers[i].longitude_deg)
                for (i,group) in enumerate(plan.groups)))
        _write_result_csv(joinpath(directory,"selection_pairs.csv"),
            ("event_i","event_j","cc","ct","total","strong"),
            (Tuple(row) for row in plan.pair_counts))
    end
    _write_result_csv(joinpath(directory, "components.csv"),
        ("component_id", "status", "event_count", "input_count", "used_count", "iterations",
         "reference_event_ids", "automatic_reference", "free_parameter_count", "jacobian_rank",
         "qr_diagonal_ratio", "weighted_sum_squares", "secant_gradient_count", "datum", "differential_count", "absolute_count",
         "initialization","initialization_coordinates","initial_x_km","initial_y_km","initial_z_km",
         "solver","row_weighting","stop_reason","attempts","successful_iterations","airquake_attempts",
         "active_event_count","solver_rank","retained_rank","condition_number","condition_operator",
         "regularized_condition_number","retained_condition_number",
         "damp","svd_relative_cutoff","airquake_policy","airquake_depth_km","diagnostic_reason",
         "scheduled","schedule_updates","final_stage_index","final_stage_damp","final_applied_damp","active_station_ids","exited_reference_ids",
         "linear_solver","linear_solver_status","linear_iterations","normal_residual_norm","relative_normal_residual",
         "condition_estimate","condition_estimate_operator","spectrum_computed","design_rows","design_columns",
         "design_nonzeros","design_storage_bytes","lsqr_rtol","lsqr_max_iterations","lsqr_condition_limit"),
        ((row.component_id, row.status, row.event_count, row.input_count, row.used_count, row.iterations,
          join(row.reference_event_ids, '|'), row.automatic_reference, row.free_parameter_count,
          row.jacobian_rank, row.qr_diagonal_ratio, row.weighted_sum_squares, row.secant_gradient_count,
          row.datum,row.differential_count,row.absolute_count,row.initialization,row.initialization_coordinates,
          (row.initial_source_km===nothing ? (nothing,nothing,nothing) : row.initial_source_km)...,
          get(row,:solver,:damped),get(row,:row_weighting,:empirical),get(row,:stop_reason,nothing),
          get(row,:attempts,row.iterations),get(row,:successful_iterations,row.iterations),
          get(row,:airquake_attempts,0),get(row,:active_event_count,row.event_count),
          get(row,:solver_rank,nothing),get(row,:retained_rank,nothing),get(row,:condition_number,nothing),
          get(row,:condition_operator,:column_scaled_data),get(row,:regularized_condition_number,nothing),
          get(row,:retained_condition_number,nothing),get(row,:damp,nothing),
          get(row,:svd_relative_cutoff,nothing),get(row,:airquake_policy,:bounds),get(row,:airquake_depth_km,nothing),
          get(row,:diagnostic_reason,nothing),get(row,:scheduled,false),get(row,:schedule_updates,nothing),
          get(row,:final_stage_index,nothing),get(row,:final_stage_damp,nothing),get(row,:final_applied_damp,nothing),
          join(get(row,:active_station_ids,String[]),'|'),join(get(row,:exited_reference_ids,String[]),'|'),
          get(row,:linear_solver,nothing),get(row,:linear_solver_status,nothing),get(row,:linear_iterations,nothing),
          get(row,:normal_residual_norm,nothing),get(row,:relative_normal_residual,nothing),
          get(row,:condition_estimate,nothing),get(row,:condition_estimate_operator,nothing),get(row,:spectrum_computed,nothing),
          get(row,:design_rows,nothing),get(row,:design_columns,nothing),get(row,:design_nonzeros,nothing),
          get(row,:design_storage_bytes,nothing),get(row,:lsqr_rtol,nothing),get(row,:lsqr_max_iterations,nothing),get(row,:lsqr_condition_limit,nothing))
          for row in result.components))
    _write_result_csv(joinpath(directory,"absolute_residuals.csv"),
        ("event_id","station_id","phase","arrival_s","sigma_s","arrival_pred_s","residual_s",
         "normalized_residual","used","unused_reason","component_id","equation_weight","weighted_residual","observation_id"),
        ((row.event_id,row.station_id,row.phase,row.arrival_s,row.sigma_s,row.arrival_pred_s,row.residual_s,
          row.normalized_residual,row.used,row.unused_reason,row.component_id,
          get(row,:equation_weight,inv(row.sigma_s)),get(row,:weighted_residual,row.normalized_residual),get(row,:observation_id,nothing)) for row in result.absolute_residuals))
    _write_relative_error_results(directory,result;diagnostics=true)
    history=nothing
    if detail==:full
        _write_result_csv(joinpath(directory, "iterations.csv"),
            ("component_id", "iteration", "weighted_sum_squares", "position_step_km", "origin_step_s",
             "attempt","successful_iterations","accepted","reason","airquake_event_ids","removed_or_disconnected_ids",
             "step_fraction","retained_rank","solver_rank","condition_number",
             "stage_index","stage_damp","applied_damp","removed_by_weighting_ids","active_station_ids",
             "linear_solver","linear_solver_status","linear_iterations","normal_residual_norm","relative_normal_residual",
             "condition_estimate","undamped_linear_solver_status","undamped_linear_converged"),
            ((component.component_id, row.iteration, row.weighted_sum_squares, row.position_step_km,
              row.origin_step_s,get(row,:attempt,row.iteration),get(row,:successful_iterations,row.iteration),
              get(row,:accepted,true),get(row,:reason,:updated),join(get(row,:airquake_event_ids,String[]),'|'),
              join(get(row,:removed_or_disconnected_ids,String[]),'|'),get(row,:step_fraction,nothing),
              get(row,:retained_rank,nothing),get(row,:solver_rank,nothing),get(row,:condition_number,nothing),
              get(row,:stage_index,nothing),get(row,:stage_damp,nothing),get(row,:applied_damp,nothing),
              join(get(row,:removed_by_weighting_ids,String[]),'|'),join(get(row,:active_station_ids,String[]),'|'),
              get(row,:linear_solver,nothing),get(row,:linear_solver_status,nothing),get(row,:linear_iterations,nothing),
              get(row,:normal_residual_norm,nothing),get(row,:relative_normal_residual,nothing),get(row,:condition_estimate,nothing),
              get(row,:undamped_linear_solver_status,nothing),get(row,:undamped_linear_converged,nothing))
              for component in result.components for row in component.trace))
        history=relative_iteration_diagnostics(result)
        _write_result_csv(joinpath(directory,"source_iterations.csv"),
            ("component_id","attempt","iteration","successful_iterations","accepted","reason","event_id",
             "x_km","y_km","z_km","origin_s","raw_dx_km","raw_dy_km","raw_dz_km","raw_dorigin_s",
             "full_dx_km","full_dy_km","full_dz_km","full_dorigin_s",
             "actual_dx_km","actual_dy_km","actual_dz_km","actual_dorigin_s","solver","stage_index","is_reference","depth_fixed","active","exit_reason","origin_profiled"),
            ((r.component_id,r.attempt,r.iteration,r.successful_iterations,r.accepted,r.reason,r.event_id,
              (getproperty(r,k) for k in (:x_km,:y_km,:z_km,:origin_s,:raw_dx_km,:raw_dy_km,:raw_dz_km,:raw_dorigin_s,
                :full_dx_km,:full_dy_km,:full_dz_km,:full_dorigin_s,:actual_dx_km,:actual_dy_km,:actual_dz_km,:actual_dorigin_s,
                :solver,:stage_index,:is_reference,:depth_fixed,:active,:exit_reason,:origin_profiled))...) for r in history.sources))
        centroid_fields=(:component_id,:attempt,:iteration,:successful_iterations,:solver,:stage_index,:accepted,:reason,:origin_profiled,
            :scope,:event_count,:event_ids,:initial_x_km,:initial_y_km,:initial_z_km,:initial_origin_s,
            :before_x_km,:before_y_km,:before_z_km,:before_origin_s,:current_x_km,:current_y_km,:current_z_km,:current_origin_s,
            :drift_dx_km,:drift_dy_km,:drift_dz_km,:drift_dorigin_s,:actual_dx_km,:actual_dy_km,:actual_dz_km,:actual_dorigin_s)
        _write_result_csv(joinpath(directory,"centroid_iterations.csv"),String.(centroid_fields),
            ((k==:event_ids ? join(r.event_ids,'|') : getproperty(r,k) for k in centroid_fields) for r in history.centroids))
        correction_fields=(:component_id,:attempt,:iteration,:successful_iterations,:solver,:stage_index,:accepted,:reason,:origin_profiled,
            :scope,:step_kind,:event_count,:mean_dx_km,:mean_dy_km,:mean_dz_km,:mean_dorigin_s,
            :rms_dx_km,:rms_dy_km,:rms_dz_km,:rms_dorigin_s,:max_abs_dx_km,:max_abs_dy_km,:max_abs_dz_km,:max_abs_dorigin_s,
            :mean_position_step_km,:rms_position_step_km,:max_position_step_km)
        _write_result_csv(joinpath(directory,"correction_iterations.csv"),String.(correction_fields),
            ((getproperty(r,k) for k in correction_fields) for r in history.corrections))
    end
    component_fields=(:component_id,:status,:stop_reason,:datum,:reference_event_ids,:automatic_reference,
        :solver,:row_weighting,:event_count,:active_event_count,:input_count,:used_count,:differential_count,:absolute_count,
        :iterations,:attempts,:successful_iterations,:airquake_attempts,:weighted_sum_squares,:free_parameter_count,
        :jacobian_rank,:solver_rank,:retained_rank,:qr_diagonal_ratio,:condition_number,:condition_operator,
        :regularized_condition_number,:retained_condition_number,:diagnostic_reason,:damp,:svd_relative_cutoff,
        :scheduled,:schedule_updates,:final_stage_index,:final_stage_damp,:final_applied_damp,:active_station_ids,:exited_reference_ids,
        :linear_solver,:linear_solver_status,:linear_iterations,:normal_residual_norm,:relative_normal_residual,
        :condition_estimate,:condition_estimate_operator,:spectrum_computed,:design_rows,:design_columns,
        :design_nonzeros,:design_storage_bytes,:lsqr_rtol,:lsqr_max_iterations,:lsqr_condition_limit)
    calculation=joinpath(dirname(directory),"calculation.toml")
    parameters=isfile(calculation) ? TOML.parsefile(calculation) : nothing
    summary=(;parameters,
        components=[Dict(String(k)=>get(c,k,nothing) for k in component_fields) for c in result.components],
        statistics=statistics.summary,statistics_assumptions=statistics.assumptions,
        unused_reasons=[r for r in statistics.unused_reasons if r.scope==:summary],
        final_corrections=detail==:full ? [r for r in history.corrections if r.scope==:all && r.step_kind==:actual &&
            r.attempt==maximum(x.attempt for x in history.corrections if x.component_id==r.component_id)] : nothing,
        history_files=detail==:full ?
            (sources="source_iterations.csv",centroids="centroid_iterations.csv",corrections="correction_iterations.csv",
                attempts="iterations.csv",weights="weight_iterations.csv",selections="selection_events.csv") :
            (sources=nothing,centroids=nothing,corrections=nothing,attempts=nothing,weights=nothing,
                selections=plan===nothing ? nothing : "selection_events.csv"),
        replay=(input_config=parameters===nothing ? nothing : get(parameters,"input_config",nothing),
            command="julia --startup-file=no --project=. bin/maoyeloc.jl <input_config>",
            input_files_present=parameters===nothing ? nothing : all(isfile,values(parameters["input_files"]))),
        uncertainty_scope=:conditional_on_recorded_model_weights_selection_and_datum)
    if parameters!==nothing && (get(parameters,"input_config",nothing)=="in-memory" ||
            any(get(source,"kind",nothing)=="memory" && get(source,"used",true)
                for source in values(get(parameters,"input_sources",Dict()))))
        summary=merge(summary,(replay=(input_config=parameters["input_config"],command=nothing,
            input_files_present=all(isfile,values(parameters["input_files"])),
            input_records_available_from_files=false,
            instruction="Run the original Julia script again with its configuration and in-memory inputs."),))
    end
    if detail==:core
        summary=merge(summary,(output_detail=:core,final_corrections_status=:not_exported,
            history_export_status=:not_exported,in_memory_history_retained=true,
            not_exported_files=not_exported_files))
    end
    open(joinpath(directory,"run_summary.toml"),"w") do io
        TOML.print(io,_result_metadata(summary))
    end
    _write_result_csv(joinpath(directory,"column_scales.csv"),
        ("component_id","iteration","event_id","parameter","column_rms","zero_column","operator"),
        ((component.component_id,component.iterations,component.parameter_keys[k]...,
          component.diagnostics.column_scales[k],k in component.diagnostics.zero_columns,component.condition_operator)
          for component in result.components if get(component,:diagnostics,nothing)!==nothing
          for k in eachindex(component.parameter_keys)))
    if detail==:full
        _write_result_csv(joinpath(directory,"solver_directions.csv"),
            ("component_id","iteration","direction","singular_value","relative_singular_value","retained","gain",
             "event_id","parameter","scaled_component","physical_component","operator"),
            ((component.component_id,component.iterations,j,component.diagnostics.singular_values[j],
              maximum(component.diagnostics.singular_values;init=0.)==0 ? 0. :
                  component.diagnostics.singular_values[j]/maximum(component.diagnostics.singular_values),
              component.diagnostics.retained[j],component.diagnostics.gains[j],component.parameter_keys[k]...,
              component.diagnostics.scaled_vectors[k,j],component.diagnostics.physical_vectors[k,j],component.condition_operator)
              for component in result.components if get(component,:diagnostics,nothing)!==nothing
              for j in eachindex(component.diagnostics.singular_values) for k in eachindex(component.parameter_keys)))
    end
    _write_output_detail(directory,detail,not_exported_files)
    return directory
end
