# Supply the actual equation multiplier directly. Statistical sigma remains
# metadata and is never divided out and multiplied back into empirical rows.
function _dd_scaled_equations(model,stations,data,arrivals,index,parameters,columns,step,weights;
        build_design=true,ray_cache=nothing,compact_design=false,row_indices=nothing,absolute_indices=nothing)
    _dd_equations(model,stations,data,index,parameters,columns,step;
        absolute_observations=arrivals,build_design,equation_weights=weights,ray_cache,compact_design,row_indices,absolute_indices)
end

# Full spectra remain available for the original small SVD problem. LSQR
# uses only sparse products, including the re-linearized final diagnostics.
function _dd_update_step(design,residual,solver,damp,cutoff,linear_solver,rtol,max_iterations,condition_limit)
    sparse_step = linear_solver==:lsqr ||
        linear_solver==:auto && (size(design,2)>128 || length(design)>1_000_000)
    if sparse_step
        controls = max_iterations==0 ? (;) : (;max_iterations)
        result = _dd_lsqr_step(design,residual;damp,rtol,condition_limit,controls...)
        return merge(result,(linear_solver=:lsqr,linear_iterations=result.iterations,
            design_rows=size(design,1),design_columns=size(design,2),design_nonzeros=nnz(design),
            design_storage_bytes=_dd_design_bytes(design),
            spectrum_computed=false))
    end
    _dd_dense_update_step(_dd_scaled_step(design,residual;damp,svd_relative_cutoff=cutoff),design)
end

function _dd_dense_update_step(result,design)
    merge(result,
        (linear_solver=:dense,converged=true,linear_solver_status=:direct,
         linear_iterations=0,normal_residual_norm=nothing,relative_normal_residual=nothing,
         condition_estimate=nothing,design_rows=size(design,1),design_columns=size(design,2),
         design_nonzeros=nnz(design),
         design_storage_bytes=_dd_design_bytes(design),
         spectrum_computed=true))
end

# The raw, zero-DAMP step still decides convergence. Dense positive-DAMP
# updates share its SVD; sparse LSQR retains its existing separate solves.
function _dd_update_steps(design,residual,solver,damp,cutoff,linear_solver,rtol,max_iterations,condition_limit,scheduled)
    dense_damped = damp>0 && !(linear_solver==:lsqr ||
        linear_solver==:auto && (size(design,2)>128 || length(design)>1_000_000))
    if dense_damped
        prepared = _dd_scaled_svd(design,residual;damp,svd_relative_cutoff=cutoff)
        diagnostics = _dd_dense_update_step(_dd_scaled_step_from_svd(prepared,damp,cutoff),design)
        raw = _dd_dense_update_step(_dd_scaled_step_from_svd(prepared,0.,cutoff),design)
        return diagnostics,raw,false
    end
    diagnostics = _dd_update_step(design,residual,solver,damp,cutoff,linear_solver,rtol,max_iterations,condition_limit)
    diagnostics.converged || return diagnostics,nothing,false
    raw_not_requested = scheduled && diagnostics.linear_solver==:lsqr && damp>0
    raw = raw_not_requested ? nothing : damp>0 ?
        _dd_update_step(design,residual,:damped,0.,cutoff,linear_solver,rtol,max_iterations,condition_limit) : diagnostics
    diagnostics,raw,raw_not_requested
end

# The original scaled algorithms update a complete XYZ/O vector. Neither a
# separately optimized clock nor an automatic hard reference belongs here.
function _relocate_scaled_component(model,stations,events,observations,fixed,bounds,
        max_iterations,location_tolerance,time_tolerance,difference_step;
        absolute_observations,solver,row_weighting,weights,damp,svd_relative_cutoff,
        airquake_policy,airquake_depth_km,linear_solver=:auto,lsqr_rtol=1e-6,
        lsqr_max_iterations=0,lsqr_condition_limit=1e12,
        weight_schedule=nothing,default_quality=nothing,min_weight=1e-5,_history::Symbol=:full)
    # Full relocation reruns with perturbed observations need the complete
    # sequence of selection, exit and solver decisions, not the
    # display corrections and per-reading weight report. The ordinary solve
    # keeps its full history; both modes feed the same path-signature code.
    full_history = _history==:full
    n,nd = length(events),length(observations)
    ids = String[event.event_id for event in events]
    index = Dict(id=>i for (i,id) in enumerate(ids))
    ray_indices = _dd_ray_indices(observations,absolute_observations,index)
    parameters::Matrix{Float64} = reduce(vcat,[reshape(Float64[event.source_km...,event.origin_s],1,4) for event in events])
    active = trues(n)
    reasons = Union{Nothing,Symbol}[nothing for _ in 1:n]
    scheduled = weight_schedule!==nothing
    schedule = scheduled ? _dd_weight_schedule(weight_schedule) : nothing
    schedule_updates = scheduled ? last(schedule).end_update : nothing
    weights = Float64.(weights)
    row_reasons = Union{Nothing,Symbol}[scheduled && i<=nd ? nothing : iszero(w) ? :zero_row_weight : nothing
        for (i,w) in enumerate(weights)]
    scheduled && (weights[1:nd] .= 0.)
    last_equation_weights = copy(weights)
    last_equation_rows = Int[]
    last_stage_index,last_stage_damp,last_applied_damp = nothing,nothing,Float64(damp)
    weight_metadata = fill((0,0),nd)
    weight_trace = _DDWeightHistory(observations,absolute_observations,full_history,default_quality)
    quality_metadata=(classes=Symbol[row.data_class for row in weight_trace.observations],
        quality=Float64[row.quality for row in weight_trace.observations],
        quality_source=Symbol[row.quality_source for row in weight_trace.observations])
    ever_used = falses(length(weights))
    trace = NamedTuple[]
    function save_trace!(row)
        push!(trace,full_history ? row : (;row.accepted,row.reason,row.iteration,row.stage_index,
            row.applied_damp,row.retained_rank,row.solver_rank,row.step_fraction,
            row.airquake_event_ids,row.removed_or_disconnected_ids,row.removed_by_weighting_ids,
            row.active_station_ids,row.parameters_after))
    end
    attempts,successful,airquake_attempts,secant_count = 0,0,0,0
    no_observations = scheduled ? nd==0 && all(iszero,weights) : isempty(weights) || all(iszero,weights)
    status = no_observations ? :no_observations : :max_iterations
    stop_reason = no_observations ? :no_positive_weight_observations : :iteration_limit
    diagnostics = nothing
    failed_linear_solve = nothing
    parameter_keys = Tuple{String,Symbol}[]

    # Zero-weight observations remain in the output but not the equations.
    function current_rows()
        local dr,ar
        dr = findall(i->row_reasons[i]===nothing,1:nd)
        ar = findall(i->row_reasons[nd+i]===nothing,1:length(absolute_observations))
        dr,ar
    end
    function columns_for_active()
        local columns
        columns = zeros(Int,n,4)
        keys = Tuple{String,Symbol}[]
        for i in 1:n,axis in 1:4
            if active[i] && !(ids[i] in fixed)
                push!(keys,(ids[i],(:x_km,:y_km,:z_km,:origin_s)[axis]))
                columns[i,axis] = length(keys)
            end
        end
        columns,keys
    end
    function expand_step(vector,columns)
        full = zeros(n,4)
        for i in 1:n,axis in 1:4
            column = columns[i,axis]
            column>0 && (full[i,axis]=vector[column])
        end
        full
    end
    function active_station_ids()
        local dr,ar
        dr,ar = current_rows()
        present=falses(length(ray_indices.station_names))
        for j in dr
            present[ray_indices.differential[j][5]]=true
        end
        for j in ar
            present[ray_indices.absolute[j][3]]=true
        end
        sort!(ray_indices.station_names[present])
    end
    function remove_unconnected!(reason)
        local dr,ar
        dr,ar = current_rows()
        remaining = falses(n)
        for j in dr
            i,k=ray_indices.differential[j][1:2]
            remaining[i]=true
            remaining[k]=true
        end
        for j in ar
            remaining[ray_indices.absolute[j][1]]=true
        end
        removed = findall(i->active[i] && !remaining[i],1:n)
        for i in removed
            active[i] = false
            reasons[i] = reason
        end
        removed
    end
    function save_unapplied_attempt(reason,stage,applied_damp;
            base=nothing,removed_ids=String[],linear_diagnostics=nothing,row_weights=Float64[])
        metadata = base===nothing ? (;attempt=attempts,iteration=successful+1,event_ids=copy(ids),solver,
            full_raw_step=nothing,full_step=nothing,weighted_sum_squares=nothing,
            row_weights=Float64[],column_scales=Float64[],retained_rank=nothing,solver_rank=nothing,
            condition_number=nothing,stage_index=scheduled ? stage.stage_index : nothing,stage,
            stage_damp=scheduled ? stage.damp : nothing,applied_damp,
            removed_by_weighting_ids=removed_ids,active_station_ids=active_station_ids()) : base
        if linear_diagnostics!==nothing
            d = linear_diagnostics
            metadata = merge(metadata,(;row_weights=full_history ? copy(row_weights) : nothing,
                column_scales=full_history ? copy(d.column_scales) : nothing,
                retained_rank=d.retained_rank,solver_rank=d.numerical_rank,condition_number=d.condition_number,
                linear_solver=d.linear_solver,linear_solver_status=d.linear_solver_status,
                linear_iterations=d.linear_iterations,normal_residual_norm=d.normal_residual_norm,
                relative_normal_residual=d.relative_normal_residual,condition_estimate=d.condition_estimate,
                undamped_linear_solver_status=:not_requested,undamped_linear_converged=false))
        end
        save_trace!(merge(metadata,(;successful_iterations=successful,accepted=false,reason,
            airquake_event_ids=String[],removed_or_disconnected_ids=String[],
            parameters_after=copy(parameters),actual_step=full_history ? zeros(n,4) : nothing,step_fraction=0.,
            position_step_km=0.,origin_step_s=0.)))
    end

    ray_cache = _dd_ray_cache(ray_indices.ray_count)
    raw_residual_cache=nothing
    while !no_observations && successful<max_iterations
        stage = scheduled ? _dd_weight_stage(schedule,successful) : nothing
        if scheduled && stage===nothing
            status,stop_reason = :completed_schedule,:schedule_completed
            break
        end
        attempts += 1
        applied_damp = scheduled ? stage.damp : Float64(damp)
        removed_by_weighting = Int[]
        # An accepted trial is the next iteration's identical source state.
        # Its forward rays can be reused while weights are recalculated.
        if scheduled
            dr,ar = current_rows()
            raw_residual = try
                if raw_residual_cache===nothing
                    first(_dd_scaled_equations(model,stations,view(observations,dr),NamedTuple[],index,parameters,
                        zeros(Int,n,4),difference_step,ones(length(dr));build_design=false,ray_cache,
                        row_indices=view(ray_indices.differential,dr)))
                else
                    raw_residual_cache
                end
            catch exception
                exception isa DomainError || rethrow()
                status,stop_reason = :coincident_source_station,:invalid_current_geometry
                save_unapplied_attempt(stop_reason,stage,applied_damp)
                break
            end
            raw_residual_cache=nothing
            separation=Float64[sqrt(sum((parameters[row[1],axis]-parameters[row[2],axis])^2 for axis in 1:3))
                for row in view(ray_indices.differential,dr)]
            selected_quality=(classes=view(quality_metadata.classes,dr),quality=view(quality_metadata.quality,dr),
                quality_source=view(quality_metadata.quality_source,dr))
            row_weights = _dd_empirical_weights(view(observations,dr),nothing,raw_residual,stage;
                min_weight,default_quality,quality_metadata=selected_quality,separation_km=separation)
            _dd_save_weight_attempt!(weight_trace,attempts,successful,dr,ar,row_weights,raw_residual,weights[nd.+ar],stage)
            effective_weights=row_weights.effective_weights::Vector{Float64}
            removal_reasons=row_weights.reasons::Vector{Union{Nothing,Symbol}}
            for (k,j) in enumerate(dr)
                weights[j] = effective_weights[k]
                row_reasons[j] = removal_reasons[k]
                weight_metadata[j] = (length(weight_trace.attempts),k)
            end
            removed_by_weighting = remove_unconnected!(:disconnected_after_weighting)
        end
        dr,ar = current_rows()
        if isempty(dr) && isempty(ar)
            status,stop_reason = :no_observations,:no_remaining_observations
            save_unapplied_attempt(stop_reason,stage,applied_damp;removed_ids=ids[removed_by_weighting])
            break
        end
        columns,parameter_keys = columns_for_active()
        if isempty(parameter_keys)
            status,stop_reason = :all_events_fixed,:all_events_fixed
            save_unapplied_attempt(stop_reason,stage,applied_damp;removed_ids=ids[removed_by_weighting])
            break
        end
        selected_weights = weights[[dr;nd.+ar]]
        equations = try
            _dd_scaled_equations(model,stations,view(observations,dr),view(absolute_observations,ar),index,
                parameters,columns,difference_step,selected_weights;ray_cache,compact_design=linear_solver!=:dense,
                row_indices=view(ray_indices.differential,dr),absolute_indices=view(ray_indices.absolute,ar))
        catch exception
            exception isa DomainError || rethrow()
            status,stop_reason = :coincident_source_station,:invalid_current_geometry
            save_unapplied_attempt(stop_reason,stage,applied_damp;removed_ids=ids[removed_by_weighting])
            break
        end
        residual,design,secants = equations
        secant_count += secants
        diagnostics,raw_diagnostics,raw_not_requested = _dd_update_steps(design,residual,solver,applied_damp,svd_relative_cutoff,
            linear_solver,lsqr_rtol,lsqr_max_iterations,lsqr_condition_limit,scheduled)
        if !diagnostics.converged
            status,stop_reason = :linear_solver_failed,diagnostics.linear_solver_status
            failed_linear_solve = (;diagnostics,weights=copy(weights),damp=applied_damp)
            save_unapplied_attempt(stop_reason,stage,applied_damp;removed_ids=ids[removed_by_weighting],
                linear_diagnostics=diagnostics,row_weights=selected_weights)
            break
        end
        last_equation_weights = copy(weights)
        last_equation_rows = [dr;nd.+ar]
        last_stage_index = scheduled ? stage.stage_index : nothing
        last_stage_damp = scheduled ? stage.damp : nothing
        last_applied_damp = applied_damp
        full_raw = raw_diagnostics===nothing ? nothing : expand_step(raw_diagnostics.step,columns)
        full_step = expand_step(diagnostics.step,columns)
        iteration = successful+1
        position_step = full_raw===nothing ? Inf : maximum(norm(full_raw[i,1:3]) for i in 1:n)
        origin_step = full_raw===nothing ? Inf : maximum(abs,full_raw[:,4])
        base = (;attempt=attempts,iteration,event_ids=full_history ? copy(ids) : ids,solver,
            full_raw_step=full_raw,full_step,weighted_sum_squares=sum(abs2,residual),
            row_weights=full_history ? copy(selected_weights) : nothing,
            column_scales=full_history ? copy(diagnostics.column_scales) : nothing,
            retained_rank=diagnostics.retained_rank,solver_rank=diagnostics.numerical_rank,
            condition_number=diagnostics.condition_number,
            linear_solver=diagnostics.linear_solver,linear_solver_status=diagnostics.linear_solver_status,
            linear_iterations=diagnostics.linear_iterations,normal_residual_norm=diagnostics.normal_residual_norm,
            relative_normal_residual=diagnostics.relative_normal_residual,
            condition_estimate=diagnostics.condition_estimate,
            undamped_linear_solver_status=raw_not_requested ? :not_requested_for_schedule : raw_diagnostics.linear_solver_status,
            undamped_linear_converged=raw_not_requested ? false : raw_diagnostics.converged,
            stage_index=scheduled ? stage.stage_index : nothing,stage,
            stage_damp=scheduled ? stage.damp : nothing,applied_damp,
            removed_by_weighting_ids=ids[removed_by_weighting],active_station_ids=active_station_ids())
        proposed = parameters+full_step
        bad = airquake_policy==:exclude_repeat ?
            findall(i->active[i] && proposed[i,3]<airquake_depth_km,1:n) : Int[]
        if !isempty(bad)
            # Reject EVERY update in this attempt before writing any XYZ/O.
            airquake_attempts += 1
            bad_mask = falses(n)
            bad_mask[bad].=true
            for i in bad
                active[i]=false
                reasons[i]=:airquake_removed
            end
            for j in eachindex(observations)
                i,k=ray_indices.differential[j][1:2]
                (bad_mask[i] || bad_mask[k]) &&
                    (!scheduled || row_reasons[j]===nothing) && (row_reasons[j]=:airquake_event)
            end
            for j in eachindex(absolute_observations)
                bad_mask[ray_indices.absolute[j][1]] && (!scheduled || row_reasons[nd+j]===nothing) &&
                    (row_reasons[nd+j]=:airquake_event)
            end
            disconnected = remove_unconnected!(:disconnected_after_airquake)
            if scheduled
                _dd_weight_history_exits!(weight_trace,row_reasons)
            end
            save_trace!(merge(base,(;successful_iterations=successful,accepted=false,
                reason=:airquake_rejected,airquake_event_ids=ids[bad],
                removed_or_disconnected_ids=ids[sort([bad;disconnected])],
                parameters_after=copy(parameters),actual_step=full_history ? zeros(n,4) : nothing,step_fraction=0.,
                position_step_km=0.,origin_step_s=0.)))
            if count(active)<2
                status,stop_reason = :insufficient_events_after_airquake,:fewer_than_two_events
                break
            end
            # No source moved: cached rays still describe the current geometry.
            # Retry the same successful iteration on the remaining rows.
            continue
        end
        if !scheduled && raw_diagnostics!==nothing && raw_diagnostics.converged && position_step<=location_tolerance && origin_step<=time_tolerance
            status = applied_damp>0 ? :converged_regularized :
                (diagnostics.numerical_rank===nothing || diagnostics.numerical_rank<length(parameter_keys)) ? :converged_minimum_norm : :converged
            stop_reason = :undamped_step
            ever_used[[dr;nd.+ar]] .= true
            save_trace!(merge(base,(;successful_iterations=successful,accepted=false,
                reason=:step_tolerance,airquake_event_ids=String[],removed_or_disconnected_ids=String[],
                parameters_after=copy(parameters),actual_step=full_history ? zeros(n,4) : nothing,step_fraction=0.,
                position_step_km=0.,origin_step_s=0.)))
            break
        end
        # Bounds select a feasible whole direction. The original airquake
        # policy above is evaluated first and cannot be mimicked by clipping.
        fraction = 1.
        depth_ratios = fill(Inf,n)
        bound_values = zeros(n)
        for i in 1:n
            active[i] || continue
            lower = max(bounds[1],_source_top(model,Tuple(parameters[i,1:3])))
            dz = full_step[i,3]
            if dz<0
                depth_ratios[i]=(lower-parameters[i,3])/dz
                bound_values[i]=lower
            elseif dz>0
                depth_ratios[i]=(bounds[2]-parameters[i,3])/dz
                bound_values[i]=bounds[2]
            end
            fraction=min(fraction,depth_ratios[i])
        end
        if fraction<=0
            status,stop_reason = :step_blocked_by_depth_bound,:nonzero_raw_step_at_bound
            save_unapplied_attempt(stop_reason,stage,applied_damp;base)
            break
        end
        candidate = parameters+fraction*full_step
        for i in 1:n
            fraction==depth_ratios[i] && (candidate[i,3]=bound_values[i])
        end
        trial_cache = _dd_ray_cache(ray_indices.ray_count)
        trial_seconds = try
            if all(i->!active[i] || candidate[i,3]>=_source_top(model,Tuple(candidate[i,1:3])),1:n)
                first(_dd_scaled_equations(model,stations,view(observations,dr),
                    view(absolute_observations,ar),index,candidate,columns,difference_step,ones(length(selected_weights));
                    build_design=false,ray_cache=trial_cache,
                    row_indices=view(ray_indices.differential,dr),absolute_indices=view(ray_indices.absolute,ar)))
            else
                nothing
            end
        catch exception
            exception isa DomainError || rethrow()
            nothing
        end
        if trial_seconds===nothing || !all(isfinite,trial_seconds)
            status,stop_reason = :invalid_trial_source,:trial_outside_model_domain
            save_unapplied_attempt(stop_reason,stage,applied_damp;base)
            break
        end
        actual_step = candidate-parameters
        parameters .= candidate
        ray_cache = trial_cache
        # The accepted source and row order are the next reweighting state.
        scheduled && (raw_residual_cache=trial_seconds[1:length(dr)])
        trial_residual=trial_seconds.*selected_weights
        successful += 1
        ever_used[[dr;nd.+ar]] .= true
        save_trace!(merge(base,(;successful_iterations=successful,accepted=true,reason=:updated,
            airquake_event_ids=String[],removed_or_disconnected_ids=String[],
            parameters_after=copy(parameters),actual_step,step_fraction=fraction,
            weighted_sum_squares=sum(abs2,trial_residual),
            position_step_km=maximum(norm(actual_step[i,1:3]) for i in 1:n),
            origin_step_s=maximum(abs,actual_step[:,4]))))
    end
    if scheduled && successful==schedule_updates && status==:max_iterations
        status,stop_reason = :completed_schedule,:schedule_completed
    end

    # Diagnostics are re-linearized at the saved final state, not reused from
    # the position before the last update. Removed events stay in the catalog.
    dr,ar = current_rows()
    columns,parameter_keys = columns_for_active()
    diagnostics = failed_linear_solve===nothing ? nothing : failed_linear_solve.diagnostics
    data_rank = nothing
    diagnostic_reason = isempty(parameter_keys) ? :no_free_parameters : :no_remaining_observations
    final_weights = failed_linear_solve!==nothing ? failed_linear_solve.weights :
        scheduled ? [row_reasons[i]===nothing ? last_equation_weights[i] : 0. for i in eachindex(weights)] : weights
    final_damp = failed_linear_solve!==nothing ? failed_linear_solve.damp : scheduled ? last_applied_damp : Float64(damp)
    if failed_linear_solve!==nothing
        diagnostic_reason = diagnostics.linear_solver_status
    elseif !isempty(parameter_keys) && (!isempty(dr) || !isempty(ar)) && (!scheduled || !isempty(last_equation_rows))
        try
            residual,design,secants = _dd_scaled_equations(model,stations,view(observations,dr),view(absolute_observations,ar),
                index,parameters,columns,difference_step,final_weights[[dr;nd.+ar]];compact_design=linear_solver!=:dense,
                ray_cache=_dd_ray_cache(ray_indices.ray_count),
                row_indices=view(ray_indices.differential,dr),absolute_indices=view(ray_indices.absolute,ar))
            secant_count += secants
            diagnostics = _dd_update_step(design,residual,solver,final_damp,svd_relative_cutoff,
                linear_solver,lsqr_rtol,lsqr_max_iterations,lsqr_condition_limit)
            if diagnostics.linear_solver==:dense
                _,data_rank = _absolute_svd(Matrix(design))
                diagnostic_reason = nothing
            else
                diagnostic_reason = diagnostics.converged ? :sparse_no_full_spectrum : diagnostics.linear_solver_status
            end
        catch exception
            exception isa DomainError || rethrow()
            diagnostic_reason = :coincident_source_station
        end
    end
    final_seconds = first(_dd_scaled_equations(model,stations,observations,absolute_observations,index,parameters,
        zeros(Int,n,4),difference_step,ones(length(weights));build_design=false))
    sigmas=[Float64[row.sigma_s for row in observations];Float64[row.sigma_s for row in absolute_observations]]
    final_normalized = final_seconds./sigmas
    saved_equation_weights=last_equation_weights::Vector{Float64}
    used = [row_reasons[i]===nothing && ever_used[i] for i in eachindex(weights)]
    unused_reason(i) = row_reasons[i]!==nothing ? row_reasons[i] : used[i] ? nothing :
        max_iterations==0 ? :not_iterated : status
    residuals = [merge(row,scheduled ? _dd_weight_history_metadata(weight_trace,weight_metadata[i]) : (;),
        (dt_pred_s=row.dt_s-final_seconds[i],
        residual_s=final_seconds[i],normalized_residual=final_normalized[i],
        equation_weight=final_weights[i],weighted_residual=final_weights[i]*final_seconds[i],
        last_used_equation_weight=saved_equation_weights[i],
        used=used[i],unused_reason=unused_reason(i))) for (i,row) in enumerate(observations)]
    absolute_residuals = [merge(row,(arrival_pred_s=row.arrival_s-final_seconds[nd+i],
        residual_s=final_seconds[nd+i],normalized_residual=final_normalized[nd+i],
        equation_weight=final_weights[nd+i],weighted_residual=final_weights[nd+i]*final_seconds[nd+i],
        used=used[nd+i],unused_reason=unused_reason(nd+i))) for (i,row) in enumerate(absolute_observations)]
    involved_rows = [Int[] for _ in ids]
    for (j,row) in enumerate(observations)
        push!(involved_rows[index[row.event_i]],j)
        push!(involved_rows[index[row.event_j]],j)
    end
    for (j,row) in enumerate(absolute_observations)
        push!(involved_rows[index[row.event_id]],nd+j)
    end
    results = NamedTuple[]
    for i in 1:n
        involved = involved_rows[i]
        event_status = reasons[i]!==nothing ? reasons[i] : isempty(involved) ? :no_observations :
            ids[i] in fixed ? :reference_fixed : status
        push!(results,(event_id=ids[i],source_km=Tuple(parameters[i,1:3]),origin_s=parameters[i,4],
            status=event_status,input_count=length(involved),used_count=count(used[involved]),
            solver,iterations=successful,depth_at_bound=parameters[i,3]==max(bounds[1],_source_top(model,Tuple(parameters[i,1:3]))) ||
                parameters[i,3]==bounds[2],exit_reason=reasons[i]))
    end
    objective = sum((final_weights[i]*final_seconds[i])^2
        for i in [dr;nd.+ar];init=0.)
    summary = (;status,stop_reason,event_count=n,active_event_count=count(active),
        input_count=length(weights),differential_count=nd,absolute_count=length(absolute_observations),
        used_count=count(used),iterations=successful,successful_iterations=successful,attempts,airquake_attempts,
        free_parameter_count=length(parameter_keys),jacobian_rank=data_rank,qr_diagonal_ratio=nothing,
        weighted_sum_squares=objective,secant_gradient_count=secant_count,trace,solver,row_weighting,
        damp=final_damp,svd_relative_cutoff,airquake_policy,airquake_depth_km,original_membership_retained=true,
        scheduled,weight_schedule=schedule,schedule_updates,weight_trace,
        final_stage_index=last_stage_index,final_stage_damp=last_stage_damp,final_applied_damp=final_damp,
        final_equation_weights=final_weights,last_equation_rows,active_station_ids=active_station_ids(),
        exited_reference_ids=String[id for id in fixed if !active[index[id]]],
        solver_rank=diagnostics===nothing ? nothing : diagnostics.numerical_rank,
        retained_rank=diagnostics===nothing ? nothing : diagnostics.retained_rank,
        condition_number=diagnostics===nothing ? nothing : diagnostics.condition_number,
        regularized_condition_number=diagnostics===nothing ? nothing : diagnostics.regularized_condition_number,
        retained_condition_number=diagnostics===nothing ? nothing : diagnostics.retained_condition_number,
        linear_solver=diagnostics===nothing ? linear_solver : diagnostics.linear_solver,
        linear_solver_status=diagnostics===nothing ? nothing : diagnostics.linear_solver_status,
        linear_iterations=diagnostics===nothing ? nothing : diagnostics.linear_iterations,
        normal_residual_norm=diagnostics===nothing ? nothing : diagnostics.normal_residual_norm,
        relative_normal_residual=diagnostics===nothing ? nothing : diagnostics.relative_normal_residual,
        condition_estimate=diagnostics===nothing ? nothing : diagnostics.condition_estimate,
        condition_estimate_operator=diagnostics===nothing ? nothing :
            get(diagnostics,:condition_estimate_operator,nothing),
        design_rows=diagnostics===nothing ? nothing : diagnostics.design_rows,
        design_columns=diagnostics===nothing ? nothing : diagnostics.design_columns,
        design_nonzeros=diagnostics===nothing ? nothing : diagnostics.design_nonzeros,
        design_storage_bytes=diagnostics===nothing ? nothing : diagnostics.design_storage_bytes,
        spectrum_computed=diagnostics===nothing ? false : diagnostics.spectrum_computed,
        lsqr_rtol,lsqr_max_iterations,lsqr_condition_limit,
        condition_operator=:column_scaled_data,
        diagnostics,diagnostic_reason,parameter_keys,diagnostic_source_parameters=copy(parameters))
    (;events=results,residuals,absolute_residuals,summary)
end
