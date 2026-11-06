# EIGTOL is an absolute singular-value cutoff in the physical km/second
# columns. Only a common ROW normalization is applied, never column scaling
# or clock elimination. design and residual already include 1/sigma.
function _absolute_spectral_step(design,residual,sigma;
        axes=(1,2,3,4),eigtol=0.,singular_value_shift=.006,row_normalization=nothing)
    n = length(sigma)
    normalization = row_normalization===nothing ? (n==0 ? 1. : n/sum(inv,sigma)) : Float64(row_normalization)
    B = normalization .* design[:,collect(axes)]
    rhs = normalization .* residual
    p = length(axes)
    decomposition = svd(B;full=n<p)
    singular_values = [decomposition.S;zeros(max(0,p-length(decomposition.S)))]
    numerical_floor = max(n,p)*eps(Float64)*maximum(singular_values;init=0.)
    retained = (singular_values .> eigtol) .& (singular_values .> numerical_floor)
    step4,raw_step4 = zeros(4),zeros(4)
    for j in findall(retained)
        amplitude = dot(decomposition.U[:,j],rhs)
        direction = decomposition.V[:,j]
        raw_step4[collect(axes)] .+= direction .* (amplitude/singular_values[j])
        step4[collect(axes)] .+= direction .* (amplitude/(singular_values[j]+singular_value_shift))
    end
    vectors4 = zeros(4,p)
    vectors4[collect(axes),:] .= decomposition.V
    (;step4,raw_step4,singular_values,retained,retained_rank=count(retained),
      numerical_rank=count(>(numerical_floor),singular_values),row_normalization=normalization,
      vectors4,axes=Tuple(axes),eigtol,singular_value_shift)
end

function _validate_absolute_controls(;eigtol,singular_value_shift,dxfix_km,dzmax_km,dzair,
        damp,late_damping_factor,dxmax_km,rback_s,bacfac,d2far_km,stopping_rule,
        dquit_km,drqt_s,stop_start_iteration)
    (eigtol===nothing || (isfinite(eigtol) && eigtol>=0)) &&
        isfinite(singular_value_shift) && singular_value_shift>=0 ||
        throw(ArgumentError("EIGTOL and the singular-value shift must be finite and nonnegative."))
    isfinite(dxfix_km) && dxfix_km>=0 && dzmax_km>0 && dxmax_km>0 && d2far_km>0 ||
        throw(ArgumentError("DXFIX must be finite nonnegative; distance limits must be positive (or Inf)."))
    isfinite(damp) && damp>0 && isfinite(late_damping_factor) && late_damping_factor>0 &&
        isfinite(dzair) && 0<dzair<=1 || throw(ArgumentError("Invalid absolute step factors."))
    (rback_s===nothing || (isfinite(rback_s) && rback_s>=0)) && isfinite(bacfac) && 0<bacfac<1 ||
        throw(ArgumentError("RBACK must be nonnegative seconds; BACFAC must lie strictly between 0 and 1."))
    stopping_rule in (:step,:legacy) && isfinite(dquit_km) && dquit_km>0 &&
        isfinite(drqt_s) && drqt_s>0 && stop_start_iteration>=1 ||
        throw(ArgumentError("Use stopping_rule=:step/:legacy and positive CON tolerances/start iteration."))
end

# A complete-vector iteration is necessary here: independently profiling the
# origin would move outside a retained mixed spatial/clock singular subspace.
function _locate_spectral(model,stations,observations,source,origin,bounds,
        constraint,constraint_source,event_id;
        max_iterations,location_tolerance_km,time_tolerance_s,finite_difference_step_km,
        eigtol,singular_value_shift,dxfix_km,dzmax_km,dzair,damp,late_damping_factor,
        dxmax_km,rback_s,bacfac,d2far_km,station_site_ids,stopping_rule,dquit_km,drqt_s,
        stop_start_iteration,legacy_iteration_limit)
    n = length(observations)
    lower,upper = bounds
    fixed_depth = lower==upper
    temporarily_locked = !fixed_depth && dxfix_km>0
    active_depth = fixed_depth
    sigma = Float64[p.sigma_s for p in observations]
    inverse_variance_sum = sum(s->inv(s^2),sigma;init=0.)
    trace = NamedTuple[]
    status = n==0 ? :no_observations : :max_iterations
    stop_reason = n==0 ? :no_observations : legacy_iteration_limit===nothing ? :iteration_limit : :hard_iteration_limit
    iterations,secant_count = 0,0
    solver_used = false
    previous_rms,previous_step = nothing,nothing
    minimum_rms = Inf
    second_distance = nothing

    for iteration in 1:(n==0 ? 0 : max_iterations)
        second_distance = _absolute_second_station_distance(stations,observations,source,station_site_ids)
        if isfinite(d2far_km) && (second_distance===nothing || second_distance>d2far_km)
            status = second_distance===nothing ? :second_station_unavailable : :second_station_too_far
            stop_reason = status
            break
        end
        residual,design,secants = _absolute_equations(model,stations,observations,source,origin,finite_difference_step_km)
        secant_count += secants
        objective = sum(abs2,residual)
        rms = sqrt(objective/inverse_variance_sum) # seconds, same W² RMS as the original
        local_lower = max(lower,_source_top(model,source))
        active_depth = fixed_depth || temporarily_locked
        axes = active_depth ? (1,2,4) : (1,2,3,4)
        spectral = _absolute_spectral_step(design,residual,sigma;axes,eigtol,singular_value_shift)
        if !active_depth && ((source[3]==local_lower && spectral.raw_step4[3]<0) ||
                             (source[3]==upper && spectral.raw_step4[3]>0))
            active_depth = true
            axes = (1,2,4)
            spectral = _absolute_spectral_step(design,residual,sigma;axes,eigtol,singular_value_shift)
        end
        if spectral.retained_rank==0
            status,stop_reason = :no_retained_directions,:no_retained_directions
            break
        end
        solver_used = true
        raw_small = norm(spectral.raw_step4[1:3])<=location_tolerance_km &&
                    abs(spectral.raw_step4[4])<=time_tolerance_s
        # A small locked-depth step releases z, it does not prove the full
        # requested problem solved. Depth is never re-locked later.
        if temporarily_locked && raw_small
            temporarily_locked = false
            minimum_rms=min(minimum_rms,rms)
            continue
        end
        legacy_reason = _absolute_legacy_stop(rms,previous_rms,previous_step;
            iteration,start_iteration=stop_start_iteration,temporarily_locked,dquit_km,drqt_s,
            legacy_iteration_limit,minimum_rms,rback_s)
        minimum_rms = min(minimum_rms,rms)
        if stopping_rule==:legacy && legacy_reason!==nothing
            status,stop_reason = :stopped_by_legacy_rule,legacy_reason
            break
        elseif stopping_rule==:step && raw_small
            status = spectral.retained_rank<length(axes) ? :converged_truncated :
                active_depth && !fixed_depth ? :converged_at_depth_bound : :converged
            stop_reason = spectral.retained_rank<length(axes) ? :retained_subspace_step : :undamped_step
            break
        end
        control = _absolute_controlled_step(spectral.step4,source;iteration,
            max_iterations=legacy_iteration_limit===nothing ? max_iterations : legacy_iteration_limit,
            damp,late_damping_factor,dzmax_km,dzair,dxmax_km,model_top_km=_source_top(model,source))
        step = control.step
        depth_fraction = step[3]<0 ? (local_lower-source[3])/step[3] :
                         step[3]>0 ? (upper-source[3])/step[3] : Inf
        bound_fraction = min(1.,depth_fraction)
        fraction = bound_fraction
        candidate_source,candidate_origin = copy(source),origin
        candidate_objective,candidate_rms = Inf,Inf
        accepted = false
        backtracks = 0
        for attempt in 0:50
            candidate_source .= source .+ fraction .* step[1:3]
            if fraction==depth_fraction
                candidate_source[3] = step[3]<0 ? local_lower : upper
            end
            candidate_origin = origin+fraction*step[4]
            candidate_objective = _absolute_objective(model,stations,observations,candidate_source,candidate_origin)
            candidate_rms = sqrt(candidate_objective/inverse_variance_sum)
            # With RBACK enabled, tolerance is measured against the historical
            # minimum weighted RMS; successive returned fractions are (1-BACFAC)^j.
            accepted = rback_s===nothing ? candidate_objective<objective ||
                (raw_small && isfinite(candidate_objective)) :
                isfinite(candidate_objective) && candidate_rms<=minimum_rms+rback_s
            accepted && break
            fraction *= rback_s===nothing ? .5 : 1-bacfac
            backtracks += 1
        end
        if !accepted
            status,stop_reason = :line_search_failed,:no_acceptable_controlled_step
            break
        end
        actual_step = [candidate_source-source;candidate_origin-origin]
        position_step = norm(actual_step[1:3])
        was_locked = temporarily_locked
        temporarily_locked && position_step<dxfix_km && (temporarily_locked=false)
        previous_rms,previous_step = rms,position_step
        source .= candidate_source
        origin = candidate_origin
        iterations = iteration
        push!(trace,(;iteration,source_km=Tuple(source),origin_s=origin,
            weighted_sum_squares=candidate_objective,position_step_km=position_step,
            origin_step_s=actual_step[4],active_depth,damping=nothing,
            weighted_rms_s=candidate_rms,retained_rank=spectral.retained_rank,
            row_normalization=spectral.row_normalization,temporarily_locked=was_locked,
            depth_released=was_locked && !temporarily_locked,
            spectral_raw_step=Tuple(spectral.raw_step4),spectral_step=Tuple(spectral.step4),
            actual_step=Tuple(actual_step),step_fraction=fraction,backtracks,
            damping_factor=control.damping_factor,air_factor=control.air_factor,
            vertical_limited=control.vertical_limited,spatial_factor=control.spatial_factor,
            second_station_distance_km=second_distance))
    end

    unused_reason = solver_used ? nothing : status==:max_iterations ? :not_iterated : status
    residuals = map(observations) do pick
        predicted = origin+_forward_time(model,Tuple(source),stations,pick.station_id,pick.phase)
        residual = pick.arrival_s-predicted
        merge(pick,(arrival_pred_s=predicted,residual_s=residual,
            normalized_residual=residual/pick.sigma_s,used=solver_used,unused_reason))
    end
    final_residual,final_design,final_secants = _absolute_equations(
        model,stations,observations,source,origin,finite_difference_step_km)
    secant_count += final_secants
    second_distance = _absolute_second_station_distance(stations,observations,source,station_site_ids)
    full_spectral = _absolute_spectral_step(final_design,final_residual,sigma;eigtol,singular_value_shift)
    final_axes = active_depth ? (1,2,4) : (1,2,3,4)
    final_spectral = _absolute_spectral_step(final_design,final_residual,sigma;
        axes=final_axes,eigtol,singular_value_shift)
    truncated = final_spectral.retained_rank<length(final_axes)
    condition_number = final_spectral.numerical_rank==length(final_axes) ?
        first(final_spectral.singular_values)/last(final_spectral.singular_values) : Inf
    unconstrained_condition_number = full_spectral.numerical_rank==4 ?
        first(full_spectral.singular_values)/last(full_spectral.singular_values) : Inf
    (;source_km=Tuple(source),origin_s=origin,status,stop_reason,iterations,event_id,
      constraint,constraint_source,free_parameter_count=fixed_depth ? 3 : 4,
      residuals,input_count=n,used_count=solver_used ? n : 0,unused_count=solver_used ? 0 : n,
      unused_reason,forward_evaluated_count=n,weighted_sum_squares=sum(abs2,final_residual),
      rms_residual_s=n==0 ? nothing : sqrt(sum(row.residual_s^2 for row in residuals)/n),
      weighted_rms_s=n==0 ? nothing : sqrt(sum(abs2,final_residual)/inverse_variance_sum),
      jacobian_rank=full_spectral.numerical_rank,free_parameter_rank=final_spectral.numerical_rank,
      retained_rank=final_spectral.retained_rank,weak_direction_truncation=truncated,
      condition_number,unconstrained_condition_number,depth_bounds=bounds,
      depth_at_bound=source[3]==max(lower,_source_top(model,source)) || source[3]==upper,
      depth_constraint_active=active_depth,temporarily_locked,second_station_distance_km=second_distance,
      secant_gradient_count=secant_count,finite_difference_evaluations=6secant_count,
      solver_method=:full_spectral,spectral=final_spectral,legacy_iteration_limit,uncertainty=nothing,trace)
end

function _absolute_legacy_stop(rms,previous_rms,previous_step;
        iteration,start_iteration,temporarily_locked,dquit_km,drqt_s,
        legacy_iteration_limit=nothing,minimum_rms=Inf,rback_s=nothing)
    iteration>=start_iteration || return nothing
    if !temporarily_locked && previous_rms!==nothing
        abs(previous_rms-rms)<drqt_s && return :rms_change_below_drqt
        previous_step<dquit_km && return :controlled_step_below_dquit
    end
    if legacy_iteration_limit!==nothing && iteration>legacy_iteration_limit &&
            isfinite(minimum_rms) && rms<minimum_rms+rback_s
        return :legacy_limit_with_rms_recovered
    end
    nothing
end
