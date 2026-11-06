# Hypoinverse's quality/station/S factors multiply equation ROWS. They are
# empirical coefficients, not measurement sigmas or inverse variances.
function _absolute_quality(pick,default_quality=nothing)
    if hasproperty(pick,:quality_code)
        value = getproperty(pick,:quality_code)
        source = :explicit
    elseif default_quality!==nothing
        value = default_quality
        source = :explicit_default
    else
        throw(ArgumentError("Absolute location requires quality_code on every arrival, or an explicit default_quality from 0 through 9. A missing code is not full weight, and sigma_s is not converted into a quality code."))
    end
    value isa Integer && !(value isa Bool) && 0<=value<=9 ||
        throw(ArgumentError("Phase quality_code must be an integer from 0 through 9."))
    Int(value),source
end

_absolute_weight_option(options,key,default) = options isa AbstractDict ?
    get(options,String(key),default) : get(options,key,default)

function _absolute_cosine_weight(value,inner,outer)
    value<=inner && return 1.
    value>=outer && return 0.
    (1+cospi((value-inner)/(outer-inner)))/2
end

function _absolute_distance_options(options;initial=false)
    options===nothing && return nothing
    start = _absolute_weight_option(options,:start_iteration,initial ? 100 : 4)
    radius = Float64(_absolute_weight_option(options,:minimum_radius_km,50.))
    inner = Float64(_absolute_weight_option(options,:inner_factor,initial ? 3. : 1.))
    outer = Float64(_absolute_weight_option(options,:outer_factor,initial ? 6. : 3.))
    start isa Integer && start>=1 && all(isfinite,(radius,inner,outer)) && radius>=0 && 0<=inner<outer ||
        throw(ArgumentError("Distance weights need a positive start iteration/radius and 0 <= inner_factor < outer_factor."))
    (;start_iteration=Int(start),minimum_radius_km=radius,inner_factor=inner,outer_factor=outer)
end

function _absolute_residual_options(options)
    options===nothing && return nothing
    start = _absolute_weight_option(options,:start_iteration,4)
    floor_s = Float64(_absolute_weight_option(options,:floor_s,.16))
    inner = Float64(_absolute_weight_option(options,:inner_factor,1.5))
    outer = Float64(_absolute_weight_option(options,:outer_factor,3.))
    passes = _absolute_weight_option(options,:passes,2)
    start isa Integer && start>=1 && passes isa Integer && passes>=1 &&
        all(isfinite,(floor_s,inner,outer)) && floor_s>=0 && 0<=inner<outer ||
        throw(ArgumentError("Residual weights need positive start/passes/floor_s and 0 <= inner_factor < outer_factor."))
    (;start_iteration=Int(start),floor_s,inner_factor=inner,outer_factor=outer,passes=Int(passes))
end

function _absolute_event_rules(options)
    options===nothing && return nothing
    minimum=_absolute_weight_option(options,:minimum_readings,4)
    policy=Symbol(_absolute_weight_option(options,:insufficient_readings,"exit"))
    minimum isa Integer && minimum>=4 && policy in (:exit,:cancel_reweighting) ||
        throw(ArgumentError("MIN requires minimum_readings >= 4; insufficient_readings is exit or cancel_reweighting."))
    (;minimum_readings=Int(minimum),insufficient_readings=policy)
end

function _validate_absolute_weight_factors(quality_weights,station_time_weights,s_weight)
    length(quality_weights)==4 && all(w->w isa Real && isfinite(w) && w>=0,quality_weights) &&
        s_weight isa Real && isfinite(s_weight) && s_weight>=0 &&
        all(w->w isa Real && isfinite(w) && w>=0,values(station_time_weights)) ||
        throw(ArgumentError("Quality, station and S equation factors must be finite and nonnegative; WET has four entries."))
    nothing
end

function _absolute_empirical_weights(observations,stations,source,residual_s;
        quality_weights=(1.,.75,.5,.25),station_time_weights=Dict(),s_weight=1.,
        iteration=1,station_site_ids=Dict(),distance_weighting=nothing,
        initial_distance_weighting=nothing,residual_weighting=nothing,
        reweighting_enabled=true,force_reweighting=false,default_quality=nothing)
    _validate_absolute_weight_factors(quality_weights,station_time_weights,s_weight)
    n = length(observations)
    length(residual_s)==n || throw(ArgumentError("Residuals must match the empirical-weight observations."))
    decoded = [_absolute_quality(p,default_quality) for p in observations]
    quality_code = first.(decoded)
    quality_source = last.(decoded)
    quality = [q<4 ? Float64(quality_weights[q+1]) : 0. for q in quality_code]
    station = [Float64(get(station_time_weights,p.station_id,1.)) for p in observations]
    phase = [p.phase=="S" ? Float64(s_weight) : 1. for p in observations]
    base = quality.*station.*phase
    all(isfinite,base) || throw(ArgumentError("The product of empirical equation factors overflowed."))
    distance = ones(n)
    # Original eligibility is phase presence and q<4, independent of WET/SWT.
    # Physical-site de-duplication fixes the old equal-distance/phase counting.
    candidates = [p for (p,q) in zip(observations,quality_code) if q<4]
    second = _absolute_second_station_distance(stations,candidates,source,station_site_ids)
    dis = _absolute_distance_options(distance_weighting)
    di1 = _absolute_distance_options(initial_distance_weighting;initial=true)
    rms = _absolute_residual_options(residual_weighting)
    dis_active = reweighting_enabled && dis!==nothing && (force_reweighting || iteration>=dis.start_iteration)
    di1_active = reweighting_enabled && di1!==nothing && (force_reweighting || iteration>=di1.start_iteration)
    status = :ok
    if (dis_active || di1_active) && second===nothing
        status = :second_station_unavailable
    elseif dis_active || di1_active
        for (i,pick) in enumerate(observations)
            d = hypot(source[1]-stations[pick.station_id][1],source[2]-stations[pick.station_id][2])
            # DIS takes precedence only where its inner-radius condition holds.
            a = dis_active ? max(second,dis.minimum_radius_km)*dis.inner_factor : Inf
            if dis_active && d>=a
                b = max(second,dis.minimum_radius_km)*dis.outer_factor
                distance[i]=_absolute_cosine_weight(d,a,b)
            elseif di1_active
                r = max(second,di1.minimum_radius_km)
                distance[i]=_absolute_cosine_weight(d,r*di1.inner_factor,r*di1.outer_factor)
            end
        end
    end
    weights = base.*distance
    residual_factors = ones(n)
    references = Float64[]
    rms_active = reweighting_enabled && rms!==nothing && (force_reweighting || iteration>=rms.start_iteration)
    if rms_active
        for _ in 1:rms.passes
            denominator = sum(abs2,weights)
            denominator>0 || break
            reference = max(rms.floor_s,sqrt(sum(abs2,weights.*residual_s)/denominator))
            push!(references,reference)
            reference==0 && continue # all current residuals are exactly zero
            factors = [_absolute_cosine_weight(abs(r),rms.inner_factor*reference,rms.outer_factor*reference)
                for r in residual_s]
            residual_factors .*= factors
            weights .*= factors
        end
    end
    positive = count(>(0),weights)
    normalization = positive==0 ? nothing : positive/sum(weights)
    effective = normalization===nothing ? zeros(n) : normalization.*weights
    positive==0 && status==:ok && (status=:all_weights_zero)
    reasons = [quality[i]==0 ? :zero_quality_weight : station[i]==0 ? :zero_station_weight :
        phase[i]==0 ? :zero_phase_weight : distance[i]==0 ? :distance_weight_zero :
        residual_factors[i]==0 ? :residual_weight_zero : nothing for i in 1:n]
    (;status,quality_code,quality_source,quality,station,phase,base,distance,residual_factors,
      residual_references_s=references,raw_weights=weights,effective_weights=effective,
      positive_count=positive,normalization,second_station_distance_km=second,
      weighted_reading_count=count(>(.1),effective),
      weighted_s_count=count(i->effective[i]>.1 && observations[i].phase=="S",1:n),
      distance_active=dis_active || di1_active,residual_active=rms_active,reasons)
end

const _AbsoluteCachedRay = NamedTuple{(:time_s,:dsource_s_km),
    Tuple{Float64,Union{Nothing,NTuple{3,Float64}}}}

# This cache belongs to one event at one source position. Observation weights
# and origin time can change without changing its geometric/phase predictions.
function _absolute_cached_ray!(cache,model,stations,pick,source)
    get!(cache,(pick.station_id,pick.phase)) do
        prediction=travel_time(model,Tuple(source),stations,pick.station_id,pick.phase)
        _AbsoluteCachedRay((prediction.time_s,prediction.dsource_s_km))
    end
end

function _absolute_cached_forward_time!(cache,model,stations,pick,source)
    try
        _absolute_cached_ray!(cache,model,stations,pick,source).time_s
    catch exception
        exception isa DomainError && Tuple(source)==Tuple(stations[pick.station_id]) || rethrow()
        _forward_time(model,Tuple(source),stations,pick.station_id,pick.phase)
    end
end
# Regional forward mixtures can be valid even when a component source gradient
# is unavailable. Retain their explicit forward-only evaluation.
_absolute_cached_forward_time!(cache,model::RegionalModels,stations,pick,source)=
    _forward_time(model,Tuple(source),stations,pick.station_id,pick.phase)

function _absolute_weighted_equations(model,stations,observations,source,origin,step,weights;ray_cache=nothing)
    # Only strictly positive rows reach this function. No derivative of an
    # unused/coincident phase is required. Sigma remains unchanged in the data.
    residual=zeros(length(observations))
    design=zeros(length(observations),4)
    secants=0
    for (i,pick) in enumerate(observations)
        prediction=ray_cache===nothing ? travel_time(model,Tuple(source),stations,pick.station_id,pick.phase) :
            _absolute_cached_ray!(ray_cache,model,stations,pick,source)
        derivative=prediction.dsource_s_km
        if derivative===nothing
            derivative=_absolute_secant(model,source,stations,pick.station_id,pick.phase,step)
            secants+=1
        end
        residual[i]=(pick.arrival_s-origin-prediction.time_s)*weights[i]
        design[i,1:3] .= derivative.*weights[i]
        design[i,4]=weights[i]
    end
    residual,design,secants
end

function _absolute_empirical_residuals(model,stations,observations,source,origin,state,solver_used,constraint;ray_cache=nothing)
    map(collect(enumerate(observations))) do (i,p)
        time=ray_cache===nothing ? _forward_time(model,Tuple(source),stations,p.station_id,p.phase) :
            _absolute_cached_forward_time!(ray_cache,model,stations,p,source)
        predicted = origin+time
        residual = p.arrival_s-predicted
        fitted = solver_used && state.raw_weights[i]>0 && constraint!=:fixed_all
        reason = state.reasons[i]!==nothing ? state.reasons[i] : constraint==:fixed_all ? :all_parameters_fixed :
            fitted ? nothing : :not_fitted
        merge(p,(arrival_pred_s=predicted,residual_s=residual,normalized_residual=residual/p.sigma_s,
            quality_code=state.quality_code[i],quality_code_source=state.quality_source[i],
            quality_weight=state.quality[i],station_weight=state.station[i],phase_weight=state.phase[i],
            base_weight=state.base[i],distance_weight=state.distance[i],residual_weight=state.residual_factors[i],
            raw_weight=state.raw_weights[i],equation_weight=state.effective_weights[i],
            weighted_residual=state.effective_weights[i]*residual,used=fitted,unused_reason=reason))
    end
end

# Empirical quality, station, phase, distance and residual factors. sigma_s
# stays on the reading and is not the equation weight.
function _locate_empirical(model,stations,observations,source,origin,bounds,constraint,constraint_source,event_id;
        quality_weights,station_time_weights,s_weight,distance_weighting,initial_distance_weighting,
        residual_weighting,station_site_ids,max_iterations,location_tolerance_km,time_tolerance_s,
        finite_difference_step_km,eigtol,singular_value_shift,dxfix_km,dzmax_km,dzair,damp,
        late_damping_factor,dxmax_km,rback_s,bacfac,d2far_km,stopping_rule,dquit_km,drqt_s,stop_start_iteration,
        event_rules,legacy_iteration_limit,default_quality)
    n = length(observations)
    lower,upper = bounds
    force = constraint in (:fixed_position,:fixed_all)
    rules=_absolute_event_rules(event_rules)
    reweighting_enabled=true
    reweighting_cancelled=false
    cancellation_iteration=nothing
    rays=Dict{Tuple{String,String},_AbsoluteCachedRay}()
    ray_source=Ref(Tuple(source))
    function rays_at(point)
        position=Tuple(point)
        if position!=ray_source[]
            empty!(rays)
            ray_source[]=position
        end
        rays
    end
    function raw_residual(point,time)
        cache=rays_at(point)
        [p.arrival_s-time-_absolute_cached_forward_time!(cache,model,stations,p,point) for p in observations]
    end
    weights_at(point,time,iteration) = _absolute_empirical_weights(observations,stations,point,raw_residual(point,time);
        quality_weights,station_time_weights,s_weight,iteration,station_site_ids,distance_weighting,
        initial_distance_weighting,residual_weighting,reweighting_enabled=(reweighting_enabled::Bool),force_reweighting=force,
        default_quality)
    state = weights_at(source,origin,1)
    trace = NamedTuple[]
    weight_trace = NamedTuple[]
    iterations,secants = 0,0
    solver_used = false
    locked = !force && lower!=upper && dxfix_km>0
    active_depth = lower==upper
    status = n==0 ? :no_observations : :max_iterations
    stop_reason = n==0 ? :no_observations : legacy_iteration_limit===nothing ? :iteration_limit : :hard_iteration_limit
    previous_rms,previous_step,minimum_rms = nothing,nothing,Inf
    previous_actual_step=nothing
    backup_step=nothing
    backup_count=0
    final_iteration = 1
    function record_weights(iteration,state;weight_pass=:selected,selected_for_equations=true)
        enabled=reweighting_enabled::Bool
        cancelled=reweighting_cancelled::Bool
        for i in eachindex(observations)
            push!(weight_trace,(iteration=iteration,observation_id=get(observations[i],:observation_id,nothing),
                station_id=observations[i].station_id,phase=observations[i].phase,
                quality_code=state.quality_code[i],quality_weight=state.quality[i],station_weight=state.station[i],
                phase_weight=state.phase[i],base_weight=state.base[i],distance_weight=state.distance[i],
                residual_weight=state.residual_factors[i],raw_weight=state.raw_weights[i],
                normalization=state.normalization,equation_weight=state.effective_weights[i],
                rms_references_s=copy(state.residual_references_s),positive_count=state.positive_count,
                weighted_reading_count=state.weighted_reading_count,unused_reason=state.reasons[i],
                weight_pass,selected_for_equations,reweighting_enabled=enabled,reweighting_cancelled=cancelled))
        end
    end
    function select_weights(iteration)
        candidate=weights_at(source,origin::Float64,iteration)
        if rules!==nothing && candidate.weighted_reading_count<4 && reweighting_enabled
            record_weights(iteration,candidate;weight_pass=:insufficient_trial,selected_for_equations=false)
            if rules.insufficient_readings==:cancel_reweighting
                # Rebuild BEFORE solving, rather than applying one invalid old
                # weighted update as the original JUN branch does.
                reweighting_enabled=false
                reweighting_cancelled=true
                cancellation_iteration=iteration
                candidate=weights_at(source,origin::Float64,iteration)
            else
                return candidate,false
            end
        end
        permitted=candidate.status==:ok
        record_weights(iteration,candidate;weight_pass=reweighting_cancelled ? :base_after_cancellation : :selected,
            selected_for_equations=permitted && constraint!=:fixed_all)
        candidate,permitted
    end
    if force
        state,permitted=select_weights(1)
        if constraint==:fixed_all
            status,stop_reason = :fully_fixed,:all_parameters_fixed
        elseif permitted
            residual = raw_residual(source,origin)
            origin += sum(state.raw_weights.^2 .* residual)/sum(abs2,state.raw_weights)
            solver_used=true
            status,stop_reason = :converged,:exact_weighted_origin
        else
            status=n==0 ? :no_observations : state.status==:ok ? :insufficient_weighted_readings : state.status
            stop_reason=status
        end
    else
        for iteration in 1:(n==0 ? 0 : max_iterations)
            final_iteration=iteration
            state,permitted=select_weights(iteration)
            if !permitted
                status=state.status==:ok ? :insufficient_weighted_readings : state.status
                stop_reason=status
                solver_used=false
                break
            end
            second=state.second_station_distance_km
            if isfinite(d2far_km) && (second===nothing || second>d2far_km)
                status=second===nothing ? :second_station_unavailable : :second_station_too_far
                stop_reason=status
                break
            end
            selected=findall(>(0),state.raw_weights)
            picks=observations[selected]
            w=state.raw_weights[selected]
            residual,design,used_secants=_absolute_weighted_equations(model,stations,picks,source,origin,finite_difference_step_km,w;ray_cache=rays_at(source))
            secants+=used_secants
            local_lower=max(lower,_source_top(model,source))
            active_depth=lower==upper || locked
            axes=active_depth ? (1,2,4) : (1,2,3,4)
            spectrum=_absolute_spectral_step(design,residual,inv.(w);axes,eigtol,singular_value_shift,row_normalization=state.normalization)
            if !active_depth && ((source[3]==local_lower && spectrum.raw_step4[3]<0) || (source[3]==upper && spectrum.raw_step4[3]>0))
                active_depth=true
                axes=(1,2,4)
                spectrum=_absolute_spectral_step(design,residual,inv.(w);axes,eigtol,singular_value_shift,row_normalization=state.normalization)
            end
            if spectrum.retained_rank==0
                status,stop_reason=:no_retained_directions,:no_retained_directions
                break
            end
            solver_used=true
            objective=sum(abs2,residual)
            rms=sqrt(objective/sum(abs2,w))
            raw_small=norm(spectrum.raw_step4[1:3])<=location_tolerance_km && abs(spectrum.raw_step4[4])<=time_tolerance_s
            if locked && raw_small
                locked=false
                minimum_rms=min(minimum_rms,rms)
                continue
            end
            # Do not stop before a configured distance/residual stage has run.
            dis_start=distance_weighting===nothing ? 0 : _absolute_distance_options(distance_weighting).start_iteration
            rms_start=residual_weighting===nothing ? 0 : _absolute_residual_options(residual_weighting).start_iteration
            stages_ready=iteration>max(dis_start,rms_start)
            legacy_reason=_absolute_legacy_stop(rms,previous_rms,previous_step;iteration,
                start_iteration=stop_start_iteration,temporarily_locked=locked,dquit_km,drqt_s,
                legacy_iteration_limit,minimum_rms,rback_s)
            minimum_rms=min(minimum_rms,rms)
            if stages_ready && ((stopping_rule==:step && raw_small) || (stopping_rule==:legacy && legacy_reason!==nothing))
                status=stopping_rule==:legacy ? :stopped_by_legacy_rule : spectrum.retained_rank<length(axes) ?
                    :converged_truncated : active_depth && lower!=upper ? :converged_at_depth_bound : :converged
                stop_reason=stopping_rule==:legacy ? legacy_reason : :undamped_frozen_weight_step
                break
            end
            # Reweighting can raise RMS at the same location, even if the
            # preceding frozen objective decreased. The original RBACK moves
            # back along that preceding XYZ/O step and reevaluates all weights
            # next iteration. It must not try a new descent step first.
            if stopping_rule==:legacy && rback_s!==nothing && rms>minimum_rms+rback_s && previous_actual_step!==nothing
                backup_step=backup_count==0 ? -bacfac.*previous_actual_step : (1-bacfac).*backup_step
                backup_count+=1
                candidate=source.+backup_step[1:3]
                candidate_origin=origin+backup_step[4]
                candidate_residual=raw_residual(candidate,candidate_origin)
                candidate_objective=sum(abs2,w.*candidate_residual[selected])
                previous_rms,previous_step=rms,norm(backup_step[1:3])
                source .= candidate
                origin=candidate_origin
                iterations=iteration
                push!(trace,(;iteration,source_km=Tuple(source),origin_s=origin,
                    weighted_sum_squares=candidate_objective,position_step_km=previous_step,
                    origin_step_s=backup_step[4],active_depth,damping=nothing,
                    weighted_rms_s=sqrt(candidate_objective/sum(abs2,w)),retained_rank=spectrum.retained_rank,
                    row_normalization=state.normalization,temporarily_locked=locked,depth_released=false,
                    actual_step=Tuple(backup_step),step_fraction=1.,backtracks=backup_count,
                    second_station_distance_km=state.second_station_distance_km,
                    reweighting_enabled,reweighting_cancelled,weighted_reading_count=state.weighted_reading_count,
                    iteration_action=:legacy_rms_backup))
                continue
            end
            control=_absolute_controlled_step(spectrum.step4,source;iteration,
                max_iterations=legacy_iteration_limit===nothing ? max_iterations : legacy_iteration_limit,damp,late_damping_factor,
                dzmax_km,dzair,dxmax_km,model_top_km=_source_top(model,source))
            step=control.step
            depth_fraction=step[3]<0 ? (local_lower-source[3])/step[3] : step[3]>0 ? (upper-source[3])/step[3] : Inf
            fraction=min(1.,depth_fraction)
            candidate=copy(source)
            candidate_origin=origin
            candidate_objective=Inf
            backtracks=0
            accepted=false
            for _ in 0:50
                candidate .= source.+fraction.*step[1:3]
                fraction==depth_fraction && (candidate[3]=step[3]<0 ? local_lower : upper)
                candidate_origin=origin+fraction*step[4]
                candidate_objective=try
                    r=raw_residual(candidate,candidate_origin)
                    sum(abs2,w.*r[selected])
                catch exception
                    exception isa DomainError || rethrow()
                    Inf
                end
                candidate_rms=sqrt(candidate_objective/sum(abs2,w))
                accepted=rback_s===nothing ? (candidate_objective<objective || raw_small && isfinite(candidate_objective)) :
                    isfinite(candidate_objective) && candidate_rms<=minimum_rms+rback_s
                accepted && break
                fraction *= rback_s===nothing ? .5 : 1-bacfac
                backtracks+=1
            end
            if !accepted
                status,stop_reason=:line_search_failed,:no_acceptable_frozen_weight_step
                break
            end
            actual=[candidate-source;candidate_origin-origin]
            previous_actual_step=copy(actual)
            backup_step=nothing
            backup_count=0
            was_locked=locked
            locked && norm(actual[1:3])<dxfix_km && (locked=false)
            previous_rms,previous_step=rms,norm(actual[1:3])
            source .= candidate
            origin=candidate_origin
            iterations=iteration
            push!(trace,(;iteration,source_km=Tuple(source),origin_s=origin,weighted_sum_squares=candidate_objective,
                position_step_km=norm(actual[1:3]),origin_step_s=actual[4],active_depth,damping=nothing,
                weighted_rms_s=sqrt(candidate_objective/sum(abs2,w)),retained_rank=spectrum.retained_rank,
                row_normalization=state.normalization,temporarily_locked=was_locked,depth_released=was_locked && !locked,
                spectral_raw_step=Tuple(spectrum.raw_step4),spectral_step=Tuple(spectrum.step4),actual_step=Tuple(actual),
                step_fraction=fraction,backtracks,damping_factor=control.damping_factor,air_factor=control.air_factor,
                vertical_limited=control.vertical_limited,spatial_factor=control.spatial_factor,
                second_station_distance_km=state.second_station_distance_km,
                reweighting_enabled,reweighting_cancelled,weighted_reading_count=state.weighted_reading_count,
                iteration_action=:solve))
        end
    end
    # Preserve the last actual frozen weights. A new IRLS update must not be
    # substituted for the weights that defined the last accepted solve.
    rows=_absolute_empirical_residuals(model,stations,observations,source,origin,state,solver_used,constraint;ray_cache=rays_at(source))
    full_spectrum=nothing
    spectrum=nothing
    diagnostic_reason=nothing
    selected=findall(>(0),state.raw_weights)
    axes=force ? constraint==:fixed_all ? () : (4,) : active_depth ? (1,2,4) : (1,2,3,4)
    if !force && !isempty(selected)
        try
            r,A,k=_absolute_weighted_equations(model,stations,observations[selected],source,origin,finite_difference_step_km,state.raw_weights[selected];ray_cache=rays_at(source))
            secants+=k
            full_spectrum=_absolute_spectral_step(A,r,inv.(state.raw_weights[selected]);eigtol,singular_value_shift,row_normalization=state.normalization)
            spectrum=axes==(1,2,3,4) ? full_spectrum :
                _absolute_spectral_step(A,r,inv.(state.raw_weights[selected]);axes,eigtol,singular_value_shift,row_normalization=state.normalization)
        catch exception
            exception isa DomainError || rethrow()
            diagnostic_reason=:coincident_source_station
        end
    end
    used_count=count(row->row.used,rows)
    objective=sum((state.raw_weights[i]*rows[i].residual_s)^2 for i in 1:n;init=0.)
    base_result=(;source_km=Tuple(source),origin_s=origin,status,stop_reason,iterations,event_id,constraint,constraint_source,
        free_parameter_count=length(force ? axes : lower==upper ? (1,2,4) : (1,2,3,4)),residuals=rows,input_count=n,
        used_count,unused_count=n-used_count,unused_reason=used_count==0 ? constraint==:fixed_all ? :all_parameters_fixed : status : nothing,
        forward_evaluated_count=n,weighted_sum_squares=objective,rms_residual_s=n==0 ? nothing : sqrt(sum(row.residual_s^2 for row in rows)/n),
        weighted_rms_s=sum(abs2,state.raw_weights)==0 ? nothing : sqrt(objective/sum(abs2,state.raw_weights)),
        jacobian_rank=full_spectrum===nothing ? nothing : full_spectrum.numerical_rank,
        free_parameter_rank=force ? constraint==:fixed_all ? 0 : solver_used ? 1 : 0 : spectrum===nothing ? nothing : spectrum.numerical_rank,
        retained_rank=spectrum===nothing ? nothing : spectrum.retained_rank,
        weak_direction_truncation=spectrum!==nothing && spectrum.retained_rank<length(axes),
        condition_number=spectrum===nothing ? force && solver_used ? 1. : nothing :
            spectrum.numerical_rank==length(axes) ? first(spectrum.singular_values)/last(spectrum.singular_values) : Inf,
        unconstrained_condition_number=full_spectrum===nothing ? nothing : full_spectrum.numerical_rank==4 ?
            first(full_spectrum.singular_values)/last(full_spectrum.singular_values) : Inf,
        depth_bounds=bounds,depth_at_bound=source[3]==max(lower,_source_top(model,source)) || source[3]==upper,
        depth_constraint_active=force || active_depth,temporarily_locked=locked,
        second_station_distance_km=state.second_station_distance_km,secant_gradient_count=secants,
        finite_difference_evaluations=6secants,solver_method=force ? constraint==:fixed_all ? :fixed_all_forward : :fixed_position_origin : :empirical_full_spectral,
        row_weighting=:empirical,weight_normalization=state.normalization,positive_weight_count=state.positive_count,
        weighted_reading_count=state.weighted_reading_count,weighted_s_count=state.weighted_s_count,
        weights_iteration=final_iteration,diagnostic_reason,reweighting_cancelled,cancellation_iteration,
        catalog_eligible=rules===nothing || state.weighted_reading_count>=rules.minimum_readings,
        catalog_exclusion_reason=rules!==nothing && state.weighted_reading_count<rules.minimum_readings ? :below_minimum_readings : nothing,
        minimum_readings=rules===nothing ? nothing : rules.minimum_readings,
        legacy_iteration_limit,uncertainty=nothing,trace,weight_trace)
    spectrum===nothing ? base_result : merge(base_result,(spectral=spectrum,))
end
