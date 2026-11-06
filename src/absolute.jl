# An explicit ID distinguishes separate readings of the same station/phase.
# Legacy API inputs without IDs keep their former station/phase key.
function _absolute_observation_key(row)
    value(key)=row isa AbstractDict ? get(row,String(key),nothing) : get(row,key,nothing)
    id=value(:observation_id)
    if id!==nothing
        id isa AbstractString && !isempty(strip(id)) ||
            throw(ArgumentError("observation_id must be a nonempty string."))
        return (:observation_id,String(id))
    end
    return (:station_phase,String(value(:station_id)),String(value(:phase)))
end

"""
    locate(model, stations, picks; initial_source_km, initial_origin_s, ...)

Estimate one event's `(x, y, z)` in km and origin time in seconds from absolute
P/S arrivals. This is the only absolute-location path: empirical quality,
station, phase, distance and residual row factors, then an SVD correction of
position and origin with step-length control and the legacy RMS/step stop.
`stations` maps station IDs to `(x, y, z)`; each pick supplies `station_id`,
`phase`, `arrival_s`, positive `sigma_s`, and `quality_code` (or an explicit
integer `default_quality` from 0 through 9). A missing code is not quality 0,
and `sigma_s` is not converted into a quality. `sigma_s` stays the recorded
arrival scatter; it is not the equation weight.
Optional `event_rules` separates the four-reading iteration rule from final
MIN catalog eligibility. `legacy_iteration_limit` with legacy stopping is
ITRLIM, while `max_iterations` remains the hard attempt cap. Zero ITRLIM
fixes XYZ and solves origin unless all four parameters were already fixed.
With `stopping_rule=:step`, `legacy_iteration_limit` defaults to `nothing`;
the retained undamped correction is checked against the step tolerances.

Position and origin enter one SVD correction. Along each retained singular
direction, the projected residual is multiplied by `1/(s+singular_value_shift)`;
that multiplier is the reported spectral gain, unrelated to instrument gain.
Depth bounds are constraints in the common z-positive-down coordinate system,
not a change to the velocity model. `:converged_truncated` refers only to the
retained subspace. `:stopped_by_legacy_rule` records the requested RMS and
controlled-step rule. This routine does not estimate uncertainties.

When a first-arrival branch has no unique analytic derivative, use local
finite-difference secants, reporting their use. These approximate slopes at a
corner, not a claim that the travel-time function is differentiable there.
"""
function locate(model, stations::AbstractDict, picks;
                initial_source_km, initial_origin_s,
                constraint::Symbol=:free, constraint_source::AbstractString="unconstrained",
                depth_bounds=(model_top(model), Inf), max_iterations::Integer=50,
                location_tolerance_km::Real=1e-6, time_tolerance_s::Real=1e-7,
                finite_difference_step_km::Real=1e-5,
                eigtol::Real=.012, singular_value_shift::Real=.006,
                dxfix_km::Real=7., dzmax_km::Real=30., dzair::Real=.5,
                damp::Real=.9, late_damping_factor::Real=.5, dxmax_km::Real=50.,
                rback_s::Real=.02, bacfac::Real=.6, d2far_km::Real=250.,
                station_site_ids=Dict{String,String}(), stopping_rule::Symbol=:legacy,
                dquit_km::Real=.04, drqt_s::Real=.001, stop_start_iteration::Integer=1,
                row_weighting::Symbol=:empirical,quality_weights=(1.,.75,.5,.25),
                station_time_weights=Dict(),s_weight::Real=.5,distance_weighting=(;),
                initial_distance_weighting=(;),residual_weighting=(;),event_rules=nothing,
                default_quality=nothing,
                legacy_iteration_limit::Union{Nothing,Integer}=stopping_rule==:legacy ? 20 : nothing)
    length(initial_source_km) == 3 ||
        throw(ArgumentError("initial_source_km must contain x, y, z in km."))
    source = Float64[initial_source_km...]
    origin = Float64(initial_origin_s)
    all(isfinite, source) && isfinite(origin) ||
        throw(ArgumentError("Initial position and origin time must be finite."))
    lower, upper = Float64.(depth_bounds)
    isfinite(lower) && lower >= model_top(model) && upper >= lower ||
        throw(ArgumentError("Depth bounds must lie within the model and satisfy lower <= upper."))
    lower <= source[3] <= upper ||
        throw(ArgumentError("Initial depth lies outside depth_bounds; it is not silently moved."))
    source[3] >= _source_top(model,source) ||
        throw(ArgumentError("Initial source is above an active model top."))
    constraint in (:free,:fixed_depth,:fixed_position,:fixed_all) ||
        throw(ArgumentError("Unknown absolute-location constraint $constraint."))
    isempty(strip(constraint_source)) && throw(ArgumentError("Constraint source must be explicit and nonempty."))
    constraint!=:free && constraint_source=="unconstrained" && (constraint_source="locate_constraint_keyword")
    if constraint==:fixed_depth
        lower = upper = source[3]
    elseif constraint==:free && lower==upper
        constraint = :fixed_depth
        constraint_source = "depth_bounds"
    end
    max_iterations >= 0 || throw(ArgumentError("max_iterations must be nonnegative."))
    if legacy_iteration_limit!==nothing
        legacy_iteration_limit isa Integer && legacy_iteration_limit>=0 && stopping_rule==:legacy ||
            throw(ArgumentError("legacy_iteration_limit is a nonnegative ITRLIM and requires stopping_rule=:legacy."))
        rback_s===nothing && (rback_s=.02)
        if legacy_iteration_limit==0 && !(constraint in (:fixed_position,:fixed_all))
            constraint=:fixed_position
            constraint_source="legacy_zero_iteration_limit"
        end
    end
    location_tolerance_km > 0 && time_tolerance_s > 0 &&
        finite_difference_step_km > 0 ||
        throw(ArgumentError("Step tolerances and finite-difference step must be positive."))

    observations = collect(picks)
    seen = Set{Tuple}()
    for pick in observations
        haskey(stations, pick.station_id) ||
            throw(ArgumentError("Unknown station ID $(pick.station_id)."))
        pick.phase in ("P", "S") || throw(ArgumentError("Phase must be P or S."))
        isfinite(pick.arrival_s) && isfinite(pick.sigma_s) && pick.sigma_s > 0 ||
            throw(ArgumentError("Arrivals must be finite and sigma_s must be finite and positive."))
        key = _absolute_observation_key(pick)
        key in seen && throw(ArgumentError("Duplicate reading identity; same-station/phase readings require distinct observation_id values."))
        push!(seen, key)
    end
    # Absolute picks for more than one event must be located separately.
    event_ids = unique(pick.event_id for pick in observations if hasproperty(pick, :event_id))
    length(event_ids) <= 1 ||
        throw(ArgumentError("locate accepts observations for one event at a time."))

    row_weighting==:empirical || throw(ArgumentError("Absolute location uses empirical quality, station, phase, distance and residual weights."))
    _validate_absolute_controls(;eigtol,singular_value_shift,dxfix_km,dzmax_km,dzair,
        damp,late_damping_factor,dxmax_km,rback_s,bacfac,d2far_km,stopping_rule,
        dquit_km,drqt_s,stop_start_iteration)
    default_quality===nothing || (default_quality isa Integer && !(default_quality isa Bool) && 0<=default_quality<=9) ||
        throw(ArgumentError("default_quality must be an integer from 0 through 9."))
    all(id->haskey(stations,id),keys(station_time_weights)) || throw(ArgumentError("Unknown station in station_time_weights."))
    return _locate_empirical(model,stations,observations,source,origin,(lower,upper),constraint,
        String(constraint_source),isempty(event_ids) ? nothing : only(event_ids);
        quality_weights,station_time_weights,s_weight,distance_weighting,initial_distance_weighting,residual_weighting,
        station_site_ids,max_iterations,location_tolerance_km,time_tolerance_s,finite_difference_step_km,
        eigtol=Float64(eigtol),singular_value_shift,dxfix_km,dzmax_km,dzair,damp,
        late_damping_factor,dxmax_km,rback_s,bacfac,d2far_km,stopping_rule,dquit_km,drqt_s,stop_start_iteration,
        event_rules,legacy_iteration_limit,default_quality)
end

function _absolute_svd(design)
    factorization = svd(design; full=false)
    threshold = isempty(factorization.S) ? 0.0 :
        max(size(design)...) * eps(Float64) * first(factorization.S)
    rank = count(value -> value > threshold, factorization.S)
    return factorization, rank
end

function _absolute_equations(model, stations, observations, source, origin, difference_step)
    n = length(observations)
    residual = zeros(n)
    design = zeros(n, 4)
    secant_count = 0
    for (row, pick) in enumerate(observations)
        prediction = travel_time(model, Tuple(source), stations, pick.station_id, pick.phase)
        derivative = prediction.dsource_s_km
        if isnothing(derivative)
            derivative = _absolute_secant(model, source, stations, pick.station_id, pick.phase, difference_step)
            secant_count += 1
        end
        residual[row] = (pick.arrival_s - origin - prediction.time_s) / pick.sigma_s
        design[row, 1:3] .= derivative ./ pick.sigma_s
        design[row, 4] = 1 / pick.sigma_s
    end
    return residual, design, secant_count
end

function _absolute_secant(model, source, station, phase, step;
                          source_top_km=model_top(model))
    return ntuple(3) do coordinate
        before = copy(source)
        after = copy(source)
        before[coordinate] -= step
        after[coordinate] += step
        coordinate == 3 && (before[3] = max(before[3], source_top_km))
        t_before = travel_time(model, Tuple(before), station, phase).time_s
        t_after = travel_time(model, Tuple(after), station, phase).time_s
        (t_after - t_before) / (after[coordinate] - before[coordinate])
    end
end

_absolute_secant(model,source,stations::AbstractDict,station_id,phase,step) =
    _absolute_secant(_station_model(model,station_id),source,stations[station_id],phase,step;
                     source_top_km=model_top(model))

function _absolute_secant(model::RegionalModels,source,stations::AbstractDict,station_id,phase,step)
    ntuple(3) do coordinate
        before,after = collect(source),collect(source)
        before[coordinate] -= step
        after[coordinate] += step
        coordinate == 3 && (before[3]=max(before[3],model_top(model)))
        times = map((before,after)) do point
            try
                _forward_time(model,Tuple(point),stations,station_id,phase)
            catch exception
                exception isa DomainError || rethrow()
                nothing
            end
        end
        t_before,t_after = times
        if t_before === nothing
            t_after === nothing && throw(DomainError(source,"No feasible regional secant on this axis."))
            center = _forward_time(model,Tuple(source),stations,station_id,phase)
            return (t_after-center)/(after[coordinate]-source[coordinate])
        elseif t_after === nothing
            center = _forward_time(model,Tuple(source),stations,station_id,phase)
            return (center-t_before)/(source[coordinate]-before[coordinate])
        end
        (t_after-t_before)/(after[coordinate]-before[coordinate])
    end
end

function _absolute_objective(model, stations, observations, source, origin)
    objective = 0.0
    for pick in observations
        prediction = try
            travel_time(model, Tuple(source), stations, pick.station_id, pick.phase)
        catch error
            # The homogeneous derivative is undefined at source == station;
            # a line-search trial there is rejected, without shifting the source.
            error isa DomainError || rethrow()
            return Inf
        end
        objective += ((pick.arrival_s - origin - prediction.time_s) / pick.sigma_s)^2
    end
    return objective
end
