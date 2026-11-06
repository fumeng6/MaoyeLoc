"""
    absolute_uncertainty(model, stations, picks, location;
                         observation_covariance_s2=nothing, assume_independent=false)

Check that `picks` are the fitted absolute readings, then report that a
covariance is not estimated for the empirical row-factor location. The fitted
weights are the saved equation factors, not `1/sigma_s`. A supplied
observation covariance is still checked for shape and positive
semidefiniteness; it is not turned into a parameter covariance here.
`absolute_formal_uncertainty` is the separate frozen-weight calculation.

Fully fixed events return `:all_parameters_fixed`. A fixed-position event
with no readings returns `:no_observations`. Neither result is filled with
zeros.
"""
function absolute_uncertainty(model, stations::AbstractDict, picks, location;
                              observation_covariance_s2=nothing,
                              assume_independent::Bool=false)
    observations = collect(picks)
    n = length(observations)
    n == location.input_count || throw(ArgumentError("Picks must match the located observation count."))
    fitted = Dict(_absolute_observation_key(row)=>row for row in location.residuals)
    length(fitted) == length(location.residuals) == n ||
        throw(ArgumentError("Fitted residuals must identify the located observations uniquely."))
    source = location.source_km
    length(source) == 3 && all(isfinite, source) && isfinite(location.origin_s) ||
        throw(ArgumentError("Location coordinates and origin time must be finite."))
    lower, upper = location.depth_bounds
    isfinite(lower) && lower >= model_top(model) && lower <= source[3] <= upper ||
        throw(ArgumentError("Location and depth bounds must lie within the model."))
    constraint = get(location,:constraint,lower == upper ? :fixed_depth : :free)
    fixed_position = constraint in (:fixed_position,:fixed_all)
    fixed_depth = fixed_position || constraint == :fixed_depth || lower == upper
    axes = constraint == :fixed_all ? () : fixed_position ? (4,) : fixed_depth ? (1, 2, 4) : (1, 2, 3, 4)
    parameter_names = Tuple((:x_km,:y_km,:z_km,:origin_s)[axis] for axis in axes)
    supplied_covariance = observation_covariance_s2 !== nothing
    get(location,:row_weighting,:empirical)==:empirical ||
        throw(ArgumentError("Absolute uncertainty describes the empirical row-factor location."))
    assumptions = (estimator=:empirical_absolute_rows,
        selected_search_minimum_conditioned=false,
        observation_covariance=supplied_covariance ? :supplied_full : assume_independent ? :independent_sigma : :unspecified,
        covariance_order=:input_picks_order, location_refitted=false,
        linearization=:first_order, residual_curvature_ignored=true,
        weight_selection_uncertainty_included=false, velocity_model_fixed=true,
        station_coordinates_fixed=true, phase_assignment_fixed=true,
        station_delays_fixed=true, station_delay_errors_included=false,
        first_arrival_branch_fixed=true, fixed_depth_conditional=fixed_depth,
        fixed_position_conditional=fixed_position,
        fixed_parameters_conditional=fixed_depth,
        residual_rescaling=false, model_bias_included=false,
        real_accuracy_claim=false)
    unavailable(reason; rank=nothing) = (status=:not_estimated, reason=reason,
        parameter_names=parameter_names, covariance=nothing,
        standard_deviations=(x_km=nothing,y_km=nothing,z_km=nothing,origin_s=nothing),
        rank=rank, observation_count=n, assumptions=assumptions)

    seen = Set{Tuple}()
    ids = unique(pick.event_id for pick in observations if hasproperty(pick,:event_id))
    length(ids) <= 1 || throw(ArgumentError("Uncertainty accepts picks for one event."))
    if !isempty(ids) && location.event_id !== nothing
        only(ids) == location.event_id || throw(ArgumentError("Pick event ID does not match location."))
    end
    for pick in observations
        haskey(stations,pick.station_id) || throw(ArgumentError("Unknown station ID $(pick.station_id)."))
        pick.phase in ("P","S") || throw(ArgumentError("Phase must be P or S."))
        isfinite(pick.arrival_s) && isfinite(pick.sigma_s) && pick.sigma_s > 0 ||
            throw(ArgumentError("Picks require finite arrivals and positive finite sigma_s."))
        key = _absolute_observation_key(pick)
        key in seen && throw(ArgumentError("Duplicate absolute observation identity."))
        original = get(fitted,key,nothing)
        original !== nothing && String(pick.station_id) == String(original.station_id) &&
            String(pick.phase) == String(original.phase) &&
            pick.arrival_s == original.arrival_s && pick.sigma_s == original.sigma_s &&
            (hasproperty(pick,:event_id) ? pick.event_id : nothing) ==
            (hasproperty(original,:event_id) ? original.event_id : nothing) ||
            throw(ArgumentError("Picks and sigma weights must match the fitted absolute location."))
        push!(seen,key)
        if fixed_position
            predicted_time = _forward_time(model,source,stations,pick.station_id,pick.phase)
            isfinite(predicted_time) || throw(ArgumentError("Predicted travel time must be finite."))
        else
            try
                travel_time(model,source,stations,pick.station_id,pick.phase)
            catch exception
                exception isa DomainError && Tuple(source)==Tuple(stations[pick.station_id]) || rethrow()
            end
        end
    end
    supplied_covariance && _observation_covariance_factor(observation_covariance_s2,n)

    constraint == :fixed_all && return unavailable(:all_parameters_fixed)
    fixed_position && n == 0 && return unavailable(:no_observations)
    return unavailable(:empirical_estimator_covariance_not_implemented)
end

# The factor is used only to propagate known observation errors, never to
# invert C. Exact dependencies (shared picks, common errors) are legitimate.
function _observation_covariance_factor(covariance_s2,n)
    covariance_s2 isa AbstractMatrix && size(covariance_s2)==(n,n) ||
        throw(ArgumentError("Observation covariance must be an n by n matrix in input observation order."))
    matrix = Matrix{Float64}(covariance_s2)
    all(isfinite,matrix) && isapprox(matrix,matrix';rtol=100eps(Float64),atol=0.) ||
        throw(ArgumentError("Observation covariance must be finite and symmetric."))
    decomposition = eigen(Symmetric((matrix+matrix')/2))
    tolerance = 100eps(Float64)*max(n,1)*maximum(abs,decomposition.values;init=0.)
    all(value->value >= -tolerance,decomposition.values) ||
        throw(ArgumentError("Observation covariance must be positive semidefinite."))
    # Tiny negative eigenvalues are roundoff; multiplying this factor keeps
    # propagated variances nonnegative even for singular C.
    return decomposition.vectors * Diagonal(sqrt.(max.(decomposition.values,0.)))
end
