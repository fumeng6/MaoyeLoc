"""
    _dd_scaled_step(design, residual; damp=0., svd_relative_cutoff=1e-7)

One column-scaled damped update. Each event contributes four columns
`(x_km, y_km, z_km, origin_s)`. `design` and the residual in seconds already
carry the equation weights. There is no extra row normalization and no
separate origin elimination.

N[j,j] is the RMS of column j. The step minimizes

    ||design*step - residual||^2 + damp^2*||N*step||^2.

damp=0 selects the minimum-norm solution in that scaled space. A column that
is entirely zero is reported and given a zero update; it is not treated as a
known parameter. `svd_relative_cutoff` is retained only so existing
diagnostics keep a field; it does not change this damped step.
"""
function _dd_scaled_step(design::AbstractMatrix, residual;
        damp::Real=0., svd_relative_cutoff::Real=1e-7)
    prepared = _dd_scaled_svd(design, residual; damp, svd_relative_cutoff)
    _dd_scaled_step_from_svd(prepared, damp, svd_relative_cutoff)
end

# A damped update and its undamped comparison use the same scaled matrix.
function _dd_scaled_svd(design, residual; damp, svd_relative_cutoff)
    m, p = size(design)
    p > 0 && p % 4 == 0 || throw(ArgumentError("DD scaled steps require four columns per event."))
    length(residual) == m && all(isfinite, design) && all(isfinite, residual) ||
        throw(ArgumentError("The design and matching residual must contain finite values."))
    isfinite(damp) && damp >= 0 && isfinite(svd_relative_cutoff) && 0 <= svd_relative_cutoff <= 1 ||
        throw(ArgumentError("DAMP must be finite nonnegative and relative cutoff must lie in [0, 1]."))
    data = Matrix{Float64}(design)
    rhs = Float64[residual...]
    rows = size(data, 1)
    column_scales = rows == 0 ? zeros(p) : vec(sqrt.(sum(abs2, data; dims=1) ./ rows))
    zero_columns = findall(iszero, column_scales)
    divisors = [iszero(scale) ? 1. : scale for scale in column_scales]
    B = data ./ permutedims(divisors)
    decomposition = svd(B; full=rows < p)
    (; data, rhs, rows, p, column_scales, zero_columns, divisors, decomposition)
end

function _dd_scaled_step_from_svd(prepared, damp, svd_relative_cutoff)
    (; data, rhs, rows, p, column_scales, zero_columns, divisors, decomposition) = prepared
    singular_values = [decomposition.S; zeros(max(0, p - length(decomposition.S)))]
    largest = maximum(singular_values; init=0.)
    numerical_floor = max(rows, p) * eps(Float64) * largest
    numerical_rank = count(>(numerical_floor), singular_values)
    retained = singular_values .> (damp > 0 ? 0. : numerical_floor)
    scaled_step = zeros(p)
    for j in findall(retained)
        value = singular_values[j]
        gain = damp > 0 ? value / (value^2 + damp^2) : inv(value)
        @views scaled_step .+= decomposition.V[:, j] .* (gain * dot(decomposition.U[:, j], rhs))
    end
    step = scaled_step ./ divisors
    step[zero_columns] .= 0.
    scaled_step[zero_columns] .= 0.
    data_residual = data * step - rhs
    penalty = Float64(damp) .* column_scales .* step
    residual_norm, penalty_norm = norm(data_residual), norm(penalty)
    condition_number = numerical_rank == p ? first(singular_values) / last(singular_values) : Inf
    regularized_condition_number = damp > 0 && isempty(zero_columns) ?
        hypot(largest, damp) / hypot(last(singular_values), damp) : condition_number
    retained_condition_number = any(retained) ? maximum(singular_values[retained]) / minimum(singular_values[retained]) : Inf
    physical_vectors = isempty(zero_columns) ? decomposition.V ./ divisors :
        Matrix{Union{Float64,Nothing}}(decomposition.V ./ divisors)
    for column in zero_columns
        physical_vectors[column, :] .= nothing
    end
    gains = [retained[j] ? (damp > 0 ? singular_values[j] / (singular_values[j]^2 + damp^2) : inv(singular_values[j])) : 0.
        for j in eachindex(singular_values)]
    return (; step, column_scales, zero_columns, singular_values, retained, gains,
        scaled_vectors=decomposition.V, physical_vectors,
        retained_rank=count(retained), numerical_rank, condition_number, regularized_condition_number,
        retained_condition_number, scaled_step,
        residual=data_residual, residual_norm, penalty_norm,
        objective=residual_norm^2 + penalty_norm^2, method=:damped, damp=Float64(damp),
        svd_relative_cutoff=Float64(svd_relative_cutoff))
end
