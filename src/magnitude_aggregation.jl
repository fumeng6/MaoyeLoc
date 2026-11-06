"""
    magnitude_observation_weight(quality_code, station_weight)

Return `(weight, status, reason)` for a magnitude observation. Quality is an
integer 0--9 and the independent station magnitude weight is in 0--10. Codes
0--3 give `NINT(2.5*(4-quality)*station_weight)`; codes 4--9 give zero. NINT
rounds halfway away from zero. Location sigma, WET and location reweighting
do not enter this calculation. Missing/invalid inputs have `weight=nothing`,
never a fabricated zero; a valid zero retains its exclusion reason.
"""
function magnitude_observation_weight(quality_code,station_weight)
    (quality_code===nothing || quality_code===missing) &&
        return (;weight=nothing,status=:unavailable,reason=:missing_quality_code)
    quality_code isa Integer && !(quality_code isa Bool) && 0<=quality_code<=9 ||
        return (;weight=nothing,status=:unavailable,reason=:invalid_quality_code)
    (station_weight===nothing || station_weight===missing) &&
        return (;weight=nothing,status=:unavailable,reason=:missing_station_weight)
    station_weight isa Real && !(station_weight isa Bool) && isfinite(station_weight) && 0<=station_weight<=10 ||
        return (;weight=nothing,status=:unavailable,reason=:invalid_station_weight)
    weight = quality_code<4 ? round(2.5*(4-quality_code)*Float64(station_weight),RoundNearestTiesAway) : 0.0
    reason = quality_code>=4 ? :quality_excluded : station_weight==0 ? :zero_station_weight :
        weight==0 ? :rounded_zero_weight : :calculated
    return (;weight,status=:valid,reason)
end

# The old SORT2 is a max-heap: equal children choose the left one, and an
# equal candidate stops sifting. The resulting order of tied values matters
# for interpolation between cumulative-weight midpoints in Hypoinverse's
# MEDWT routine. Work on indices only.
function _magnitude_heap_sift!(order,values,root,last)
    candidate=order[root]
    while 2*root<=last
        child=2*root
        child<last && values[order[child]]<values[order[child+1]] && (child+=1)
        values[candidate]>=values[order[child]] && break
        order[root]=order[child]
        root=child
    end
    order[root]=candidate
end

function _magnitude_heap_order(values)
    order=collect(eachindex(values))
    for root in fld(length(order),2):-1:1
        _magnitude_heap_sift!(order,values,root,length(order))
    end
    for last in length(order):-1:2
        order[1],order[last]=order[last],order[1]
        _magnitude_heap_sift!(order,values,1,last-1)
    end
    return order
end

function _magnitude_midpoint_median(values,weights,legacy)
    order = legacy ? _magnitude_heap_order(values) : sortperm(values)
    length(order)==1 && return (only(values),order)
    target=sum(weights)/2
    cumulative=weights[order[1]]
    for j in 2:length(order)
        previous,current=order[j-1],order[j]
        upper=cumulative+weights[current]/2
        if upper>target
            lower=cumulative-weights[previous]/2
            fraction=(target-lower)/(upper-lower)
            return (muladd(fraction,values[current]-values[previous],values[previous]),order)
        end
        cumulative+=weights[current]
    end
    return (values[order[end]],order)
end

# A weighted median minimizes sum(w*abs(value-center)). At an exact half
# total-weight boundary the whole interval minimizes it; choose its midpoint.
function _magnitude_weighted_median(values,weights)
    order=sortperm(values)
    target=sum(weights)/2
    cumulative=0.0
    for (j,index) in enumerate(order)
        cumulative+=weights[index]
        if cumulative>=target
            center=values[index]
            cumulative==target && j<length(order) && (center=center/2+values[order[j+1]]/2)
            return (center,order)
        end
    end
    return (values[order[end]],order)
end

"""
    aggregate_magnitudes(values, weights; method=:median, quantization=:native)

Summarize finite station magnitudes with nonnegative observation weights.
Zero weights do not enter; negative magnitudes remain valid. Invalid vector
elements throw ArgumentError, so callers must distinguish unavailable station
calculations before calling. Empty/all-zero input returns an unavailable
result with `magnitude=nothing` and `scatter_mag=nothing`.

Native `:median` minimizes sum(w*abs(m-center)). If exactly half the weight
lies on either side of an empty interval, its midpoint selects one minimizer.
Its MAD is the same weighted median of absolute deviations from that center.
Tied values can be reordered or split into rows with the same total weight
without changing the weighted empirical distribution or these estimates.
MAD is a descriptive scatter, without a Gaussian scaling or confidence claim.

`:mean` uses sum(w*m)/sum(w); its scatter is
sqrt(sum(w^2*(m-center)^2)/sum(w^2)), the original squared-weight RMS, without
a degrees-of-freedom correction or a statistical confidence interpretation.

`:native` retains Float64 station values, centers and scatter. The explicit
`:legacy_centimag` mode instead retains Hypoinverse MEDWT's interpolation
between cumulative-weight midpoints and its original heap tie ordering.
It also uses NINT rounding: station values to 0.01 mag, center to 0.01 mag, then MAD/RMS
about that rounded center to 0.01 mag. It does not copy Float32 arithmetic or
the old missing-magnitude sentinel/clamping. All modes count positive-weight
measurement rows, including repeated station IDs supplied by the caller.
"""
function aggregate_magnitudes(values,weights;method=:median,quantization=:native)
    method in (:median,:mean) || throw(ArgumentError("Magnitude method must be :median or :mean."))
    quantization in (:native,:legacy_centimag) ||
        throw(ArgumentError("Magnitude quantization must be :native or :legacy_centimag."))
    length(values)==length(weights) || throw(ArgumentError("Magnitudes and weights must have equal lengths."))
    magnitudes,positive_weights=Float64[],Float64[]
    for (magnitude,weight) in zip(values,weights)
        magnitude isa Real && !(magnitude isa Bool) && isfinite(magnitude) && isfinite(Float64(magnitude)) ||
            throw(ArgumentError("Station magnitudes must be finite numeric values."))
        weight isa Real && !(weight isa Bool) && isfinite(weight) && isfinite(Float64(weight)) && weight>=0 ||
            throw(ArgumentError("Magnitude weights must be finite and nonnegative."))
        weight==0 && continue
        push!(magnitudes,Float64(magnitude))
        push!(positive_weights,Float64(weight))
    end
    n_measurements=length(magnitudes)
    total_weight=sum(positive_weights)
    isfinite(total_weight) || throw(ArgumentError("Magnitude total weight must be finite."))
    scatter_kind = method==:median ? (quantization==:native ? :weighted_mad : :weighted_midpoint_mad) : :squared_weight_rms
    n_measurements==0 && return (;magnitude=nothing,scatter_mag=nothing,method,quantization,
        scatter_kind,n_measurements,total_weight,status=:unavailable,reason=:no_positive_weight)
    legacy=quantization==:legacy_centimag
    if legacy
        magnitudes=round.(100 .* magnitudes,RoundNearestTiesAway)
        all(isfinite,magnitudes) || throw(ArgumentError("Centimagnitude values must remain finite."))
    end
    center = if method==:median
        middle,order=legacy ? _magnitude_midpoint_median(magnitudes,positive_weights,true) :
            _magnitude_weighted_median(magnitudes,positive_weights)
        magnitudes,positive_weights=magnitudes[order],positive_weights[order]
        middle
    else
        sum(positive_weights .* magnitudes)/total_weight
    end
    legacy && (center=round(center,RoundNearestTiesAway))
    scatter = if n_measurements==1
        0.0
    elseif method==:median
        first(legacy ? _magnitude_midpoint_median(abs.(magnitudes .- center),positive_weights,true) :
            _magnitude_weighted_median(abs.(magnitudes .- center),positive_weights))
    else
        sqrt(sum((positive_weights .* (magnitudes .- center)).^2)/sum(abs2,positive_weights))
    end
    legacy && (scatter=round(scatter,RoundNearestTiesAway))
    isfinite(center) && isfinite(scatter) || throw(ArgumentError("Magnitude aggregation produced a nonfinite result."))
    scale=legacy ? .01 : 1.0
    return (;magnitude=center*scale,scatter_mag=scatter*scale,method,quantization,
        scatter_kind,n_measurements,total_weight,status=:valid,reason=:calculated)
end
