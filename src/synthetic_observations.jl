"""
    synthetic_differences(model, stations, events; phases=("P","S"),
        noise_half_width_s=.002, seed=0, station_ids=sort(collect(keys(stations))))

Generate complete arrival differences O(i)-O(j)+T(i)-T(j)+Uniform(-a,a)
for every original-order event pair, selected station and phase. Independent
uniform noise has standard deviation a/sqrt(3). Zero noise uses a labeled
1 s numerical row scale, not a measurement uncertainty. Rows have source CC
as simulated independent differential measurements, not shared-pick CT errors.
"""
function synthetic_differences(model,stations,events;phases=("P","S"),noise_half_width_s=.002,seed=0,
        station_ids=sort(collect(keys(stations))))
    !isempty(phases) && all(p in ("P","S") for p in phases) && length(unique(phases))==length(phases) ||
        throw(ArgumentError("phases must contain distinct P and/or S entries."))
    noise_half_width_s isa Real && isfinite(noise_half_width_s) && noise_half_width_s>=0 ||
        throw(ArgumentError("Synthetic noise half-width must be finite nonnegative seconds."))
    seed isa Integer || throw(ArgumentError("Synthetic random seed must be an integer."))
    catalog=collect(events)
    length(unique(event.event_id for event in catalog))==length(catalog) || throw(ArgumentError("Synthetic event IDs must be unique."))
    length(catalog)>=2 || throw(ArgumentError("Synthetic differential generation requires at least two events."))
    ids=String.(collect(station_ids))
    length(unique(ids))==length(ids) && all(haskey(stations,id) for id in ids) || throw(ArgumentError("Synthetic station order must contain unique known IDs."))
    rng=MersenneTwister(seed)
    half_width=Float64(noise_half_width_s)
    sigma=half_width==0 ? 1. : half_width/sqrt(3.)
    rows,times=NamedTuple[],NamedTuple[]
    for phase in phases,station_id in ids
        values=Float64[_forward_time(model,event.source_km,stations,station_id,phase) for event in catalog]
        append!(times,[(;event_id=catalog[j].event_id,station_id,phase,travel_time_s=values[j]) for j in eachindex(catalog)])
        for j in 1:length(catalog)-1,k in j+1:length(catalog)
            exact=values[j]-values[k]
            noise=half_width==0 ? 0. : half_width*(2rand(rng)-1)
            origin_difference=get(catalog[j],:origin_s,0.)-get(catalog[k],:origin_s,0.)
            push!(rows,(;event_i=catalog[j].event_id,event_j=catalog[k].event_id,station_id,phase,
                dt_s=origin_difference+exact+noise,sigma_s=sigma,source="CC",
                truth_dt_s=origin_difference+exact,travel_difference_s=exact,synthetic_dt_s=exact+noise,
                reference_origin_difference_s=origin_difference,noise_s=noise))
        end
    end
    (;events=[merge(event,(origin_s=Float64(get(event,:origin_s,0.)),)) for event in catalog],observations=rows,travel_times=times,
      original_events=catalog,metadata=(;mode=:independent_synthetic_differences,phases,seed,
        noise_distribution=:independent_uniform,noise_half_width_s=half_width,
        sigma_s=sigma,sigma_role=half_width==0 ? :deterministic_reader_placeholder : :uniform_second_moment,
        rng=:MersenneTwister,julia_version=string(VERSION),station_order=ids,
        pair_order=:original_catalog_j_less_than_k,origin_definition=:native_relative_utc_with_zero_initial_corrections,
        synthetic_difference=:travel_time_difference_plus_noise,native_difference=:reference_origin_difference_plus_synthetic_difference,
        row_classes=:independent_differential_measurements))
end
