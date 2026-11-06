# H44: source quantities remain separate from preference and display choices.
const _PREFERRED_MAGNITUDE_SLOTS = (:coda1,:amplitude1,:external,:amplitude2,:coda2)

function _preferred_slot(value)
    value isa AbstractString && Symbol(value) in _PREFERRED_MAGNITUDE_SLOTS ||
        throw(ArgumentError("Magnitude slot must be coda1, amplitude1, external, amplitude2, or coda2."))
    return Symbol(value)
end

function _preferred_lookup(summaries,slot)
    if summaries isa NamedTuple || summaries isa AbstractDict
        return _magcal_value(summaries,slot)
    end
    found = nothing
    for row in summaries
        candidate = _magcal_value(row,:slot)
        (candidate == slot || candidate == String(slot)) || continue
        found === nothing || throw(ArgumentError("Preferred magnitude requires one event and one summary per slot."))
        found = row
    end
    return found
end

function _preferred_record(record)
    if record === nothing
        return (;magnitude=nothing,label="",n_readings=nothing,reading_status=:not_supplied,
            scatter_mag=nothing,scatter_status=:not_supplied,status=:unavailable,reason=:missing_magnitude)
    end
    value = _magcal_value(record,:magnitude)
    numeric = value isa Real && !(value isa Bool) && isfinite(value)
    input_status = _magcal_value(record,:status,:valid)
    valid_status = input_status == :valid || input_status == "valid"
    status,reason = if !valid_status
        (:unavailable,_magcal_value(record,:reason,:source_status_not_valid))
    elseif value === nothing || value === missing
        (:unavailable,:missing_magnitude)
    elseif !numeric
        (:unavailable,:invalid_magnitude)
    else
        (:valid,:available)
    end
    label = _magcal_value(record,:label,"")
    label === nothing && (label = "")
    label isa AbstractString || throw(ArgumentError("Magnitude label must be a string."))
    count = _magcal_value(record,:n_readings,_magcal_value(record,:n_measurements))
    count_valid = count isa Real && !(count isa Bool) && isfinite(count) && count >= 0 && isinteger(count)
    reading_status = count === nothing || count === missing ? :not_supplied : count_valid ? :available : :invalid
    count_valid || (count = nothing)
    scatter = _magcal_value(record,:scatter_mag)
    scatter_valid = scatter isa Real && !(scatter isa Bool) && isfinite(scatter) && scatter >= 0
    scatter_status = scatter === nothing || scatter === missing ? :not_supplied : scatter_valid ? :available : :invalid
    scatter_valid || (scatter = nothing)
    return (;magnitude=numeric ? value : nothing,label=String(label),n_readings=count,reading_status,
        scatter_mag=scatter,scatter_status,status,reason)
end

"""
    preferred_magnitude(summaries; external=nothing, rules=[])

Choose the first eligible original PRE rule for one event, without changing
any of its magnitudes. `summaries` is the event's four event_magnitudes.events
rows (a vector/tuple with `slot` fields), or a NamedTuple keyed by slot.
Internal reading counts are `n_measurements`; external records use
`magnitude`, `label`, `n_readings`, `status`, and optional `scatter_mag`.
Missing status defaults to valid only when a finite numeric magnitude exists;
an explicit non-valid status cannot be made available by a numeric value.
Zero and negative magnitudes are valid; the old VNOMAG test is not used.

Each rule has a string `slot`, optional nonnegative integer `min_readings`
(default 0), and inclusive `min_magnitude`/`max_magnitude` (default -Inf/Inf).
Slots correspond to original PRE codes 1=coda1, 2=amplitude1, 3=external,
4=amplitude2, 5=coda2. Repeated slots are allowed and evaluated in order.
An empty rules list selects nothing. A zero minimum disables the count gate:
an unknown external count remains nothing but may qualify. A positive minimum
requires a known valid count meeting the threshold. No count is fabricated.

Return magnitude, label, n_readings, scatter_mag/scatter_status, source slot,
rule_index, status/reason and short rule_results for the attempted rules.
Unknown external scatter remains nothing, rather than the original forced 0.
The rule results distinguish unavailable quantities, unknown/insufficient
counts and out-of-range magnitudes. Selection stops at the first eligible
rule. Neither MAR's readable-display swap nor the legacy external-A output
convention enters this preference calculation.
"""
function preferred_magnitude(summaries;external=nothing,rules=[])
    rules isa AbstractVector || rules isa Tuple ||
        throw(ArgumentError("Preferred magnitude rules must be an ordered list."))
    checked = NamedTuple[]
    for (index,rule) in enumerate(rules)
        slot = _preferred_slot(_magcal_value(rule,:slot))
        minimum = _magcal_value(rule,:min_readings,0)
        minimum isa Real && !(minimum isa Bool) && isfinite(minimum) && minimum >= 0 && isinteger(minimum) ||
            throw(ArgumentError("Preferred magnitude min_readings must be a nonnegative integer."))
        low,high = _magcal_value(rule,:min_magnitude,-Inf),_magcal_value(rule,:max_magnitude,Inf)
        low isa Real && !(low isa Bool) && !isnan(low) && high isa Real && !(high isa Bool) && !isnan(high) && low <= high ||
            throw(ArgumentError("Preferred magnitude bounds must be ordered real values, allowing infinite bounds."))
        source = _preferred_record(slot == :external ? external : _preferred_lookup(summaries,slot))
        eligible,reason = if source.status != :valid
            (false,source.reason)
        elseif minimum > 0 && source.n_readings === nothing
            (false,:reading_count_unavailable)
        elseif minimum > 0 && source.n_readings < minimum
            (false,:insufficient_readings)
        elseif source.magnitude < low || source.magnitude > high
            (false,:magnitude_outside_range)
        else
            (true,:eligible)
        end
        push!(checked,(;rule_index=index,slot,eligible,reason,magnitude=source.magnitude,
            n_readings=source.n_readings,reading_status=source.reading_status))
        if eligible
            return (;magnitude=source.magnitude,label=source.label,n_readings=source.n_readings,
                scatter_mag=source.scatter_mag,scatter_status=source.scatter_status,slot,rule_index=index,
                status=:valid,reason=:first_eligible_rule,rule_results=checked)
        end
    end
    return (;magnitude=nothing,label="",n_readings=nothing,scatter_mag=nothing,scatter_status=:not_supplied,
        slot=nothing,rule_index=nothing,status=:unavailable,
        reason=isempty(rules) ? :no_rules_configured : :no_eligible_rule,rule_results=checked)
end

"""
    magnitude_display(summaries; external=nothing, settings=Dict(),
        in_fmag2_region=false, legacy_external_amplitude=nothing)

Prepare readable magnitude columns for one event, without mutating sources.
`settings.order` is an ordered list of 0--5 native slot strings, allowing
repeated slots. Its default is ["coda1","amplitude1","external"], the original
RDM defaults 1,3,5. RDM and PRE numeric codes differ: RDM is 1=coda1, 2=coda2,
3=amplitude1, 4=amplitude2, 5=external; this API uses unambiguous slot names.

`blank_missing=true` (default) produces empty display quantities for missing
magnitudes. If false, only `display_magnitude` uses `no_magnitude_value`
(default 0); `magnitude` remains nothing and status remains unavailable.
`swap_coda_in_region=true` swaps the actual displayed coda sources only when
`in_fmag2_region=true`, as SWITCH12 does for a MAR readable summary. Values
and labels move together; preferred_magnitude and calculation slots do not.

By default external-A is kept separate. Only an explicit
`use_legacy_external_amplitude=true` permits an available
`legacy_external_amplitude` to be shown in the amplitude2 column. Each entry
retains `slot` (requested), `source_slot` (actually displayed), and the original
requested `calculated_magnitude`/`calculated_label`; source summaries remain
unchanged. This removes HYSUM's mutation of XMAG2 while preserving its useful
explicit output convention. Unavailable external-A does not replace a valid
computed amplitude2. No display placeholder can become a preferred magnitude.
"""
function magnitude_display(summaries;external=nothing,settings=Dict(),in_fmag2_region=false,
        legacy_external_amplitude=nothing)
    order = _magcal_value(settings,:order,["coda1","amplitude1","external"])
    (order isa AbstractVector || order isa Tuple) && length(order) <= 5 ||
        throw(ArgumentError("Magnitude display order must contain zero to five slot names."))
    blank = _magcal_value(settings,:blank_missing,true)
    swap = _magcal_value(settings,:swap_coda_in_region,false)
    legacy = _magcal_value(settings,:use_legacy_external_amplitude,false)
    blank isa Bool && swap isa Bool && legacy isa Bool && in_fmag2_region isa Bool ||
        throw(ArgumentError("Magnitude display switches must be booleans."))
    no_magnitude = _magcal_value(settings,:no_magnitude_value,0.0)
    no_magnitude isa Real && !(no_magnitude isa Bool) && isfinite(no_magnitude) ||
        throw(ArgumentError("The missing-magnitude display value must be finite."))
    entries = NamedTuple[]
    for (position,name) in enumerate(order)
        slot = _preferred_slot(name)
        source_slot = swap && in_fmag2_region && slot == :coda1 ? :coda2 :
            swap && in_fmag2_region && slot == :coda2 ? :coda1 : slot
        original = _preferred_record(slot == :external ? external : _preferred_lookup(summaries,slot))
        source = _preferred_record(source_slot == :external ? external : _preferred_lookup(summaries,source_slot))
        reason = source_slot == slot ? :source_quantity : :mar_coda_display_swap
        if legacy && slot == :amplitude2
            alternative = _preferred_record(legacy_external_amplitude)
            if alternative.status == :valid
                source,source_slot,reason = alternative,:legacy_external_amplitude,:legacy_external_amplitude
            end
        end
        valid = source.status == :valid
        push!(entries,(;position,slot,source_slot,magnitude=valid ? source.magnitude : nothing,
            label=source.label,n_readings=source.n_readings,
            scatter_mag=source.scatter_mag,scatter_status=source.scatter_status,
            status=source.status,reason=valid ? reason : source.reason,
            display_magnitude=valid ? source.magnitude : blank ? nothing : no_magnitude,
            display_label=valid || !blank ? source.label : "",
            calculated_magnitude=original.status == :valid ? original.magnitude : nothing,
            calculated_label=original.label))
    end
    return entries
end
