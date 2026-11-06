_magcal_value(row, key, default=nothing) = row isa AbstractDict ? get(row, String(key), default) : get(row, key, default)

function _magcal_utc(value; require_hour=false)
    if value isa DateTime
        stamp = value
    else
        value isa AbstractString || throw(ArgumentError("Magnitude times require a UTC DateTime or an ISO timestamp ending in Z."))
        matched = match(r"^(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})(?:\.(\d+))?Z$",value)
        matched === nothing && throw(ArgumentError("Magnitude times require an ISO timestamp ending in Z."))
        fraction = matched.captures[2]
        # Preserve the original event timestamp in the public result. Only this
        # hourly lookup helper truncates sub-millisecond digits, without rounding.
        millis = fraction === nothing ? 0 : parse(Int,rpad(first(fraction,min(3,length(fraction))),3,'0'))
        stamp = DateTime(matched.captures[1])+Millisecond(millis)
        require_hour && fraction !== nothing && any(!=('0'),fraction) &&
            throw(ArgumentError("Magnitude history expiration is specified to the UTC hour."))
    end
    require_hour && stamp != floor(stamp,Hour) &&
        throw(ArgumentError("Magnitude history expiration is specified to the UTC hour."))
    return stamp
end

# Active Hypoinverse 1.40 CALSV values, in order for 0:6:78 dB.
# The old approximate logarithmic attenuation formula is not this table.
const _MAGCAL_ATTENUATION = (25.2,11.65,5.576,2.795,1.4,.702,.352,.176,
                            .0885,.044,.0222,.0111,.00557,.00279)
const _MAGCAL_MOTOR = (nothing,1.0,nothing,1.0,.448,8.0,24.0,15.0)

function _magcal_history_value(event_hour, station_id, history)
    matched = [(i,row) for (i,row) in enumerate(history)
               if _magcal_value(row,:station_id) == station_id]
    isempty(matched) && return nothing
    previous = nothing
    for (j,(_,row)) in enumerate(matched)
        raw_expiry = _magcal_value(row,:expires_utc)
        expiry = raw_expiry === nothing ? nothing : _magcal_utc(raw_expiry;require_hour=true)
        if expiry === nothing
            j == length(matched) || throw(ArgumentError("An open-ended magnitude history record must be last for $station_id."))
        else
            previous === nothing || expiry > previous || throw(ArgumentError("Magnitude history expirations must increase for $station_id."))
            previous = expiry
        end
    end
    for (i,row) in matched
        raw_expiry = _magcal_value(row,:expires_utc)
        expiry = raw_expiry === nothing ? nothing : _magcal_utc(raw_expiry;require_hour=true)
        if expiry === nothing || event_hour < expiry
            return (; value=_magcal_value(row,:value), history_row=i, expires_utc=expiry, covered=true)
        end
    end
    return (; value=nothing,history_row=nothing,expires_utc=previous,covered=false)
end

"""
    magnitude_calibration(time_utc, channel; phase_cal=nothing,
                          calibration_history=[], attenuation_history=[])

Resolve Hypoinverse 1.40 magnitude calibration for one complete channel.
`channel` is a NamedTuple or string-keyed dictionary with `station_id` and
`instrument_type` (0:7). Optional station gain sources are exactly one of
`calibration_factor`, `attenuation_db`, `digital_system_factor_uv_per_count`.
For the digital source also supply `amplitude_units=2` (Menlo/Earthworm),
3 (HVO), or 4 (Menlo measurements encoded in hundreds of counts). CAL uses
the true-count system factor for both 2 and 4. It is dimensionless; digital F is microvolts/count, and the sensor
motor constant is separately fixed by the original instrument type in V/(cm/s).
Types 0/2 are Wood-Anderson/Sprengnether and have no velocity motor constant.

Each history row has the exact complete `station_id`, numeric `value` (CAL or
dB, according to its list), and optional `expires_utc`. Rows for each channel
are chronological; the first covers the past and the final omitted expiration
covers the future. No prefix matching, station-name broadcasting, or inferred
components are used. Type 1 selects attenuation history; other types select
CAL history. A matching history supersedes a station gain source; uncovered
history does not silently fall back to the station value. Nonzero phase CAL
supersedes both; phase zero means no override, as in the traditional input.

UTC times are DateTime values (caller-declared UTC) or ISO strings ending in Z.
Both event time and expiration are resolved to the original hourly granularity.
Intervals are half open: at the expiration hour the next record applies. This
is stateless and deliberately removes the original initial-loader/runtime
disagreement at exact expiration and the requirement for ordered events.

Return the selected value/source/epoch and distinct valid, zero/unspecified,
invalid, or expired-history states. CAL=0 permits a zero coda gain correction
but cannot calibrate an amplitude magnitude. Invalid calibration is not a
reference gain and returns no coda correction. Positive CAL gives
`coda_gain_correction=-log10(CAL/3.95)` before the separate H41/DUG suppression.
This function does not calculate a station or event magnitude.
"""
function magnitude_calibration(time_utc, channel; phase_cal=nothing,
                               calibration_history=[], attenuation_history=[])
    event_hour = floor(_magcal_utc(time_utc),Hour)
    station_id = _magcal_value(channel,:station_id)
    station_id isa AbstractString && !isempty(station_id) ||
        throw(ArgumentError("A complete station_id identifying the channel is required."))
    raw_type = _magcal_value(channel,:instrument_type)
    raw_type isa Integer && !(raw_type isa Bool) && 0 <= raw_type <= 7 ||
        throw(ArgumentError("Magnitude instrument_type must be an integer from 0 to 7."))
    instrument_type = Int(raw_type)
    motor = _MAGCAL_MOTOR[instrument_type+1]
    source, value, kind = :unspecified, 0.0, :calibration
    history_row, expiry, system_factor, units = nothing,nothing,nothing,nothing
    history_covered = true
    if phase_cal !== nothing && phase_cal != 0
        source, value = :phase_calibration, phase_cal
    else
        history = instrument_type == 1 ? attenuation_history : calibration_history
        record = _magcal_history_value(event_hour,station_id,history)
        if record !== nothing
            source = instrument_type == 1 ? :attenuation_history : :calibration_history
            kind = instrument_type == 1 ? :attenuation : :calibration
            value,history_row,expiry,history_covered = record.value,record.history_row,record.expires_utc,record.covered
        else
            sources = [(key,_magcal_value(channel,key)) for key in
                       (:calibration_factor,:attenuation_db,:digital_system_factor_uv_per_count)
                       if _magcal_value(channel,key) !== nothing]
            length(sources) <= 1 || throw(ArgumentError("Specify only one station gain source for $station_id."))
            if !isempty(sources)
                key,value = only(sources)
                if key == :calibration_factor
                    source = :station_calibration
                elseif key == :attenuation_db
                    source,kind = :station_attenuation,:attenuation
                else
                    source,kind = :digital_system_factor,:digital
                    system_factor = value
                    units = _magcal_value(channel,:amplitude_units)
                end
            end
        end
    end
    valid_number = value isa Real && !(value isa Bool) && isfinite(value)
    cal, status, reason = nothing,:invalid_calibration,:invalid_calibration_factor
    if !history_covered
        status,reason = :calibration_history_expired,:no_history_record_for_event
    elseif kind == :attenuation
        if instrument_type != 1
            reason = :attenuation_requires_instrument_type_1
        elseif valid_number && 0 <= value <= 78 && isinteger(value/6)
            cal,status,reason = _MAGCAL_ATTENUATION[Int(value/6)+1],:valid,:valid
        else
            reason = :attenuation_must_be_0_to_78_in_steps_of_6_db
        end
    elseif kind == :digital
        if motor === nothing
            reason = :digital_system_factor_requires_velocity_instrument
        elseif !(units isa Integer) || units isa Bool || !(units in (2,3,4))
            reason = :digital_system_factor_requires_amplitude_units_2_3_or_4
        elseif valid_number && value > 0
            cal = (units == 3 ? .3456 : 1.382)/value
            status,reason = isfinite(cal) && cal > 0 ? (:valid,:valid) : (:invalid_calibration,:invalid_digital_system_factor)
            status == :valid || (cal = nothing)
        else
            reason = :digital_system_factor_must_be_positive_finite
        end
    elseif valid_number && value >= 0
        cal = Float64(value)
        status,reason = cal > 0 ? (:valid,:valid) :
            (:unknown_calibration,source == :unspecified ? :calibration_not_supplied : :zero_calibration)
    end
    coda_gain = status == :valid ? log10(3.95)-log10(cal) : status == :unknown_calibration ? 0.0 : nothing
    digital_sensitivity = kind == :digital && status == :valid ? motor*1e8/system_factor : nothing
    return (; station_id=String(station_id),event_time_utc=time_utc,event_hour_utc=event_hour,
            instrument_type,motor_constant_v_per_cm_s=motor,
            calibration_factor=cal,calibration_source=source,history_row,expires_utc=expiry,
            system_factor_uv_per_count=kind == :digital && valid_number ? Float64(system_factor) : nothing,
            digital_sensitivity_counts_per_m_s=digital_sensitivity,
            status,reason,amplitude_calibrated=status == :valid,coda_gain_correction=coda_gain)
end
