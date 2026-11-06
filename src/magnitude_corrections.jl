# Station corrections are additive magnitude increments; weights and gain use
# separate explicit fields, independent of the correction value.
function _magcorr_station(event_hour, station_id, channel, history, coda, use_unlisted)
    weight_key = coda ? :coda_weight : :amplitude_weight
    correction_key = coda ? :coda_correction_mag : :amplitude_correction_mag
    weight = _magcal_value(channel,weight_key,10.0)
    weight isa Real && !(weight isa Bool) && isfinite(weight) && 0 <= weight <= 10 ||
        throw(ArgumentError("Station magnitude weights must lie between 0 and 10."))
    gain_enabled = coda ? _magcal_value(channel,:coda_gain_enabled,true) : true
    gain_enabled isa Bool || throw(ArgumentError("coda_gain_enabled must be true or false."))
    record = _magcal_history_value(event_hour,station_id,history)
    if record !== nothing && !record.covered
        return (; correction=nothing,event_weight=0.0,gain_enabled=nothing,
                source=:history,status=:history_expired,history_row=nothing,
                expires_utc=record.expires_utc,source_label="")
    end
    if record === nothing
        value = _magcal_value(channel,correction_key,0.0)
        weight = use_unlisted ? Float64(weight) : 0.0
        source,status,row,expiry,label = :station,:valid,nothing,nothing,""
    else
        value = record.value
        row,expiry = record.history_row,record.expires_utc
        source,status = :history,:valid
        label = string(_magcal_value(history[row],:source,""))
        weight = _magcal_value(history[row],:event_weight,weight)
        coda && (gain_enabled = _magcal_value(history[row],:gain_enabled,gain_enabled))
    end
    value isa Real && !(value isa Bool) && isfinite(value) ||
        throw(ArgumentError("A station magnitude correction must be finite."))
    weight isa Real && !(weight isa Bool) && isfinite(weight) && 0 <= weight <= 10 ||
        throw(ArgumentError("Station magnitude weights must lie between 0 and 10."))
    gain_enabled isa Bool || throw(ArgumentError("Coda gain_enabled must be true or false."))
    return (; correction=Float64(value),event_weight=Float64(weight),
            gain_enabled,source,status,history_row=row,
            expires_utc=expiry,source_label=label)
end

function _magcorr_components(component, corrections)
    total,rows,sources = 0.0,Int[],String[]
    for (i,row) in enumerate(corrections)
        requested = _magcal_value(row,:component)
        requested isa AbstractString && !isempty(requested) ||
            throw(ArgumentError("A component magnitude correction requires its complete component code."))
        requested == component || continue
        value = _magcal_value(row,:value)
        value isa Real && !(value isa Bool) && isfinite(value) ||
            throw(ArgumentError("A component magnitude correction must be finite."))
        total += Float64(value)
        push!(rows,i)
        push!(sources,string(_magcal_value(row,:source,"")))
    end
    return (;total,rows,sources)
end

function _magcorr_instrument(station_id, amplitude_history)
    override = nothing
    for row in amplitude_history
        _magcal_value(row,:station_id) == station_id || continue
        value = _magcal_value(row,:instrument_type)
        value === nothing && continue
        value isa Integer && !(value isa Bool) && 0 <= value <= 7 ||
            throw(ArgumentError("An amplitude-history instrument_type override must be an integer from 0 to 7."))
        override === nothing || override == value ||
            throw(ArgumentError("Amplitude-history instrument type is static for a channel, not a time history."))
        override = Int(value)
    end
    return override
end

"""
    magnitude_corrections(time_utc, channel; duration_history=[],
        amplitude_history=[], component_coda=[], component_amplitude=[],
        use_unlisted_coda=true, use_unlisted_amplitude=true)

Resolve additive station and component magnitude corrections for one complete
channel. `channel` is a NamedTuple or string-keyed dictionary with `station_id`
and the complete `component` code when component corrections are configured.
Static `coda_correction_mag` and `amplitude_correction_mag` default to zero and
are direct magnitude increments, without encoded control flags. Separate
`coda_weight` and `amplitude_weight` are event-aggregation weights from 0 to 10
(default 10); `coda_gain_enabled` is a Boolean (default true).

History rows contain exact `station_id`, direct magnitude increment `value`,
optional `expires_utc` and `source`. Optional `event_weight` overrides the
corresponding channel weight; otherwise that weight is retained. A duration
history row may also set Boolean `gain_enabled`; omitted values retain the
channel's `coda_gain_enabled`. Histories do not force a full station weight.

The first history interval extends into the past; omitted final expiry extends
into the future. Expiries are increasing UTC hours per channel. Intervals are
half open and events may arrive in any order. Expired history returns unknown
correction, zero event weight and `history_expired`, without reusing a static
correction. Unlisted channels use their static correction and weight, or zero
weight when the corresponding `use_unlisted_*` is false. Zero event weight
does not remove an otherwise calculable station magnitude or its correction.

Amplitude-history rows may specify `instrument_type` (0:7). These overrides
must agree across a channel's rows: the override is static for that channel,
not an instrument time history. Apply it before selecting the calibration.

Component rows contain the complete `component`, direct magnitude increment
`value` and optional `source`. Exact matches add cumulatively. Coda component
corrections apply to Md1 and Mt, while station corrections apply to Md1, Md2
and Mt. Amplitude component and station corrections apply to every amplitude
magnitude type. No identity-prefix matching or magnitude quantization is used.
This function returns additive terms, weights and gain selection; it does not
calculate station or event magnitudes.
"""
function magnitude_corrections(time_utc, channel; duration_history=[],
        amplitude_history=[], component_coda=[], component_amplitude=[],
        use_unlisted_coda=true, use_unlisted_amplitude=true)
    event_time = _magcal_utc(time_utc)
    event_hour = floor(event_time,Hour)
    station_id = _magcal_value(channel,:station_id)
    station_id isa AbstractString && !isempty(station_id) ||
        throw(ArgumentError("A complete station_id identifying the channel is required."))
    component = _magcal_value(channel,:component,"")
    component isa AbstractString || throw(ArgumentError("Component must be a string."))
    isempty(component) && (!isempty(component_coda) || !isempty(component_amplitude)) &&
        throw(ArgumentError("Specify the actual component when component corrections are configured."))
    coda = _magcorr_station(event_hour,station_id,channel,duration_history,true,use_unlisted_coda)
    amplitude = _magcorr_station(event_hour,station_id,channel,amplitude_history,false,use_unlisted_amplitude)
    fcm = _magcorr_components(component,component_coda)
    xcm = _magcorr_components(component,component_amplitude)
    instrument_type_override = _magcorr_instrument(station_id,amplitude_history)
    md1 = coda.correction === nothing ? nothing : coda.correction+fcm.total
    mx = amplitude.correction === nothing ? nothing : amplitude.correction+xcm.total
    return (;station_id=String(station_id),component=String(component),
            event_time_utc=time_utc,event_hour_utc=event_hour,
            md1_correction_mag=md1,md2_correction_mag=coda.correction,
            mt_correction_mag=md1,amplitude_correction_mag=mx,
            coda_station_correction_mag=coda.correction,
            amplitude_station_correction_mag=amplitude.correction,
            coda_component_correction_mag=fcm.total,
            amplitude_component_correction_mag=xcm.total,
            coda_event_weight=coda.event_weight,amplitude_event_weight=amplitude.event_weight,
            coda_gain_enabled=coda.gain_enabled,instrument_type_override,
            coda_status=coda.status,amplitude_status=amplitude.status,
            coda_source=coda.source,amplitude_source=amplitude.source,
            duration_history_row=coda.history_row,amplitude_history_row=amplitude.history_row,
            duration_expires_utc=coda.expires_utc,amplitude_expires_utc=amplitude.expires_utc,
            duration_source_label=coda.source_label,amplitude_source_label=amplitude.source_label,
            coda_component_rows=fcm.rows,amplitude_component_rows=xcm.rows,
            coda_component_sources=fcm.sources,amplitude_component_sources=xcm.sources)
end
