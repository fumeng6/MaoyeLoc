# H42: original independent coda/amplitude event sets. Calculation and
# selection are separate; a zero event weight does not erase a station value.
function _magsel_components(selection,component)
    selection == "all" && return true
    selection == "none" && return false
    selection isa AbstractVector && all(c -> c isa AbstractString && !isempty(c),selection) ||
        throw(ArgumentError("Magnitude components must be all, none, or a list of complete component codes."))
    return component in selection
end

function _magsel_instruments(selection,instrument)
    selection == "all" && return true
    selection == "none" && return false
    selection isa AbstractVector && all(i -> i isa Integer && !(i isa Bool) && 0 <= i <= 7,selection) ||
        throw(ArgumentError("Magnitude instrument_types must be all, none, or a list of integers from 0 to 7."))
    return instrument in selection
end

function _magsel_region(settings,latitude_deg,longitude_deg)
    region = _magcal_value(settings,:fmag2_region)
    region === nothing && return :disabled
    bounds = map((:lat_min_deg,:lat_max_deg,:lon_min_deg,:lon_max_deg)) do key
        value = _magcal_value(region,key)
        value isa Real && !(value isa Bool) && isfinite(value) ||
            throw(ArgumentError("FMAG2 region requires finite $key in degrees, with east-positive longitude."))
        Float64(value)
    end
    lat_min,lat_max,lon_min,lon_max = bounds
    -90 <= lat_min <= lat_max <= 90 && -180 <= lon_min <= lon_max <= 180 ||
        throw(ArgumentError("FMAG2 region requires ordered latitude [-90,90] and east-longitude [-180,180] bounds without a date-line crossing."))
    latitude_deg isa Real && !(latitude_deg isa Bool) && isfinite(latitude_deg) && -90 <= latitude_deg <= 90 &&
        longitude_deg isa Real && !(longitude_deg isa Bool) && isfinite(longitude_deg) && -180 <= longitude_deg <= 180 ||
        throw(ArgumentError("FMAG2 region requires the event's actual latitude_deg and east-positive longitude_deg."))
    return lat_min <= latitude_deg <= lat_max && lon_min <= longitude_deg <= lon_max ? :inside : :outside
end

function _magsel_slot(slot,label,magnitude_type,value,matched,selection_reason,quality_code,station_weight)
    label isa AbstractString && length(label) <= 1 ||
        throw(ArgumentError("Magnitude labels must be empty or one character."))
    magnitude = value === nothing ? nothing : _magcal_value(value,:magnitude)
    numeric = magnitude isa Real && !(magnitude isa Bool) && isfinite(magnitude)
    calculated = value !== nothing && _magcal_value(value,:status,:valid) == :valid
    available = calculated && numeric
    selected = available && matched
    status,reason = if !available
        (:not_calculated,value === nothing ? :station_magnitude_not_calculated :
            !calculated ? _magcal_value(value,:reason,:unavailable_station_magnitude) : :invalid_station_magnitude)
    elseif !matched
        (:not_selected,selection_reason)
    else
        (:selected,:selected)
    end
    return (;slot,label=String(label),magnitude_type,magnitude=numeric ? magnitude : nothing,
        selected,status,reason,quality_code,station_weight)
end

"""
    magnitude_selection(row, component, settings;
                        latitude_deg=nothing, longitude_deg=nothing)

Choose which of four separately configured Hypoinverse event-magnitude sets
may use one prepared measurement. These are separate selection rules, not a
claim of statistical independence. `row` is a magnitude_preparation row, with `measurement`,
`calibration`, `corrections`, `coda` and `amplitude`. `settings` is exactly the
contents of [magnitudes.selection], not the complete magnitudes table.
Identity/component matching is exact on the complete native component code.

Tables `coda1`/`coda2` choose `magnitude_type="md1"|"mt"|"md2"`,
`components="all"|"none"|[complete codes]`, and `label=""` or one character.
Defaults are Md1/all for coda1, Md1/none for coda2, and empty labels. The two
sets independently select from the calculated Md1/Mt/Md2 results. They may
use different values for the same measurement; no shared station-display
slot or switch of the printed order replaces either calculation array.

Tables `amplitude1`/`amplitude2` have the same component and label selection,
defaulting to all/none respectively. `amplitude_selection_mode="component"`
(default) chooses by component; "instrument" instead chooses each table's
`instrument_types` (default [], or "all"/"none"/a list of 0:7). These two
selectors are alternatives, not an AND. The effective instrument comes from
the calibration after H41's static XMC override. Each amplitude table also
applies its independent `magnitude_type="all"|"ml"|"mx"` filter (default all).

An optional `fmag2_region` table contains `lat_min_deg`, `lat_max_deg`,
`lon_min_deg`, `lon_max_deg`. Longitude is explicitly east-positive and the
box cannot cross the date line. It requires the actual event latitude and
longitude. Inside, including every boundary, all components are enabled for
coda2; outside all are disabled. This replaces coda2's component selection,
but never its magnitude formula, quality code, or weight. This is the source
MAR behavior, not an extra AND with the component list.

Return `coda1`, `coda2`, `amplitude1`, `amplitude2` plus `region_status`.
Each returned entry retains its calculated magnitude, selection status/reason, label,
actual magnitude type, original measurement quality code (missing field
defaults to 0), and H41 station event weight (default 10). Availability and
selector membership determine `selected`. Zero station weight or quality
codes >=4 do not erase the station value or selector membership; numerical
weighting/aggregation is H43. Valid negative or zero magnitudes are available.
Unavailable station calculations remain distinct from valid unselected ones.
"""
function magnitude_selection(row,component,settings;latitude_deg=nothing,longitude_deg=nothing)
    component isa AbstractString || throw(ArgumentError("Magnitude selection requires a component string."))
    mode = _magcal_value(settings,:amplitude_selection_mode,"component")
    mode in ("component","instrument") ||
        throw(ArgumentError("amplitude_selection_mode must be component or instrument."))
    region_status = _magsel_region(settings,latitude_deg,longitude_deg)
    measurement = _magcal_value(row,:measurement,(;))
    corrections = _magcal_value(row,:corrections,(;))
    corrections === nothing && (corrections = (;))
    coda_quality = _magcal_value(measurement,:coda_quality_code,0)
    amplitude_quality = _magcal_value(measurement,:amplitude_quality_code,0)
    coda_weight = _magcal_value(corrections,:coda_event_weight,10.0)
    amplitude_weight = _magcal_value(corrections,:amplitude_event_weight,10.0)
    coda = _magcal_value(row,:coda)
    amplitude = _magcal_value(row,:amplitude)
    calibration = _magcal_value(row,:calibration)
    instrument = calibration === nothing ? nothing : _magcal_value(calibration,:instrument_type)
    coda_slots = map((:coda1,:coda2),("all","none")) do slot,default_components
        table = _magcal_value(settings,slot,(;))
        kind = _magcal_value(table,:magnitude_type,"md1")
        kind in ("md1","mt","md2") ||
            throw(ArgumentError("Coda selection magnitude_type must be md1, mt, or md2."))
        matched,reason = if slot == :coda2 && region_status != :disabled
            (region_status == :inside,:outside_fmag2_region)
        else
            (_magsel_components(_magcal_value(table,:components,default_components),component),:component_not_selected)
        end
        value = coda === nothing ? nothing : _magcal_value(coda,Symbol(kind))
        _magsel_slot(slot,_magcal_value(table,:label,""),Symbol(kind),value,matched,reason,coda_quality,coda_weight)
    end
    amplitude_slots = map((:amplitude1,:amplitude2),("all","none")) do slot,default_components
        table = _magcal_value(settings,slot,(;))
        filter = _magcal_value(table,:magnitude_type,"all")
        filter in ("all","ml","mx") ||
            throw(ArgumentError("Amplitude selection magnitude_type must be all, ml, or mx."))
        matched = mode == "component" ?
            _magsel_components(_magcal_value(table,:components,default_components),component) :
            _magsel_instruments(_magcal_value(table,:instrument_types,[]),instrument)
        reason = mode == "component" ? :component_not_selected : :instrument_not_selected
        kind = amplitude === nothing ? nothing : _magcal_value(amplitude,:magnitude_type)
        if matched && filter != "all" && kind != Symbol(filter)
            matched,reason = false,:magnitude_type_not_selected
        end
        _magsel_slot(slot,_magcal_value(table,:label,""),kind,amplitude,matched,reason,amplitude_quality,amplitude_weight)
    end
    return (;coda1=coda_slots[1],coda2=coda_slots[2],
        amplitude1=amplitude_slots[1],amplitude2=amplitude_slots[2],region_status)
end
