# H35--H37: station coda magnitudes. Time is seconds, R is epicentral km,
# and z is the calculation depth, before any reporting-datum conversion.
function _coda_number(value; positive=false, nonnegative=false)
    (value === nothing || value === missing) && return (nothing,:missing)
    value isa Real && !(value isa Bool) || return (nothing,:nonnumeric)
    isfinite(value) || return (nothing,:nonfinite)
    number = Float64(value)
    isfinite(number) || return (nothing,:nonfinite)
    positive && number <= 0 && return (nothing,:nonpositive)
    nonnegative && number < 0 && return (nothing,:negative)
    return (number,:valid)
end

function _coda_coefficient(parameters,key,default)
    value,status = _coda_number(_magcal_value(parameters,key,default))
    status == :valid || throw(ArgumentError("Coda parameter $key must be finite."))
    return value
end

function _coda_duration_parameters(parameters)
    break_s = _coda_coefficient(parameters,:duration_break_s,210.0)
    gain = _coda_coefficient(parameters,:gain_coefficient,0.0)
    branches = map((:lower,:upper),(-5.2,-.905),(3.89,2.026)) do branch,constant,log_coefficient
        row = _magcal_value(parameters,branch,(;))
        (; constant_mag=_coda_coefficient(row,:constant_mag,constant),
           log10_coefficient=_coda_coefficient(row,:log10_coefficient,log_coefficient),
           depth_coefficient_mag_per_km=_coda_coefficient(row,:depth_coefficient_mag_per_km,.013),
           distance_coefficient_mag_per_km=_coda_coefficient(row,:distance_coefficient_mag_per_km,.0037),
           linear_coefficient_mag_per_s=_coda_coefficient(row,:linear_coefficient_mag_per_s,0.0))
    end
    return (; duration_break_s=break_s,gain_coefficient=gain,lower=branches[1],upper=branches[2])
end

function _coda_parameters(parameters,component)
    md1 = _coda_duration_parameters(_magcal_value(parameters,:md1,(;)))
    md2 = _coda_duration_parameters(_magcal_value(parameters,:md2,(;)))
    row = _magcal_value(parameters,:mt,(;))
    mt = (; constant_mag=_coda_coefficient(row,:constant_mag,-1.03),
        log10_coefficient=_coda_coefficient(row,:log10_coefficient,2.10),
        log10_squared_coefficient=_coda_coefficient(row,:log10_squared_coefficient,0.0),
        linear_coefficient_mag_per_s=_coda_coefficient(row,:linear_coefficient_mag_per_s,.00268),
        depth_coefficient_mag_per_km=_coda_coefficient(row,:depth_coefficient_mag_per_km,0.0),
        gain_coefficient=_coda_coefficient(row,:gain_coefficient,1.0))
    row = _magcal_value(parameters,:additional_corrections,(;))
    additional = (; near_distance_break_km=_coda_coefficient(row,:near_distance_break_km,0.0),
        near_distance_coefficient_mag_per_km=_coda_coefficient(row,:near_distance_coefficient_mag_per_km,0.0),
        far_distance_break_km=_coda_coefficient(row,:far_distance_break_km,0.0),
        far_distance_coefficient_mag_per_km=_coda_coefficient(row,:far_distance_coefficient_mag_per_km,0.0),
        depth_break_km=_coda_coefficient(row,:depth_break_km,0.0),
        depth_coefficient_mag_per_km=_coda_coefficient(row,:depth_coefficient_mag_per_km,0.0))
    selection = _magcal_value(parameters,:gain_components,"all")
    allowed = if selection == "all"
        true
    elseif selection == "none"
        false
    elseif selection isa AbstractVector && all(s -> s isa AbstractString && !isempty(s),selection)
        component isa AbstractString && !isempty(component) ||
            throw(ArgumentError("Coda gain component selection requires the actual complete component."))
        component in selection
    else
        throw(ArgumentError("Coda gain_components must be all, none, or a list of complete component codes."))
    end
    return (;md1,md2,mt,additional,gain_component_allowed=allowed)
end

_coda_unavailable(reason,branch=nothing) =
    (;magnitude=nothing,status=:not_calculated,reason,duration_branch=branch,terms=nothing)

function _coda_result(terms,branch=nothing)
    magnitude = sum(values(terms))
    isfinite(magnitude) || return _coda_unavailable(:nonfinite_magnitude,branch)
    return (;magnitude,status=:valid,reason=:calculated,duration_branch=branch,terms)
end

function _coda_duration_magnitude(duration,distance,depth,parameters,gain,station,component,additional)
    branch = duration < parameters.duration_break_s ? :lower : :upper
    row = getproperty(parameters,branch)
    gain === nothing && parameters.gain_coefficient != 0 &&
        return _coda_unavailable(:unavailable_gain_correction,branch)
    terms = (; constant_mag=row.constant_mag,
        log_mag=row.log10_coefficient*log10(duration),log_squared_mag=0.0,
        linear_time_mag=row.linear_coefficient_mag_per_s*duration,
        depth_mag=row.depth_coefficient_mag_per_km*depth,
        distance_mag=row.distance_coefficient_mag_per_km*distance,
        gain_mag=parameters.gain_coefficient == 0 ? 0.0 : parameters.gain_coefficient*gain,
        station_mag=station,component_mag=component,
        near_distance_mag=additional.near_distance_mag,
        far_distance_mag=additional.far_distance_mag,
        additional_depth_mag=additional.additional_depth_mag)
    return _coda_result(terms,branch)
end

"""
    coda_magnitudes(duration_s, calculated_p_time_s, distance_km, calculation_depth_km;
                    parameters=Dict(), gain_correction=nothing,
                    corrections=nothing, component="")

Estimate station magnitudes from the measured time between the P arrival and
the end of the coda (F-P). Hypoinverse 1.40 calls the two duration relations
DUR (Md1) and DUB (Md2); TAU (Mt) uses that duration plus the model P travel time.
`duration_s` is the externally measured F-P duration; `calculated_p_time_s` is
the model P travel time, excluding origin time, known station delay, and any
empirical S intercept. It is computed as P even for an S-only station record.
`distance_km` is horizontal source-station distance and depth is calculation
z, not the optional reporting-datum depth. No observed P arrival is substituted.

Parameters accept string-keyed TOML dictionaries or NamedTuples. `md1` and
`md2` each contain `duration_break_s` (210), `gain_coefficient` (0), and `lower`
and `upper` tables. Each branch has `constant_mag` (-5.2/-.905),
`log10_coefficient` (3.89/2.026), `depth_coefficient_mag_per_km` (.013),
`distance_coefficient_mag_per_km` (.0037), and `linear_coefficient_mag_per_s`
(0). Duration equal to its break uses the upper branch. A zero or negative
break selects the upper branch for all valid positive durations. The native Md2 linear
and gain defaults are explicitly zero: the reference BLOCK DATA does not
initialize those three DUB coefficients.

The `mt` table contains `constant_mag` (-1.03), `log10_coefficient` (2.10),
`log10_squared_coefficient` (0), `linear_coefficient_mag_per_s` (.00268),
`depth_coefficient_mag_per_km` (0), and `gain_coefficient` (1). Its lapse time
is tau=calculated_p_time_s+duration_s, required to be strictly greater than .5 s.

`additional_corrections` contains `near_distance_break_km`,
`near_distance_coefficient_mag_per_km`, `far_distance_break_km`,
`far_distance_coefficient_mag_per_km`, `depth_break_km`, and
`depth_coefficient_mag_per_km`, all default zero. DU2 adds the near term
`c_near*(R-R_near)` for R<R_near, the far term `c_far*(R-R_far)` for R>R_far,
and `c_depth*(z-z_break)` for z>z_break, independently, to Md1 and Mt only.
The three corrections can overlap; empirical breakpoints need only be finite.

Pass H40's raw `coda_gain_correction` as `gain_correction` and H41's result as
`corrections`. Its FMC station correction applies to all three magnitudes,
FCM component correction to Md1/Mt only, and `coda_gain_enabled=false`
suppresses gain for all three. With no H41 result the station/component
corrections default to zero and the event weight to 10. Manually supplied
corrections can use the same `coda_station_correction_mag`,
`coda_component_correction_mag`, `coda_gain_enabled`, and `coda_event_weight`
fields. `parameters.gain_components` is "all" (default), "none", or a list of
complete component codes. This DUG selection gates only Md1/Md2, not Mt.
An unavailable calibration is acceptable only when its actual gain term is
suppressed or its coefficient is zero; it is never silently replaced by zero.

Return each station magnitude with additive terms, status and reason. Missing
or invalid measurements produce `nothing`, distinct from real zero or negative
magnitudes. There is no VNOMAG clipping, observation-quality selection,
multi-set selection, or event averaging. The H41 event weight is retained as
metadata; weight zero does not remove an otherwise valid station magnitude.
"""
function coda_magnitudes(duration_s,calculated_p_time_s,distance_km,calculation_depth_km;
        parameters=Dict(),gain_correction=nothing,corrections=nothing,component="")
    configured = _coda_parameters(parameters,component)
    corrections === nothing && (corrections = (;))
    duration,duration_status = _coda_number(duration_s;positive=true)
    ptime,calculated_p_time_status = _coda_number(calculated_p_time_s;nonnegative=true)
    distance,distance_status = _coda_number(distance_km;nonnegative=true)
    depth,depth_status = _coda_number(calculation_depth_km)
    raw_gain,gain_status = _coda_number(gain_correction)
    station,station_status = _coda_number(_magcal_value(corrections,:coda_station_correction_mag,0.0))
    component_term,component_status = _coda_number(_magcal_value(corrections,:coda_component_correction_mag,0.0))
    gain_enabled = _magcal_value(corrections,:coda_gain_enabled,true)
    gain_enabled isa Bool || gain_enabled === nothing ||
        throw(ArgumentError("Coda gain-enabled flag must be true, false, or unknown."))
    event_weight = _magcal_value(corrections,:coda_event_weight,10.0)
    effective_gain = gain_enabled === false ? 0.0 : gain_enabled === true ? raw_gain : nothing
    duration_gain = configured.gain_component_allowed ? effective_gain : 0.0
    tau = duration === nothing || ptime === nothing ? nothing : duration+ptime
    unavailable = if duration_status != :valid
        Symbol(duration_status,"_duration")
    elseif distance_status != :valid
        Symbol(distance_status,"_distance")
    elseif depth_status != :valid
        Symbol(depth_status,"_calculation_depth")
    elseif station_status != :valid
        :unavailable_station_correction
    else
        nothing
    end
    if unavailable !== nothing
        md1 = md2 = mt = _coda_unavailable(unavailable)
    else
        row = configured.additional
        additional = (;
            near_distance_mag=distance < row.near_distance_break_km ?
                row.near_distance_coefficient_mag_per_km*(distance-row.near_distance_break_km) : 0.0,
            far_distance_mag=distance > row.far_distance_break_km ?
                row.far_distance_coefficient_mag_per_km*(distance-row.far_distance_break_km) : 0.0,
            additional_depth_mag=depth > row.depth_break_km ?
                row.depth_coefficient_mag_per_km*(depth-row.depth_break_km) : 0.0)
        no_additional = (;near_distance_mag=0.0,far_distance_mag=0.0,additional_depth_mag=0.0)
        md1 = component_status == :valid ?
            _coda_duration_magnitude(duration,distance,depth,configured.md1,duration_gain,station,component_term,additional) :
            _coda_unavailable(:unavailable_component_correction)
        md2 = _coda_duration_magnitude(duration,distance,depth,configured.md2,duration_gain,station,0.0,no_additional)
        mt = if calculated_p_time_status != :valid
            _coda_unavailable(Symbol(calculated_p_time_status,"_calculated_p_time"))
        elseif !isfinite(tau)
            _coda_unavailable(:nonfinite_tau)
        elseif tau <= .5
            _coda_unavailable(:tau_not_above_half_second)
        elseif component_status != :valid
            _coda_unavailable(:unavailable_component_correction)
        elseif effective_gain === nothing && configured.mt.gain_coefficient != 0
            _coda_unavailable(:unavailable_gain_correction)
        else
            row = configured.mt
            tlog = log10(tau)
            _coda_result((; constant_mag=row.constant_mag,
                log_mag=row.log10_coefficient*tlog,log_squared_mag=row.log10_squared_coefficient*tlog^2,
                linear_time_mag=row.linear_coefficient_mag_per_s*tau,
                depth_mag=row.depth_coefficient_mag_per_km*depth,distance_mag=0.0,
                gain_mag=row.gain_coefficient == 0 ? 0.0 : row.gain_coefficient*effective_gain,
                station_mag=station,component_mag=component_term,
                near_distance_mag=additional.near_distance_mag,
                far_distance_mag=additional.far_distance_mag,
                additional_depth_mag=additional.additional_depth_mag))
        end
    end
    return (;md1,md2,mt,tau_s=tau,duration_status,calculated_p_time_status,distance_status,
        calculation_depth_status=depth_status,gain_status,
        distance_km=distance,calculation_depth_km=depth,gain_component_allowed=configured.gain_component_allowed,
        raw_coda_gain_correction=raw_gain,effective_coda_gain_correction=effective_gain,
        duration_gain_correction=duration_gain,coda_event_weight=event_weight)
end
