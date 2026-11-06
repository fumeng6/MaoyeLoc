# H38/H39. These are the empirical amplitude relations actually used by
# Hypoinverse 1.40, not an arbitrary-instrument waveform response calculation.
const _AMPLITUDE_LOG_RESPONSE = (.288,.432,.561,.680,.786,.891,.983,1.066,1.138,
    1.205,1.276,1.355,1.443,1.535,1.630,1.726,1.822,1.916,
    2.007,2.090,2.145,2.099,1.878,1.546,1.172,.771)
const _AMPLITUDE_NORDQUIST_DISTANCE = [0.,5.,10.,20.,35.,55.,73.,96.,220.,334.,
    445.,600.,750.,1010.,1340.,1600.,1750.,1960.,2230.]
const _AMPLITUDE_NORDQUIST_CORRECTION = (1.4,1.4,1.5,1.7,2.3,2.7,2.85,3.,
    3.65,4.2,4.6,5.,5.3,5.7,6.,6.2,6.3,6.4,6.5)
const _AMPLITUDE_RICHTER_CORRECTION = (1.4,1.4,1.5,1.6,1.7,1.9,2.1,2.3,2.4,
    2.5,2.6,2.7,2.8,2.8,2.8,2.85,2.9,2.9,3.,3.,3.,3.1,3.1,3.2,3.2,
    3.3,3.3,3.4,3.4,3.5,3.5,3.6,3.65,3.7,3.7,3.8,3.8,3.9,3.9,4.,4.,
    4.1,4.1,4.2,4.2,4.3,4.3,4.3,4.4,4.4,4.5,4.5,4.5,4.6,4.6,4.6,4.6,
    4.7,4.7,4.7,4.7,4.8,4.8,4.8,4.8,4.8,4.9,4.9,4.9,4.9,4.9)

function _amplitude_number(value; positive=false, nonnegative=false)
    (value === nothing || value === missing) && return (nothing,:missing)
    value isa Real && !(value isa Bool) || return (nothing,:nonnumeric)
    isfinite(value) || return (nothing,:nonfinite)
    number = Float64(value)
    isfinite(number) || return (nothing,:nonfinite)
    positive && number <= 0 && return (nothing,:nonpositive)
    nonnegative && number < 0 && return (nothing,:negative)
    return (number,:valid)
end

"""
    amplitude_distance_correction(relation, epicentral_distance_km, calculation_depth_km)

Return the distance-dependent magnitude correction -log10(A0) used by
Hypoinverse 1.40, where A0 is the reference amplitude for the chosen empirical
magnitude relation. The correction is in magnitude units; distances are km.
Choices are 1 Eaton, 2 Bakun/Joyner, 3 Richter approximation, 4 Nordquist,
5 Richter's discrete table, 6 CISN/Uhrhammer; 0 explicitly omits this term.
The original slant distance is s=hypot(R,z)+0.01 km, using calculation depth,
not the optional report-depth datum. Relations 4/5 use epicentral distance.
Richter's table uses nearest knots with positive half-ties rounded upwards,
and its final 4.9 value also applies beyond 600 km. Nordquist uses linear
interpolation on its original table for 0<=R<=2224 km. CISN requires
0.1<s<=500 km. Unknown distance terms return no correction, not the old
zero/-9 sentinels. Native coordinates are not rounded to legacy file precision.

Eaton's far branch includes -2.21, missing in the manual's displayed formula.
The Richter far branch is -3.38+1.5*log10(s^2), the expanded source expression.
The response/distance functions are empirical magnitude scales, not a model
of three-dimensional attenuation or a location-accuracy estimate.
"""
function amplitude_distance_correction(relation,epicentral_distance_km,calculation_depth_km)
    relation isa Integer && !(relation isa Bool) && 0 <= relation <= 6 ||
        throw(ArgumentError("Amplitude distance relation must be an integer from 0 to 6."))
    distance,distance_status = _amplitude_number(epicentral_distance_km;nonnegative=true)
    depth,depth_status = _amplitude_number(calculation_depth_km)
    slant = distance === nothing || depth === nothing ? nothing : hypot(distance,depth)+.01
    metadata = (;relation=Int(relation),epicentral_distance_km=distance,
        calculation_depth_km=depth,slant_distance_km=slant)
    unavailable(reason) = (;metadata...,correction_mag=nothing,status=:not_calculated,
        reason,branch=nothing,table_index=nothing,terms=nothing)
    distance_status == :valid || return unavailable(Symbol(distance_status,"_distance"))
    depth_status == :valid || return unavailable(Symbol(depth_status,"_calculation_depth"))
    isfinite(slant) || return unavailable(:nonfinite_slant_distance)
    base,near,polynomial = 0.0,0.0,0.0
    branch,table_index = :none,nothing
    if relation == 1
        branch = slant < 185.3 ? :near : :far
        base = branch == :near ? .821*log10(slant)+.00405*slant+.955 : 2.55*log10(slant)-2.21
        distance < 70 && (near = -.09*sin(.07*(distance-25)))
    elseif relation == 2
        branch,base = :bakun_joyner,log10(slant)+.003*slant+.7
    elseif relation == 3
        l = 2*log10(slant)
        branch = l <= 4.6 ? :near : :far
        base = branch == :near ? -.15+.8*l : -3.38+1.5*l
    elseif relation == 4
        distance <= 2224 || return unavailable(:outside_nordquist_distance_range)
        i = searchsortedlast(_AMPLITUDE_NORDQUIST_DISTANCE,distance)
        lo,hi = _AMPLITUDE_NORDQUIST_DISTANCE[i],_AMPLITUDE_NORDQUIST_DISTANCE[i+1]
        a,b = _AMPLITUDE_NORDQUIST_CORRECTION[i],_AMPLITUDE_NORDQUIST_CORRECTION[i+1]
        base = a+(b-a)*(distance-lo)/(hi-lo)
        branch,table_index = :linear_table,i
    elseif relation == 5
        table_index = distance >= 600 ? 71 : distance < 100 ?
            floor(Int,distance/5+.5)+1 : floor(Int,distance/10+.5)+11
        base = _AMPLITUDE_RICHTER_CORRECTION[table_index]
        branch = :nearest_table
    elseif relation == 6
        .1 < slant <= 500 || return unavailable(:outside_cisn_slant_distance_range)
        if slant <= 8
            b = (2.6182-1.5429)/(log10(60.)-log10(8.))
            base = 1.5429+b*(log10(slant)-log10(8.))
            branch = :short_distance_extrapolation
        else
            x = -1+2*(log10(slant)-log10(8.))/(log10(500.)-log10(8.))
            previous,current = 1.0,x
            for coefficient in (.056,-.031,-.053,-.080,-.028,.015)
                polynomial += coefficient*current
                previous,current = current,2*x*current-previous
            end
            base = 1.11*log10(slant)+.00189*slant+.591+.0054
            branch = :chebyshev
        end
    end
    terms = (;base_mag=base,near_epicentral_mag=near,chebyshev_mag=polynomial)
    correction_mag = sum(values(terms))
    isfinite(correction_mag) || return unavailable(:nonfinite_distance_correction)
    return (;metadata...,correction_mag,status=:valid,reason=:calculated,branch,table_index,terms)
end

"""
    amplitude_magnitude(amplitude, period_s, epicentral_distance_km, calculation_depth_km;
        amplitude_units=0, amplitude_type=0, default_period_s=.2,
        calibration, corrections=nothing, distance_relation=1)

Estimate one station magnitude from measured amplitude and period, instrument
calibration, and an empirical distance correction, using Hypoinverse 1.40's
ML/MX definitions. Return the magnitude and its additive correction terms.
Unit codes 0/1 are peak-to-peak/zero-to-peak mm. Codes 2/3 are peak-to-peak
Menlo/Earthworm and HVO counts, converted by .04883/.012207 mm/count; code 4
is Menlo counts expressed in hundreds, converted by 4.883 mm/input unit.
Code 4 uses the same electronic CAL as code 2, not a different sensor gain.
For H40 digital-F calibration, F remains microvolts per actual count; choose
the Menlo CAL convention (1.382/F), independent of the input's hundred-count
encoding. The old literal conversion factors are retained, not replaced by
their slightly different ideal rational approximations.

Amplitude type 0 infers ML for instruments 0/2 and MX for 1/3:7. Type 1 means
amplitude already expressed as a Wood-Anderson (WA) recording (ML), type 2
velocity (MX), type 3 acceleration and type
4 explicit no-magnitude. A velocity amplitude requires an actual velocity
instrument motor constant; the old type-0 invalid array access is not copied.
The H40 `calibration` result provides CAL and instrument type. Apply H41's
static XMC instrument override before obtaining that calibration. Its zero,
invalid and expired states never supply a fictitious reference gain.

Positive phase period takes priority. Missing/zero/negative phase period
uses the supplied station default. Non-numeric/nonfinite values are invalid.
The original global period range .02:6.3 seconds applies to ML and MX.
The source's log-response table is interpolated in log frequency. Its last
eight values differ from the older table printed on manual p61. This is the
NCSN analog L4C relative-WA approximation used by the original for all
velocity instruments; the physically supported common-response band is
3--10 Hz, not arbitrary response deconvolution. For Sprengnether type 2,
the source adds .41+.56*log10(.2/T), whereas manual p64 prints a minus sign;
the fit's stated period range is .1--1.9 s. Both physical-range flags are
returned without discarding the original allowed-period calculations.

H41 station XMC and cumulative component XCM are added to all types. A zero
event weight does not remove the station result. Valid negative magnitudes
are preserved. Missing or invalid measurements return magnitude=nothing,
with a distinct reason, and are never confused with a true zero magnitude.
Event-set selection, observation weights and event aggregation are separate.
"""
function amplitude_magnitude(amplitude,period_s,epicentral_distance_km,calculation_depth_km;
        amplitude_units=0,amplitude_type=0,default_period_s=.2,
        calibration,corrections=nothing,distance_relation=1)
    corrections === nothing && (corrections = (;))
    distance = amplitude_distance_correction(distance_relation,epicentral_distance_km,calculation_depth_km)
    amp,amplitude_status = _amplitude_number(amplitude;positive=true)
    period,period_status = _amplitude_number(period_s)
    period_source = :phase
    if period_status == :missing || (period_status == :valid && period <= 0)
        period,period_status = _amplitude_number(default_period_s;positive=true)
        period_source = :station_default
    end
    instrument = calibration === nothing ? nothing : _magcal_value(calibration,:instrument_type)
    cal,cal_status = _amplitude_number(calibration === nothing ? nothing :
        _magcal_value(calibration,:calibration_factor))
    calibration_status = calibration === nothing ? :missing_calibration :
        _magcal_value(calibration,:status,:unspecified)
    calibration_reason = calibration === nothing ? :calibration_not_supplied :
        _magcal_value(calibration,:reason,:unspecified)
    station,station_status = _amplitude_number(_magcal_value(corrections,:amplitude_station_correction_mag,0.0))
    component,component_status = _amplitude_number(_magcal_value(corrections,:amplitude_component_correction_mag,0.0))
    event_weight = _magcal_value(corrections,:amplitude_event_weight,10.0)
    unit_valid = amplitude_units isa Integer && !(amplitude_units isa Bool) && 0 <= amplitude_units <= 4
    type_valid = amplitude_type isa Integer && !(amplitude_type isa Bool) && 0 <= amplitude_type <= 4
    instrument_valid = instrument isa Integer && !(instrument isa Bool) && 0 <= instrument <= 7
    magnitude_type = !type_valid || !instrument_valid || amplitude_type >= 3 ? nothing :
        amplitude_type == 1 ? :ml : amplitude_type == 2 ? :mx : instrument in (0,2) ? :ml : :mx
    unit_factor = unit_valid ? (1.,2.,.04883,.012207,4.883)[amplitude_units+1] : nothing
    peak_to_peak = amp === nothing || unit_factor === nothing ? nothing : amp*unit_factor
    metadata = (;magnitude_type,instrument_type=instrument,amplitude_units,amplitude_type,
        amplitude_status,period_status,period_s=period,period_source,
        unit_factor_mm=unit_factor,amplitude_peak_to_peak_mm=peak_to_peak,
        calibration_factor=cal,calibration_status,calibration_reason,
        amplitude_event_weight=event_weight,distance)
    unavailable(reason) = (;metadata...,magnitude=nothing,status=:not_calculated,reason,
        terms=nothing,log10_response=nothing,motor_constant_v_per_cm_s=nothing,
        velocity_response_in_3_to_10_hz=nothing,sprengnether_period_in_0p1_to_1p9_s=nothing)
    type_valid || return unavailable(:invalid_amplitude_type)
    amplitude_type == 3 && return unavailable(:acceleration_amplitude_not_supported)
    amplitude_type == 4 && return unavailable(:magnitude_explicitly_disabled)
    amplitude_status == :valid || return unavailable(Symbol(amplitude_status,"_amplitude"))
    unit_valid || return unavailable(:invalid_amplitude_units)
    isfinite(peak_to_peak) || return unavailable(:nonfinite_converted_amplitude)
    calibration === nothing && return unavailable(:calibration_not_supplied)
    instrument_valid || return unavailable(:unavailable_instrument_type)
    override = _magcal_value(corrections,:instrument_type_override)
    override === nothing || override == instrument ||
        return unavailable(:instrument_override_not_applied_to_calibration)
    _magcal_value(calibration,:amplitude_calibrated,true) === true ||
        return unavailable(calibration_reason)
    cal_status == :valid || return unavailable(Symbol(cal_status,"_calibration"))
    cal == 0 && return unavailable(:zero_calibration)
    cal > 0 || return unavailable(:negative_calibration)
    period_status == :valid || return unavailable(Symbol(period_status,"_period"))
    .02 <= period <= 6.3 || return unavailable(:outside_allowed_period_range)
    distance.status == :valid || return unavailable(distance.reason)
    station_status == :valid || return unavailable(:unavailable_station_correction)
    component_status == :valid || return unavailable(:unavailable_component_correction)
    motor,log_response = nothing,0.0
    velocity_band,sprengnether_band = nothing,nothing
    if magnitude_type == :mx
        motor = _MAGCAL_MOTOR[instrument+1]
        motor === nothing && return unavailable(:velocity_amplitude_requires_velocity_instrument)
        x = 10*log10(1/period)+9
        i = floor(Int,x)
        log_response = _AMPLITUDE_LOG_RESPONSE[i]+(x-i)*
            (_AMPLITUDE_LOG_RESPONSE[i+1]-_AMPLITUDE_LOG_RESPONSE[i])
        velocity_band = 3 <= 1/period <= 10
    end
    sprengnether = instrument == 2 ? .41+.56*log10(.2/period) : 0.0
    instrument == 2 && (sprengnether_band = .1 <= period <= 1.9)
    terms = (;log_amplitude_calibration_mag=log10(amp)+log10(unit_factor)-log10(2.)-log10(cal),
        response_mag=-log_response,motor_mag=motor === nothing ? 0.0 : -log10(motor),
        sprengnether_mag=sprengnether,distance_mag=distance.correction_mag,
        station_mag=station,component_mag=component)
    magnitude = sum(values(terms))
    isfinite(magnitude) || return unavailable(:nonfinite_magnitude)
    return (;metadata...,magnitude,status=:valid,reason=:calculated,terms,
        log10_response=log_response,motor_constant_v_per_cm_s=motor,
        velocity_response_in_3_to_10_hz=velocity_band,
        sprengnether_period_in_0p1_to_1p9_s=sprengnether_band)
end
