"Resolve calibration, corrections and requested station coda/amplitude magnitudes."
function magnitude_preparation(events,stations,settings;time_reference_utc,model=nothing,frame=nothing)
    endswith(time_reference_utc,"Z") || throw(ArgumentError("Magnitude reference time must be declared UTC with Z."))
    channels=get(settings,"channels",Dict())
    measurements=get(settings,"measurements",[])
    channels isa AbstractDict && measurements isa AbstractVector ||
        throw(ArgumentError("magnitudes requires channels tables and a measurements array."))
    calculate_coda=haskey(settings,"coda_magnitudes") || any(haskey(row,"duration_s") for row in measurements)
    calculate_amplitude=haskey(settings,"amplitude_magnitudes") || any(haskey(row,"amplitude") for row in measurements)
    calculate_events=any(haskey(settings,key) for key in ("selection","aggregation","preferred","display","external","legacy_external_amplitude"))
    calculate_coda && model===nothing && throw(ArgumentError("Station coda magnitudes require the travel-time model."))
    by_id=Dict(event.event_id=>event for event in events)
    reference=DateTime(time_reference_utc[1:end-1])
    rows=NamedTuple[]
    for (input_row,measurement) in enumerate(measurements)
        id=get(measurement,"event_id",nothing)
        station=get(measurement,"station_id",nothing)
        haskey(by_id,id) && haskey(stations,station) && haskey(channels,station) ||
            throw(ArgumentError("Magnitude measurement needs an existing event, station and explicit channel definition."))
        channel=merge(Dict("station_id"=>station),channels[station])
        channel["station_id"]==station || throw(ArgumentError("Magnitude channel table key must match station_id."))
        time=_origin_utc(reference,by_id[id].origin_s)
        corrections=magnitude_corrections(time,channel;
            duration_history=get(settings,"duration_history",[]),amplitude_history=get(settings,"amplitude_history",[]),
            component_coda=get(settings,"component_coda",[]),component_amplitude=get(settings,"component_amplitude",[]),
            use_unlisted_coda=get(settings,"use_unlisted_coda",true),use_unlisted_amplitude=get(settings,"use_unlisted_amplitude",true))
        selected=corrections.instrument_type_override===nothing ? channel :
            merge(channel,Dict("instrument_type"=>corrections.instrument_type_override))
        units=get(measurement,"amplitude_units",get(channel,"amplitude_units",0))
        if haskey(selected,"digital_system_factor_uv_per_count") && units in (2,3,4)
            selected=merge(selected,Dict("amplitude_units"=>units))
        end
        calibration=magnitude_calibration(time,selected;phase_cal=get(measurement,"phase_cal",nothing),
            calibration_history=get(settings,"calibration_history",[]),attenuation_history=get(settings,"attenuation_history",[]))
        gain=corrections.coda_gain_enabled===true ? calibration.coda_gain_correction :
            corrections.coda_gain_enabled===false ? 0. : nothing
        coda=nothing
        source=by_id[id].source_km
        station_xyz=stations[station]
        distance=hypot(source[1]-station_xyz[1],source[2]-station_xyz[2])
        if calculate_coda
            # Mt uses the model P travel time without origin time or station
            # delay, even for a record containing only S or coda measurements.
            # A missing P ray does not invalidate a duration-based Md.
            ptime,pstatus=try
                ray=travel_time(model,source,stations,station,"P")
                (get(ray,:ray_time_s,ray.time_s),:calculated)
            catch exception
                if exception isa DomainError && Tuple(source)==Tuple(station_xyz)
                    (0.0,:coincident_endpoints)
                elseif exception isa DomainError || exception isa ArgumentError
                    (nothing,:unavailable_model_p_time)
                else
                    rethrow()
                end
            end
            values=coda_magnitudes(get(measurement,"duration_s",nothing),ptime,distance,source[3];
                parameters=get(settings,"coda_magnitudes",Dict()),gain_correction=calibration.coda_gain_correction,
                corrections,component=get(channel,"component",""))
            coda=merge(values,(calculated_p_time_s=ptime,model_p_time_status=pstatus,
                coda_quality_code=get(measurement,"coda_quality_code",nothing)))
        end
        amplitude_settings=get(settings,"amplitude_magnitudes",Dict())
        distance_relation=get(get(amplitude_settings,"component_distance_relations",Dict()),get(channel,"component",""),
            get(amplitude_settings,"distance_relation",1))
        amplitude=calculate_amplitude ? amplitude_magnitude(get(measurement,"amplitude",nothing),
            get(measurement,"period_s",nothing),distance,source[3];amplitude_units=units,
            amplitude_type=get(measurement,"amplitude_type",0),default_period_s=get(channel,"default_period_s",.2),
            calibration,corrections,distance_relation) : nothing
        push!(rows,(;input_row,event_id=String(id),station_id=String(station),event_time_utc=time,
            location_status=get(by_id[id],:status,:supplied_catalog),measurement,calibration,corrections,
            effective_coda_gain_correction=gain,coda,amplitude))
    end
    aggregated=calculate_events ? event_magnitudes(events,rows,settings;frame) : nothing
    (;rows,event_magnitudes=aggregated,assumptions=(quantity=calculate_events ? :event_magnitudes :
        calculate_coda || calculate_amplitude ? :station_magnitudes : :calibration_and_station_corrections,
        calibration_time=:located_event_origin_utc,history_boundary=:hourly_end_exclusive,
        channel_identity=:exact_complete_station_id,measurement_values_external=true,
        coda_p_time=:bare_model_p_ray_without_station_delay,
        coda_distance=:horizontal_epicentral_km,coda_depth=:calculation_z_km,
        coda_magnitudes_computed=calculate_coda,amplitude_magnitudes_computed=calculate_amplitude,
        station_magnitudes_computed=calculate_coda || calculate_amplitude,event_magnitudes_computed=calculate_events))
end

function _write_magnitude_preparation(directory,result)
    metadata=_result_metadata(result)
    _write_result_csv(joinpath(directory,"magnitude_calibration.csv"),
        ("input_row","event_id","station_id","event_time_utc","location_status","instrument_type","calibration_factor",
         "status","reason","calibration_source","history_row","amplitude_calibrated","raw_coda_gain_correction","effective_coda_gain_correction",
         "md1_correction_mag","md2_correction_mag","mt_correction_mag","amplitude_correction_mag","coda_event_weight","amplitude_event_weight","coda_gain_enabled"),
        ((r.input_row,r.event_id,r.station_id,r.event_time_utc,r.location_status,r.calibration.instrument_type,
          r.calibration.calibration_factor,r.calibration.status,r.calibration.reason,r.calibration.calibration_source,
          r.calibration.history_row,r.calibration.amplitude_calibrated,r.calibration.coda_gain_correction,r.effective_coda_gain_correction,
          r.corrections.md1_correction_mag,r.corrections.md2_correction_mag,r.corrections.mt_correction_mag,r.corrections.amplitude_correction_mag,
          r.corrections.coda_event_weight,r.corrections.amplitude_event_weight,r.corrections.coda_gain_enabled) for r in result.rows))
    open(joinpath(directory,"magnitude_calibration.toml"),"w") do io
        TOML.print(io,metadata)
    end
    if result.assumptions.coda_magnitudes_computed
        terms=(:constant_mag,:log_mag,:log_squared_mag,:linear_time_mag,:depth_mag,:distance_mag,
            :gain_mag,:station_mag,:component_mag,:near_distance_mag,:far_distance_mag,:additional_depth_mag)
        _write_result_csv(joinpath(directory,"coda_station_magnitudes.csv"),
            ("input_row","event_id","station_id","location_status","magnitude_type","magnitude","status","reason",
             "duration_s","calculated_p_time_s","model_p_time_status","tau_s","distance_km","calculation_depth_km",
             "duration_branch","coda_quality_code","coda_event_weight",String.(terms)...),
            ((r.input_row,r.event_id,r.station_id,r.location_status,kind,m.magnitude,m.status,m.reason,
              get(r.measurement,"duration_s",nothing),r.coda.calculated_p_time_s,r.coda.model_p_time_status,r.coda.tau_s,
              r.coda.distance_km,r.coda.calculation_depth_km,m.duration_branch,r.coda.coda_quality_code,r.coda.coda_event_weight,
              (m.terms===nothing ? ntuple(_->nothing,length(terms)) : Tuple(m.terms))...)
             for r in result.rows for kind in (:md1,:md2,:mt) for m in (getproperty(r.coda,kind),)))
        open(joinpath(directory,"coda_station_magnitudes.toml"),"w") do io
            TOML.print(io,metadata)
        end
    end
    if result.assumptions.amplitude_magnitudes_computed
        terms=(:log_amplitude_calibration_mag,:response_mag,:motor_mag,:sprengnether_mag,:distance_mag,:station_mag,:component_mag)
        _write_result_csv(joinpath(directory,"amplitude_station_magnitudes.csv"),
            ("input_row","event_id","station_id","location_status","magnitude_type","magnitude","status","reason",
             "amplitude","amplitude_units","amplitude_type","period_s","period_source","unit_factor_mm",
             "amplitude_peak_to_peak_mm","calibration_factor","amplitude_quality_code","amplitude_event_weight",
             "distance_relation","epicentral_distance_km","calculation_depth_km","slant_distance_km",String.(terms)...),
            ((r.input_row,r.event_id,r.station_id,r.location_status,m.magnitude_type,m.magnitude,m.status,m.reason,
              get(r.measurement,"amplitude",nothing),m.amplitude_units,m.amplitude_type,m.period_s,m.period_source,
              m.unit_factor_mm,m.amplitude_peak_to_peak_mm,m.calibration_factor,get(r.measurement,"amplitude_quality_code",nothing),
              m.amplitude_event_weight,m.distance.relation,m.distance.epicentral_distance_km,m.distance.calculation_depth_km,
              m.distance.slant_distance_km,(m.terms===nothing ? ntuple(_->nothing,length(terms)) : Tuple(m.terms))...)
             for r in result.rows for m in (r.amplitude,)))
        open(joinpath(directory,"amplitude_station_magnitudes.toml"),"w") do io
            TOML.print(io,metadata)
        end
    end
    result.event_magnitudes===nothing || _write_event_magnitudes(directory,result.event_magnitudes)
end
