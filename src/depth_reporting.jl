# Reporting depth has an explicit datum; a table uses its receiver plane.
_depth_model_zero(model) = 0.
_depth_model_zero(model::TravelTimeTableModel) = model.station_plane_km
_depth_model_zero(model::PhaseModels) = _depth_model_zero(model.p_model)
_depth_model_zero(model::StationModels) = _depth_model_zero(model.models[model.default_model_id])

"""
    depth_report(model,source,stations,station_ids; model_datum=:nearest5,
                 coordinate_datum=:geoid,report_datum=:geoid,
                 station_selection=:nearest5,...)

Express the same computed source depth relative to either a model reference
plane or the caller-declared geoid (the vertical reference for physical
elevations). All depths and elevations are in km, positive down and up,
respectively. This routine does not compute geoid undulation or convert
ellipsoidal heights. `model_datum`, `coordinate_datum` and `report_datum`
select the reporting convention, the input z reference and the output depth
reference, respectively; they do not change the propagation model.

`station_ids` contains ALL matched event station
records, including unweighted or unused records. Select up to five nearest
horizontal sites with known elevations, deduplicating components using
`station_site_ids`. Missing elevation-map entries mean unknown; a known zero
is valid. If no map is given, elevations are -station_z only for geoid
coordinates. Converting model coordinates with model datum :nearest5
requires explicit physical elevations.

The default :nearest5 selection has no distance/elevation cutoff.
Optional :nearby5 uses distance <= 100 km and elevation >= -0.120 km.
If no site qualifies, use the nearest record with known elevation, even
outside the radius or below the elevation cutoff. Unknown elevations never
become zero. The selection is applied only for model datum :nearest5;
datum :geoid needs no station selection. This option adds
`datum_selection` with the strategy, fallback status and per-record reasons;
the default return fields remain unchanged. Local horizontal distances,
ID tie-breaking and unrounded elevations retain the native conventions,
and does not alter source/station coordinates. Result writers
export the optional records to `datum_stations.csv`; this API does not write.

With model datum :nearest5, H is the unrounded mean elevation (km). Geoid
coordinates give z_geoid=z and z_model=z+H. Model coordinates give
z_model=z-model_zero_km and z_geoid=z_model-H. The default numerical model zero
is the travel-time table's station plane, or zero for other models. No source, station,
ray, residual or covariance is changed. Datum :geoid uses H=0 and reports G
regardless of report_datum. Missing H leaves the input-datum depth available.
`status` tracks whether H is known; read the individual depth fields as well.
Downstream Hypo71 quality classification uses reported depth and can change
when the selected datum elevation changes. Datum elevation uncertainty and
station-selection changes are not propagated into the location covariance.

`regional_model_datums` may override the convention by the ID of the region
with the largest combined model weight. The Hypoinverse ALT station-model
choice and SAL independent S-wave structure do not separately choose a depth
reference plane. Equal weights keep node
definition order. Unknown datum is reported, never replaced by zero.
"""
function depth_report(model,source,stations::AbstractDict,station_ids;
                      model_datum::Symbol=:nearest5,coordinate_datum::Symbol=:geoid,
                      report_datum::Symbol=:geoid,model_zero_km=nothing,
                      station_selection::Symbol=:nearest5,
                      station_elevations_km=nothing,station_site_ids::AbstractDict=Dict(),
                      regional_model_datums::AbstractDict=Dict())
    length(source)==3 && all(isfinite,source) || throw(ArgumentError("Depth reporting needs finite source xyz in km."))
    coordinate_datum in (:model,:geoid) && report_datum in (:model,:geoid) ||
        throw(ArgumentError("coordinate_datum and report_datum must be model (model reference plane) or geoid (caller-declared geoid)."))
    station_selection in (:nearest5,:nearby5) ||
        throw(ArgumentError("Depth-reference station_selection must be nearest5 or nearby5."))
    branch,dominant = model,nothing
    if model isa RegionalModels
        weights = regional_weights(model,source)
        dominant = weights[argmax([w.weight for w in weights])].model_id
        branch = model.models[dominant]
        all(id -> haskey(model.models,id),keys(regional_model_datums)) ||
            throw(ArgumentError("regional_model_datums contains an unknown regional model ID."))
        model_datum = Symbol(get(regional_model_datums,dominant,model_datum))
    elseif !isempty(regional_model_datums)
        throw(ArgumentError("Regional reporting conventions need RegionalModels."))
    end
    model_datum in (:nearest5,:geoid) || throw(ArgumentError("model_datum must be nearest5 (selected stations' mean-elevation plane) or geoid (caller-declared geoid)."))
    zero_km = model_zero_km === nothing ? _depth_model_zero(branch) : Float64(model_zero_km)
    isfinite(zero_km) || throw(ArgumentError("model_zero_km must be finite."))
    station_elevations_km === nothing && coordinate_datum==:model && model_datum==:nearest5 &&
        throw(ArgumentError("Model-coordinate reporting needs explicit physical station_elevations_km."))
    ids = unique(String.(collect(station_ids)))
    all(id -> haskey(stations,id),ids) || throw(ArgumentError("Report station list contains an unknown station ID."))
    chosen = String[]
    heights = Float64[]
    selection_records = NamedTuple[]
    if model_datum==:nearest5
        sort!(ids;by=id -> (hypot(source[1]-stations[id][1],source[2]-stations[id][2]),id))
        sites = Set{String}()
        for id in ids
            elevation = station_elevations_km === nothing ? -stations[id][3] : get(station_elevations_km,id,nothing)
            if station_selection==:nearby5
                elevation === nothing || elevation isa Real && isfinite(elevation) ||
                    throw(ArgumentError("Station report elevations must be finite km or unknown."))
                distance = hypot(source[1]-stations[id][1],source[2]-stations[id][2])
                site = String(get(station_site_ids,id,id))
                isempty(site) && throw(ArgumentError("Station site identities must be nonempty."))
                # Filter before deduplication: a rejected component must not hide a usable one.
                reason = elevation === nothing ? :unknown_elevation : distance>100. ? :outside_radius :
                    elevation<-.120 ? :below_minimum_elevation : site in sites ? :duplicate_site :
                    length(chosen)>=5 ? :beyond_nearest_five : :selected
                selected = reason==:selected
                if selected
                    push!(sites,site)
                    push!(chosen,id)
                    push!(heights,Float64(elevation))
                end
                push!(selection_records,(station_id=id,site_id=site,horizontal_distance_km=distance,
                    elevation_km=elevation === nothing ? nothing : Float64(elevation),selected,reason))
                continue
            end
            elevation === nothing && continue
            elevation isa Real && isfinite(elevation) || throw(ArgumentError("Station report elevations must be finite km or unknown."))
            site = String(get(station_site_ids,id,id))
            isempty(site) && throw(ArgumentError("Station site identities must be nonempty."))
            site in sites && continue
            push!(sites,site)
            push!(chosen,id)
            push!(heights,Float64(elevation))
            length(chosen)==5 && break
        end
        if station_selection==:nearby5 && isempty(chosen)
            # Bypass both cutoffs on fallback, but keep native unknown elevations missing.
            fallback = findfirst(row -> row.elevation_km!==nothing,selection_records)
            if fallback!==nothing
                row = selection_records[fallback]
                push!(chosen,row.station_id)
                push!(heights,row.elevation_km)
                selection_records[fallback] = merge(row,(selected=true,reason=:nearest_known_station_fallback))
            end
        end
    end
    height = model_datum==:geoid ? 0. : isempty(heights) ? nothing : sum(heights)/length(heights)
    z_model,z_geoid = if coordinate_datum==:model
        z = Float64(source[3])-zero_km
        (z,height === nothing ? nothing : z-height)
    else
        z = Float64(source[3])
        (height === nothing ? nothing : z+height,z)
    end
    target = model_datum==:geoid ? :geoid : report_datum
    reported = target==:geoid ? z_geoid : z_model
    report = (;model_depth_km=z_model,geoid_depth_km=z_geoid,reported_depth_km=reported,
      report_datum=target,model_datum,coordinate_depth_datum=coordinate_datum,
      model_zero_km=zero_km,datum_elevation_km=height,datum_station_count=length(chosen),
      datum_station_ids=chosen,dominant_model_id=dominant,
      depth_datum_flag=target==:geoid ? :G : :M,
      status=height === nothing ? :no_known_elevations : :available)
    station_selection==:nearest5 && return report
    selection_status = model_datum==:geoid ? :not_required : isempty(chosen) ? :no_known_elevations :
        any(row -> row.reason==:nearest_known_station_fallback,selection_records) ? :nearest_known_station_fallback : :nearby_stations
    merge(report,(datum_selection=(strategy=station_selection,status=selection_status,
        radius_km=100.,minimum_elevation_km=-.120,stations=selection_records),))
end

# Only explicit nearby5 produces this table; default exports keep their existing layout.
function _write_depth_stations(directory,events)
    reports = [event for event in events if hasproperty(event,:depth_report) &&
        hasproperty(event.depth_report,:datum_selection)]
    isempty(reports) && return
    _write_result_csv(joinpath(directory,"datum_stations.csv"),
        ("event_id","station_selection","selection_status","station_id","site_id",
         "horizontal_distance_km","elevation_km","selected","reason"),
        ((event.event_id,event.depth_report.datum_selection.strategy,event.depth_report.datum_selection.status,
          row.station_id,row.site_id,row.horizontal_distance_km,row.elevation_km,row.selected,row.reason)
         for event in reports for row in event.depth_report.datum_selection.stations))
end

const _DEPTH_RESULT_HEADER = ("model_depth_km","geoid_depth_km","reported_depth_km",
    "report_datum","model_datum","coordinate_depth_datum","model_zero_km",
    "datum_elevation_km","datum_station_count","datum_station_ids","dominant_model_id",
    "depth_datum_flag","depth_report_status")

_depth_result_header(events) = any(event -> hasproperty(event,:depth_report),events) ? _DEPTH_RESULT_HEADER : ()
function _depth_result_values(event)
    hasproperty(event,:depth_report) || return ()
    d = event.depth_report
    (d.model_depth_km,d.geoid_depth_km,d.reported_depth_km,d.report_datum,d.model_datum,
     d.coordinate_depth_datum,d.model_zero_km,d.datum_elevation_km,d.datum_station_count,
     join(d.datum_station_ids,'|'),d.dominant_model_id,d.depth_datum_flag,d.status)
end

function _with_depth_reports(result,model,stations,event_stations,config)
    isempty(config) && return result
    all(haskey(config,key) for key in ("model_datum","coordinate_datum","report_datum")) ||
        throw(ArgumentError("depth_reporting requires explicit model_datum, coordinate_datum and report_datum."))
    options = (model_datum=Symbol(config["model_datum"]),coordinate_datum=Symbol(config["coordinate_datum"]),
        report_datum=Symbol(config["report_datum"]),model_zero_km=get(config,"model_zero_km",nothing),
        station_selection=Symbol(get(config,"station_selection","nearest5")),
        station_elevations_km=get(config,"station_elevations_km",nothing),
        station_site_ids=get(config,"station_site_ids",Dict()),
        regional_model_datums=get(config,"regional_model_datums",Dict()))
    extra = get(config,"event_stations",Dict())
    ids = Set(event.event_id for event in result.events)
    all(id -> id in ids,keys(extra)) || throw(ArgumentError("Depth report station inventory refers to an unknown event ID."))
    rows = [merge(event,(depth_report=depth_report(model,event.source_km,stations,
        [get(event_stations,event.event_id,String[]);get(extra,event.event_id,String[])];options...),))
        for event in result.events]
    merge(result,(events=rows,))
end
