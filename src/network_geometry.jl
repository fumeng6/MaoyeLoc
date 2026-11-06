# Hypo71 labels are descriptive combinations of the original one-standard-
# axis errors and coverage criteria; they are not confidence probabilities.
function _hypo71_quality(rms,erh,erz,nwr,gap,depth,dmin)
    geometry = gap===nothing || depth===nothing || dmin===nothing ? nothing :
        nwr>=6 && gap<=90 && (dmin<=depth || dmin<=5) ? 1 :
        nwr>=6 && gap<=135 && (dmin<=2depth || dmin<=10) ? 2 :
        nwr>=6 && gap<=180 && dmin<=50 ? 3 : 4
    error = rms===nothing || erh===nothing || erz===nothing ? nothing :
        rms<.15 && erh<=1 && erz<=2 ? 1 :
        rms<.30 && erh<=2.5 && erz<=5 ? 2 :
        rms<.50 && erh<=5 ? 3 : 4
    letter(q)=q===nothing ? nothing : string(Char(64+q))
    # All ranks are positive; Fortran NINT rounds exact halves UP, unlike
    # Julia's default ties-to-even round.
    combined=geometry===nothing || error===nothing ? nothing : fld(geometry+error+1,2)
    (;quality_code=letter(combined),error_quality_code=letter(error),
      geometry_quality_code=letter(geometry),
      quality_reason=error===nothing ? :formal_errors_unavailable : geometry===nothing ? :coverage_unavailable : nothing)
end

function _geometry_source_angle(ray)
    if hasproperty(ray,:components)
        components=[(;model_id=c.model_id,weight=c.weight,
            _geometry_source_angle(c.prediction)...) for c in ray.components]
        if length(components)==1
            c=only(components)
            return (;source_angle_deg=c.source_angle_deg,angle_definition=c.angle_definition,
                angle_reason=c.angle_reason,source_angle_components=components)
        end
        return (;source_angle_deg=nothing,angle_definition=:regional_mixture,
            angle_reason=:no_single_ray_for_model_mixture,source_angle_components=components)
    end
    if hasproperty(ray,:source_angle_deg)
        # A table supplies the normal to its interpolated equal-time surface,
        # rather than an independently traced physical ray direction.
        smooth=ray.dsource_s_km!==nothing
        return (;source_angle_deg=smooth ? ray.source_angle_deg : nothing,angle_definition=:table_isochrone_normal,
            angle_reason=smooth ? nothing : :nondifferentiable_ray_without_unique_direction,source_angle_components=NamedTuple[])
    end
    gradient=ray.dsource_s_km
    angle=gradient===nothing ? nothing : atand(hypot(gradient[1],gradient[2]),-gradient[3])
    (;source_angle_deg=angle,angle_definition=:isotropic_ray_from_nadir,
      angle_reason=angle===nothing ? :nondifferentiable_ray_without_unique_direction : nothing,
      source_angle_components=NamedTuple[])
end

"""
    absolute_geometry(model,stations,location; reported_depth_km=location.source_km[3],
                      ellipsoid=nothing,station_site_ids=Dict())

Final per-reading epicentral DIST, clockwise-from-north AZM and source angle
from nadir (0 downward, 180 upward), plus GAP/DMIN and Hypo71 A-D labels.
Smooth isotropic-ray direction is -grad_source(T), by Fermat's endpoint
relation, independent of speed. At a nonunique first arrival no direction is
invented. Table angles describe the normal to the interpolated equal-time
surface. A regional mixture has component directions but no single physical ray.

Empirical absolute rows use saved final equation_weight, without treating
`used=false` on fixed events as zero weight. Sigma-WLS uses actual 1/sigma
normalized to mean one, explicitly identifying this native convention.
NWR/NWS count weights >0.1; GAP/DMIN use weights >=0.1, matching the original
distinct inequalities. The valid-reading count uses quality_code <4 with a
finite P/S arrival, independently of station or residual reweighting. Rows
without a quality code explicitly default to zero for this count only; that
default is not an empirical location weight. GAP/DMIN deduplicate
physical sites using station_site_ids (default exact station ID).
Each reading and its regional ray components retain observation_id; separate
IDs at one station do not create extra physical sites or azimuthal coverage.

Hypo71 error class uses legacy_erh_km/legacy_erz_km from the supplied
ellipsoid, NEVER its confidence-scaled extents. Missing errors give no overall
quality. reported_depth_km is the same explicitly reported depth used by the
original quality function. All labels are descriptive, not field accuracy.
"""
function absolute_geometry(model,stations::AbstractDict,location;
        reported_depth_km=location.source_km[3],ellipsoid=nothing,station_site_ids::AbstractDict=Dict())
    source=location.source_km
    length(source)==3 && all(isfinite,source) || throw(ArgumentError("Geometry requires finite source xyz km."))
    reported_depth_km===nothing || reported_depth_km isa Real && isfinite(reported_depth_km) ||
        throw(ArgumentError("Reported quality depth must be finite km or missing."))
    readings=collect(location.residuals)
    empirical=get(location,:row_weighting,:empirical)==:empirical
    empirical || throw(ArgumentError("Geometry uses the saved empirical equation weights."))
    weights=Float64[row.equation_weight for row in readings]
    all(w->isfinite(w)&&w>=0,weights) || throw(ArgumentError("Geometry row weights must be finite and nonnegative."))
    rows=NamedTuple[]
    for (i,row) in enumerate(readings)
        haskey(stations,row.station_id) && row.phase in ("P","S") || throw(ArgumentError("Geometry needs known P/S station identities."))
        station=stations[row.station_id]
        dx,dy=station[1]-source[1],station[2]-source[2]
        distance=hypot(dx,dy)
        azimuth=distance==0 ? nothing : mod(atand(dx,dy),360.)
        site=String(get(station_site_ids,row.station_id,row.station_id))
        isempty(site) && throw(ArgumentError("Physical site identities must be nonempty."))
        if hasproperty(row,:quality_code)
            value=getproperty(row,:quality_code)
            value isa Integer && !(value isa Bool) && 0<=value<=9 ||
                throw(ArgumentError("Geometry quality_code must be an integer from 0 through 9."))
            quality=Int(value)
        else
            quality=0
        end
        phase_present=isfinite(row.arrival_s)
        angle=try
            _geometry_source_angle(travel_time(model,source,stations,row.station_id,row.phase))
        catch error
            error isa DomainError || rethrow()
            (;source_angle_deg=nothing,angle_definition=:unavailable,angle_reason=:undefined_ray,
              source_angle_components=NamedTuple[])
        end
        observation_id=get(row,:observation_id,nothing)
        angle=merge(angle,(source_angle_components=[merge(c,(;observation_id)) for c in angle.source_angle_components],))
        push!(rows,(;observation_id,station_id=row.station_id,phase=row.phase,site_id=site,distance_km=distance,
            azimuth_deg=azimuth,equation_weight=weights[i],used_for_fit=get(row,:used,false),
            weighted_reading=weights[i]>.1,geometry_eligible=weights[i]>=.1,
            valid_assigned_reading=quality<4 && phase_present,
            quality_code=quality,quality_source=get(row,:quality_code_source,
                hasproperty(row,:quality_code) ? :explicit : :default_zero),angle...))
    end
    # Closest component represents a supplied physical site, with an explicit
    # stable tie rule. P/S duplicates do not create extra angular coverage.
    sites=NamedTuple[]
    seen=Set{String}()
    for row in sort([r for r in rows if r.geometry_eligible];by=r->(r.distance_km,r.station_id,r.phase))
        row.site_id in seen && continue
        push!(seen,row.site_id)
        push!(sites,(;site_id=row.site_id,station_id=row.station_id,distance_km=row.distance_km,azimuth_deg=row.azimuth_deg))
    end
    azimuths=sort(unique(Float64[r.azimuth_deg for r in sites if r.azimuth_deg!==nothing]))
    gap=isempty(azimuths) ? nothing : maximum(diff(vcat(azimuths,first(azimuths)+360)))
    dmin=isempty(sites) ? nothing : minimum(r.distance_km for r in sites)
    nwr=count(r->r.weighted_reading,rows)
    nws=count(r->r.weighted_reading && r.phase=="S",rows)
    nvr=count(r->r.valid_assigned_reading,rows)
    denominator=sum(abs2,weights)
    rms=denominator==0 ? nothing : sqrt(sum((weights[i]*readings[i].residual_s)^2 for i in eachindex(weights))/denominator)
    erh=ellipsoid===nothing ? nothing : get(ellipsoid,:legacy_erh_km,nothing)
    erz=ellipsoid===nothing ? nothing : get(ellipsoid,:legacy_erz_km,nothing)
    labels=_hypo71_quality(rms,erh,erz,nwr,gap,reported_depth_km,dmin)
    (;rows,sites,gap_deg=gap,nearest_station_distance_km=dmin,weighted_reading_count=nwr,
      weighted_p_count=nwr-nws,weighted_s_count=nws,valid_reading_count=nvr,
      effective_station_count=length(sites),azimuth_station_count=count(r->r.azimuth_deg!==nothing,sites),
      weighted_rms_s=rms,legacy_erh_km=erh,legacy_erz_km=erz,reported_depth_km,
      row_weight_definition=:saved_empirical_absolute_equation_weight,
      row_normalization=1.,geometry_reason=gap===nothing ? :no_nonvertical_eligible_station : nothing,
      labels...)
end
