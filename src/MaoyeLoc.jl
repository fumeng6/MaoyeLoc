"""
Local-earthquake absolute location and double-difference relocation.

P/S first arrivals in homogeneous, layered and linear-gradient media, and
explicit empirical or independent S phase relations.
"""
module MaoyeLoc

using LinearAlgebra
using SparseArrays
using TOML
using Dates
using Printf
using Random

const _APP_VERSION = v"1.0.0"

include("execution.jl")
include("result_metadata.jl")

export HomogeneousModel, LayeredModel, LinearGradientModel, FiniteGradientModel, PiecewiseGradientModel, StationModels, travel_time, model_top, locate, relocate
export read_model, read_stations, read_events, read_picks, read_differentials
export read_location_results, differential_geometry
export locate_catalog, write_absolute_results, write_relative_results, run_case
export prepare_case, check_case, effective_config
export initialize_catalog
export LocalFrame, geographic_to_local, local_to_geographic
export catalog_differences, synthetic_differences
export write_geographic_results, write_location_report, write_magnitude_details
export review_catalog, run_catalog_review
export save_physical_state, load_physical_state
export prepare_location_plot, plot_location_results
export calculation_defaults, run_batch
export pair_catalog, pair_neighbors
export absolute_uncertainty,absolute_formal_uncertainty,error_ellipsoid,relative_uncertainty
export absolute_diagnostics,absolute_geometry
export geographic_region
export processing_metadata
export relative_statistics,relative_observation_diagnostics
export relative_iteration_diagnostics
export relative_algorithm_uncertainty
export magnitude_calibration,magnitude_corrections,magnitude_preparation,coda_magnitudes
export amplitude_distance_correction,amplitude_magnitude
export magnitude_observation_weight,aggregate_magnitudes,event_magnitudes
export magnitude_selection
export preferred_magnitude,magnitude_display
export TravelTimeTable, TravelTimeTableModel, read_travel_time_table, write_travel_time_table, generate_table, evaluate_table
export PhaseModels
export RayTransformModel
export RegionalNode, RegionalModels, regional_weights
export depth_report, trace_ray, ray_geometric_density

"""
    HomogeneousModel(vp_km_s, vs_km_s; top_km)

Constant P/S speeds (km/s) in a half-space extending from `top_km` downward.
Coordinates are in km, with z positive downward. The model top and both
endpoints must use the same vertical reference; station elevation h means z=-h.
"""
struct HomogeneousModel
    vp_km_s::Float64
    vs_km_s::Float64
    top_km::Float64

    function HomogeneousModel(vp_km_s::Real, vs_km_s::Real; top_km::Real)
        vp, vs, top = Float64(vp_km_s), Float64(vs_km_s), Float64(top_km)
        all(isfinite, (vp, vs, top)) ||
            throw(ArgumentError("Model speeds and top must be finite."))
        vp > vs > 0 || throw(ArgumentError("Require Vp > Vs > 0 km/s."))
        new(vp, vs, top)
    end
end

"""
    travel_time(model::HomogeneousModel, source_km, station_km, phase)

Return `(time_s, dsource_s_km)` for phase `"P"` or `"S"`, where each endpoint
contains (x, y, z) in km. `dsource_s_km` is the tuple
(dT/dx_source, dT/dy_source, dT/dz_source), in s/km, holding the station fixed.
These are derivatives of predicted travel time, not of observed-minus-predicted
residual. The homogeneous first arrival follows the straight line:
T = r/v and dT/dsource = (source - station)/(v*r).

Endpoints above the model top are rejected, without clipping or extrapolation.
Coincident endpoints are rejected because the source gradient is undefined,
although the travel time alone would be zero. No reflections, conversions,
velocity gradients, layered refractions, or station corrections are included.
"""
function travel_time(model::HomogeneousModel, source_km, station_km,
                     phase::AbstractString)
    length(source_km) == 3 && length(station_km) == 3 ||
        throw(ArgumentError("Source and station must each contain (x, y, z) in km."))
    source = ntuple(i -> Float64(source_km[i]), 3)
    station = ntuple(i -> Float64(station_km[i]), 3)
    all(isfinite, source) && all(isfinite, station) ||
        throw(ArgumentError("Source and station coordinates must be finite."))
    source[3] >= model.top_km && station[3] >= model.top_km ||
        throw(ArgumentError("Source and station must be at or below the model top."))

    velocity = if phase == "P"
        model.vp_km_s
    elseif phase == "S"
        model.vs_km_s
    else
        throw(ArgumentError("Phase must be P or S."))
    end

    delta = ntuple(i -> source[i] - station[i], 3)
    distance = hypot(delta...)
    distance > 0 ||
        throw(DomainError(distance, "Coincident endpoints have no unique source gradient."))
    time_s = distance / velocity
    dsource_s_km = ntuple(i -> (delta[i] / distance) / velocity, 3)
    return (; time_s, dsource_s_km)
end

include("layered.jl")
include("gradient.jl")
include("finite_gradient.jl")
include("piecewise_gradient.jl")
include("ray_diagnostics.jl")
include("travel_time_tables.jl")
include("travel_time_generation.jl")
include("ray_transform.jl")
include("phase_models.jl")
include("station_models.jl")
include("regional_models.jl")

# Coincident endpoints have zero geometric time and no unique gradient.
# An explicit empirical S intercept can remain in the phase prediction.
# Fixed references and failure residuals need only that forward time.
function _forward_time(model,source,station,phase)
    try
        return travel_time(model,source,station,phase).time_s
    catch exception
        exception isa DomainError && Tuple(source)==Tuple(station) || rethrow()
        return _coincident_time_s(model,phase) # Endpoints and phase have been validated.
    end
end

function _forward_time(model,source,stations::AbstractDict,station_id,phase)
    try
        return travel_time(model,source,stations,station_id,phase).time_s
    catch exception
        exception isa DomainError && Tuple(source)==Tuple(stations[station_id]) || rethrow()
        # The selected ray/table model has validated both endpoints and phase.
        # Coincidence removes geometric time; an empirical phase intercept and
        # known station delay still apply to a fixed reference prediction.
        return _coincident_time_s(_station_model(model,station_id),phase)+_station_delay(model,station_id,phase)
    end
end

include("coordinates.jl")
include("absolute.jl")
include("absolute_step_controls.jl")
include("absolute_controls.jl")
include("absolute_weights.jl")
include("absolute_uncertainty.jl")
include("absolute_formal_uncertainty.jl")
include("error_ellipsoid.jl")
include("absolute_diagnostics.jl")
include("network_geometry.jl")
include("geographic_regions.jl")
include("magnitude_calibration.jl")
include("magnitude_corrections.jl")
include("coda_magnitudes.jl")
include("amplitude_magnitudes.jl")
include("magnitude_aggregation.jl")
include("magnitude_selection.jl")
include("preferred_magnitudes.jl")
include("processing_metadata.jl")
include("relative_scaled_steps.jl")
include("relative_design.jl")
include("relative_lsqr.jl")
include("relative_selection.jl")
include("relative_weights.jl")
include("relative.jl")
include("relative_weight_history.jl")
include("relative_scaled_location.jl")
include("relative_uncertainty.jl")
include("relative_statistics.jl")
include("relative_observation_diagnostics.jl")
include("relative_iteration_diagnostics.jl")
include("relative_algorithm_uncertainty.jl")
include("input.jl")
include("differential_time_input.jl")
include("pairing.jl")
include("neighbor_pairing.jl")
include("depth_reporting.jl")
include("case_validation.jl")
include("case.jl")
include("workflow.jl")
include("result_output.jl")
include("geographic_output.jl")
include("magnitude_workflow.jl")
include("event_magnitudes.jl")
include("location_report.jl")
include("catalog_review.jl")
include("synthetic_observations.jl")
include("physical_state.jl")
include("location_plot.jl")
include("batch_commands.jl")
include("command_line.jl")

end
