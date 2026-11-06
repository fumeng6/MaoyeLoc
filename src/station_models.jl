const StationModelLeaf = Union{PhaseModelLeaf,PhaseModels}

"""
    StationModels(models, station_model_ids; default_model_id, station_delays_s=Dict())

Select one complete, fixed one-dimensional P/S model by station ID. `models`
maps model IDs to supported ray/lookup-table models or explicit P/S model
pairs; `station_model_ids`
maps station IDs to model IDs. Unmapped stations use `default_model_id`.
The choice does not depend on source position. Optional known `station_delays_s`
maps `(model_id, station_id, phase)` to finite seconds, independently for P/S.
The selected model's delay is added to its predicted phase time; missing
entries mean zero. For an affine `PhaseModels` P/S relation, only P delay is supplied
and S delay is s_scale*P_delay. Independent phase pairs use their own P/S
delays. Delays do not change source derivatives and are not estimated.

The source must lie at or below the deepest model top among the default and
explicitly referenced models. Stored but unreferenced models do not limit it.
Each station itself only needs to lie within its selected model.
"""
struct StationModels
    models::Dict{String,StationModelLeaf}
    station_model_ids::Dict{String,String}
    default_model_id::String
    station_delays_s::Dict{Tuple{String,String,String},Float64}

    function StationModels(models::AbstractDict, station_model_ids::AbstractDict;
                           default_model_id::AbstractString,
                           station_delays_s::AbstractDict=Dict())
        leaves = Dict{String,StationModelLeaf}()
        for (id, leaf) in models
            id isa AbstractString && !isempty(id) ||
                throw(ArgumentError("Model IDs must be nonempty strings."))
            leaf isa StationModelLeaf ||
                throw(ArgumentError("Station models must be supported one-dimensional ray/table models or P/S model pairs; nested StationModels are not supported."))
            leaves[String(id)] = leaf
        end
        default = String(default_model_id)
        haskey(leaves,default) || throw(ArgumentError("Unknown default model ID $default."))
        assignments = Dict{String,String}()
        for (station_id, model_id) in station_model_ids
            station_id isa AbstractString && !isempty(station_id) && model_id isa AbstractString ||
                throw(ArgumentError("Station assignments must map nonempty station IDs to model ID strings."))
            haskey(leaves,model_id) || throw(ArgumentError("Unknown model ID $model_id for station $station_id."))
            assignments[String(station_id)] = String(model_id)
        end
        delays = Dict{Tuple{String,String,String},Float64}()
        for (key,value) in station_delays_s
            key isa Tuple && length(key)==3 && all(id -> id isa AbstractString && !isempty(id),key) ||
                throw(ArgumentError("Station delay keys must be (model_id, station_id, phase) strings."))
            model_id,station_id,phase = String.(key)
            haskey(leaves,model_id) || throw(ArgumentError("Unknown delay model ID $model_id."))
            phase in ("P","S") || throw(ArgumentError("Station delay phase must be P or S."))
            leaf = leaves[model_id]
            phase == "S" && leaf isa PhaseModels && leaf.s_model === nothing &&
                throw(ArgumentError("Affine S models derive their S delay from s_scale * P delay; do not specify a separate S delay."))
            value isa Real && isfinite(value) || throw(ArgumentError("Station delays must be finite seconds."))
            delay = Float64(value)
            isfinite(delay) || throw(ArgumentError("Station delays must be finite Float64 seconds."))
            delays[(model_id,station_id,phase)] = delay
        end
        new(leaves,assignments,default,delays)
    end
end

model_top(model::StationModels) = max(model_top(model.models[model.default_model_id]),
    maximum((model_top(model.models[id]) for id in values(model.station_model_ids));init=-Inf))

_station_model(model, station_id) = model
_station_model(model::StationModels, station_id) =
    model.models[get(model.station_model_ids,station_id,model.default_model_id)]

_station_delay(model, station_id, phase) = 0.0
function _station_delay(model::StationModels,station_id,phase)
    id = get(model.station_model_ids,station_id,model.default_model_id)
    leaf = model.models[id]
    if phase == "S" && leaf isa PhaseModels && leaf.s_model === nothing
        return leaf.s_scale*get(model.station_delays_s,(id,station_id,"P"),0.)
    end
    get(model.station_delays_s,(id,station_id,phase),0.0)
end

"""
    travel_time(model, source_km, stations, station_id, phase)

Evaluate the model assigned to this station. For a single model this is the
same as `travel_time(model, source_km, stations[station_id], phase)`.
Station IDs remain explicit so every inversion and covariance path uses the
same fixed choice. A `StationModels` value requires this five-argument form.
When the set contains configured delays, `time_s=ray_time_s+phase_offset_s+station_delay_s`
and the selected `model_id` are returned explicitly; gradients remain ray
derivatives. phase_offset_s is absent (and zero) on ordinary ray/table models.
A set without delays preserves the selected model's original return value.
"""
travel_time(model, source_km, stations::AbstractDict, station_id::AbstractString,
            phase::AbstractString) = travel_time(model,source_km,stations[station_id],phase)

function travel_time(model::StationModels, source_km, stations::AbstractDict,
                     station_id::AbstractString, phase::AbstractString)
    length(source_km) == 3 || throw(ArgumentError("Source must contain (x, y, z) in km."))
    source_km[3] >= model_top(model) ||
        throw(ArgumentError("Source lies above the deepest top of the default and assigned station models."))
    ray = travel_time(_station_model(model,station_id),source_km,stations[station_id],phase)
    isempty(model.station_delays_s) && return ray
    delay = _station_delay(model,station_id,phase)
    return merge(ray,(time_s=ray.time_s+delay,ray_time_s=get(ray,:ray_time_s,ray.time_s),
        station_delay_s=delay,model_id=get(model.station_model_ids,station_id,model.default_model_id)))
end
