"""A circular regional node in the same local east/north km frame as sources."""
struct RegionalNode
    model_id::String
    x_km::Float64
    y_km::Float64
    inner_radius_km::Float64
    outer_radius_km::Float64

    function RegionalNode(model_id::AbstractString,x::Real,y::Real,inner::Real,outer::Real)
        values = Float64.((x,y,inner,outer))
        !isempty(model_id) && all(isfinite,values) && 0 <= values[3] < values[4] ||
            throw(ArgumentError("Regional nodes need a model ID, finite xy and 0 <= inner < outer radii in km."))
        new(String(model_id),values...)
    end
end

const RegionalModelBranch = Union{StationModelLeaf,StationModels}

"""
    RegionalModels(models,nodes;default_model_id)

Predict arrival times by weighting one-dimensional model predictions according
to source position. This empirical mixture does not trace a ray through a
three-dimensional velocity structure. The weights follow Hypoinverse 1.40
HYTRA node order:
first inner circle wins; otherwise use the first three containing outer
circles. One/two nodes receive a default complement if their sum is below
one; three nodes are always normalized. Duplicate model IDs are combined.
Each regional model can be a ray/table model, `PhaseModels` P/S relation or
fixed `StationModels` selection.
The full phase prediction, including intercept and known delay, is mixed.
All branches use one local Cartesian frame and vertical datum.
"""
struct RegionalModels
    models::Dict{String,RegionalModelBranch}
    nodes::Vector{RegionalNode}
    default_model_id::String

    function RegionalModels(models::AbstractDict,nodes::AbstractVector{RegionalNode};
                            default_model_id::AbstractString)
        branches = Dict{String,RegionalModelBranch}()
        for (id,branch) in models
            id isa AbstractString && !isempty(id) && branch isa RegionalModelBranch ||
                throw(ArgumentError("Regional branches must map model IDs to ray/table, phase or fixed station models."))
            branches[String(id)] = branch
        end
        haskey(branches,default_model_id) || throw(ArgumentError("Unknown default regional model."))
        all(node -> haskey(branches,node.model_id),nodes) ||
            throw(ArgumentError("Regional node refers to an unknown model."))
        new(branches,collect(nodes),String(default_model_id))
    end
end

# This is the shallowest possible model top. At a given source position, only
# models with nonzero regional weight restrict its depth; an inactive region must not exclude
# an otherwise valid negative-depth event. Unreferenced models are ignored.
model_top(model::RegionalModels) = minimum(model_top(model.models[id]) for id in
    unique([model.default_model_id; [node.model_id for node in model.nodes]]))

_source_top(model,source) = model_top(model)
_source_top(model::RegionalModels,source) = maximum(model_top(model.models[w.model_id])
    for w in regional_weights(model,source))

function _check_regional_domain(model,source,weights)
    top = maximum(model_top(model.models[w.model_id]) for w in weights)
    source[3] >= top || throw(DomainError(source,"Source above an active regional model top ($top km)."))
end

function _regional_weights(model::RegionalModels,source)
    length(source)==3 && all(isfinite,source) ||
        throw(ArgumentError("Source must contain finite (x,y,z) in km."))
    distances = [hypot(source[1]-node.x_km,source[2]-node.y_km) for node in model.nodes]
    # Circle switches, the three-node cap and the default-complement switch
    # can be nonsmooth or discontinuous for overlapping nodes. Do not claim a
    # unique derivative at them. The solver can explicitly use its secant.
    boundary = false
    weights = NamedTuple{(:model_id,:weight,:dweight_dsource_km),
                         Tuple{String,Float64,NTuple{3,Float64}}}[]
    for (node,d) in zip(model.nodes,distances)
        if d <= node.inner_radius_km
            push!(weights,(model_id=node.model_id,weight=1.,dweight_dsource_km=(0.,0.,0.)))
            return weights,d == node.inner_radius_km ? :boundary : :smooth
        end
    end
    last_selected = length(model.nodes)
    for (i,(node,d)) in enumerate(zip(model.nodes,distances))
        d < node.outer_radius_km || continue
        width = node.outer_radius_km-node.inner_radius_km
        # cospi near the outer circle can round 1+cos(angle) to zero;
        # this equivalent half-angle form retains its small positive weight.
        weight = sinpi((node.outer_radius_km-d)/(2width))^2
        slope = -pi/width*sinpi((node.outer_radius_km-d)/(2width))*
                          sinpi((d-node.inner_radius_km)/(2width))
        gradient = (slope*(source[1]-node.x_km)/d,
                    slope*(source[2]-node.y_km)/d,0.)
        push!(weights,(model_id=node.model_id,weight=weight,dweight_dsource_km=gradient))
        if length(weights)==3
            last_selected = i
            break
        end
    end
    boundary = any(i -> distances[i]==model.nodes[i].outer_radius_km,1:last_selected)
    if isempty(weights)
        push!(weights,(model_id=model.default_model_id,weight=1.,dweight_dsource_km=(0.,0.,0.)))
    else
        total = sum(w.weight for w in weights)
        gradient = ntuple(k -> sum(w.dweight_dsource_km[k] for w in weights),3)
        at_normalization_boundary = length(weights)==2 && abs(total-1)<=8eps(Float64)
        boundary |= at_normalization_boundary
        if length(weights)<3 && total<1 && !at_normalization_boundary
            push!(weights,(model_id=model.default_model_id,weight=1-total,
                           dweight_dsource_km=ntuple(k -> -gradient[k],3)))
        else
            weights = [(model_id=w.model_id,weight=w.weight/total,
                dweight_dsource_km=ntuple(k -> (w.dweight_dsource_km[k]-w.weight/total*gradient[k])/total,3))
                for w in weights]
        end
    end
    combined = typeof(weights)()
    for w in weights
        i = findfirst(other -> other.model_id==w.model_id,combined)
        if i === nothing
            push!(combined,w)
        else
            previous = combined[i]
            combined[i] = (model_id=w.model_id,weight=previous.weight+w.weight,
                dweight_dsource_km=previous.dweight_dsource_km.+w.dweight_dsource_km)
        end
    end
    combined,boundary ? :boundary : :smooth
end

"""Normalized regional model weights and their source-coordinate derivatives."""
regional_weights(model::RegionalModels,source) = first(_regional_weights(model,source))

function _regional_prediction(model::RegionalModels,source,predict)
    weights,status = _regional_weights(model,source)
    _check_regional_domain(model,source,weights)
    components = [merge(w,(prediction=predict(model.models[w.model_id]),)) for w in weights]
    time = sum(c.weight*c.prediction.time_s for c in components)
    derivative = status==:boundary || any(c -> c.prediction.dsource_s_km === nothing,components) ? nothing :
        ntuple(k -> sum(c.weight*c.prediction.dsource_s_km[k]+
                        c.prediction.time_s*c.dweight_dsource_km[k] for c in components),3)
    ray_time = sum(c.weight*get(c.prediction,:ray_time_s,c.prediction.time_s) for c in components)
    offset = sum(c.weight*get(c.prediction,:phase_offset_s,0.) for c in components)
    delay = sum(c.weight*get(c.prediction,:station_delay_s,0.) for c in components)
    # A mixture is an empirical prediction surface, not a single ray. Preserve
    # each component's branch/table diagnostics without inventing a ray angle.
    (;time_s=time,dsource_s_km=derivative,ray_time_s=ray_time,phase_offset_s=offset,
      station_delay_s=delay,components,region_status=status)
end

travel_time(model::RegionalModels,source,station,phase::AbstractString) =
    _regional_prediction(model,source,branch -> travel_time(branch,source,station,phase))

travel_time(model::RegionalModels,source,stations::AbstractDict,station_id::AbstractString,phase::AbstractString) =
    _regional_prediction(model,source,branch -> travel_time(branch,source,stations,station_id,phase))

# Fixed coincident references still retain the weighted b and station delays.
function _forward_time(model::RegionalModels,source,stations::AbstractDict,station_id,phase)
    weights = regional_weights(model,source)
    _check_regional_domain(model,source,weights)
    sum(w.weight*_forward_time(model.models[w.model_id],source,stations,station_id,phase)
        for w in weights)
end

function _forward_time(model::RegionalModels,source,station,phase)
    weights = regional_weights(model,source)
    _check_regional_domain(model,source,weights)
    sum(w.weight*_forward_time(model.models[w.model_id],source,station,phase)
        for w in weights)
end
