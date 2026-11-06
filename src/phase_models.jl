# Complete existing ray/table definitions used as the P and S components below.
const PhaseModelLeaf = Union{HomogeneousModel,LayeredModel,LinearGradientModel,
                            FiniteGradientModel,PiecewiseGradientModel,TravelTimeTableModel,RayTransformModel}

"""
    PhaseModels(p_model; s_scale, s_offset_s=0)
    PhaseModels(p_model, s_model)

An explicit P/S relation. P uses the P channel of p_model. In affine mode,
S=a*T_P+b and its source gradient is a*grad(T_P). Here a=`s_scale` is
dimensionless and b=`s_offset_s` is an empirical time offset in seconds,
not a station delay or a velocity layer. A StationModels wrapper derives its
S delay as a*d_P. In independent mode S uses the S channel of s_model, with its
own structure and explicit S delay. The unused channels do not contribute.
No regional weights, station selection, or observation errors are inferred.
"""
struct PhaseModels
    p_model::PhaseModelLeaf
    s_model::Union{Nothing,PhaseModelLeaf}
    s_scale::Float64
    s_offset_s::Float64

    function PhaseModels(p_model::PhaseModelLeaf; s_scale::Real, s_offset_s::Real=0)
        a,b = Float64(s_scale),Float64(s_offset_s)
        isfinite(a) && a > 0 && isfinite(b) ||
            throw(ArgumentError("An affine S relation requires positive finite s_scale and finite s_offset_s."))
        new(p_model,nothing,a,b)
    end
    PhaseModels(p_model::PhaseModelLeaf,s_model::PhaseModelLeaf) = new(p_model,s_model,1.,0.)
end

model_top(model::PhaseModels) = model.s_model === nothing ? model_top(model.p_model) :
    max(model_top(model.p_model),model_top(model.s_model))

function travel_time(model::PhaseModels,source_km,station_km,phase::AbstractString)
    phase in ("P","S") || throw(ArgumentError("Phase must be P or S."))
    length(source_km) == 3 || throw(ArgumentError("Source must contain (x,y,z) in km."))
    source_km[3] >= model_top(model) || throw(ArgumentError("Source lies above the deeper of the P- and S-model tops."))
    if phase == "P"
        ray = travel_time(model.p_model,source_km,station_km,"P")
        return merge(ray,(ray_time_s=get(ray,:ray_time_s,ray.time_s),phase_offset_s=get(ray,:phase_offset_s,0.),phase_relation=:p))
    elseif model.s_model !== nothing
        ray = travel_time(model.s_model,source_km,station_km,"S")
        return merge(ray,(ray_time_s=get(ray,:ray_time_s,ray.time_s),phase_offset_s=get(ray,:phase_offset_s,0.),phase_relation=:independent_s))
    end
    ray = travel_time(model.p_model,source_km,station_km,"P")
    a,b = model.s_scale,model.s_offset_s
    scaled_time = a*ray.time_s
    result = merge(ray,(time_s=scaled_time+b,
        dsource_s_km=ray.dsource_s_km === nothing ? nothing : a.*ray.dsource_s_km,
        ray_time_s=scaled_time,phase_offset_s=b,phase_relation=:affine_s))
    # Table derivatives and node-error bounds describe the selected phase.
    # Multiplying every time by a positive scale preserves angular direction.
    if hasproperty(ray,:dtdr_s_km)
        bounds = ray.node_error_bounds
        result = merge(result,(dtdr_s_km=a*ray.dtdr_s_km,dtdz_s_km=a*ray.dtdz_s_km,
            node_error_bounds=(time_s=a*bounds.time_s,dtdr_s_km=a*bounds.dtdr_s_km,
                               dtdz_s_km=a*bounds.dtdz_s_km)))
    end
    result
end

# Coincidence removes geometric time, but not an empirical affine intercept.
_coincident_time_s(model,phase) = 0.0
_coincident_time_s(model::PhaseModels,phase) =
    phase == "S" && model.s_model === nothing ? model.s_offset_s : 0.0
