# Explicit coordinate shifts and affine time relation for one physical ray.
# The transform does not introduce a different ray solver.
const UntransformedRayModel = Union{HomogeneousModel,LayeredModel,LinearGradientModel,
    FiniteGradientModel,PiecewiseGradientModel,TravelTimeTableModel}

struct RayTransformModel
    ray_model::UntransformedRayModel
    source_depth_shift_km::Float64
    receiver_depth_shift_km::Float64
    receiver_plane_km::Union{Nothing,Float64}
    time_scale::Float64
    time_offset_s::Float64
    function RayTransformModel(model::UntransformedRayModel; source_depth_shift_km=0.,
            receiver_depth_shift_km=0.,receiver_plane_km=nothing,time_scale=1.,time_offset_s=0.)
        all(isfinite,(source_depth_shift_km,receiver_depth_shift_km,time_scale,time_offset_s)) && time_scale>0 ||
            throw(ArgumentError("Ray transforms need finite shifts/offset and positive time_scale."))
        receiver_plane_km===nothing || isfinite(receiver_plane_km) || throw(ArgumentError("Receiver plane must be finite."))
        new(model,source_depth_shift_km,receiver_depth_shift_km,receiver_plane_km,time_scale,time_offset_s)
    end
end

model_top(model::RayTransformModel)=model_top(model.ray_model)-model.source_depth_shift_km

function travel_time(model::RayTransformModel,source_km,station_km,phase::AbstractString)
    source=(source_km[1],source_km[2],source_km[3]+model.source_depth_shift_km)
    depth=model.receiver_plane_km===nothing ? station_km[3] : model.receiver_plane_km
    station=(station_km[1],station_km[2],depth+model.receiver_depth_shift_km)
    # Coincidence is defined in ray coordinates, which may differ from the
    # external source/station coordinates. A fixed event still has a time.
    ray=try
        travel_time(model.ray_model,source,station,phase)
    catch exception
        exception isa DomainError && source==station || rethrow()
        (time_s=0.,dsource_s_km=nothing,branch=:coincident)
    end
    a,b=model.time_scale,model.time_offset_s
    result=merge(ray,(time_s=a*ray.time_s+b,
        dsource_s_km=ray.dsource_s_km===nothing ? nothing : a.*ray.dsource_s_km,
        ray_time_s=a*ray.time_s,phase_offset_s=b,
        ray_source_depth_km=source[3],ray_receiver_depth_km=station[3]))
    if hasproperty(ray,:dtdr_s_km)
        bounds=ray.node_error_bounds
        result=merge(result,(dtdr_s_km=a*ray.dtdr_s_km,dtdz_s_km=a*ray.dtdz_s_km,
            node_error_bounds=(time_s=a*bounds.time_s,dtdr_s_km=a*bounds.dtdr_s_km,dtdz_s_km=a*bounds.dtdz_s_km)))
    end
    result
end
