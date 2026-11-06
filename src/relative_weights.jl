# HypoDD's dimensionless quality and stage factors multiply equation rows.
# These helpers neither infer a statistical sigma nor mutate the input rows.
_dd_weight_option(options,key,default) = options isa AbstractDict ?
    get(options,String(key),default) : get(options,key,default)

# Default stages for the one relative path. Class weights, then residual and
# distance downweighting, with damping on the column-scaled correction.
# These numbers are this program's default, not a copy of one hypoDD control file.
function _dd_default_stages()
    [Dict{String,Any}("iterations"=>4,"wtccp"=>1.,"wtccs"=>1.,"wtctp"=>1.,"wtcts"=>1.,
            "wrcc"=>-9.,"wrct"=>-9.,"wdcc_km"=>-9.,"wdct_km"=>-9.,"damp"=>5.),
        Dict{String,Any}("iterations"=>4,"wtccp"=>1.,"wtccs"=>1.,"wtctp"=>1.,"wtcts"=>0.6,
            "wrcc"=>6.,"wrct"=>6.,"wdcc_km"=>-9.,"wdct_km"=>10.,"damp"=>5.),
        Dict{String,Any}("iterations"=>4,"wtccp"=>1.,"wtccs"=>0.6,"wtctp"=>1.,"wtcts"=>0.6,
            "wrcc"=>6.,"wrct"=>6.,"wdcc_km"=>8.,"wdct_km"=>8.,"damp"=>5.)]
end

const _DDWeightStage=NamedTuple{
    (:stage_index,:iterations,:start_update,:end_update,:wtccp,:wtccs,:wtctp,:wtcts,
     :wrcc,:wrct,:wdcc_km,:wdct_km,:damp),
    Tuple{Int,Int,Int,Int,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64}}

function _dd_weight_schedule(options)
    options isa AbstractVector && !isempty(options) ||
        throw(ArgumentError("DD weighting stages must be a nonempty array of explicit stage tables."))
    stages = _DDWeightStage[]
    completed = 0
    for (index,options_stage) in enumerate(options)
        iterations = _dd_weight_option(options_stage,:iterations,nothing)
        iterations isa Integer && iterations>0 || throw(ArgumentError("Each DD stage needs positive integer iterations."))
        settings = (;(key=>Float64(_dd_weight_option(options_stage,key,default)) for (key,default) in
            ((:wtccp,1.),(:wtccs,1.),(:wtctp,1.),(:wtcts,1.),(:wrcc,-9.),(:wrct,-9.),
             (:wdcc_km,-9.),(:wdct_km,-9.),(:damp,0.)))...)
        all(isfinite,values(settings)) || throw(ArgumentError("DD stage values must be finite."))
        all(value==-9 || value>=0 for value in (settings.wtccp,settings.wtccs,settings.wtctp,settings.wtcts,
            settings.wrcc,settings.wrct)) || throw(ArgumentError("DD weights and WR need nonnegative values or -9."))
        settings.damp>=0 || throw(ArgumentError("The DD damping magnitude must be nonnegative."))
        all(value==-9 || value>0 for value in (settings.wdcc_km,settings.wdct_km)) ||
            throw(ArgumentError("DD distance cutoffs are positive km or -9 (disabled)."))
        push!(stages,(;stage_index=index,iterations=Int(iterations),start_update=completed,
            end_update=completed+iterations,settings...))
        completed += iterations
    end
    stages
end

# The counter is the number of accepted whole-group updates. An airquake
# rejection/retry does not advance this counter or consume a stage iteration.
function _dd_weight_stage(schedule,successful_updates)
    successful_updates isa Integer && successful_updates>=0 || throw(ArgumentError("DD successful_updates must be nonnegative."))
    index=findfirst(stage->successful_updates<stage.end_update,schedule)
    index===nothing ? nothing : schedule[index]
end

function _dd_observation_class(row)
    row.source in ("CC","CT") && row.phase in ("P","S") ||
        throw(ArgumentError("DD empirical rows require source CC/CT and phase P/S."))
    row.source=="CC" ? (row.phase=="P" ? :ccp : :ccs) : (row.phase=="P" ? :ctp : :cts)
end

_dd_class_weight(stage,class) = class==:ccp ? stage.wtccp : class==:ccs ? stage.wtccs :
    class==:ctp ? stage.wtctp : stage.wtcts

function _dd_row_quality(row,default_quality)
    for key in (:original_quality,:quality)
        if hasproperty(row,key)
            value=getproperty(row,key)
            value isa Real && isfinite(value) || throw(ArgumentError("DD quality must be a finite numerical row factor."))
            return Float64(value),key
        end
    end
    default_quality===nothing && throw(ArgumentError("DD empirical weighting requires explicit quality; use default_quality explicitly for native unit weights."))
    default_quality isa Real && isfinite(default_quality) && default_quality>=0 ||
        throw(ArgumentError("The explicit default DD quality must be finite and nonnegative."))
    Float64(default_quality),:explicit_default
end

function _dd_source_separation(sources,row)
    # Source values are XYZ tuples/vectors in the existing local km frame.
    # Catalog and current source dictionaries are deliberately separate inputs.
    a,b=sources[row.event_i],sources[row.event_j]
    sqrt(sum((a[k]-b[k])^2 for k in 1:3))
end

function _dd_initial_distance_selection(observations,catalog_sources,first_stage)
    distances=Float64[]
    reasons=Union{Nothing,Symbol}[]
    for row in observations
        _dd_observation_class(row)
        d=Float64(_dd_source_separation(catalog_sources,row))
        cutoff=row.source=="CC" ? first_stage.wdcc_km : first_stage.wdct_km
        push!(distances,d)
        push!(reasons,cutoff>0 && d>cutoff ? :initial_pair_distance_cutoff : nothing)
    end
    (;retained_mask=isnothing.(reasons),reasons,separation_km=distances)
end

function _dd_weight_median(values)
    n=length(values)
    n==0 && return nothing
    ordered=copy(values)
    isodd(n) && return partialsort!(ordered,(n+1)÷2)
    middle=partialsort!(ordered,n÷2:n÷2+1)
    (middle[1]+middle[2])/2
end

function _dd_residual_scale(residual,wr)
    isempty(residual) && return (median_s=nothing,mad_s=nothing,scale_s=nothing,cutoff_s=nothing,status=:no_observations)
    wr<=0 && return (median_s=nothing,mad_s=nothing,scale_s=nothing,cutoff_s=nothing,status=:disabled)
    wr<1 && return (median_s=nothing,mad_s=nothing,scale_s=nothing,cutoff_s=Float64(wr),status=:absolute_seconds)
    center=_dd_weight_median(residual)
    mad=_dd_weight_median(abs.(residual.-center))
    scale=mad/.67449
    (;median_s=center,mad_s=mad,scale_s=scale,cutoff_s=wr*scale,status=mad==0 ? :zero_mad : :mad)
end

function _dd_empirical_weights(observations,sources,residual_s,stage;min_weight=1e-5,default_quality=nothing,
        quality_metadata=nothing,separation_km=nothing)
    length(observations)==length(residual_s) && all(isfinite,residual_s) ||
        throw(ArgumentError("DD weights require one finite residual in seconds per observation."))
    min_weight isa Real && isfinite(min_weight) && min_weight>0 || throw(ArgumentError("DD min_weight must be positive."))
    prepared_quality=if quality_metadata===nothing
        quality_pairs=[_dd_row_quality(row,default_quality) for row in observations]
        (classes=_dd_observation_class.(observations),quality=first.(quality_pairs),
            quality_source=last.(quality_pairs))
    else
        quality_metadata
    end
    classes=prepared_quality.classes
    quality=prepared_quality.quality
    quality_source=prepared_quality.quality_source
    length(classes)==length(quality)==length(quality_source)==length(observations) ||
        throw(DimensionMismatch("DD quality metadata must match the observation rows."))
    class_weight=[_dd_class_weight(stage,class) for class in classes]
    # Negative quality is explicitly excluded; it is never converted to a
    # positive row multiplier or used as a pairing preference.
    base=[class_weight[i]==-9 || quality[i]<0 ? 0. : class_weight[i]*quality[i] for i in eachindex(observations)]
    separation=separation_km===nothing ? Float64[_dd_source_separation(sources,row) for row in observations] : separation_km
    length(separation)==length(observations) || throw(DimensionMismatch("DD distances must match the observation rows."))
    scale_cc=_dd_residual_scale(Float64[residual_s[i] for i in eachindex(observations) if classes[i] in (:ccp,:ccs)],stage.wrcc)
    scale_ct=_dd_residual_scale(Float64[residual_s[i] for i in eachindex(observations) if classes[i] in (:ctp,:cts)],stage.wrct)
    n=length(observations)
    distance_factor,residual_factor=ones(n),ones(n)
    reasons=Union{Nothing,Symbol}[nothing for _ in observations]
    raw=copy(base)
    @inbounds for i in eachindex(observations)
        cc=classes[i] in (:ccp,:ccs)
        wd=cc ? stage.wdcc_km : stage.wdct_km
        if wd!=-9
            exponent=cc ? 5 : 3
            distance_factor[i]=separation[i]>=wd ? 0. : (1-(separation[i]/wd)^exponent)^exponent
        end
        raw[i]*=distance_factor[i]
        scale=cc ? scale_cc : scale_ct
        # Match the original application threshold after distance weighting.
        if scale.cutoff_s!==nothing && raw[i]>1e-6
            cutoff=scale.cutoff_s
            r=abs(residual_s[i])
            residual_factor[i]=cutoff==0 ? (r==0 ? 1. : 0.) : r>=cutoff ? 0. : (1-(r/cutoff)^3)^3
        end
        raw[i]*=residual_factor[i]
        reasons[i]=class_weight[i]==-9 ? :data_class_disabled : quality[i]<0 ? :negative_quality_weight :
            quality[i]==0 ? :zero_quality_weight : class_weight[i]==0 ? :zero_class_weight :
            distance_factor[i]==0 ? :pair_distance_cutoff :
            residual_factor[i]==0 ? (scale.status==:zero_mad ? :zero_mad_nonzero_residual : :residual_cutoff) :
            raw[i]<min_weight ? :below_min_weight : nothing
    end
    retained=raw.>=min_weight
    effective=[retained[i] ? raw[i] : 0. for i in eachindex(raw)]
    (;classes,quality,quality_source,class_weight,base,distance_factor,residual_factor,
      raw_weights=raw,effective_weights=effective,retained_mask=retained,reasons,
      separation_km=separation,residual_scale=(cc=scale_cc,ct=scale_ct),min_weight=Float64(min_weight))
end
