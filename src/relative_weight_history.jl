# Keep each attempt's numerical arrays once, and share its stage/scales and
# the original observation identities. Rows are reconstructed on access, so
# complete history does not require a boxed NamedTuple per reading per update.
struct _DDWeightAttempt{W,S,C,T}
    attempt::Int
    successful_updates::Int
    differential_rows::Vector{Int}
    absolute_rows::Vector{Int}
    weights::W
    residual_s::Vector{Float64}
    absolute_weights::Vector{Float64}
    absolute_reasons::Vector{Union{Nothing,Symbol}}
    stage::S
    scale_cc::C
    scale_ct::T
end

struct _DDWeightHistory{D,A} <: AbstractVector{NamedTuple}
    observations::D
    arrivals::A
    attempts::Vector{_DDWeightAttempt}
    ends::Vector{Int}
    full::Bool
end

function _DDWeightHistory(observations,arrivals,full,default_quality)
    identities = [begin
        quality,quality_source = _dd_row_quality(row,default_quality)
        (;row.event_i,row.event_j,row.station_id,row.phase,row.source,quality,quality_source,
            data_class=_dd_observation_class(row))
    end for row in observations]
    absolute_identities = [(;row.event_id,row.station_id,row.phase) for row in arrivals]
    _DDWeightHistory(identities,absolute_identities,_DDWeightAttempt[],Int[],full)
end
Base.size(history::_DDWeightHistory) = (isempty(history.ends) ? 0 : last(history.ends),)
Base.IndexStyle(::Type{<:_DDWeightHistory}) = IndexLinear()

function _dd_save_weight_attempt!(history,attempt,successful,dr,ar,weights,residual,absolute_weights,stage)
    # Only retain columns actually present in the history. Selection masks and
    # temporary quality pairs are not part of a saved observation record.
    columns = (;weights.distance_factor,weights.residual_factor,weights.raw_weights,
        weights.effective_weights,weights.reasons,weights.separation_km)
    block = _DDWeightAttempt(attempt,successful,dr,ar,columns,residual,absolute_weights,
        Union{Nothing,Symbol}[nothing for _ in ar],stage,weights.residual_scale.cc,weights.residual_scale.ct)
    push!(history.attempts,block)
    push!(history.ends,length(history)+length(dr)+length(ar))
    nothing
end

function _dd_weight_history_row(history,block::_DDWeightAttempt,k)
    nd = length(block.differential_rows)
    w = block.weights
    if k<=nd
        j = block.differential_rows[k]
        observation = history.observations[j]
        metadata = _dd_weight_history_metadata(history,block,k)
        row = (;metadata...,
            attempt=block.attempt,iteration=block.successful_updates+1,successful_updates=block.successful_updates,
            input_row=j,row_kind=:differential,event_i=observation.event_i,event_j=observation.event_j,
            station_id=observation.station_id,phase=observation.phase,source=observation.source,
            residual_s=block.residual_s[k],separation_km=w.separation_km[k],
            equation_weight=w.effective_weights[k],removal_reason=w.reasons[k],stage=block.stage,
            residual_scale_cc=block.scale_cc,residual_scale_ct=block.scale_ct)
    else
        a = k-nd
        j = block.absolute_rows[a]
        observation = history.arrivals[j]
        weight = block.absolute_weights[a]
        row = (;attempt=block.attempt,iteration=block.successful_updates+1,successful_updates=block.successful_updates,
            input_row=length(history.observations)+j,row_kind=:absolute,event_i=observation.event_id,event_j=nothing,
            station_id=observation.station_id,phase=observation.phase,source=:absolute,
            quality=nothing,quality_source=nothing,data_class=nothing,class_weight=nothing,
            empirical_base_weight=weight,distance_weight=1.,residual_weight=1.,raw_weight=weight,
            stage_index=block.stage.stage_index,equation_weight=weight,removal_reason=block.absolute_reasons[a],
            residual_s=nothing,separation_km=nothing,stage=block.stage,
            residual_scale_cc=block.scale_cc,residual_scale_ct=block.scale_ct)
    end
    history.full ? row : (;row.attempt,row.input_row,row.row_kind,row.source,
        row.removal_reason,row.empirical_base_weight,row.distance_weight,row.residual_s,
        residual_scale_cc=(status=row.residual_scale_cc.status,),
        residual_scale_ct=(status=row.residual_scale_ct.status,))
end

function _dd_weight_history_metadata(history,position)
    b,k = position
    b==0 && return (;)
    _dd_weight_history_metadata(history,history.attempts[b],k)
end

function _dd_weight_history_metadata(history,block::_DDWeightAttempt,k)
    w = block.weights
    observation = history.observations[block.differential_rows[k]]
    quality = observation.quality
    class_weight = _dd_class_weight(block.stage,observation.data_class)
    base = class_weight==-9 || quality<0 ? 0. : class_weight*quality
    (;quality,quality_source=observation.quality_source,
        data_class=observation.data_class,class_weight,
        empirical_base_weight=base,distance_weight=w.distance_factor[k],
        residual_weight=w.residual_factor[k],raw_weight=w.raw_weights[k],stage_index=block.stage.stage_index)
end

function Base.getindex(history::_DDWeightHistory,k::Int)
    @boundscheck checkbounds(history,k)
    b = searchsortedfirst(history.ends,k)
    offset = b==1 ? 0 : history.ends[b-1]
    _dd_weight_history_row(history,history.attempts[b],k-offset)
end

function Base.iterate(history::_DDWeightHistory,state=(1,1))
    b,k = state
    while b<=length(history.attempts)
        block = history.attempts[b]
        if k<=length(block.differential_rows)+length(block.absolute_rows)
            return _dd_weight_history_row(history,block,k),(b,k+1)
        end
        b,k = b+1,1
    end
    nothing
end

# Rejected airquake attempts keep their original weights but record the exit
# reason on the affected rows. Earlier attempts must remain unchanged.
function _dd_weight_history_exits!(history,row_reasons)
    block = last(history.attempts)
    for (k,j) in enumerate(block.differential_rows)
        row_reasons[j]!==nothing && (block.weights.reasons[k]=row_reasons[j])
    end
    nd = length(history.observations)
    for (k,j) in enumerate(block.absolute_rows)
        row_reasons[nd+j]!==nothing && (block.absolute_reasons[k]=row_reasons[nd+j])
    end
end
