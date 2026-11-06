mutable struct _RelativeStatisticsGroup
    input_count::Int
    absolute_input_count::Int
    counts::Vector{Int}
    squares::Vector{Float64}
    absolute_counts::Vector{Int}
    absolute_square::Float64
    unused::Dict{Tuple{String,String,Symbol},Int}
end
_RelativeStatisticsGroup()=_RelativeStatisticsGroup(0,0,zeros(Int,4),zeros(2),zeros(Int,2),0.,
    Dict{Tuple{String,String,Symbol},Int}())

function _relative_unused_key(row,absolute)
    original=get(row,:unused_reason,nothing)
    reason=original===nothing || ismissing(original) || isempty(String(original)) ? :unspecified : Symbol(original)
    (absolute ? "ABS" : String(row.source),String(row.phase),reason)
end

function _relative_record_dd!(event_groups,component_groups,station_groups,summary,row)
    groups=(event_groups[row.event_i],event_groups[row.event_j],
        get!(_RelativeStatisticsGroup,component_groups,row.component_id),
        get!(_RelativeStatisticsGroup,station_groups,(row.component_id,String(row.station_id))),summary)
    if row.used
        row.residual_s isa Real && isfinite(row.residual_s) ||
            throw(ArgumentError("A used DD row requires a finite residual in seconds."))
        source_index=row.source=="CC" ? 1 : 2
        class_index=2(source_index-1)+(row.phase=="P" ? 1 : 2)
        squared=row.residual_s^2
        for group in groups
            group.input_count+=1
            group.counts[class_index]+=1
            group.squares[source_index]+=squared
        end
    else
        key=_relative_unused_key(row,false)
        for group in groups
            group.input_count+=1
            group.unused[key]=get(group.unused,key,0)+1
        end
    end
    nothing
end

function _relative_record_absolute!(event_groups,component_groups,station_groups,summary,row)
    groups=(event_groups[row.event_id],
        get!(_RelativeStatisticsGroup,component_groups,row.component_id),
        get!(_RelativeStatisticsGroup,station_groups,(row.component_id,String(row.station_id))),summary)
    if row.used
        row.residual_s isa Real && isfinite(row.residual_s) ||
            throw(ArgumentError("A used ABS row requires a finite residual in seconds."))
        phase_index=row.phase=="P" ? 1 : 2
        squared=row.residual_s^2
        for group in groups
            group.absolute_input_count+=1
            group.absolute_counts[phase_index]+=1
            group.absolute_square+=squared
        end
    else
        key=_relative_unused_key(row,true)
        for group in groups
            group.absolute_input_count+=1
            group.unused[key]=get(group.unused,key,0)+1
        end
    end
    nothing
end

"""
    relative_statistics(result)

Count the final `used=true` DD rows by CC-P, CC-S, CT-P and CT-S for every
event, (component, station) and component. A DD row counts once for its station
and component, and once for EACH of its two events. CC/CT RMS is
`sqrt(sum(residual_s^2)/n)` in seconds, without centering or equation weights.
ABS arrivals are counted separately and never enter the DD counts or RMS.
An absent class has `nothing` RMS; unused rows retain their separate reasons.

The station and event definitions follow HypoDD 1.3 `hypoDD.f` station and
source statistics (lines 685–732 and 899–949), using the original residual
seconds rather than the legacy millisecond print conversion. Counts are row
counts, not numbers of statistically independent measurements.
"""
function relative_statistics(result)
    new_group = _RelativeStatisticsGroup
    event_groups = Dict(String(row.event_id)=>new_group() for row in result.events)
    component_groups = Dict{Union{Nothing,Int},typeof(new_group())}()
    station_groups = Dict{Tuple{Union{Nothing,Int},String},typeof(new_group())}()
    summary_group = new_group()
    for row in result.events
        get!(new_group,component_groups,row.component_id)
    end
    for row in result.components
        get!(new_group,component_groups,row.component_id)
    end
    for row in result.residuals
        _relative_record_dd!(event_groups,component_groups,station_groups,summary_group,row)
    end
    for row in result.absolute_residuals
        _relative_record_absolute!(event_groups,component_groups,station_groups,summary_group,row)
    end
    unused_reasons = NamedTuple[]
    function record_unused!(group,scope,component_id,event_id,station_id)
        for ((source,phase,reason),count) in sort!(collect(group.unused);by=first)
            push!(unused_reasons,(;scope,component_id,event_id,station_id,source,phase,reason,count))
        end
    end
    events = NamedTuple[]
    for row in result.events
        group = event_groups[row.event_id]
        push!(events,merge((event_id=String(row.event_id),component_id=row.component_id,status=row.status,
            selection_status=get(row,:selection_status,:selected),exit_reason=get(row,:exit_reason,nothing)),
            _relative_statistics_counts(group)))
        record_unused!(group,:event,row.component_id,String(row.event_id),nothing)
    end
    component_order(id) = id===nothing ? typemax(Int) : id
    stations = NamedTuple[]
    for ((component_id,station_id),group) in sort!(collect(station_groups);by=p->(component_order(first(p)[1]),first(p)[2]))
        push!(stations,merge((;component_id,station_id),_relative_statistics_counts(group)))
        record_unused!(group,:station,component_id,nothing,station_id)
    end
    components = NamedTuple[]
    original_components = Dict(row.component_id=>row for row in result.components)
    for (component_id,group) in sort!(collect(component_groups);by=p->component_order(first(p)))
        original = get(original_components,component_id,(;))
        station_rows = filter(row->row.component_id===component_id,stations)
        push!(components,merge((;component_id,status=get(original,:status,nothing),
            event_count=count(row->row.component_id===component_id,result.events)),
            _relative_statistics_counts(group),_relative_station_maxima(station_rows)))
        record_unused!(group,:component,component_id,nothing,nothing)
    end
    summary = merge((event_count=length(events),station_component_count=length(stations),
        component_count=count(row->row.component_id!==nothing,components),
        event_status_counts=_relative_status_counts(result.events,:status),
        event_selection_status_counts=_relative_status_counts(result.events,:selection_status),
        event_exit_reason_counts=_relative_status_counts(result.events,:exit_reason)),
        _relative_statistics_counts(summary_group),_relative_station_maxima(stations))
    record_unused!(summary_group,:summary,nothing,nothing,nothing)
    assumptions = (row_selection=:final_used_flag,residual_units=:seconds,
        rms_definition=:sqrt_sum_squared_residual_over_count,rms_weighting=:unweighted,
        dd_event_counting=:one_per_endpoint,dd_station_counting=:one_per_row,
        station_grouping=:component_id_and_station_id,absolute_rows=:separate,
        absent_rms=nothing,independent_observations_assumed=false,
        component_identity=:original_result_membership,
        unused_reason_counts=:separate_at_each_scope,
        maximum_tie_rule=:first_in_component_then_station_order)
    (;events,stations,components,summary,unused_reasons,assumptions)
end

function _relative_statistics_counts(group::_RelativeStatisticsGroup)
    counts,squares=group.counts,group.squares
    absolute_counts,absolute_square=group.absolute_counts,group.absolute_square
    cc_count,ct_count = sum(counts[1:2]),sum(counts[3:4])
    used_count = cc_count+ct_count
    absolute_used_count = sum(absolute_counts)
    (input_count=group.input_count,used_count=used_count,
        unused_count=group.input_count-used_count,
        cc_p_count=counts[1],cc_s_count=counts[2],ct_p_count=counts[3],ct_s_count=counts[4],
        cc_count,ct_count,cc_rms_s=cc_count==0 ? nothing : sqrt(squares[1]/cc_count),
        ct_rms_s=ct_count==0 ? nothing : sqrt(squares[2]/ct_count),
        absolute_input_count=group.absolute_input_count,absolute_used_count,
        absolute_unused_count=group.absolute_input_count-absolute_used_count,
        absolute_p_count=absolute_counts[1],absolute_s_count=absolute_counts[2],
        absolute_rms_s=absolute_used_count==0 ? nothing : sqrt(absolute_square/absolute_used_count))
end

function _relative_station_maxima(stations)
    cc = filter(row->row.cc_count>0,stations)
    ct = filter(row->row.ct_count>0,stations)
    cc_max = isempty(cc) ? nothing : cc[argmax(getproperty.(cc,:cc_rms_s))]
    ct_max = isempty(ct) ? nothing : ct[argmax(getproperty.(ct,:ct_rms_s))]
    (max_station_cc_rms_s=cc_max===nothing ? nothing : cc_max.cc_rms_s,
        max_station_cc_station_id=cc_max===nothing ? nothing : cc_max.station_id,
        max_station_cc_component_id=cc_max===nothing ? nothing : cc_max.component_id,
        max_station_ct_rms_s=ct_max===nothing ? nothing : ct_max.ct_rms_s,
        max_station_ct_station_id=ct_max===nothing ? nothing : ct_max.station_id,
        max_station_ct_component_id=ct_max===nothing ? nothing : ct_max.component_id)
end

function _relative_status_counts(events,key)
    counts = Dict{Symbol,Int}()
    for row in events
        value = get(row,key,nothing)
        value===nothing || ismissing(value) || isempty(String(value)) ||
            (counts[Symbol(value)] = get(counts,Symbol(value),0)+1)
    end
    [(reason=reason,count=count) for (reason,count) in sort!(collect(counts);by=first)]
end
