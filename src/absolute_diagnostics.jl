"""
    absolute_diagnostics(model, stations, picks, location)

Quantify each arrival's contribution to the locally linearized absolute-location
fit, with the final fitted weights held fixed. For its free columns `A = W*J`,
retain the same singular directions
as the fitted result and report `INFO[i] = sum(U[i,k]^2 for retained k)`.
These are dimensionless projection leverages: their sum is the retained rank,
including the fitted origin direction. A station's contribution is the sum of
its P and S rows. `info_milli` is the original printed 1000*INFO convention.
Distinct `observation_id` values keep repeated same-station/phase readings
separate; station contributions still sum every reading at that station.

Saved empirical equation weights and the final saved spectral cutoff are
used without another reweighting. Fixed parameters do not contribute columns. An active depth bound or temporary depth lock conditions
the diagnostics on z; a fully fixed event has rank and contributions zero.
Iteration damping and singular-value shifts do not change this projection:
they control updates, rather than the frozen local least-squares estimator.

The result describes local information even when the nonlinear solver has
not converged; solution_status is retained explicitly. It is not a timing
variance, covariance, probability, or a statement of location accuracy.
Nonsmooth fitted rays return not_estimated rather than silently inventing a
secant step or assigning a physical derivative at a branch/interface corner.
The supplied model and stations must be those used for location.
"""
function absolute_diagnostics(model,stations::AbstractDict,picks,location)
    observations=collect(picks)
    n=length(observations)
    n==location.input_count==length(location.residuals) ||
        throw(ArgumentError("Diagnostic picks must match the fitted absolute observation count."))
    fitted=Dict(_absolute_observation_key(row)=>row for row in location.residuals)
    length(fitted)==n || throw(ArgumentError("Fitted diagnostic observation identities must be unique."))
    event_ids=unique(pick.event_id for pick in observations if hasproperty(pick,:event_id))
    length(event_ids)<=1 && (isempty(event_ids) || location.event_id===nothing || only(event_ids)==location.event_id) ||
        throw(ArgumentError("Diagnostic picks must belong to the fitted event."))
    constraint=get(location,:constraint,:free)
    fixed_position=constraint in (:fixed_position,:fixed_all)
    fixed_depth=constraint==:fixed_depth || location.depth_bounds[1]==location.depth_bounds[2]
    saved=get(location,:spectral,nothing)
    axes=constraint==:fixed_all ? () : fixed_position ? (4,) : saved!==nothing ? saved.axes :
        fixed_depth || get(location,:depth_constraint_active,false) ? (1,2,4) : (1,2,3,4)
    p=length(axes)
    parameter_names=Tuple((:x_km,:y_km,:z_km,:origin_s)[axis] for axis in axes)
    empirical=get(location,:row_weighting,:empirical)==:empirical
    empirical || throw(ArgumentError("Absolute diagnostics describe the empirical row-factor location."))
    weights=zeros(n)
    seen=Set{Tuple}()
    for (i,pick) in enumerate(observations)
        key=_absolute_observation_key(pick)
        key in seen && throw(ArgumentError("Duplicate absolute observation identity in diagnostic input."))
        push!(seen,key)
        original=get(fitted,key,nothing)
        original!==nothing && String(pick.station_id)==String(original.station_id) &&
            String(pick.phase)==String(original.phase) &&
            pick.arrival_s==original.arrival_s && pick.sigma_s==original.sigma_s &&
            get(pick,:event_id,nothing)==get(original,:event_id,nothing) ||
            throw(ArgumentError("Observation identity, arrival and sigma must match the fitted absolute solution."))
        haskey(stations,pick.station_id) && pick.phase in ("P","S") && isfinite(pick.sigma_s) && pick.sigma_s>0 ||
            throw(ArgumentError("Diagnostic phases need known stations and positive finite sigma metadata."))
        weight=original.equation_weight
        isfinite(weight) && weight>=0 || throw(ArgumentError("Saved diagnostic equation weights must be finite and nonnegative."))
        weights[i]=weight
    end
    selected=findall(>(0),weights)
    assumptions=(quantity=:frozen_weight_retained_projection,weights_fixed=true,
        fitted_estimator=:empirical_absolute_rows,
        branch_fixed=true,location_refitted=false,solution_convergence_required=false,
        fixed_parameters_conditional=p<4,iteration_damping_included=false,
        singular_value_shift_included=false,station_aggregation=:station_id,
        directions_coordinate_order=(:x_km,:y_km,:z_km,:origin_s),
        directions_column_scaling=:physical_km_and_seconds,
        information_scale=:dimensionless,information_milli_factor=1000.)

    function result(info,projections,directions,rank,numerical_rank;status=:computed,reason=nothing)
        rows=[(event_id=get(pick,:event_id,location.event_id),observation_id=get(pick,:observation_id,nothing),
            station_id=String(pick.station_id),
            phase=String(pick.phase),equation_weight=weights[i],included=p>0 && weights[i]>0,
            used=get(fitted[_absolute_observation_key(pick)],:used,false),
            info=info===nothing ? nothing : info[i],
            info_milli=info===nothing ? nothing : 1000info[i],
            retained_projections=projections===nothing ? nothing : Tuple(projections[i,:]))
            for (i,pick) in enumerate(observations)]
        station_ids=unique(String(pick.station_id) for pick in observations)
        contributions=map(station_ids) do id
            ids=findall(row->row.station_id==id,rows)
            phase_info(phase)=info===nothing ? nothing : sum(rows[i].info for i in ids if rows[i].phase==phase;init=0.)
            (station_id=id,phase_count=length(ids),included_count=count(i->rows[i].included,ids),
                p_info=phase_info("P"),s_info=phase_info("S"),
                info=info===nothing ? nothing : sum(rows[i].info for i in ids;init=0.))
        end
        information_sum=info===nothing ? nothing : sum(info)
        (;status,reason,solution_status=location.status,parameter_names,
            conditioned_depth=!(3 in axes),equation_weights=weights,
            row_normalization=get(location,:weight_normalization,1.),
            retained_rank=rank,numerical_rank,information_sum,
            information_sum_minus_rank=information_sum===nothing || rank===nothing ? nothing : information_sum-rank,
            rows,station_contributions=contributions,directions,assumptions)
    end

    p==0 && return result(zeros(n),zeros(n,0),NamedTuple[],0,0;reason=:all_parameters_fixed)
    (get(location,:horizontal_constraint_active,false) || get(location,:horizontal_at_bound,false)) &&
        return result(nothing,nothing,NamedTuple[],nothing,nothing;
            status=:not_estimated,reason=:active_search_boundary)
    isempty(selected) && return result(zeros(n),zeros(n,p),NamedTuple[],0,0;reason=:no_positive_weight_observations)
    A=zeros(length(selected),p)
    for (row,i) in enumerate(selected)
        pick=observations[i]
        if fixed_position
            A[row,1]=weights[i]
        else
            ray=try
                travel_time(model,location.source_km,stations,pick.station_id,pick.phase)
            catch error
                error isa DomainError || rethrow()
                nothing
            end
            (ray===nothing || ray.dsource_s_km===nothing) &&
                return result(nothing,nothing,NamedTuple[],get(location,:retained_rank,nothing),nothing;
                    status=:not_estimated,reason=:nondifferentiable_travel_time)
            gradient=[ray.dsource_s_km...;1.]
            A[row,:] .= weights[i].*gradient[collect(axes)]
        end
    end
    if saved===nothing
        decomposition=svd(A;full=length(selected)<p)
        singular_values=[decomposition.S;zeros(max(0,p-length(decomposition.S)))]
        floor=max(size(A)...)*eps(Float64)*maximum(singular_values;init=0.)
        retained=singular_values.>floor
        numerical_rank=count(retained)
        vectors4=zeros(4,p)
        vectors4[collect(axes),:] .= decomposition.V
    else
        singular_values=saved.singular_values
        retained=saved.retained
        numerical_rank=saved.numerical_rank
        vectors4=saved.vectors4
    end
    projections=zeros(n,p)
    for k in findall(retained)
        projections[selected,k] .= (A*vectors4[collect(axes),k])./singular_values[k]
    end
    info=vec(sum(abs2,projections;dims=2))
    directions=[(direction=k,singular_value=singular_values[k],retained=retained[k],
        vector4=(x_km=vectors4[1,k],y_km=vectors4[2,k],z_km=vectors4[3,k],origin_s=vectors4[4,k])) for k in 1:p]
    result(info,projections,directions,count(retained),numerical_rank)
end
