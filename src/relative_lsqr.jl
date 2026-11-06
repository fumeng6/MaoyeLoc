function _dd_lsqr_residuals!(step,data_residual,normal,design,rhs,
        scaled_step,divisors,zero_columns,damping)
    step.=scaled_step./divisors
    step[zero_columns].=0.
    mul!(data_residual,design,step)
    data_residual.-=rhs
    mul!(normal,transpose(design),data_residual)
    normal .= normal./divisors .+ (damping^2).*scaled_step
    normal_norm=norm(normal)
    augmented_norm=hypot(norm(data_residual),damping*norm(scaled_step))
    normal_norm,augmented_norm
end

"""
    _dd_lsqr_step(design, residual; damp=0., rtol=1e-6,
                    max_iterations=0, condition_limit=1e12)

Solve the column-RMS-scaled DD increment problem without dense matrices:
`min ||design*step-residual||² + damp²*||N*step||²`. There are no coordinate
sum rows or hard references. Golub--Kahan LSQR is applied to the augmented
operator `[design/N; damp*I]`, using only products with the supplied design
and its transpose. Starting at zero selects the scaled minimum-norm solution
when damp=0. A zero data column retains its actual zero scale and zero update.

LSQR recurrence estimates trigger explicit residual checks near convergence,
with an additional check every 32 iterations. Every accepted result satisfies
the freshly evaluated augmented and scaled normal residual criteria.
Iteration/condition limits and a nonstationary breakdown are failures. The
LSQR condition estimate describes the augmented scaled operator, not the
exact data condition number or rank. No full spectrum/directions are formed.
`max_iterations=0` chooses max(20,2p), where p is the free-parameter count.
"""
function _dd_lsqr_step(design::AbstractMatrix,residual;
        damp::Real=0.,rtol::Real=1e-6,max_iterations::Integer=0,
        condition_limit::Real=1e12)
    m,p=size(design)
    p>0 && p%4==0 || throw(ArgumentError("DD scaled steps require four columns per event."))
    values=design isa _DDDesign ? design.coefficients : design isa SparseMatrixCSC ? nonzeros(design) : design
    length(residual)==m && all(isfinite,values) && all(isfinite,residual) ||
        throw(ArgumentError("The DD design and matching residual must be finite."))
    isfinite(damp) && damp>=0 && isfinite(rtol) && rtol>0 && max_iterations>=0 &&
        !isnan(condition_limit) && condition_limit>1 ||
        throw(ArgumentError("Require nonnegative DAMP/iteration limit, positive rtol and condition_limit>1."))
    limit=max_iterations==0 ? max(20,2p) : Int(max_iterations)
    damping=Float64(damp)
    rhs=Float64[value for value in residual]
    column_scales=zeros(p)
    if m>0
        if design isa _DDDesign
            _dd_column_scales!(column_scales,design)
        elseif design isa SparseMatrixCSC
            entries=nonzeros(design)
            for j in 1:p
                total=0.
                for k in nzrange(design,j)
                    total+=abs2(entries[k])
                end
                column_scales[j]=sqrt(total/m)
            end
        else
            for j in 1:p
                column_scales[j]=norm(view(design,:,j))/sqrt(m)
            end
        end
    end
    zero_columns=findall(iszero,column_scales)
    divisors=[iszero(scale) ? 1. : scale for scale in column_scales]

    # u belongs to the augmented data space; no (m+p)-by-p matrix is stored.
    u=zeros(m+p)
    data_u,penalty_u=view(u,1:m),view(u,m+1:m+p)
    data_u.=rhs
    bnorm=norm(rhs)
    bnorm>0 && (u./=bnorm)
    v=zeros(p)
    mul!(v,transpose(design),data_u)
    v./=divisors
    alpha=norm(v)
    alpha>0 && (v./=alpha)
    w=copy(v)
    scaled_step=zeros(p)
    step=zeros(p)
    data_residual=-rhs
    forward=zeros(design isa _DDDesign ? 0 : m)
    physical_work=zeros(p)
    transpose_work=zeros(p)
    normal=zeros(p)
    mul!(normal,transpose(design),data_residual)
    normal./=divisors
    normal_norm=norm(normal)
    augmented_norm=bnorm
    operator_norm_estimate=0.
    condition_estimate=0.
    direction_norm_sum=0.
    rhobar,phibar=alpha,bnorm
    iterations=0
    last_checked_iteration=0
    converged=bnorm==0 || alpha==0
    linear_solver_status=bnorm==0 ? :zero_rhs : alpha==0 ? :zero_normal_rhs : :iteration_limit

    if !converged
        for iteration in 1:limit
            iterations=iteration
            previous_alpha=alpha
            physical_work.=v./divisors
            _dd_forward_update!(data_u,forward,design,physical_work,previous_alpha)
            penalty_u.=damping.*v.-previous_alpha.*penalty_u
            beta=norm(u)
            beta>0 && (u./=beta)
            mul!(transpose_work,transpose(design),data_u)
            transpose_work.=transpose_work./divisors.+damping.*penalty_u.-beta.*v
            alpha=norm(transpose_work)
            if alpha>0
                v.=transpose_work./alpha
            else
                fill!(v,0.)
            end
            # Frobenius estimate from the completed projected bidiagonal B_k.
            operator_norm_estimate=hypot(operator_norm_estimate,previous_alpha,beta)
            rho=hypot(rhobar,beta)
            if rho==0 || !isfinite(rho)
                linear_solver_status=:bidiagonal_breakdown
                break
            end
            cosine,sine=rhobar/rho,beta/rho
            theta=sine*alpha
            rhobar=-cosine*alpha
            phi=cosine*phibar
            phibar=sine*phibar
            direction_norm_sum+=sum(abs2,w)/rho^2
            scaled_step.+=(phi/rho).*w
            w.=v.-(theta/rho).*w
            condition_estimate=operator_norm_estimate*sqrt(direction_norm_sum)

            # Recurrence estimates trigger exact checks near convergence.
            # Every accepted result uses the explicit augmented residuals.
            step_norm=norm(scaled_step)
            # LSQR estimates the augmented residual by |phibar| and its
            # normal residual by alpha*|sine*phi| in the same scaled problem.
            estimated_residual=abs(phibar)
            estimated_normal=alpha*abs(sine*phi)
            near_compatible=estimated_residual<=rtol*(bnorm+operator_norm_estimate*step_norm)
            near_stationary=estimated_normal<=rtol*operator_norm_estimate*estimated_residual
            terminal=iteration==limit || condition_estimate>=condition_limit || alpha==0 || beta==0
            finite_state=all(isfinite,scaled_step) && isfinite(step_norm) &&
                isfinite(operator_norm_estimate) && isfinite(estimated_normal)
            check_residual=iteration==1 || iteration%32==0 || near_compatible || near_stationary ||
                terminal || !finite_state
            check_residual || continue
            normal_norm,augmented_norm=_dd_lsqr_residuals!(step,data_residual,normal,design,rhs,
                scaled_step,divisors,zero_columns,damping)
            last_checked_iteration=iteration
            if !(all(isfinite,scaled_step) && isfinite(normal_norm) && isfinite(augmented_norm))
                linear_solver_status=:nonfinite_iteration
                break
            end
            compatible=augmented_norm<=rtol*(bnorm+operator_norm_estimate*norm(scaled_step))
            stationary=normal_norm<=rtol*operator_norm_estimate*augmented_norm
            if compatible || stationary
                converged=true
                linear_solver_status=compatible ? :compatible_residual : :normal_residual
                break
            elseif condition_estimate>=condition_limit
                linear_solver_status=:condition_limit
                break
            elseif alpha==0 || beta==0
                linear_solver_status=:bidiagonal_breakdown
                break
            end
        end
    end
    if last_checked_iteration!=iterations
        normal_norm,augmented_norm=_dd_lsqr_residuals!(step,data_residual,normal,design,rhs,
            scaled_step,divisors,zero_columns,damping)
    end
    scaled_step[zero_columns].=0.
    step.=scaled_step./divisors
    residual_norm=norm(data_residual)
    penalty_norm=damping*norm(scaled_step)
    relative_normal_residual=normal_norm/max(operator_norm_estimate*augmented_norm,floatmin(Float64))
    return (;step,column_scales,zero_columns,scaled_step,
        singular_values=Float64[],retained=falses(0),gains=Float64[],
        scaled_vectors=zeros(p,0),physical_vectors=zeros(p,0),
        retained_rank=nothing,numerical_rank=nothing,
        condition_number=nothing,regularized_condition_number=nothing,retained_condition_number=nothing,
        residual=data_residual,residual_norm,penalty_norm,objective=residual_norm^2+penalty_norm^2,
        method=:damped,damp=damping,svd_relative_cutoff=nothing,backend=:lsqr,
        converged,linear_solver_converged=converged,linear_solver_status,
        iterations,linear_solver_iterations=iterations,linear_solver_iteration_limit=limit,
        normal_residual_norm=normal_norm,relative_normal_residual,
        augmented_residual_norm=augmented_norm,relative_residual=bnorm==0 ? 0. : augmented_norm/bnorm,
        operator_norm_estimate,condition_estimate,linear_solver_condition_estimate=condition_estimate,
        condition_estimate_operator=:augmented_column_scaled_design,rtol=Float64(rtol))
end
