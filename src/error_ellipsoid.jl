# Regularized lower gamma for chi-square ranks 1..3. The half-integer gamma
# constants are analytic, so confidence scaling needs no statistics package.
function _spatial_chisq_cdf(value,rank)
    value<=0 && return 0.
    rank==2 && return -expm1(-value/2)
    a=rank/2
    x=value/2
    term,total=inv(a),inv(a)
    for k in 1:2000
        term*=x/(a+k)
        total+=term
        term<=eps(Float64)*total && break
    end
    gamma=rank==1 ? sqrt(pi) : sqrt(pi)/2
    clamp(exp(-x+a*log(x))*total/gamma,0.,1.)
end

function _spatial_chisq_quantile(probability,rank)
    rank==0 && return 0.
    rank==2 && return -2log1p(-probability)
    lower,upper=0.,1.
    while _spatial_chisq_cdf(upper,rank)<probability
        upper*=2
    end
    for _ in 1:100
        midpoint=(lower+upper)/2
        (midpoint==lower || midpoint==upper) && break
        if _spatial_chisq_cdf(midpoint,rank)<probability
            lower=midpoint
        else
            upper=midpoint
        end
    end
    (lower+upper)/2
end

"""
    error_ellipsoid(uncertainty; probability=.95)

Spatial marginal Gaussian ellipsoid/ellipse of an explicitly supplied
covariance with named local east/north/down parameter axes. Origin is
marginalized by taking the spatial block. Probability uses chi-square of the
spatial dimension, not a one-dimensional normal factor. Singular covariance
gives a conservative region, with its enclosed probability reported. Fixed-depth output is
a two-dimensional conditional ellipse; missing spatial parameters stay missing.
Legacy ERH/ERZ (largest projection of any principal standard axis) and true
ellipsoid projection extents are reported separately.
"""
function error_ellipsoid(uncertainty;probability=.95)
    probability isa Real && isfinite(probability) && 0<probability<1 ||
        throw(ArgumentError("Ellipsoid probability must lie strictly between zero and one."))
    unavailable(reason)=(;status=:not_estimated,reason,probability=Float64(probability),axes=NamedTuple[],
        parameter_names=(),spatial_dimension=0,stochastic_rank=nothing,confidence_scale=nothing,
        horizontal_radius_km=nothing,vertical_extent_km=nothing,legacy_erh_km=nothing,legacy_erz_km=nothing,
        source_assumptions=uncertainty===nothing ? nothing : get(uncertainty,:assumptions,nothing))
    uncertainty===nothing && return unavailable(:covariance_not_requested)
    uncertainty.covariance===nothing && return unavailable(get(uncertainty,:reason,:covariance_unavailable))
    names=uncertainty.parameter_names
    indices=findall(name->name in (:x_km,:y_km,:z_km),collect(names))
    isempty(indices) && return unavailable(:no_spatial_parameters)
    spatial_names=Tuple(names[i] for i in indices)
    length(unique(spatial_names))==length(spatial_names) || throw(ArgumentError("Duplicate spatial covariance axis."))
    C=Matrix{Float64}(uncertainty.covariance[indices,indices])
    all(isfinite,C) && isapprox(C,C';rtol=100eps(Float64),atol=0.) ||
        throw(ArgumentError("Spatial covariance must be finite and symmetric."))
    spectrum=eigen(Symmetric((C+C')/2))
    tolerance=100eps(Float64)*maximum(abs,spectrum.values;init=0.)
    minimum(spectrum.values)>=-tolerance || throw(ArgumentError("Spatial covariance must be positive semidefinite."))
    values=max.(spectrum.values,0.)
    rank=count(>(0),values)
    scale=sqrt(_spatial_chisq_quantile(Float64(probability),length(indices)))
    axes=NamedTuple[]
    legacy_horizontal,legacy_vertical=0.,0.
    for j in reverse(eachindex(values))
        vector=zeros(3)
        for (i,name) in enumerate(spatial_names)
            vector[findfirst(==(name),(:x_km,:y_km,:z_km))]=spectrum.vectors[i,j]
        end
        (vector[3]<0 || vector[3]==0 && (vector[2]<0 || vector[2]==0 && vector[1]<0)) && (vector .*= -1)
        horizontal=hypot(vector[1],vector[2])
        standard=sqrt(values[j])
        azimuth=horizontal<100eps(Float64) ? nothing : mod(rad2deg(atan(vector[1],vector[2])),360.)
        dip=rad2deg(atan(abs(vector[3]),horizontal))
        identifiable=all(k==j || abs(values[k]-values[j])>tolerance for k in eachindex(values))
        push!(axes,(axis=length(axes)+1,standard_error_km=standard,semiaxis_km=scale*standard,
            east_component=vector[1],north_component=vector[2],down_component=vector[3],
            azimuth_deg=azimuth,dip_deg=dip,orientation_identifiable=identifiable))
        legacy_horizontal=max(legacy_horizontal,standard*horizontal)
        legacy_vertical=max(legacy_vertical,standard*abs(vector[3]))
    end
    horizontal_indices=findall(name->name in (:x_km,:y_km),collect(spatial_names))
    horizontal_radius=length(horizontal_indices)==2 ? scale*sqrt(max(maximum(eigvals(Symmetric(C[horizontal_indices,horizontal_indices]))),0.)) : nothing
    zindex=findfirst(==(:z_km),spatial_names)
    vertical=zindex===nothing ? nothing : scale*sqrt(max(C[zindex,zindex],0.))
    (;status=:estimated_conditional,reason=nothing,probability=Float64(probability),axes,
      parameter_names=spatial_names,spatial_dimension=length(indices),stochastic_rank=rank,confidence_scale=scale,
      enclosed_probability=rank==0 ? 1. : _spatial_chisq_cdf(scale^2,rank),
      horizontal_radius_km=horizontal_radius,vertical_extent_km=vertical,
      legacy_erh_km=length(horizontal_indices)==2 ? legacy_horizontal : nothing,
      legacy_erz_km=zindex===nothing ? nothing : legacy_vertical,
      source_assumptions=get(uncertainty,:assumptions,nothing))
end
