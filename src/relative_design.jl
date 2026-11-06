# Each row stores the two four-parameter event blocks in column order.
struct _DDDesign <: AbstractMatrix{Float64}
    first_columns::Matrix{Int32}
    coefficients::Matrix{Float64}
    column_count::Int
end
Base.size(A::_DDDesign)=(size(A.coefficients,2),A.column_count)
function Base.getindex(A::_DDDesign,row::Int,column::Int)
    @boundscheck checkbounds(A,row,column)
    @inbounds for block in 1:2
        first=Int(A.first_columns[block,row])
        if first>0 && first<=column<first+4
            return A.coefficients[4(block-1)+column-first+1,row]
        end
    end
    0.
end
@inline function _dd_block_product(values,row,offset,x,first,total)
    @inbounds for axis in 0:3
        total=muladd(values[offset+axis,row],x[first+axis],total)
    end
    total
end
@inline function _dd_row_product(A::_DDDesign,row,x)
    @inbounds begin
        first=Int(A.first_columns[1,row]);second=Int(A.first_columns[2,row])
        total=0.
        first>0 && (total=_dd_block_product(A.coefficients,row,1,x,first,total))
        second>0 && (total=_dd_block_product(A.coefficients,row,5,x,second,total))
        total
    end
end
function LinearAlgebra.mul!(y::AbstractVector,A::_DDDesign,x::AbstractVector)
    length(y)==size(A,1) && length(x)==size(A,2) || throw(DimensionMismatch("DD matrix product dimensions differ."))
    @inbounds for row in eachindex(y)
        y[row]=_dd_row_product(A,row,x)
    end
    y
end

function _dd_forward_update!(data_u,forward,design,x,alpha)
    mul!(forward,design,x)
    data_u.=forward.-alpha.*data_u
end
function _dd_forward_update!(data_u,forward,design::_DDDesign,x,alpha)
    length(data_u)==size(design,1) && length(x)==size(design,2) ||
        throw(DimensionMismatch("DD matrix product dimensions differ."))
    @inbounds for row in eachindex(data_u)
        product=_dd_row_product(design,row,x)
        data_u[row]=product-alpha*data_u[row]
    end
    data_u
end
function LinearAlgebra.mul!(y::AbstractVector,At::Transpose{Float64,_DDDesign},x::AbstractVector)
    A=parent(At)
    length(y)==size(A,2) && length(x)==size(A,1) || throw(DimensionMismatch("Transposed DD product dimensions differ."))
    fill!(y,0.)
    @inbounds for row in eachindex(x)
        value=x[row]
        for block in 1:2
            first=Int(A.first_columns[block,row])
            if first>0
                for axis in 0:3
                    column=first+axis
                    y[column]=muladd(A.coefficients[4(block-1)+axis+1,row],value,y[column])
                end
            end
        end
    end
    y
end
LinearAlgebra.mul!(y::AbstractVector,At::Adjoint{Float64,_DDDesign},x::AbstractVector)=mul!(y,transpose(parent(At)),x)
function _dd_column_scales!(scales,A::_DDDesign)
    @inbounds for row in axes(A,1),block in 1:2
        first=Int(A.first_columns[block,row])
        if first>0
            for axis in 0:3
                scales[first+axis]+=abs2(A.coefficients[4(block-1)+axis+1,row])
            end
        end
    end
    for column in eachindex(scales)
        scales[column]=sqrt(scales[column]/size(A,1))
    end
    scales
end
SparseArrays.nnz(A::_DDDesign)=4count(!iszero,A.first_columns)
_dd_design_bytes(A::_DDDesign)=sizeof(A.first_columns)+sizeof(A.coefficients)
_dd_design_bytes(A::SparseMatrixCSC)=sizeof(nonzeros(A))+sizeof(rowvals(A))+sizeof(SparseArrays.getcolptr(A))

function _dd_compact_columns(columns)
    for event in axes(columns,1),axis in 1:4
        first=columns[event,1]
        columns[event,axis]==(first==0 ? 0 : first+axis-1) || return false
    end
    true
end

function _dd_set_design_row!(A,row,column_i,gradient_i,column_j,gradient_j,weight)
    @inbounds if column_i<=column_j
        A.first_columns[1,row]=column_i
        A.first_columns[2,row]=column_j
        offset_i,offset_j=1,5
    else
        A.first_columns[1,row]=column_j
        A.first_columns[2,row]=column_i
        offset_i,offset_j=5,1
    end
    @inbounds for (column,gradient,offset,sign) in ((column_i,gradient_i,offset_i,1.),(column_j,gradient_j,offset_j,-1.))
        if column>0
            gradient=gradient::NTuple{3,Float64}
            for axis in 1:3
                A.coefficients[offset+axis-1,row]=sign*gradient[axis]*weight
            end
            A.coefficients[offset+3,row]=sign*weight
        end
    end
end
