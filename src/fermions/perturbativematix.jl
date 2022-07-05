module Perturbativematrix
using LinearAlgebra

"""
    PDmatrix{T} <: AbstractMatrix{T}
    A = A0 + t A1
    A0:Diagonal matrix
    A1:Offdiagonal matrix
    t:parameter

    A0::AbstractVector{T}
    A1::AbstractMatrix{T}
    t::T
"""
struct PDmatrix <: AbstractMatrix
    A0::AbstractVector
    A1::AbstractMatrix
    t::T

    function PDmatrix(A0, A1, t)
        return new(A0, A1, t)
    end


end

function Base.getindex(A::PDmatrix, i::Int, j::Int)
    if i == j
        diag = A.A0[i, i]
    else
        diag = 0
    end
    return diag + A.t * A.A1[i, j]
end

function Base.:*(A::PDmatrix{T}, B::PDmatrix{T})
    error("not supported yet")
end

function Base.:*(A::PDmatrix, v::T) where {T<:AbstractVector}
    x = A.A0 .* v
    x += A.t * A.A1 * v
    return x
end

function LinearAlgebra.mul!(y, A::PDmatrix, x)
    N = length(y)
    for i = 1:N
        y[i] = A.A0[i] * x[i]
    end
    mul!(y, A.A1, x, A.t, 1)
end

end
