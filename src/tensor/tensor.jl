module Tensor
using Einsum
export ⊕, ⊗, ⋆
function ⊕(A::T, B::T) where {T<:AbstractMatrix}
    @einsum C[i, j, k, l] := A[k, j] * B[i, l]
end

function oplus!(C::Array{T,4}, A::Array{T,2}, B::Array{T,2}) where {T<:Number}
    K, J = size(A)
    I, L = size(B)
    for l = 1:L
        for k = 1:K
            for j = 1:J
                for i = 1:I
                    C[i, j, k, l] = A[k, j] * B[i, l]
                end
            end
        end
    end
end

function otimes!(C::Array{T,4}, A::Array{T,2}, B::Array{T,2}) where {T<:Number}
    I, J = size(A)
    K, L = size(B)
    for l = 1:L
        for k = 1:K
            for j = 1:J
                for i = 1:I
                    C[i, j, k, l] = A[i, j] * B[k, l]
                end
            end
        end
    end
end

function ⊗(A::T, B::T) where {T<:AbstractMatrix}
    @einsum C[i, j, k, l] := A[i, j] * B[k, l]
end

function ⋆(A::Array{T1,2}, B::Array{T2,4}) where {T1<:Number,T2<:Number}
    @einsum C[i, l] := A[j, k] * B[i, j, k, l]
end

function star!(
    C::Array{T3,2},
    A::Array{T1,2},
    B::Array{T2,4},
) where {T1<:Number,T2<:Number,T3<:Number}
    J, K = size(A)
    I, JB, KB, L = size(B)
    @assert J == JB "dimension mismatch! J of A = $J and J of B = $JB"
    @assert K == KB "dimension mismatch! K of A = $K and K of B = $KB"
    @. C = 0
    for l = 1:L
        for k = 1:K
            for j = 1:J
                for i = 1:I
                    C[i, l] += A[j, k] * B[i, j, k, l]
                end
            end
        end
    end
end

function star!(
    C::Array{T3,4},
    A::Array{T1,4},
    B::Array{T2,4},
) where {T1<:Number,T2<:Number,T3<:Number}
    I, M, N, L = size(A)
    MB, J, K, NB = size(B)
    @assert N == NB "dimension mismatch! N of A = $N and N of B = $NB"
    @assert M == MB "dimension mismatch! M of A = $M and M of B = $MB"
    @. C = 0
    for l = 1:L
        for k = 1:K
            for j = 1:J
                for i = 1:I
                    for n = 1:N
                        for m = 1:M
                            C[i, j, k, l] += A[i, m, n, l] * B[m, j, k, n]
                        end
                    end
                end
            end
        end
    end
end

function ⋆(A::Array{T1,4}, B::Array{T2,4}) where {T1<:Number,T2<:Number}
    @einsum C[i, j, k, l] := A[i, m, n, l] * B[m, j, k, n]
end

const I0 = [1 0; 0 1]
const IoplusI = I0 ⊕ I0
const IotimesI = I0 ⊗ I0
end

#=
using ..Tensor
A = rand(2,2)
B = rand(2,2)
C1 = A ⊕ B
display(C1)
println("")
Tensor.oplus!(C1,A,B)
display(C1)
println("")
#exit()


A = rand(2,2)
B = rand(2,2)
C = A ⊗ B
display(C)
println("")
C1 = copy(C)
Tensor.otimes!(C1,A,B)
display(C1)
println("C1")
#exit()

println("")
B = rand(2,2)
C = A ⊗ B
D = B ⋆ C
display(D)
println("\tD")
D1 = copy(D)
Tensor.star!(D1,B,C)
display(D1)
println("")


println("")
E = C1 ⋆ C
display(E)

println("\tD")
E1 = copy(E)
Tensor.star!(E1,C1,C)
display(E1)
println("")

Tensor.star!(E1,Tensor.IoplusI,C)
display(E1)
println("")

=#
