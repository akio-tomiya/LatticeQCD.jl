abstract type AbstractFermionfields_4D{NC} <: AbstractFermionfields{NC,4}
end

function get_origin(a::T1) where T1 <: AbstractFermionfields_4D
    return (1,1,1,1)
end

function clear_fermion!(a::AbstractFermionfields_4D{NC}) where NC 
    n1,n2,n3,n4,n5,n6 = size(a.f)
    @inbounds for i6=1:n6
        for i5=1:n5
            for i4=1:n4
                for i3=1:n3
                    for i2=1:n2
                        @simd for i1=1:NC
                            a.f[i1,i2,i3,i4,i5,i6]= 0
                        end
                    end
                end
            end
        end
    end
end
    
function Base.setindex!(x::T,v,i1,i2,i3,i4,i5,i6)  where T <: AbstractFermionfields_4D
    x.f[i1,i2 + x.NDW,i3 + x.NDW,i4 + x.NDW,i5 + x.NDW,i6] = v
end

function Base.getindex(x::T,i1,i2,i3,i4,i5,i6) where T <: AbstractFermionfields_4D
    return x.f[i1,i2 .+ x.NDW,i3 .+ x.NDW,i4 .+ x.NDW,i5 .+ x.NDW,i6]
end

function Base.getindex(F::Adjoint_fermionfields{T},i1,i2,i3,i4,i5,i6) where T <:AbstractFermionfields_4D #F'
    return conj(F.parent[i1,i2,i3,i4,i5,i6])
end

function Base.setindex!(F::Adjoint_fermionfields{T},v,i1,i2,i3,i4,i5,i6,μ)  where T <: AbstractFermionfields_4D
    error("type $(typeof(F)) has no setindex method. This type is read only.")
end

struct Shifted_fermionfields_4D{NC,T} <: Shifted_fermionfields{NC,4}
    parent::T
    #parent::T
    shift::NTuple{4,Int8}
    NC::Int64

    #function Shifted_Gaugefields(U::T,shift,Dim) where {T <: AbstractGaugefields}
    function Shifted_fermionfields_4D(F::AbstractFermionfields_4D{NC},shift) where NC
        return new{NC,typeof(F)}(F,shift,NC)
    end
end

        #lattice shift
function shift_fermion(F::AbstractFermionfields_4D{NC},ν::T) where {T <: Integer,NC}
    if ν == 1
        shift = (1,0,0,0)
    elseif ν == 2
        shift = (0,1,0,0)
    elseif ν == 3
        shift = (0,0,1,0)
    elseif ν == 4
        shift = (0,0,0,1)
    elseif ν == -1
            shift = (-1,0,0,0)
    elseif ν == -2
            shift = (0,-1,0,0)
    elseif ν == -3
            shift = (0,0,-1,0)
    elseif ν == -4
            shift = (0,0,0,-1)
    end

    return Shifted_fermionfields_4D(F,shift)
end

function shift_fermion(F::TF,shift::NTuple{Dim,T}) where {Dim,T <: Integer,TF <: AbstractFermionfields_4D}
    return shift_fermion(F,shift)
end

function Base.setindex!(F::T,v,i1,i2,i3,i4,i5,i6)  where T <: Shifted_fermionfields_4D
    error("type $(typeof(F)) has no setindex method. This type is read only.")
end

function Base.getindex(F::T,i1,i2,i3,i4,i5,i6)  where T <: Shifted_fermionfields_4D
    return F.parent[i1,i2.+ F.shift[1],i3.+ F.shift[2],i4.+ F.shift[3],i5.+ F.shift[4],i6]
end

function LinearAlgebra.mul!(y::AbstractFermionfields_4D{NC},A::T,x::T3) where {NC,T<:Abstractfields,T3 <:Abstractfermion}
    #@assert NC == x.NC "dimension mismatch! NC in y is $NC but NC in x is $(x.NC)"
    NX = y.NX
    NY = y.NY
    NZ = y.NZ
    NT = y.NT
    NG = y.NG

    @inbounds for ialpha=1:NG
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k1=1:NC
                            y[k1,ix,iy,iz,it,ialpha] = 0
                            @simd for k2=1:NC
                                y[k1,ix,iy,iz,it,ialpha] += A[k1,k2,ix,iy,iz,it]*x[k2,ix,iy,iz,it,ialpha]
                            end
                        end
                    end
                end
            end
        end
    end
end

function LinearAlgebra.mul!(y::AbstractFermionfields_4D{3},A::T,x::T3) where {T<:Abstractfields,T3 <:Abstractfermion}
    #@assert 3 == x.NC "dimension mismatch! NC in y is 3 but NC in x is $(x.NC)"
    NX = y.NX
    NY = y.NY
    NZ = y.NZ
    NT = y.NT
    NG = y.NG

    @inbounds for ialpha=1:NG
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        x1 = x[1,ix,iy,iz,it,ialpha]
                        x2 = x[2,ix,iy,iz,it,ialpha]
                        x3 = x[3,ix,iy,iz,it,ialpha]
                        y[1,ix,iy,iz,it,ialpha] = A[1,1,ix,iy,iz,it]*x1 + 
                                                    A[1,2,ix,iy,iz,it]*x2+ 
                                                    A[1,3,ix,iy,iz,it]*x3
                        y[2,ix,iy,iz,it,ialpha] = A[2,1,ix,iy,iz,it]*x1+ 
                                                    A[2,2,ix,iy,iz,it]*x2 + 
                                                    A[2,3,ix,iy,iz,it]*x3
                        y[3,ix,iy,iz,it,ialpha] = A[3,1,ix,iy,iz,it]*x1+ 
                                                    A[3,2,ix,iy,iz,it]*x2 + 
                                                    A[3,3,ix,iy,iz,it]*x3
                    end
                end
            end
        end
    end
end

function LinearAlgebra.mul!(y::AbstractFermionfields_4D{2},A::T,x::T3) where {T<:Abstractfields,T3 <:Abstractfermion}
    #@assert 2 == x.NC "dimension mismatch! NC in y is 2 but NC in x is $(x.NC)"
    NX = y.NX
    NY = y.NY
    NZ = y.NZ
    NT = y.NT
    NG = y.NG

    @inbounds for ialpha=1:NG
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        x1 = x[1,ix,iy,iz,it,ialpha]
                        x2 = x[2,ix,iy,iz,it,ialpha]
                        y[1,ix,iy,iz,it,ialpha] = A[1,1,ix,iy,iz,it]*x1 + 
                                                    A[1,2,ix,iy,iz,it]*x2
                        y[2,ix,iy,iz,it,ialpha] = A[2,1,ix,iy,iz,it]*x1+ 
                                                    A[2,2,ix,iy,iz,it]*x2

                    end
                end
            end
        end
    end
end

function LinearAlgebra.mul!(y::AbstractFermionfields_4D{NC},A::T,x::T3) where {NC,T<:Number,T3 <:Abstractfermion}
    @assert NC == x.NC "dimension mismatch! NC in y is $NC but NC in x is $(x.NC)"
    NX = y.NX
    NY = y.NY
    NZ = y.NZ
    NT = y.NT
    NG = y.NG

    @inbounds for ialpha=1:NG
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k1=1:NC
                            y[k1,ix,iy,iz,it,ialpha] = A*x[k1,ix,iy,iz,it,ialpha]
                        end
                    end
                end
            end
        end
    end
end

"""
mul!(y,A,x,α,β) -> α*A*x+β*y -> y
"""
function LinearAlgebra.mul!(y::AbstractFermionfields_4D{NC},A::T,x::T3,α::TA,β::TB) where {NC,T<:Number,T3 <:Abstractfermion,TA <: Number,TB <: Number}
    @assert NC == x.NC "dimension mismatch! NC in y is $NC but NC in x is $(x.NC)"
    NX = y.NX
    NY = y.NY
    NZ = y.NZ
    NT = y.NT
    NG = y.NG
    if A == one(A)
        @inbounds for ialpha=1:NG
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            for k1=1:NC
                                y[k1,ix,iy,iz,it,ialpha] = α*x[k1,ix,iy,iz,it,ialpha] + β*y[k1,ix,iy,iz,it,ialpha] 
                            end
                        end
                    end
                end
            end
        end
    else
        @inbounds for ialpha=1:NG
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            for k1=1:NC
                                y[k1,ix,iy,iz,it,ialpha] = A*α*x[k1,ix,iy,iz,it,ialpha] + β*y[k1,ix,iy,iz,it,ialpha] 
                            end
                        end
                    end
                end
            end
        end
    end
end

#Overwrite Y with X*a + Y*b, where a and b are scalars. Return Y.
function LinearAlgebra.axpby!(a::Number, X::T, b::Number, Y::AbstractFermionfields_4D{NC}) where {NC,T <: AbstractFermionfields_4D}
    n1,n2,n3,n4,n5,n6 = size(Y.f)

    @inbounds for i6=1:n6
        for i5=1:n5
            for i4=1:n4
                for i3=1:n3
                    for i2=1:n2
                        @simd for i1=1:NC
                            Y.f[i1,i2,i3,i4,i5,i6] = a*X.f[i1,i2,i3,i4,i5,i6] + b*Y.f[i1,i2,i3,i4,i5,i6]
                        end
                    end
                end
            end
        end
    end
end

function add_fermion!(c::AbstractFermionfields_4D{NC},α::Number,a::T1,β::Number,b::T2) where {NC,T1 <: Abstractfermion,T2 <: Abstractfermion}#c = alpha*a + beta*b
    n1,n2,n3,n4,n5,n6 = size(c.f)

    @inbounds for i6=1:n6
        for i5=1:n5
            for i4=1:n4
                for i3=1:n3
                    for i2=1:n2
                        @simd for i1=1:NC
                            #println(a.f[i1,i2,i3,i4,i5,i6],"\t",b.f[i1,i2,i3,i4,i5,i6] )
                            c.f[i1,i2,i3,i4,i5,i6] += α*a.f[i1,i2,i3,i4,i5,i6] + β*b.f[i1,i2,i3,i4,i5,i6] 
                        end
                    end
                end
            end
        end
    end
    return
end

function LinearAlgebra.dot(a::AbstractFermionfields_4D{NC},b::T2) where {NC, T2 <: AbstractFermionfields_4D}
    NT = a.NT
    NZ = a.NZ
    NY = a.NY
    NX = a.NX
    NG = a.NG

    α = 1
    c = 0.0im
    @inbounds for α=1:NG
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        @simd for ic=1:NC
                            c+= conj(a[ic,ix,iy,iz,it,α])*b[ic,ix,iy,iz,it,α]
                        end
                    end
                end
            end
        end
    end  
    return c
end
