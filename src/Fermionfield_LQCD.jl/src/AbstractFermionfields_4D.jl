abstract type AbstractFermionfields_4D{NC} <: AbstractFermionfields{NC,4}
end

function clear_fermion!(a::AbstractFermionfields_4D{NC}) where NC 
    n1,n2,n3,n4,n5,n6 = size(a.f)
    for i6=1:n6
        for i5=1:n5
            for i4=1:n4
                for i3=1:n3
                    for i2=1:n2
                        @simd for i1=1:n1
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
    @assert NC == x.NC "dimension mismatch! NC in y is $NC but NC in x is $(x.NC)"
    NX = y.NX
    NY = y.NY
    NZ = y.NZ
    NT = y.NT
    NG = y.NG

    for ialpha=1:NG
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