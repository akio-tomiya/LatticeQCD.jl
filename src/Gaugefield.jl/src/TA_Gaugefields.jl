import ..SUN_generator:Generator,lie2matrix!,matrix2lie!

abstract type TA_Gaugefields{NC,Dim}  <: AbstractGaugefields{NC,Dim} #Traceless antihermitian matrix
end

include("./4D/TA_gaugefields_4D.jl")

function Base.:*(x::Array{<: TA_Gaugefields{NC,Dim},1},y::Array{<: TA_Gaugefields{NC,Dim},1})  where {NC,Dim}
    s = 0
    for μ=1:Dim
        s += x[μ]*y[μ]
    end

    return s
end


function initialize_TA_Gaugefields(U::Array{<:AbstractGaugefields{NC,Dim},1}) where {NC,Dim}
    F = Array{TA_Gaugefields{NC,Dim},1}(undef,Dim)
    for μ=1:Dim
        F[μ] = initialize_TA_Gaugefields(U[μ])
    end
    return F
end


function initialize_TA_Gaugefields(u::AbstractGaugefields{NC,Dim}) where {NC,Dim}
    if Dim==4
        return TA_Gaugefields(NC,u.NX,u.NY,u.NZ,u.NT)
    else
        error("Dim = $Dim is not supoorted")
    end
end



function Base.setindex!(x::T,v,i...)  where T<: TA_Gaugefields
    error("setindex! is not implemented in type $(typeof(U)) ")
    x.a[i...] = v 
end

function Base.getindex(x::T,i...) where T<: TA_Gaugefields
    error("setindex! is not implemented in type $(typeof(U)) ")
    return x.a[i...]
end

function Base.similar(U::T) where T<: TA_Gaugefields
    error("similar! is not implemented in type $(typeof(U)) ")
end

function Traceless_antihermitian_add!(U::T,factor,temp1) where T <: TA_Gaugefields
    error("Traceless_antihermitian_add! is not implemented in type $(typeof(U)) ")
end

function Traceless_antihermitian!(vout::T,vin::T) where T <: TA_Gaugefields
    error("Traceless_antihermitian! is not implemented in type $(typeof(vout)) ")
end


