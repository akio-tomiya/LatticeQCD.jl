module AbstractMD_module
using Gaugefields
using LatticeDiracOperators
using LinearAlgebra
import ..Universe_module:Univ,get_gauge_action,is_quenched
import ..System_parameters:Params

abstract type AbstractMD{Dim,TG} end

include("./standardMD.jl")

function MD(U,gauge_action,quench,Δτ,MDsteps,
                fermi_action = nothing;
                SextonWeingargten=false)
    if SextonWeingargten 
        if quench == true
            error("The quench update does not need the SextonWeingargten method. Put SextonWeingargten = false")
        end
    else
        md = StandardMD(U,gauge_action,quench,Δτ,MDsteps,fermi_action)
    end

    return md
end

function MD(p::Params,univ::Univ)
    gauge_action = get_gauge_action(univ)
    quench = is_quenched(univ)
    md = MD(univ.U,gauge_action,quench,p.Δτ,p.MDsteps,
                SextonWeingargten = p.SextonWeingargten)

    return md
end


function runMD!(U,md::AbstractMD{Dim,TG}) where {Dim,TG}
    error("runMD! with type $(typeof(md)) is not supported")
end

function initialize_MD!(U,md::AbstractMD{Dim,TG}) where {Dim,TG}
    error("initialize_MD! with type $(typeof(md)) is not supported")
end

function U_update!(U,p,ϵ,md::AbstractMD{Dim,TG}) where {Dim,TG}
    temps = get_temporary_gaugefields(md.gauge_action)
    temp1 = temps[1]
    temp2 = temps[2]
    expU = temps[3]
    W = temps[4]

    for μ=1:Dim
        exptU!(expU,ϵ*md.Δτ,p[μ],[temp1,temp2])
        mul!(W,expU,U[μ])
        substitute_U!(U[μ],W)
    end
end

function P_update!(U,p,ϵ,md::AbstractMD{Dim,TG}) where {Dim,TG} # p -> p +factor*U*dSdUμ
    NC = U[1].NC
    temps = get_temporary_gaugefields(md.gauge_action)
    dSdUμ = temps[end]
    factor =  -ϵ*md.Δτ/(NC)
    #factor =  ϵ*md.Δτ/(NC)

    for μ=1:Dim
        calc_dSdUμ!(dSdUμ,md.gauge_action,μ,U)
        mul!(temps[1],U[μ],dSdUμ) # U*dSdUμ
        Traceless_antihermitian_add!(p[μ],factor,temps[1])
    end
    
end

function P_update_fermion!(U,p,ϵ,md::AbstractMD{Dim,TG}) where {Dim,TG}  # p -> p +factor*U*dSdUμ
    #NC = U[1].NC
    temps = get_temporary_gaugefields(md.gauge_action)
    UdSfdUμ = temps[1:Dim]
    factor =  -ϵ*md.Δτ

    calc_UdSfdU!(UdSfdUμ,md.fermi_action,U,md.η)

    for μ=1:Dim
        Traceless_antihermitian_add!(p[μ],factor,UdSfdUμ[μ])
    end
end


end