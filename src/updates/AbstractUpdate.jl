module AbstractUpdate_module
using Gaugefields
using LinearAlgebra
using LatticeDiracOperators
import ..Universe_module:Univ,get_gauge_action,is_quenched
import ..System_parameters:Params
import ..AbstractMD_module:StandardMD,MD,initialize_MD!,runMD!

abstract type AbstractUpdate end

include("./standardHMC.jl")
include("./givenconfigurations.jl")

function Updatemethod(p::Params,univ::Univ)
    gauge_action = get_gauge_action(univ)
    quench = is_quenched(univ)
    U = univ.U
    updatemethod = Updatemethod(U,gauge_action,update_method,quench,
        p.Δτ,p.MDsteps,SextonWeingargten = SextonWeingargten)
    return updatemethod
end

function Updatemethod(U,gauge_action,update_method,quench,
            Δτ=nothing,MDsteps=nothing
            ;fermi_action = nothing,
            SextonWeingargten=false,
            loadU_dir = nothing,
            loadU_format = nothing
            )
    if update_method == "HMC"
        updatemethod = StandardHMC(U,gauge_action,quench,Δτ,MDsteps,
                fermi_action,
                SextonWeingargten=SextonWeingargten)
    if update_method == "Fileloading"
        updatemethod = GivenConfigurations(U,loadU_dir,loadU_format)
    else
        error("update method $(update_method) is not supported!")
    end

    return updatemethod
end




function update!(updatemethod::T,U) where T <: AbstractUpdate
    error("updatemethod type $(updatemethod) is not supported!!")
end

end