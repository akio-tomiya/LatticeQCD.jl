module AbstractUpdate_module
using Gaugefields
using LinearAlgebra
using LatticeDiracOperators
import ..Universe_module: Univ, get_gauge_action, is_quenched
import ..System_parameters: Params
import ..AbstractMD_module: StandardMD, MD, initialize_MD!, runMD!

abstract type AbstractUpdate end

include("./standardHMC.jl")
include("./givenconfigurations.jl")
include("./heatbath.jl")

#=
function Updatemethod(p::Params, univ::Univ)
    gauge_action = get_gauge_action(univ)
    quench = is_quenched(univ)
    U = univ.U
    updatemethod = Updatemethod(
        U,
        gauge_action,
        update_method,
        quench,
        p.Δτ,
        p.MDsteps,
        SextonWeingargten = SextonWeingargten,
    )
    return updatemethod
end
=#

function Updatemethod(parameters::Params, univ::Univ)
    gauge_action = get_gauge_action(univ)
    quench = is_quenched(univ)

    updatemethod = Updatemethod(
        univ.U,
        gauge_action,
        parameters.update_method,
        quench,
        parameters.Δτ,
        parameters.MDsteps,
        fermi_action=univ.fermi_action,
        cov_neural_net=univ.cov_neural_net,
        SextonWeingargten=parameters.SextonWeingargten,
        loadU_dir=parameters.loadU_dir,
        loadU_format=parameters.loadU_format,
        isevenodd=parameters.isevenodd,
        β=parameters.β,
        ITERATION_MAX=parameters.ITERATION_MAX,
        numOR=parameters.numOR,
        useOR=parameters.useOR,
        QPQ=parameters.QPQ
    )
    return updatemethod
end

function Updatemethod(
    U,
    gauge_action,
    update_method,
    quench,
    Δτ=nothing,
    MDsteps=nothing;
    fermi_action=nothing,
    cov_neural_net=nothing,
    SextonWeingargten=false,
    loadU_dir=nothing,
    loadU_format=nothing,
    isevenodd=true,
    β=5.7,
    ITERATION_MAX=10^5,
    numOR=3,
    useOR=false,
    QPQ=true,
)
    if update_method == "HMC"
        updatemethod = StandardHMC(
            U,
            gauge_action,
            quench,
            Δτ,
            MDsteps,
            fermi_action,
            cov_neural_net,
            SextonWeingargten=SextonWeingargten,
            QPQ=QPQ,
        )
    elseif update_method == "Fileloading"
        updatemethod = GivenConfigurations(U, loadU_dir, loadU_format)
    elseif update_method == "Heatbath"
        updatemethod = Heatbathupdate(
            U,
            gauge_action,
            quench,
            isevenodd=isevenodd,
            β=β,
            ITERATION_MAX=ITERATION_MAX,
            numOR=numOR,
            useOR=useOR,
        )
    else
        error("update method $(update_method) is not supported!")
    end

    return updatemethod
end


function update!(updatemethod::T, U) where {T<:AbstractUpdate}
    error("updatemethod type $(typeof(updatemethod)) is not supported!!")
end

end
