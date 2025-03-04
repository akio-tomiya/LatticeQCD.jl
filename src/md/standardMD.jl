using InteractiveUtils
import Gaugefields.Temporalfields_module: Temporalfields, get_temp, unused!


struct StandardMD{Dim,TG,TA,quench,T_FA,TF,TC} <: AbstractMD{Dim,TG}
    gauge_action::GaugeAction{Dim,TG}
    quench::Bool
    Δτ::Float64
    MDsteps::Float64
    p::Vector{TA}
    QPQ::Bool
    fermi_action::T_FA
    η::TF
    ξ::TF
    SextonWeingargten::Bool
    Nsw::Int64
    cov_neural_net::TC
    dSdU::Union{Nothing,Vector{TG}}


    function StandardMD(
        U,
        gauge_action::GaugeAction{Dim,TG},
        quench,
        Δτ,
        MDsteps,
        fermi_action=nothing,
        cov_neural_net=nothing;
        QPQ=true,
        SextonWeingargten=false,
        Nsw=2,
    ) where {Dim,TG}
        p = initialize_TA_Gaugefields(U) #This is a traceless-antihermitian gauge fields. This has NC^2-1 real coefficients. 
        TA = eltype(p)
        T_FA = typeof(fermi_action)

        if quench
            η = nothing
            ξ = nothing
            if SextonWeingargten
                error(
                    "The quench update does not need the SextonWeingargten method. Put SextonWeingargten = false",
                )
            end
        else
            if fermi_action == nothing
                η = nothing
                ξ = nothing
            else
                η = similar(fermi_action._temporary_fermionfields[1])
                ξ = similar(η)
            end
        end
        TF = typeof(η)

        @assert Nsw % 2 == 0 "Nsw should be even number! now Nsw = $Nsw"
        TC = typeof(cov_neural_net)
        if TC != Nothing
            dSdU = similar(U)
        else
            dSdU = nothing
        end

        return new{Dim,TG,TA,quench,T_FA,TF,TC}(
            gauge_action,
            quench,
            Δτ,
            MDsteps,
            p,
            QPQ,
            fermi_action,
            η,
            ξ,
            SextonWeingargten,
            Nsw,
            cov_neural_net,
            dSdU,
        )
    end
end

function initialize_MD!(
    U,
    md::StandardMD{Dim,TG,TA,quench,T_FA,TF,TC},
) where {Dim,TG,TA,quench,T_FA,TF,TC}
    gauss_distribution!(md.p) #initial momentum


    if quench == false
        if TC != Nothing
            Uout, Uout_multi, _ = calc_smearedU(U, md.cov_neural_net)
            gauss_sampling_in_action!(md.ξ, Uout, md.fermi_action)
            sample_pseudofermions!(md.η, Uout, md.fermi_action, md.ξ)
        else
            gauss_sampling_in_action!(md.ξ, U, md.fermi_action)
            sample_pseudofermions!(md.η, U, md.fermi_action, md.ξ)
        end

        #error("not supported yet")
    end
end

function runMD!(
    U,
    md::StandardMD{Dim,TG,TA,quench,T_FA,TF,TC},
) where {Dim,TG,TA,quench,T_FA,TF,TC}
    #p = md.p

    if md.QPQ
        if md.SextonWeingargten
            runMD_QPQ_sw!(U, md)
        else
            runMD_QPQ!(U, md)
        end
    else
        if md.SextonWeingargten
            error("PQP update with SextonWeingargten is not supported")
        else
            runMD_PQP!(U, md)
        end
    end

    #error("type $(typeof(md)) is not supported")
end

function runMD_QPQ!(
    U,
    md::StandardMD{Dim,TG,TA,quench,T_FA,TF,TC},
) where {Dim,TG,TA,quench,T_FA,TF,TC}
    p = md.p

    for itrj = 1:md.MDsteps
        U_update!(U, p, 0.5, md)
        P_update!(U, p, 1.0, md)
        if quench == false
            P_update_fermion!(U, p, 1.0, md)
        end
        U_update!(U, p, 0.5, md)
    end

    #error("type $(typeof(md)) is not supported")
end

function runMD_QPQ_sw!(
    U,
    md::StandardMD{Dim,TG,TA,quench,T_FA,TF,TC},
) where {Dim,TG,TA,quench,T_FA,TF,TC}
    p = md.p

    for itrj = 1:md.MDsteps
        for isw = 1:div(md.Nsw, 2)
            U_update!(U, p, 0.5 / md.Nsw, md)
            P_update!(U, p, 1.0 / md.Nsw, md)
            U_update!(U, p, 0.5 / md.Nsw, md)
        end
        if quench == false
            P_update_fermion!(U, p, 1.0, md)
        end
        for isw = 1:div(md.Nsw, 2)
            U_update!(U, p, 0.5 / md.Nsw, md)
            P_update!(U, p, 1.0 / md.Nsw, md)
            U_update!(U, p, 0.5 / md.Nsw, md)
        end
    end

    #error("type $(typeof(md)) is not supported")
end


function runMD_PQP!(
    U,
    md::StandardMD{Dim,TG,TA,quench,T_FA,TF,TC},
) where {Dim,TG,TA,quench,T_FA,TF,TC}
    p = md.p

    for itrj = 1:md.MDsteps
        P_update!(U, p, 0.5, md)
        if quench == false
            P_update_fermion!(U, p, 0.5, md)
        end

        U_update!(U, p, 1.0, md)
        P_update!(U, p, 0.5, md)
        if quench == false
            P_update_fermion!(U, p, 0.5, md)
        end
    end

    #error("type $(typeof(md)) is not supported")
end

function P_update_fermion!(
    U,
    p,
    ϵ,
    md::StandardMD{Dim,TG,TA,quench,T_FA,TF,TC},
) where {Dim,TG,TA,quench,T_FA,TF,TC<:CovNeuralnet{Dim}}  # p -> p +factor*U*dSdUμ
    #NC = U[1].NC

    temps = get_temporary_gaugefields(md.gauge_action)

    UdSfdUμ, its_UdSfdUμ = get_temp(temps, Dim)

    #UdSfdUμ = temps[1:Dim]
    factor = -ϵ * md.Δτ

    Uout, Uout_multi, _ = calc_smearedU(U, md.cov_neural_net)

    for μ = 1:Dim
        calc_UdSfdU!(UdSfdUμ, md.fermi_action, Uout, md.η)
        mul!(md.dSdU[μ], Uout[μ]', UdSfdUμ[μ])
    end
    unused!(temps, its_UdSfdUμ)
    #calc_UdSfdU!(UdSfdUμ, md.fermi_action, U, md.η)

    dSdUbare = back_prop(md.dSdU, md.cov_neural_net, Uout_multi, U)

    temp1, it_temp1 = get_temp(temps)

    for μ = 1:Dim
        #Traceless_antihermitian_add!(p[μ], factor, UdSfdUμ[μ])
        mul!(temp1, U[μ], dSdUbare[μ]) # U*dSdUμ
        Traceless_antihermitian_add!(p[μ], factor, temp1)
    end

    unused!(temps, it_temp1)
end
