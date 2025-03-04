struct StandardHMC{Tmd,TG} <: AbstractUpdate
    md::Tmd
    Uold::Vector{TG}
end



function StandardHMC(
    U,
    gauge_action,
    quench,
    Δτ,
    MDsteps,
    fermi_action,
    cov_neural_net;
    SextonWeingargten=false,
    QPQ=true,
)

    md = MD(
        U,
        gauge_action,
        quench,
        Δτ,
        MDsteps,
        fermi_action,
        cov_neural_net;
        SextonWeingargten=SextonWeingargten,
        QPQ=QPQ,
    )
    Tmd = typeof(md)
    Uold = similar(U)
    TG = eltype(Uold)

    return StandardHMC{Tmd,TG}(md, Uold)

    #error("in StandardHMC!!")
end


function update!(updatemethod::T, U) where {T<:StandardHMC}
    NC = U[1].NC
    md = updatemethod.md
    Uold = updatemethod.Uold
    substitute_U!(Uold, U) #previous configuration

    initialize_MD!(U, md)

    Sp = md.p * md.p / 2
    Sg = -evaluate_GaugeAction(md.gauge_action, U) / NC
    println_verbose_level3(U[1], "Sp_old = $Sp, Sg_old = $Sg")
    Sold = real(Sp + Sg)
    if md.quench == false
        Sfold = real(dot(md.ξ, md.ξ))
        println_verbose_level3(U[1], "Sfold = $Sfold")
        Sold += Sfold
    end

    runMD!(U, md)

    Sp = md.p * md.p / 2
    Sg = -evaluate_GaugeAction(md.gauge_action, U) / NC
    println_verbose_level3(U[1], "Sp_new = $Sp, Sg_new = $Sg")
    Snew = real(Sp + Sg)
    if md.quench == false

        if md.cov_neural_net != Nothing
            Uout, Uout_multi, _ = calc_smearedU(U, md.cov_neural_net)
            Sfnew = evaluate_FermiAction(md.fermi_action, Uout, md.η)
        else
            Sfnew = evaluate_FermiAction(md.fermi_action, U, md.η)
        end


        println_verbose_level3(U[1], "Sfnew = $Sfnew")
        Snew += Sfnew
    end
    println_verbose_level2(U[1], "Sold = $Sold, Snew = $Snew")
    println_verbose_level2(U[1], "Snew - Sold = $(Snew-Sold)")

    accept = exp(Sold - Snew) >= rand()
    if accept
        println_verbose_level2(U[1], "Accepted")
    else
        substitute_U!(U, Uold) #back to previous configuration
        println_verbose_level2(U[1], "Rejected")
    end

    return accept
    #error("updatemethod type $(typeof(updatemethod)) is not supported!!")
end
