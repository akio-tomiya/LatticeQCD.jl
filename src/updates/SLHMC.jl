struct SLHMC{Tmd,Tmd_eff,TG} <: AbstractUpdate
    md::Tmd
    md_eff::Tmd_eff
    Uold::Vector{TG}
end



function SLHMC(U,gauge_action,quench,Δτ,MDsteps,fermi_action,gauge_eff_action,fermi_eff_action;
                SextonWeingargten=false)

    md = MD(U,gauge_action,quench,Δτ,MDsteps,fermi_action;
                SextonWeingargten=SextonWeingargten)
    Tmd = typeof(md)
    md_eff = MD(U,gauge_eff_action,quench,Δτ,MDsteps,fermi_eff_action;
        SextonWeingargten=SextonWeingargten)
    Tmd_eff = typeof(md_eff)

    Uold = similar(U)
    TG = eltype(Uold)

    return SLHMC{Tmd,Tmd_eff,TG}(md,md_eff,Uold)

    #error("in StandardHMC!!")
end


function update!(updatemethod::T,U) where T <: SLHMC
    NC = U[1].NC
    md = updatemethod.md
    md_eff = updatemethod.md_eff
    Uold = updatemethod.Uold
    substitute_U!(Uold,U) #previous configuration

    initialize_MD!(U,md)

    Sp = md.p*md.p/2
    Sg = -evaluate_GaugeAction(md.gauge_action,U)/NC
    println_verbose_level3(U[1],"Sp_old = $Sp, Sg_old = $Sg")
    Sold = real(Sp +  Sg)
    if md.quench == false
        Sfold = real(dot(md.ξ,md.ξ))
        println_verbose_level3(U[1],"Sfold = $Sfold")
        Sold += Sfold
        md_eff.ξ = md.ξ
        md_eff.η = md.η
    end

    md_eff.p = md.p
    
    runMD!(U,md_eff)
    md.p =md_eff.p

    Sp = md.p*md.p/2
    Sg = -evaluate_GaugeAction(md.gauge_action,U)/NC
    println_verbose_level3(U[1],"Sp_new = $Sp, Sg_new = $Sg")
    Snew = real(Sp +  Sg)
    if md.quench == false
        Sfnew = evaluate_FermiAction(md.fermi_action,U,md.η)
        println_verbose_level3(U[1],"Sfnew = $Sfnew")
        Snew += Sfnew
    end
    println_verbose_level2(U[1],"Sold = $Sold, Snew = $Snew")
    println_verbose_level2(U[1],"Snew - Sold = $(Snew-Sold)")

    accept = exp(Sold - Snew) >= rand()
    if accept
        println_verbose_level2(U[1],"Accepted")
    else
        substitute_U!(U,Uold) #back to previous configuration
        println_verbose_level2(U[1],"Rejected")
    end

    return accept
    #error("updatemethod type $(typeof(updatemethod)) is not supported!!")
end