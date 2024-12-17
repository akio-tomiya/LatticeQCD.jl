module Universe_module
import ..System_parameters: Params
using Gaugefields
using LatticeDiracOperators
import Gaugefields

struct Univ{Dim,TG,T_FA}
    L::NTuple{Dim,Int64}
    NC::Int64
    Nwing::Int8
    gauge_action::GaugeAction{Dim,TG}
    U::Vector{TG}
    quench::Bool
    fermi_action::T_FA
    cov_neural_net::Union{Nothing,CovNeuralnet{Dim}}
    #Uold::Vector{TG}
    verbose_print::Verbose_print

end

function get_gauge_action(univ::Univ)
    return univ.gauge_action
end

function is_quenched(univ::Univ)
    return univ.quench
end


function Univ(p::Params; MPIparallel=false, PEs=nothing)
    Dim = length(p.L)
    L = Tuple(p.L)
    NC = p.NC
    Nwing = p.Nwing




    if p.initial == "cold" || p.initial == "hot" || p.initial == "one instanton"
        if Dim == 2
            U = Initialize_Gaugefields(NC, Nwing, L..., condition=p.initial)
        else
            U = Initialize_Gaugefields(NC, Nwing, L..., condition=p.initial)
        end
    else
        if Dim == 2
            U = Initialize_Gaugefields(NC, Nwing, L..., condition="cold")
        else
            U = Initialize_Gaugefields(NC, Nwing, L..., condition="cold")
        end

    end
    close(p.load_fp)
    logfilename = pwd() * "/" * p.log_dir * "/" * p.logfile
    verbose_print =
        Verbose_print(p.verboselevel, myid=get_myrank(U[1]), filename=logfilename)

    if p.initial == "cold" || p.initial == "hot" || p.initial == "one instanton"
    else
        println_verbose_level2(verbose_print, ".....  File start")
        println_verbose_level1(verbose_print, "File name is $(p.initial)")
        if p.loadU_format == "ILDG"
            ildg = ILDG(p.initial)
            i = 1
            load_gaugefield!(U, i, ildg, L, NC)
        elseif p.loadU_format == "BridgeText"
            filename = p.initial
            load_BridgeText!(filename, U, L, NC)
        elseif p.loadU_format == "JLD"
            filename = p.initial
            U = loadU(filename)
        elseif p.loadU_format == nothing
            error("loadU_format is not specified. Please add loadU_format in System Control. loadU_format should be ILDG, BridgeText or JLD")
        else
            error("loadU_format should be ILDG, BridgeText or JLD. Now $(p.loadU_format). Please check loadU_format in System Control")
        end
    end
    #println_verbose_level1(verbose_print, ".....  test mode")




    #Uold = similar(U)

    gauge_action = GaugeAction(U)

    if p.use_autogeneratedstaples
        error("p.use_autogeneratedstaples = true is not supported yet!")
    else
        plaqloop = make_loops_fromname("plaquette", Dim=Dim)
        append!(plaqloop, plaqloop')
        βinp = p.β / 2
        push!(gauge_action, βinp, plaqloop)
    end

    TG = eltype(U)
    #println(TG)
    #show(gauge_action)

    if p.Dirac_operator == nothing
        fermi_action = nothing
    else
        params = Dict()
        parameters_action = Dict()

        if p.Dirac_operator == "Staggered"
            x = Initialize_pseudofermion_fields(U[1], "staggered")
            params["Dirac_operator"] = "staggered"
            params["mass"] = p.mass
            parameters_action["Nf"] = p.Nf
        elseif p.Dirac_operator == "Wilson"
            x = Initialize_pseudofermion_fields(U[1], "Wilson", nowing=true)
            params["Dirac_operator"] = "Wilson"
            params["κ"] = p.hop
            params["r"] = p.r
            params["faster version"] = true

        elseif p.Dirac_operator == "Domainwall"
            L5 = p.Domainwall_L5
            #println(L5)
            #error("L5")
            M = p.Domainwall_M
            mass = p.Domainwall_m
            params["Dirac_operator"] = "Domainwall"
            params["mass"] = mass
            params["L5"] = L5
            params["M"] = M
            x = Initialize_pseudofermion_fields(U[1], "Domainwall", L5=L5, nowing=true)
        else
            error("not supported")
        end
        params["eps_CG"] = p.eps
        params["verbose_level"] = p.verboselevel #2
        params["MaxCGstep"] = p.MaxCGstep
        params["boundarycondition"] = p.BoundaryCondition

        D = Dirac_operator(U, x, params)
        fermi_action = FermiAction(D, parameters_action)




    end

    T_FA = typeof(fermi_action)

    if p.smearing_for_fermion == "nothing"
        cov_neural_net = nothing
    elseif p.smearing_for_fermion == "stout"
        cov_neural_net = CovNeuralnet()
        @assert p.stout_numlayers == 1 "stout_numlayers should be one. But now $(p.stout_numlayers)"
        st = STOUT_Layer(p.stout_loops, p.stout_ρ, L)
        push!(cov_neural_net, st)
    else
        error("p.smearing_for_fermion = $(p.smearing_for_fermion) is not supported")
    end

    return Univ{Dim,TG,T_FA}(
        L,
        NC,
        Nwing,
        gauge_action,
        U,
        p.quench,
        fermi_action,
        cov_neural_net,
        verbose_print,
    )

end

function Gaugefields.println_verbose_level1(univ::T, val...) where {T<:Univ}
    println_verbose_level1(univ.verbose_print, val...)
end

function Gaugefields.println_verbose_level2(univ::T, val...) where {T<:Univ}
    println_verbose_level2(univ.verbose_print, val...)
end

function Gaugefields.println_verbose_level3(univ::T, val...) where {T<:Univ}
    println_verbose_level3(univ.verbose_print, val...)
end

end
