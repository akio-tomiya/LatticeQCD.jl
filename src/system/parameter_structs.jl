module Parameter_structs
using REPL.TerminalMenus
import ..Simpleprint: println_rank0
@enum SmearingMethod Nosmearing = 1 STOUT = 2

import QCDMeasurements:
    Fermion_parameters,
    Quench_parameters,
    Wilson_parameters,
    Staggered_parameters,
    Domainwall_parameters,
    initialize_fermion_parameters,
    Measurement_parameters,
    Plaq_parameters,
    Poly_parameters,
    Wilson_loop_parameters,
    Pion_parameters,
    ChiralCondensate_parameters,
    Energy_density_parameters,
    TopologicalCharge_parameters,
    initialize_measurement_parameters,
    prepare_measurement_from_dict,
    construct_Measurement_parameters_from_dict


#=
Base.@kwdef mutable struct System
    BoundaryCondition::Vector{Int64} = [1, 1, 1, -1]
    Nwing::Int64 = 1
    verboselevel::Int64 = 1
    randomseed::Int64 = 111
    L::Vector{Int64} = [4, 4, 4, 4]
    NC::Int64 = 3
    β::Float64 = 5.7
    initialtrj::Int64 = 1
    loadU_format::String = "nothing"
    update_method::String = "HMC"
    loadU_dir::String = ""
    loadU_fromfile::Bool = false
    loadU_filename::String = ""
    initial::String = "cold"
    Dirac_operator::String = "nothing"
    quench::Bool = true
    smearing_for_fermion::String = "nothing"
    stout_numlayers::Union{Nothing,Int64} = nothing
    stout_ρ::Union{Nothing,Array{Float64,1}} = nothing
    stout_loops::Union{Nothing,Array{String,1}} = nothing
    useOR::Bool = false #Over relaxation method for heatbath updates
    numOR::Int64 = 0 #number of Over relaxation method
    βeff::Float64 = 5.7 #effective beta for SLHMC
    firstlearn::Int64 = 1
    Nthermalization::Int64 = 0 #number of updates without measuring
    #smearing::Smearing_parameters = NoSmearing_parameters()
    log_dir::String = ""
    logfile::String = ""
    saveU_format::String = "nothing"
    saveU_dir::String = ""
    saveU_every::Int64 = 1
    isevenodd = true
end
=#



const important_parameters = [
    "L",
    "β",
    "update_method",
    "MDsteps",
    "Δτ",
    "Dirac_operator",
    "fermion_parameters",
    "methodname",
    "measurement_basedir",
    "hasgradientflow",
    "measurement_dir",
    "kinds_of_topological_charge",
    "measurements_for_flow",
    "gradientflow_measurements",
]

function check_important_parameters(key)
    findornot = findfirst(x -> x == key, important_parameters)
    if findornot == nothing
        return false
    else
        return true
    end
end

struct2dict(x::T) where {T} =
    Dict{String,Any}(string(fn) => getfield(x, fn) for fn ∈ fieldnames(T))


function generate_printlist(x::Type)
    pnames = fieldnames(x)
    plist = String[]
    for i = 1:length(pnames)
        push!(plist, String(pnames[i]))
    end
    return plist
end


# Physical setting 
Base.@kwdef mutable struct Print_Physical_parameters
    L::Vector{Int64} = [4, 4, 4, 4]
    β::Float64 = 5.7
    NC::Int64 = 3
    Nthermalization::Int64 = 0
    Nsteps::Int64 = 100
    initial::String = "cold"
    initialtrj::Int64 = 1
    update_method::String = "HMC"
    useOR::Bool = false
    numOR::Int64 = 0
    Nwing::Int64 = 0
end


# Physical setting(fermions)
Base.@kwdef mutable struct Print_Fermions_parameters
    quench::Bool = true
    Dirac_operator::Union{Nothing,String} = nothing
    Clover_coefficient::Float64 = 1.5612
    r::Float64 = 1
    hop::Float64 = 0.141139
    Nf::Int64 = 4
    mass::Float64 = 0.5
    Domainwall_M::Union{Nothing,Float64} = nothing
    Domainwall_m::Union{Nothing,Float64} = nothing
    Domainwall_L5::Union{Nothing,Int64} = nothing
    BoundaryCondition::Array{Int64,1} = [1, 1, 1, -1]
    smearing_for_fermion::String = "nothing"
    stout_numlayers::Union{Nothing,Int64} = nothing
    stout_ρ::Union{Nothing,Array{Float64,1}} = nothing
    stout_loops::Union{Nothing,Array{String,1}} = nothing
    # These three params should be revised.
    M::Union{Nothing,Float64} = nothing
    m::Union{Nothing,Float64} = nothing
    N5::Union{Nothing,Int64} = nothing
end

# System Control
Base.@kwdef mutable struct Print_System_control_parameters
    log_dir::String = ""
    logfile::String = ""
    loadU_format::Union{Nothing,String} = nothing
    loadU_dir::String = ""
    loadU_fromfile::Bool = false
    loadU_filename::String = ""
    saveU_dir::String = ""
    saveU_format::Union{String,Nothing} = nothing
    saveU_every::Int64 = 1
    verboselevel::Int64 = 1
    randomseed::Int64 = 111
    measurement_basedir::String = ""
    measurement_dir::String = ""
    julian_random_number::Bool = false
    isevenodd = true
    hasgradientflow::Bool = false
    eps_flow::Float64 = 0.01
    numflow::Int64 = 100
    Nflow::Int64 = 1
end


# HMC related
Base.@kwdef mutable struct Print_HMCrelated_parameters
    Δτ::Float64 = 0.05
    SextonWeingargten::Bool = false
    N_SextonWeingargten::Int64 = 2
    MDsteps::Int64 = 20
    eps::Float64 = 1e-19
    MaxCGstep::Int64 = 3000
    QPQ::Bool = true
end

# Gradient Flow
Base.@kwdef mutable struct Print_Gradientflow_parameters
    hasgradientflow::Bool = false
    eps_flow::Float64 = 0.01
    numflow::Int64 = 100
    Nflow::Int64 = 1
    #gradientflow_measurements::Vector{Dict} = []
end




# Action parameter for SLMC
struct Print_SLMC_parameters
    βeff::Union{Float64,Array{Float64,1}}
    firstlearn::Int64
    use_autogeneratedstaples::Bool
    couplingcoeff::Array{Float64,1}
    couplinglist::Array{String,1}
    coupling_loops::Union{Nothing,Array{Array{Array{Tuple{Int,Int},1},1},1}}
end


# Measurement set
Base.@kwdef mutable struct Print_Measurement_parameters
    measurement_methods::Array{Dict,1} = Dict[]
end

const printlist_physical = generate_printlist(Print_Physical_parameters)
const printlist_fermions = generate_printlist(Print_Fermions_parameters)
const printlist_systemcontrol = generate_printlist(Print_System_control_parameters)
const printlist_HMCrelated = generate_printlist(Print_HMCrelated_parameters)
const printlist_parameters_SLMC = generate_printlist(Print_SLMC_parameters)
const printlist_measurement = generate_printlist(Print_Measurement_parameters)


Base.@kwdef mutable struct Action
    use_autogeneratedstaples::Bool = false
    couplinglist::Vector{String} = []
    couplingcoeff::Vector{ComplexF64} = []
end

Base.@kwdef mutable struct MD
    MDsteps::Int64 = 20
    Δτ::Float64 = 1 / 20
    SextonWeingargten::Bool = false
    N_SextonWeingargten::Int64 = 2
end

abstract type Smearing_parameters end

Base.@kwdef mutable struct NoSmearing_parameters <: Smearing_parameters end

const kindsof_loops = [
    "plaquette",
    "rectangular",
    "chair",
    "polyakov_x",
    "polyakov_y",
    "polyakov_z",
    "polyakov_t",
]

Base.@kwdef mutable struct Stout_parameters <: Smearing_parameters
    numlayers::Int64 = 1
    ρ::Vector{Float64} = []
    stout_loops::Union{Nothing,Array{String,1}} = nothing
end


#=
abstract type Fermion_parameters end

Base.@kwdef mutable struct Quench_parameters <: Fermion_parameters
    Dirac_operator::String = "nothing"
end

Base.@kwdef mutable struct Wilson_parameters <: Fermion_parameters
    Dirac_operator::String = "Wilson"
    hop::Float64 = 0.141139
    r::Float64 = 1
    hasclover::Bool = false
    Clover_coefficient::Float64 = 1.5612
end

Base.@kwdef mutable struct Staggered_parameters <: Fermion_parameters
    Dirac_operator::String = "Staggered"
    mass::Float64 = 0.5 #mass
    Nf::Int64 = 2 #flavor 
end

Base.@kwdef mutable struct Domainwall_parameters <: Fermion_parameters
    Dirac_operator::String = "Domainwall"
    Domainwall_L5::Int64 = 4
    Domainwall_M::Float64 = -1 #mass for Wilson operator which should be negative
    Domainwall_m::Float64 = 0.1 #physical mass
end

function initialize_fermion_parameters(fermion_type)
    if fermion_type == "nothing"
        fermion_parameter = Quench_parameters()
    elseif fermion_type == "Wilson" || fermion_type == "WilsonClover"
        fermion_parameter = Wilson_parameters()
    elseif fermion_type == "Staggered"
        fermion_parameter = Staggered_parameters()
    elseif fermion_type == "Domainwall"
        fermion_parameter = Domainwall_parameters()
    else
        @error "$fermion_type is not implemented in parameter_structs.jl"
    end
    return fermion_parameter
end

=#

Base.@kwdef mutable struct ConjugateGradient
    eps::Float64 = 1e-19
    MaxCGstep::Int64 = 3000
end


#=
abstract type Measurement_parameters end


Base.@kwdef mutable struct Plaq_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Plaquette"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
end

Base.@kwdef mutable struct Poly_parameters <: Measurement_parameters
    methodname::String = "Polyakov_loop"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    #common::Measurement_common_parameters = Measurement_common_parameters()
end

Base.@kwdef mutable struct Energy_density_parameters <: Measurement_parameters
    methodname::String = "Energy_density"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    #common::Measurement_common_parameters = Measurement_common_parameters()
end

Base.@kwdef mutable struct TopologicalCharge_parameters <: Measurement_parameters
    methodname::String = "Topological_charge"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    #common::Measurement_common_parameters = Measurement_common_parameters()
    #numflow::Int64 = 1 #number of flows
    #Nflowsteps::Int64 = 1
    #eps_flow::Float64 = 0.01
    verbose_level::Int64 = 2
    printvalues::Bool = true
    kinds_of_topological_charge::Vector{String} = ["plaquette","clover"]
end

Base.@kwdef mutable struct ChiralCondensate_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Chiral_condensate"
    measure_every::Int64 = 10
    fermiontype::String = "Staggered"
    Nf::Int64 = 4
    eps::Float64 = 1e-19
    mass::Float64 = 0.5
    MaxCGstep::Int64 = 3000
    smearing_for_fermion::String = "nothing"
    stout_numlayers::Union{Nothing,Int64} = nothing
    stout_ρ::Union{Nothing,Array{Float64,1}} = nothing
    stout_loops::Union{Nothing,Array{String,1}} = nothing
    verbose_level::Int64 = 2
    printvalues::Bool = true
    Nr = 10
    #smearing::Smearing_parameters = Stout_parameters()
end

Base.@kwdef mutable struct Pion_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Pion_correlator"
    measure_every::Int64 = 10
    fermiontype::String = "Wilson"
    eps::Float64 = 1e-19
    MaxCGstep::Int64 = 3000
    smearing_for_fermion::String = "nothing"
    stout_numlayers::Union{Nothing,Int64} = nothing
    stout_ρ::Union{Nothing,Array{Float64,1}} = nothing
    stout_loops::Union{Nothing,Array{String,1}} = nothing
    #smearing::Smearing_parameters = NoSmearing_parameters()
    fermion_parameters::Fermion_parameters = Wilson_parameters()
    verbose_level::Int64 = 2
    printvalues::Bool = true
end

Base.@kwdef mutable struct Wilson_loop_parameters <: Measurement_parameters
    #common::Measurement_common_parameters = Measurement_common_parameters()
    methodname::String = "Wilson_loop"
    measure_every::Int64 = 10
    fermiontype::String = "nothing"
    verbose_level::Int64 = 2
    printvalues::Bool = true
    Tmax::Int64 = 4
    Rmax::Int64 = 4
end

=#

#=


function initialize_measurement_parameters(methodname)
    if methodname == "Plaquette"
        method = Plaq_parameters()
    elseif methodname == "Polyakov_loop"
        method = Poly_parameters()
    elseif methodname == "Topological_charge"
        method = TopologicalCharge_parameters()
    elseif methodname == "Chiral_condensate"
        method = ChiralCondensate_parameters()
    elseif methodname == "Pion_correlator"
        method = Pion_parameters()
    elseif methodname == "Energy_density"
        method = Energy_density_parameters()
    elseif methodname == "Wilson_loop"
        method = Wilson_loop_parameters()
    else
        @error "$methodname is not implemented in parameter_structs.jl"
    end
    return method
end

=#


Base.@kwdef mutable struct Measurement_parameterset
    measurement_methods::Vector{Measurement_parameters} = []
    #measurement_basedir::String = ""
    #measurement_dir::String = ""
end

function transform_measurement_dictvec(value)
    #println(value)
    flow_dict = Dict()
    nummeasure = length(value)
    value_out = Vector{Measurement_parameters}(undef, nummeasure)
    hasgradientflow = false
    for i = 1:nummeasure
        if haskey(value[i], "methodname")
            if value[i]["methodname"] == "Topological_charge"
                hasgradientflow = true
                value_out[i] =
                    transform_topological_charge_measurement!(flow_dict, value[i])
            else
                value_out[i] = construct_Measurement_parameters_from_dict(value[i])
            end
        else
            @error("method name in measurement should be set")
        end
    end
    return value_out, flow_dict, hasgradientflow
end

#=

function construct_Measurement_parameters_from_dict(value_i::Dict)
    #println(value)
    @assert haskey(value_i, "methodname") "methodname should be set in measurement."
    methodname = value_i["methodname"]
    method = initialize_measurement_parameters(methodname)
    method_dict = struct2dict(method)
    #println("value_i ",value_i)
    if haskey(value_i, "Dirac_operator")
        fermiontype = value_i["Dirac_operator"]
    else
        if haskey(value_i, "fermiontype")
            if value_i["fermiontype"] == nothing
                fermiontype = "nothing"
            else
                fermiontype = value_i["fermiontype"]
            end
        else
            fermiontype = "nothing"
        end
    end
    #println("fermiontype $fermiontype")
    fermion_parameters = initialize_fermion_parameters(fermiontype)
    fermion_parameters_dict = struct2dict(fermion_parameters)
    #println("femriontype ",fermiontype)

    for (key_ii, value_ii) in value_i
        #println("$key_ii $value_ii")
        if haskey(method_dict, key_ii)
            if typeof(value_ii) != Nothing
                keytype = typeof(getfield(method, Symbol(key_ii)))
                setfield!(method, Symbol(key_ii), keytype(value_ii))
            end
        else
            if haskey(fermion_parameters_dict, key_ii)
                #println("fermion $key_ii $value_ii")
                keytype = typeof(getfield(fermion_parameters, Symbol(key_ii)))
                setfield!(fermion_parameters, Symbol(key_ii), keytype(value_ii))
            else
                @warn "$key_ii is not found! in $(typeof(method))"
            end
        end
    end

    if haskey(method_dict, "fermion_parameters")
        setfield!(method, Symbol("fermion_parameters"), fermion_parameters)
    end
    value_out = deepcopy(method)

    return value_out
end

=#

function transform_topological_charge_measurement!(flow_dict, measurement)
    @assert haskey(measurement, "methodname") "method name in measurement should be set $(measurement)"
    @assert measurement["methodname"] == "Topological_charge" "this function is for topological charge measurement"

    measurement_revised = Dict()

    for (key, value) in measurement
        #println((key,value))
        if key == "Nflowsteps"
            flow_dict["Nflow"] = value
        elseif key == "numflow"
            flow_dict["numflow"] = value
        elseif key == "eps_flow"
            flow_dict["eps_flow"] = value
        else
            measurement_revised[key] = value
        end
    end
    measurement_revised["fermiontype"] = "nothing"
    #println(measurement_revised)

    valuem = construct_Measurement_parameters_from_dict(measurement_revised)
    flow_dict["measurements_for_flow"] = Dict()
    flow_dict["measurements_for_flow"]["Topological_charge"] = measurement_revised


    return valuem


end





function Plaq_parameters_interactive()
    method = Plaq_parameters()
    println_rank0("You measure Plaquette loops")
    method.methodname = "Plaquette"
    method.measure_every =
        parse(Int64, Base.prompt("How often measure Plaquette loops?", default="1"))
    return method
end



function Poly_parameters_interactive()
    method = Poly_parameters()
    println_rank0("You measure Polyakov loops")
    method.methodname = "Polyakov_loop"
    method.measure_every =
        parse(Int64, Base.prompt("How often measure Polyakov loops?", default="1"))
    return method
end

function Wilson_loop_parameters_interactive(L)
    method = Wilson_loop_parameters()
    println_rank0("You measure RxT Wilson loops")
    method.methodname = "Wilson_loop"
    Rmax0 = min(L[1], L[2], L[3]) ÷ 2
    Tmax0 = L[4] ÷ 2
    method.measure_every =
        parse(Int64, Base.prompt("How often measure Wilson loops?", default="1"))
    method.Rmax =
        parse(Int64, Base.prompt("maximum R for RxT Wilson loop?", default="$Rmax0"))
    method.Tmax =
        parse(Int64, Base.prompt("maximum T for RxT Wilson loop??", default="$Tmax0"))
    return method
end


function Energy_density_parameters_interactive()
    method = Energy_density_parameters()
    println_rank0("You measure Energy density")
    method.methodname = "Energy_density"
    method.measure_every =
        parse(Int64, Base.prompt("How often measure Plaquette loops?", default="1"))
    return method
end



function TopologicalCharge_parameters_interactive()
    method = TopologicalCharge_parameters()
    println_rank0("You measure a topological charge")
    method.methodname = "Topological_charge"
    method.measure_every =
        parse(Int64, Base.prompt("How often measure a topological charge?", default="1"))

    #=
    method.numflow = parse(
    Int64,
    Base.prompt(
        "How many times do you want to flow gauge fields to measure the topological charge?",
        default = "10",
    ),
    )
    =#

    #=
    Nflowsteps = 1#L[1]
    eps_flow = 0.01

    method.Nflowsteps = parse(Int64, Base.prompt("Nflowsteps?", default = "$Nflowsteps"))
    method.eps_flow = parse(Float64, Base.prompt("eps_flow?", default = "$eps_flow"))
    =#

    return method
end









function MD_interactive(; Dirac_operator=nothing)
    md = MD()
    println_rank0("Choose parameters for MD")
    MDsteps = parse(Int64, Base.prompt("Input MD steps", default="20"))
    Δτ = parse(Float64, Base.prompt("Input Δτ", default="$(1/MDsteps)"))
    md.MDsteps = MDsteps
    md.Δτ = Δτ

    if Dirac_operator != nothing
        SW = request(
            "Use SextonWeingargten method? multi-time scale",
            RadioMenu(["false", "true"]),
        )
        SextonWeingargten = ifelse(SW == 1, false, true)

        if SextonWeingargten
            N_SextonWeingargten = parse(
                Int64,
                Base.prompt("Input number of SextonWeingargten steps", default="2"),
            )
        else
            N_SextonWeingargten = 2
        end
        md.SextonWeingargten = SextonWeingargten
        md.N_SextonWeingargten = N_SextonWeingargten
    end
    return md

end





function Stout_parameters_interactive()
    stout = Stout_parameters()

    stout_menu = MultiSelectMenu(kindsof_loops)
    choices =
        request("Select the kinds of loops you want to add in stout smearing:", stout_menu)
    count = 0
    ρs = Float64[]
    loops = String[]
    for i in choices
        count += 1
        ρ = parse(
            Float64,
            Base.prompt("coefficient ρ for $(kindsof_loops[i]) loop?", default="0.1"),
        )
        push!(ρs, ρ)
        push!(loops, kindsof_loops[i])
    end
    stout.ρ = ρs
    stout.stout_loops = loops

    return stout
end





function CG_params_interactive()
    cg = ConjugateGradient()
    eps = parse(Float64, Base.prompt("relative error in CG loops", default="1e-19"))
    MaxCGstep =
        parse(Int64, Base.prompt("Maximum iteration steps in CG loops", default="3000"))
    if eps <= 0
        error("Invalid value for eps=$eps. This has to be strictly positive.")
    end
    if MaxCGstep <= 0
        error("Invalid value for MaxCGstep=$MaxCGstep. This has to be strictly positive.")
    end
    cg.eps = eps
    cg.MaxCGstep = MaxCGstep
    return cg
end




function ChiralCondensate_parameters_interactive(; mass=0.5)
    method = ChiralCondensate_parameters()
    println_rank0("You measure chiral condensates with the statteggred fermion")
    method.methodname = "Chiral_condensate"
    method.measure_every =
        parse(Int64, Base.prompt("How often measure chiral condensates?", default="1"))
    method.mass = parse(
        Float64,
        Base.prompt(
            "Input mass for the measurement of chiral condensates",
            default="$mass",
        ),
    )

    method.Nf = 4

    println_rank0(
        "Number of flavors (tastes) for the measurement of chiral condensates is $(method.Nf)",
    )

    eps = parse(Float64, Base.prompt("relative error in CG loops", default="1e-19"))
    MaxCGstep =
        parse(Int64, Base.prompt("Maximum iteration steps in CG loops", default="3000"))
    if eps <= 0
        error("Invalid value for eps=$eps. This has to be strictly positive.")
    end
    if MaxCGstep <= 0
        error("Invalid value for MaxCGstep=$MaxCGstep. This has to be strictly positive.")
    end
    method.eps = eps
    method.MaxCGstep = MaxCGstep

    smearingmethod =
        request(
            "Choose a configuration format for loading",
            RadioMenu(["No smearing", "stout smearing"]),
        ) |> SmearingMethod

    if smearingmethod == Nosmearing
        method.smearing_for_fermion = "nothing"
        smearing = NoSmearing_parameters()
    elseif smearingmethod == STOUT
        method.smearing_for_fermion = "stout"
        smearing = Stout_parameters_interactive()
        method.stout_ρ = smearing.ρ
        method.stout_loops = smearing.stout_loops
        method.stout_numlayers = smearing.numlayers
    end

    return method
end




function Pion_parameters_interactive()
    method = Pion_parameters()
    println_rank0("You measure Pion_correlator")
    method.methodname = "Pion_correlator"
    method.measure_every =
        parse(Int64, Base.prompt("How often measure Pion_correlator?", default="1"))

    wtype = request(
        "Choose fermion type for the measurement of Pion_correlator",
        RadioMenu([
            "Standard Wilson fermion action",
            "Staggered fermion action",
            #"Wilson+Clover fermion action",
        ]),
    )

    if wtype == 1
        println_rank0("Standard Wilson fermion action will be used for the measurement")
        method.fermiontype = "Wilson"
        fermion_parameters, cg = wilson_wizard()
    elseif wtype == 3
        println_rank0("Wilson+Clover fermion action will be used for the measurement")
        method.fermiontype = "WilsonClover"
        fermion_parameters, cg = wilson_wizard()
    elseif wtype == 2
        println_rank0("Staggered fermion action will be used for the measurement")
        method.fermiontype = "Staggered"
        fermion_parameters, cg = staggered_wizard()
    end

    method.eps = cg.eps
    method.MaxCGstep = cg.MaxCGstep

    smearingmethod =
        request(
            "Choose a configuration format for loading",
            RadioMenu(["No smearing", "stout smearing"]),
        ) |> SmearingMethod

    if smearingmethod == Nosmearing
        method.smearing_for_fermion = "nothing"
        smearing = NoSmearing_parameters()
    elseif smearingmethod == STOUT
        method.smearing_for_fermion = "stout"
        smearing = Stout_parameters_interactive()
        method.stout_ρ = smearing.ρ
        method.stout_loops = smearing.stout_loops
        method.stout_numlayers = smearing.numlayers
    end


    return method
end



function wilson_wizard()
    fermion_parameters = Wilson_parameters()

    hop = parse(Float64, Base.prompt("Input the hopping parameter κ", default="0.141139"))
    #hop = parse(Float64,readline(stdin))
    if hop <= 0
        error("Invalid value for κ=$hop. This has to be strictly positive.")
    end
    println_rank0("κ = $hop")
    fermion_parameters.hop = hop

    cg = CG_params_interactive()

    return fermion_parameters, cg
end

function staggered_wizard()
    staggered = Staggered_parameters()
    mass = parse(Float64, Base.prompt("Input mass", default="0.5"))
    if mass <= 0
        error("Invalid value for mass=$mass. This has to be strictly positive.")
    end
    staggered.mass = mass

    Nftype = request(
        "Choose the number of flavors(tastes)",
        RadioMenu([
            "2 (RHMC will be used)",
            "3 (RHMC will be used)",
            "4 (HMC will be used)",
            "8 (HMC will be used)",
            "1 (RHMC will be used)",
        ]),
    )

    if Nftype == 1
        Nf = 2
    elseif Nftype == 2
        Nf = 3
    elseif Nftype == 3
        Nf = 4
    elseif Nftype == 4
        Nf = 8
    elseif Nftype == 5
        Nf = 1
    end

    staggered.Nf = Nf

    cg = CG_params_interactive()

    return staggered, cg

end

function Domainwall_wizard()
    fermion_parameters = Domainwall_parameters()
    N5 =
        parse(Int64, Base.prompt("Input the size of the extra dimension L5", default="4"))
    #fermion_parameters.Domainwall_L5 = N5
    fermion_parameters.N5 = N5
    println_rank0("Standard Domainwall fermion action is used")

    M = parse(Float64, Base.prompt("Input M", default="-1"))
    while M >= 0
        println_rank0("M should be M < 0. ")
        M = parse(Float64, Base.prompt("Input M", default="-1"))
    end
    #fermion_parameters.Domainwall_M = M
    fermion_parameters.M = M

    m = parse(Float64, Base.prompt("Input mass", default="0.25"))
    #fermion_parameters.Domainwall_m = m
    fermion_parameters.m = m

    cg = CG_params_interactive()
    return fermion_parameters, cg
end


#=
function generate_printable_parameters(p::System)
    pnames = fieldnames(System)
    physical = Print_Physical_parameters()
    names_physical = fieldnames(typeof(physical))
    fermions = Print_Fermions_parameters()
    names_fermions = fieldnames(typeof(fermions))
    control = Print_System_control_parameters()
    names_control = fieldnames(typeof(control))
    hmc = Print_HMCrelated_parameters()
    names_hmc = fieldnames(typeof(hmc))
    #measure = Print_Measurement_parameters()
    #names_measure = fieldnames(typeof(measure))

    for pname_i in pnames
        value = getfield(p, pname_i)
        #println_rank0(value)
        hasvalue = false

        physical_index = findfirst(x -> x == pname_i, names_physical)
        if physical_index != nothing
            setfield!(physical, pname_i, value)
            #println_rank0(value, "\t", pname_i)
            hasvalue = true
        end

        fermions_index = findfirst(x -> x == pname_i, names_fermions)
        if fermions_index != nothing
            setfield!(fermions, pname_i, value)
            #println_rank0(value, "\t", pname_i)
            hasvalue = true
        end

        control_index = findfirst(x -> x == pname_i, names_control)
        if control_index != nothing
            setfield!(control, pname_i, value)
            #println_rank0(value, "\t", pname_i)
            hasvalue = true
        end

        hmc_index = findfirst(x -> x == pname_i, names_hmc)
        if hmc_index != nothing
            setfield!(hmc, pname_i, value)
            #println_rank0(value, "\t", pname_i)
            hasvalue = true
        end

        #=
        measure_index = findfirst(x -> x == pname_i, names_measure)
        if measure_index != nothing
            setfield!(measure, pname_i, value)
            #println_rank0(value, "\t", pname_i)
            hasvalue = true
        end
        =#

        #=
        if hasvalue == false
            @warn "$(pname_i) is not set!"
        end
        =#

        #@assert hasvalue "$(pname_i) is not set!"
    end

    return physical, fermions, control, hmc#, measure
end

=#

function construct_printable_parameters_fromdict!(
    key,
    value,
    physical,
    fermions,
    control,
    hmc,
)
    if key == "L"
        value = collect(value)
    elseif key == "r"
        value = Float64(value)
    end

    #println_rank0("$key $value")
    hasvalue = false
    pname_i = Symbol(key)
    physical_index = findfirst(x -> x == key, printlist_physical)
    if physical_index != nothing
        setfield!(physical, pname_i, value)
        hasvalue = true
    end

    fermions_index = findfirst(x -> x == key, printlist_fermions)
    if fermions_index != nothing
        setfield!(fermions, pname_i, value)
        hasvalue = true
    end

    control_index = findfirst(x -> x == key, printlist_systemcontrol)
    if control_index != nothing
        setfield!(control, pname_i, value)
        hasvalue = true
    end

    hmc_index = findfirst(x -> x == key, printlist_HMCrelated)
    if hmc_index != nothing
        setfield!(hmc, pname_i, value)
        hasvalue = true
    end

    #if hasvalue == false
    #    @warn "$(pname_i) is not set!"
    #end

    if hasvalue == false
        @warn "$(key) is not used!"
    end

    return hasvalue
end

function construct_printable_parameters_fromdict!(x::Dict, physical, fermions, control, hmc)


    for (key, value) in x
        hasvalue = false
        pname_i = Symbol(key)
        physical_index = findfirst(x -> x == pname_i, names_physical)
        if physical_index != nothing
            setfield!(physical, pname_i, value)
            hasvalue = true
        end

        fermions_index = findfirst(x -> x == pname_i, names_fermions)
        if fermions_index != nothing
            setfield!(fermions, pname_i, value)
            hasvalue = true
        end

        control_index = findfirst(x -> x == pname_i, names_control)
        if control_index != nothing
            setfield!(control, pname_i, value)
            hasvalue = true
        end

        hmc_index = findfirst(x -> x == pname_i, names_hmc)
        if hmc_index != nothing
            setfield!(hmc, pname_i, value)
            hasvalue = true
        end

        if hasvalue == false
            @warn "$(pname_i) is not set!"
        end
    end
end


function remove_default_values!(x::Dict, defaultsystem)
    for (key, value) in x
        if hasfield(typeof(defaultsystem), Symbol(key))
            default_value = getfield(defaultsystem, Symbol(key))
            #println_rank0(key, "\t", value, "\t", default_value)
            if value == default_value || string(value) == string(default_value)
                if check_important_parameters(key) == false
                    delete!(x, key)
                end
            else
                if value == nothing
                    x[key] = "nothing"
                end
            end
        else
            if value == nothing
                x[key] = "nothing"
            end
        end

        if typeof(value) == Vector{Measurement_parameters}
            construct_dict_from_measurement!(x, value)
        end
        if typeof(value) <: Fermion_parameters
            construct_dict_from_fermion!(x, value)
        end
    end
end

function construct_dict_from_fermion!(x, value)
    fermiondic = struct2dict(value)
    #println_rank0("fermiondic",fermiondic)
    fermiondic_default = typeof(value)()
    remove_default_values!(fermiondic, fermiondic_default)
    x["fermion_parameters"] = fermiondic

end


function construct_dict_from_measurement!(x, value)

    measuredic = Dict()
    #println_rank0(value)
    for measure in value
        methoddic = struct2dict(measure)
        measure_struct_default = typeof(measure)()
        remove_default_values!(methoddic, measure_struct_default)
        measuredic[methoddic["methodname"]] = methoddic
    end
    #println_rank0("x",x)
    #println_rank0(measuredic)
    x["measurement_methods"] = measuredic
end

function remove_default_values!(x::Dict)
    physical = Print_Physical_parameters()
    fermions = Print_Fermions_parameters()
    control = Print_System_control_parameters()
    hmc = Print_HMCrelated_parameters()
    #measure = Print_Measurement_parameters()


    #defaultsystem = System()
    for (params, paramsname) in x
        remove_default_values!(x[params], physical)
        remove_default_values!(x[params], fermions)
        remove_default_values!(x[params], control)
        remove_default_values!(x[params], hmc)
        #remove_default_values!(x[params], measure)
    end
end




end
