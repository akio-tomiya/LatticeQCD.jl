module System_parameters
using Random
export Params
import ..Transform_oldinputfile:
    default_system,
    default_md,
    default_defaultmeasures,
    default_actions,
    default_cg,
    default_wilson,
    default_staggered,
    default_measurement


import ..Parameter_structs:
    printlist_physical,
    printlist_fermions,
    printlist_systemcontrol,
    printlist_HMCrelated,
    printlist_parameters_SLMC,
    printlist_measurement

const printlists = [
    printlist_physical,
    printlist_fermions,
    printlist_systemcontrol,
    printlist_HMCrelated,
    printlist_parameters_SLMC,
    printlist_measurement,
]


system = default_system()
md = default_md()
defaultmeasures = default_defaultmeasures()
actions = default_actions()
cg = default_cg()
wilson = default_wilson()
staggered = default_staggered()
measurement = default_measurement()




struct Params
    L::Tuple  # = (2,2,2,2) # Mandatory
    β::Float64 # = 6 # Mandatory
    NC::Int64 # = 3 
    Dirac_operator::Union{Nothing,String} # = "Wilson"
    quench::Bool # = false
    initial::String # = "cold"
    BoundaryCondition::Array{Int64,1} # =[1,1,1,-1]
    Nwing::Int64  # = 1 
    randomseed::Int64 # = 111

    #For actions
    use_autogeneratedstaples::Bool # = false
    couplinglist::Array{String,1} # = []
    couplingcoeff::Array{Float64,1} # = []
    coupling_loops::Union{Nothing,Array{Array{Array{Tuple{Int,Int},1},1},1}}

    #For MD
    Δτ::Float64 # = 0.05
    MDsteps::Int64 # = 20
    SextonWeingargten::Bool # = false
    Nsteps::Int64 # = 100
    Nthermalization::Int64 # = 10

    #For Sexton-Weingargten
    N_SextonWeingargten::Int64 # = 10 # Mandatory

    #For IntegratedHMC
    #IntegratedHMC::Bool # = false

    #For Heatbath
    #Heatbath::Bool # = false

    #For CG
    eps::Float64 # = 1e-19
    MaxCGstep::Int64 # = 3000

    #For Wilson

    hop::Float64 # = 0.141139#Hopping parameter # Mandatory
    r::Float64 # = 1 #Wilson term
    #For WilsonClover
    Clover_coefficient::Float64  # = 1.5612

    #For Staggered
    mass::Float64 # = 0.5 # Mandatory
    Nf::Int64 # = 4

    #verbose
    verboselevel::Int64 # = 2

    #For measurement
    measurement_methods::Vector{Dict} # = defaultmeasures

    saveU_format::Union{String,Nothing}
    saveU_every::Int64
    saveU_dir::String

    loadU_format::Union{String,Nothing}
    loadU_dir::String

    update_method::String
    βeff::Union{Float64,Array{Float64,1}}
    firstlearn::Int64
    log_dir::String
    logfile::String
    load_fp::IOStream
    measuredir::String
    integratedFermionAction::Bool
    training_data_name::String
    useOR::Bool
    numOR::Int64
    initialtrj::Int64
    loadU_fromfile::Bool
    loadU_filename::String

    smearing_for_fermion::String
    stout_numlayers::Union{Nothing,Int64}
    stout_ρ::Union{Nothing,Array{Float64,1}}
    stout_loops::Union{Nothing,Array{String,1}}

    Domainwall_M::Union{Nothing,Float64}
    Domainwall_m::Union{Nothing,Float64}
    Domainwall_L5::Union{Nothing,Int64}
    #Domainwall_b::Union{Nothing,Float64}
    ##Domainwall_c::Union{Nothing,Float64}
    #Domainwall_ωs::Union{Nothing,Array{Float64,1}}
    #Domainwall_r::Union{Nothing,Float64}

    julian_random_number::Bool
    ITERATION_MAX::Int64
    isevenodd::Bool
    hasgradientflow::Bool
    eps_flow::Float64
    numflow::Int64
    Nflow::Int64
    measurements_for_flow::Vector{Dict} #measurement in gradientflow


end



end

#using .System_parameters
#System_parameters.parameterloading(ARGS[1])
