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




mutable struct Params
    const L::Tuple  # = (2,2,2,2) # Mandatory
    const β::Float64 # = 6 # Mandatory
    const NC::Int64 # = 3 
    const Dirac_operator::Union{Nothing,String} # = "Wilson"
    quench::Bool # = false
    const initial::String # = "cold"
    const BoundaryCondition::Array{Int64,1} # =[1,1,1,-1]
    const Nwing::Int64  # = 1 
    const randomseed::Int64 # = 111

    #For actions
    const use_autogeneratedstaples::Bool # = false
    const couplinglist::Array{String,1} # = []
    const couplingcoeff::Array{Float64,1} # = []
    const coupling_loops::Union{Nothing,Array{Array{Array{Tuple{Int,Int},1},1},1}}

    #For MD
    const Δτ::Float64 # = 0.05
    const MDsteps::Int64 # = 20
    const SextonWeingargten::Bool # = false
    const Nsteps::Int64 # = 100
    const Nthermalization::Int64 # = 10

    #For Sexton-Weingargten
    const N_SextonWeingargten::Int64 # = 10 # Mandatory

    #For IntegratedHMC
    #IntegratedHMC::Bool # = false

    #For Heatbath
    #Heatbath::Bool # = false

    #For CG
    const eps::Float64 # = 1e-19
    const MaxCGstep::Int64 # = 3000

    #For Wilson

    const hop::Float64 # = 0.141139#Hopping parameter # Mandatory
    const r::Float64 # = 1 #Wilson term
    #For WilsonClover
    const Clover_coefficient::Float64  # = 1.5612

    #For Staggered
    const mass::Float64 # = 0.5 # Mandatory
    const Nf::Int64 # = 4

    #verbose
    const verboselevel::Int64 # = 2

    #For measurement
    const measurement_methods::Vector{Dict} # = defaultmeasures

    const saveU_format::Union{String,Nothing}
    const saveU_every::Int64
    const saveU_dir::String

    const loadU_format::Union{String,Nothing}
    const loadU_dir::String

    const update_method::String
    const βeff::Union{Float64,Array{Float64,1}}
    const firstlearn::Int64
    const log_dir::String
    const logfile::String
    const load_fp::IOStream
    const measuredir::String
    const integratedFermionAction::Bool
    const training_data_name::String
    const useOR::Bool
    const numOR::Int64
    const initialtrj::Int64
    const loadU_fromfile::Bool
    const loadU_filename::String

    const smearing_for_fermion::String
    const stout_numlayers::Union{Nothing,Int64}
    const stout_ρ::Union{Nothing,Array{Float64,1}}
    const stout_loops::Union{Nothing,Array{String,1}}

    const Domainwall_M::Union{Nothing,Float64}
    const Domainwall_m::Union{Nothing,Float64}
    const Domainwall_L5::Union{Nothing,Int64}
    #Domainwall_b::Union{Nothing,Float64}
    ##Domainwall_c::Union{Nothing,Float64}
    #Domainwall_ωs::Union{Nothing,Array{Float64,1}}
    #Domainwall_r::Union{Nothing,Float64}

    const julian_random_number::Bool
    const ITERATION_MAX::Int64
    const isevenodd::Bool
    const hasgradientflow::Bool
    const eps_flow::Float64
    const numflow::Int64
    const Nflow::Int64
    const measurements_for_flow::Vector{Dict} #measurement in gradientflow
    const QPQ::Bool

end



end

#using .System_parameters
#System_parameters.parameterloading(ARGS[1])
