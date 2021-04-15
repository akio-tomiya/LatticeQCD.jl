# - - parameters - - - - - - - - - - - 
# Physical setting 
system["L"] = (4, 4, 4, 4)
system["β"] = 5.7
system["NC"] = 3
system["Nthermalization"] = 0
system["Nsteps"] = 101
system["initial"] = "cold"
system["initialtrj"] = 1
system["update_method"] = "HMC"
system["Nwing"] = 1
	
# Physical setting(fermions)
system["quench"] = false
system["Dirac_operator"] = "Staggered"
staggered["Nf"] = 2
staggered["mass"] = 0.5
system["BoundaryCondition"] = [1, 1, 1, -1]
	
# System Control
system["log_dir"] = "./logs"
system["logfile"] = "HMC_L04040404_beta5.7_Staggered_mass0.5_Nf2.txt"
system["saveU_dir"] = ""
system["saveU_format"] = nothing
system["verboselevel"] = 2
system["randomseed"] = 111
measurement["measurement_basedir"] = "./measurements"
measurement["measurement_dir"] = "HMC_L04040404_beta5.7_Staggered_mass0.5_Nf2"
	
# HMC related
md["Δτ"] = 0.05
md["SextonWeingargten"] = false
md["N_SextonWeingargten"] = 2
md["MDsteps"] = 20
cg["eps"] = 1.0e-19
cg["MaxCGstep"] = 3000
	
# Action parameter for SLMC
actions["use_autogeneratedstaples"] = false
actions["couplingcoeff"] = Any[]
actions["couplinglist"] = Any[]
	
# Measurement set
measurement["measurement_methods"] = Dict[ 
  Dict{Any,Any}("methodname" => "Plaquette",
    "measure_every" => 1
  )
]
	
# - - - - - - - - - - - - - - - - - - -
