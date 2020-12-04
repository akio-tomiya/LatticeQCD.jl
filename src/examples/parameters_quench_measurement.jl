system["L"] = (8,8,8,8)
system["β"] = 5.7
system["NC"] = 3
system["quench"] = true

system["SextonWeingargten"] = false
system["N_SextonWeingargten"] = 10


md["MDsteps"] = 20
md["Δτ"] = 1/md["MDsteps"]

system["Dirac_operator"] = nothing


md["Nthermalization"] = 10

md["Nsteps"] = 100+md["Nthermalization"]

system["update_method"] = "Fileloading"
system["loadU_format"] = "JLD"
system["loadU_dir"] = "confs"

measurement["measurement_methods"] = Array{Dict,1}(undef,2)
for i=1:length(measurement["measurement_methods"])
    measurement["measurement_methods"][i] = Dict()
end
measurement["measurement_methods"][1]["methodname"] = "Plaquette"
measurement["measurement_methods"][1]["measure_every"] = 1
measurement["measurement_methods"][1]["fermiontype"] = nothing
measurement["measurement_methods"][2]["methodname"] = "Polyakov_loop"
measurement["measurement_methods"][2]["measure_every"] = 1
measurement["measurement_methods"][2]["fermiontype"] = nothing
#=
measurement["measurement_methods"][3]["methodname"] = "Topological_charge"
measurement["measurement_methods"][3]["measure_every"] = 10
measurement["measurement_methods"][3]["fermiontype"] = nothing
measurement["measurement_methods"][3]["numflow"]  = 10
measurement["measurement_methods"][4]["methodname"] = "Chiral_condensate" 
measurement["measurement_methods"][4]["measure_every"] = 20
measurement["measurement_methods"][4]["fermiontype"] = "Staggered"
measurement["measurement_methods"][4]["Nf"] = 4
measurement["measurement_methods"][5]["methodname"] = "Pion_correlator" 
measurement["measurement_methods"][5]["measure_every"] = 20
measurement["measurement_methods"][5]["fermiontype"] = "Wilson"
=#