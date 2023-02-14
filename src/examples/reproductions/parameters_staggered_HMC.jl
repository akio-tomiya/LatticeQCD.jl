# https://inspirehep.net/literature/283285
# Lattice setup
system["L"] = (12, 12, 12, 8)
#system["β"] = 5.4;system["initial"] = "cold";md["MDsteps"] = 35; md["Δτ"] = 1/md["MDsteps"] 
system["β"] = 5.3;
system["initial"] = "cold";
md["MDsteps"] = 35;
md["Δτ"] = 1 / md["MDsteps"];
#system["β"] = 5.175;initial = "cold";MDsteps = 40; Δτ = 1/MDsteps
#β = 5.1;initial = "hot";MDsteps = 80; Δτ = 0.5/MDsteps
system["Dirac_operator"] = "Staggered"
system["Nf"] = 4
system["NC"] = 3

staggered["mass"] = 0.025

# HMC
system["quench"] = false
md["SextonWeingargten"] = false
#N_SextonWeingargten = 25
md["Nsteps"] = 5000

system["update_method"] = "HMC"
system["saveU_format"] = "JLD"
system["saveU_every"] = 1
system["saveU_dir"] = "./beta53_confs"



# Measurements
measurement["measurement_methods"] = Array{Dict,1}(undef, 3)
for i = 1:length(measurement["measurement_methods"])
    measurement["measurement_methods"][i] = Dict()
end
measurement["measurement_methods"][1]["methodname"] = "Plaquette"
measurement["measurement_methods"][1]["measure_every"] = 1
measurement["measurement_methods"][1]["fermiontype"] = nothing
#
measurement["measurement_methods"][2]["methodname"] = "Polyakov_loop"
measurement["measurement_methods"][2]["measure_every"] = 1
measurement["measurement_methods"][2]["fermiontype"] = nothing
#
measurement["measurement_methods"][3]["methodname"] = "Chiral_condensate"
measurement["measurement_methods"][3]["fermiontype"] = "Staggered"
measurement["measurement_methods"][3]["measure_every"] = staggered["mass"]
measurement["measurement_methods"][3]["mass"] = 1
measurement["measurement_methods"][3]["Nf"] = 4
#mass_measurement = mass
#Nr = 10
#
#=
Title:
    The Finite Temperature Phase Transition in Four Flavor {QCD} on an 8 X 12^{3} Lattice
Authors
    MT(c) Collaboration:
    R.V. Gavai(Tata Inst. and Bielefeld U. and CERN and Kaiserslautern U. and Illinois U., Urbana), 
    Sourendu Gupta(Tata Inst. and Bielefeld U. and CERN and Kaiserslautern U. and Illinois U., Urbana), 
    A. Irback(Tata Inst. and Bielefeld U. and CERN and Kaiserslautern U. and Illinois U., Urbana), 
    F. Karsch(Tata Inst. and Bielefeld U. and CERN and Kaiserslautern U. and Illinois U., Urbana), 
    S. Meyer(Tata Inst. and Bielefeld U. and CERN and Kaiserslautern U. and Illinois U., Urbana), 
    B. Petersson(Tata Inst. and Bielefeld U. and CERN and Kaiserslautern U. and Illinois U., Urbana), 
    H. Satz(Tata Inst. and Bielefeld U. and CERN and Kaiserslautern U. and Illinois U., Urbana), 
    H.W. Wyld(Tata Inst. and Bielefeld U. and CERN and Kaiserslautern U. and Illinois U., Urbana)
Sep, 1989
Abstract: (Elsevier)
    We present results of a numerical study of lattice QCD with four dynamical flavours of staggered fermions, 
    performed by using a hybrid Monte Carlo algorithm on an 8×12 3 lattice. We find a rapid change in 
    the average value of the Polyakov loop at β c =5.25±0.025 for a quark mass ma =0.025; at this mass value, 
    the behaviour of the chiral order parameter, 〈 Ψ Ψ〉 does not yet allow an independent determination of 
    the transition point. Using existing hadron mass calculations, the value of β c we have obtained here would 
    lead to a transition temperature T ∼100 MeV .
DOI: 
    10.1016/0370-2693(89)90447-4
=#
# Results Fig 2b (By plot digitizer)
#=     β         PbP
    5.100, 0.662943264422041
    5.175, 0.3819303478073508
    5.200, 0.2874088522470827
    5.250, 0.24762703983851386
    5.300, 0.2109864936891166
    5.400, 0.1686361007473265
    5.600, 0.13783581497147912
=#
#= Plot Digitizer json for fig 2 b
{"version":[4,2],"axesColl":[{"name":"XY","type":"XYAxes","isLogX":false,"isLogY":false,"calibrationPoints":[{"px":122.93685756240822,"py":475.712187958884,"dx":"5.0","dy":"0.0","dz":null},{"px":391.7180616740088,"py":477.2393538913363,"dx":"5.7","dy":"0.0","dz":null},{"px":77.12187958883995,"py":476.4757709251101,"dx":"5.0","dy":"0.0","dz":null},{"px":80.93979441997062,"py":93.92070484581498,"dx":"5.7","dy":"1.0","dz":null}]}],"datasetColl":[{"name":"Default Dataset","axesName":"XY","metadataKeys":[],"colorRGB":[200,0,0,255],"data":[{"x":353.92070484581495,"y":425.31571218795887,"value":[5.600218584502331,0.13783581497147912]},{"x":277.56240822320115,"y":413.09838472834065,"value":[5.401048695297173,0.1686361007473265]},{"x":239.76505139500736,"y":396.6813509544787,"value":[5.302190097640081,0.2109864936891166]},{"x":230.22026431718064,"y":389.8091042584435,"value":[5.277154942674726,0.2288077930619974]},{"x":220.67547723935388,"y":382.55506607929516,"value":[5.252109864936891,0.24762703983851386]},{"x":211.51248164464025,"y":367.2834067547724,"value":[5.22785067077942,0.2874088522470827]},{"x":192.42290748898682,"y":331.0132158590308,"value":[5.177194917272429,0.3819303478073508]},{"x":164.9339207048458,"y":223.34801762114537,"value":[5.102809845658362,0.662943264422041]}],"autoDetectionData":{"fgColor":[0,0,255],"bgColor":[255,255,255],"mask":[],"colorDetectionMode":"fg","colorDistance":120,"algorithm":{"algoType":"AveragingWindowAlgo","xStep":10,"yStep":10},"name":0,"imageWidth":442,"imageHeight":520}}],"measurementColl":[]}
=#
