using LatticeQCD
using Test

@testset "LatticeQCD.jl" begin
    # Write your tests here.
    # Write your package code here.
    #univ = Universe()

    β = 6
    gparam = Setup_Gauge_action(β)
    fparam = Setup_Fermi_action()
    L = (4,4,4,4)

#    L = (12,12,12,12)
    initial = "cold"
    univ = Universe(L,gparam,initial,fparam)

    plaq = calc_plaquette(univ)
    println("plaq = ",plaq)

    mdparams = MD_parameters_standard(gparam)

    for i=1:10
        Sold = md_initialize!(univ)
        Snew = md!(univ,mdparams)

        metropolis_update!(univ,Sold,Snew)
        plaq = calc_plaquette(univ)
        println("-------------------------------------")
        println("$i-th plaq = ",plaq)
        println("-------------------------------------")
    end
end
