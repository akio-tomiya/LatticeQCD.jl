include("./LatticeQCD.jl")
using Plots
using .LatticeQCD
using Random
using Dates
#using JLD
#using GFlops

if length(ARGS) == 0
    error("""
    Use input file: 
    Like, 
    julia analyze.jl parameters.jl
    """)
end

function runtest()
    plot_plaquette(ARGS[1])
    plot_polyakov(ARGS[1])
    plot_plaq_and_poly(ARGS[1])
end



@time runtest()
