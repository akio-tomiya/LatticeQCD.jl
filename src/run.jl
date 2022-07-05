include("./LatticeQCD.jl")
using .LatticeQCD
using Random
using Dates
#using JLD


#using GFlops

if length(ARGS) == 0
    error("""
    Use input file: 
    Like, 
    julia run.jl parameters.jl
    """)
end

function runtest()
    run_LQCD(ARGS[1])
    #parameters = parameterloading()
    #run_LQCD(parameters)
end



@time runtest()
