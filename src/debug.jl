include("./LatticeQCD.jl")
using .LatticeQCD
using Random
using Dates
#using JLD

function test()
    Random.seed!(111)
    A = rand(ComplexF64, 4, 4) * 4
    A = A' * A
    n = 4
    ϕ = ComplexF64[1, 2, 3, 4]
    ϕr = LatticeQCD.calc_exactvalue(n, A, ϕ)
    println(ϕ' * ϕr)

    ϕr2 = LatticeQCD.calc_Anϕ(n, A, ϕ)
    println(ϕ' * ϕr2)

    #LatticeQCD.calc_det(n,A,ϕ)
end

test()
