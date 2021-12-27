using Test
include("../Wilsonloops.jl")
using .Wilsonloops
function test()
    loop = [(1,+1),(2,+1),(1,-1),(2,-1)]
    println(loop)
    w = Gaugeline(loop)
    println("P: ")
    display(w)
    println("P^+: ")
    display(w')
    println("staple")
    for μ=1:4
        make_staple(w,μ)
        make_staple(w',μ)
    end

    println("CUdag")
    for μ=1:4
        CUdag = make_staple_and_loop(w,μ)
        display(CUdag)
    end
    #=
    w = Wilsonline(loop)
    println("P: ")
    display(w)
    println("P^+: ")
    display(w')
    println("staple")
    make_staples(w)
    =#
end
test()