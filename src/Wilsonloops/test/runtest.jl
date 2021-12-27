using Test
include("../src/Wilsonloops.jl")
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
        println("μ = $μ")
        V1 = make_staple(w,μ)
        V2 = make_staple(w',μ)
        display(V1)
        display(V2)
    end


    println("derive w")
    for μ=1:4
        dU = derive_U(w,μ)
        for i=1:length(dU)
            display(dU[i])
        end
    end

    println("-------------------------------------------------------")
    println("C and dC/dU")
    for μ=1:4
        C = make_Cμ(w,μ)
        #=
        V1 = make_staple(w,μ)
        V2 = make_staple(w',μ)
        C = eltype(V1)[]
        for i=1:length(V1)
            push!(C,V1[i]')
        end
        for i=1:length(V2)
            push!(C,V2[i]')
        end
        =#
        println("-------------------------------------------")
        println("μ = $μ")
        for i=1:length(C)
            println("---------------------------------------")
            println("C[$i]: ")
            display(C[i])
            for ν=1:4
                println("-----------------------------")
                println("ν = $ν")
                dCdU = derive_U(C[i],ν)
                println("dC_{$μ}/dU_{$ν}: ")
                for j=1:length(dCdU)
                    display(dCdU[j])
                end
            end
        end
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