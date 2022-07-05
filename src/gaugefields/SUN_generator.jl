module SUN_generator
using LinearAlgebra

struct Generator
    NC::Int64
    generator::Array{Array{ComplexF64,2},1}
    function Generator(NC)
        return new(NC, make_generators(NC))
    end
end

function Base.length(a::Generator)
    return length(a.generator)
end

function lie2matrix!(matrix, g::Generator, a)
    matrix .= 0
    NC = g.NC
    for (i, genmatrix) in enumerate(g.generator)
        matrix .+= a[i] * genmatrix
    end
    return
end

function matrix2lie!(a, g::Generator, A)
    for i = 1:length(a)
        #println("i = $i")
        #display(g.generator[i]*g.generator[i])
        #println("\t")
        #println(tr(g.generator[i]*A)/2)
        a[i] = tr(g.generator[i] * A) / 2
    end
    return
end

sigma = []
s = ComplexF64[
    0 1
    1 0
]
push!(sigma, s)
s = ComplexF64[
    0 -im
    im 0
]
push!(sigma, s)
s = ComplexF64[
    1 0
    0 -1
]
push!(sigma, s)

function make_largematrix(σ, i, j, Nc)
    λ = zeros(ComplexF64, Nc, Nc)
    λ[i, i] = σ[1, 1]
    λ[i, j] = σ[1, 2]
    λ[j, i] = σ[2, 1]
    λ[j, j] = σ[2, 2]
    return λ
end

function normalization(A)
    Norm = sqrt(real(tr(A * A)) / 2)
    #print("    Normalization = √",real(tr(A*A))/2,"= ")
    #println("$Norm")
    return Norm
end

function make_generators(Nc)
    d = ["x", "y", "z"]
    #
    lams = Array{Complex,2}[]
    #
    if Nc < 2
        error("Invalid Nc=$Nc")
    end
    #
    # off-diagonal part
    #println("Off diagonals")
    for i = 1:Nc
        for j = i+1:Nc
            for a = 1:2
                #println("$(d[a]): i j = $i $j")
                A = make_largematrix(sigma[a], i, j, Nc)
                #showmat(A)
                N = normalization(A)
                A = A / N
                push!(lams, A)
                #println(" ")
            end
        end
    end
    #
    # diagonal part
    #println("diagonals")
    for a = 1:Nc-1
        A = zeros(ComplexF64, Nc, Nc)
        for i = 1:a
            A[i, i] = 1
        end
        A[a+1, a+1] = -tr(A)
        #showmat(A)
        N = normalization(A)
        A = A / N
        push!(lams, A)
        #println(" ")
    end
    return lams
end
end
