"""
    AlgRemez.jl is a wrapper for AlgRemez written in c++. 
    Please see
    https://github.com/maddyscientist/AlgRemez
"""
module AlgRemez
    using AlgRemez_jll
# Write your package code here.
    const exe = algremez(x -> x)

    struct AlgRemez_coeffs
        α0::Float64
        α::Array{Float64,1}
        β::Array{Float64,1}
        n::Int64
    end

    function Base.display(x::AlgRemez_coeffs)
        println("""
        f(x) = α0 + sum_i^n α[i]/(x + β[i])""")
        println("The order: $(x.n)")
        println("α0: $(x.α0)")
        println("α: $(x.α)")
        println("β: $(x.β)")
    end

    function fittedfunction(coeff::AlgRemez_coeffs)
        function func(x)
            value = coeff.α0
            for i=1:coeff.n
                value += coeff.α[i]/(x + coeff.β[i])
            end
            return value
        end
        x -> func(x)
    end

    """
Output: coeff_plus::AlgRemez_coeffs,coeff_minus::AlgRemez_coeffs
struct AlgRemez_coeffs
    α0::Float64
    α::Array{Float64,1}
    β::Array{Float64,1}
    n::Int64
end
f(x) = x^(y/z) = coeff_plus.α0 + sum_i^n coeff_plus.α[i]/(x + coeff_plus.β[i])
f(x) = x^(-y/z) = coeff_minus.α0 + sum_i^n coeff_minus.α[i]/(x + coeff_minus.β[i])
The function to be approximated f(x) = x^(y/z) and f(x) = x^(-y/z), 
with degree (n,n) over the spectral range [lambda_low,lambda_high],
using precision digits of precision in the arithmetic. 
The parameters y and z must be positive, the approximation to f(x) = x^(-y/z) is simply
the inverse of the approximation to f(x) = x^(y/z).
The default value of precision is 42. 
    """
    function calc_coefficients(y,z,n,lambda_low,lambda_high;precision=42)
        run(`$exe $y $z $n $n $lambda_low $lambda_high $precision`)
        datas = readlines("approx.dat")
        icount = 3
        αplus0 = parse(Float64,split(datas[icount],"=")[2])
        #println(αplus0)
        αplus = zeros(Float64,n)
        βplus = zeros(Float64,n)
        for i=1:n
            icount += 1
            u = split(datas[icount],",")
            αplus[i] = parse(Float64,split(u[1],"=")[2])
            βplus[i] = parse(Float64,split(u[2],"=")[2])
        end
        coeff_plus = AlgRemez_coeffs(αplus0,αplus,βplus,n)
        #println(αplus)
        #println(βplus)


        icount += 4
        αminus0 = parse(Float64,split(datas[icount],"=")[2])
        #println(αminus0)
        αminus = zeros(Float64,n)
        βminus = zeros(Float64,n)
        for i=1:n
            icount += 1
            u = split(datas[icount],",")
            αminus[i] = parse(Float64,split(u[1],"=")[2])
            βminus[i] = parse(Float64,split(u[2],"=")[2])
        end
        coeff_minus = AlgRemez_coeffs(αminus0,αminus,βminus,n)
        #println(αminus)
        #println(βminus)
        #for data in datas
        #    println(data)
        #end

        rm("approx.dat")
        rm("error.dat")

        return coeff_plus,coeff_minus
    end

end