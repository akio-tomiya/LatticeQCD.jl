mutable struct Wilson_loop_measurement{Dim,TG} <: AbstractMeasurement
    filename::Union{Nothing,String}
    _temporary_gaugefields::Vector{TG}
    Dim::Int8
    #factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool


    function Wilson_loop_measurement(
        U::Vector{T};
        filename = nothing,
        verbose_level = 2,
        printvalues = true,
    ) where {T}
        myrank = get_myrank(U)
        #=
        if U[1].mpi == false
            myrank = 0
        else
            myrank = U[1].myrank
        end
        =#
        if printvalues
            verbose_print = Verbose_print(verbose_level, myid = myrank, filename = filename)
        else
            verbose_print = nothing
        end
        Dim = length(U)


        numg = 2
        _temporary_gaugefields = Vector{T}(undef, numg)
        _temporary_gaugefields[1] = similar(U[1])
        for i = 2:numg
            _temporary_gaugefields[i] = similar(U[1])
        end

        return new{Dim,T}(filename, _temporary_gaugefields, Dim, verbose_print, printvalues)

    end
end

function Wilson_loop_measurement(U::Vector{T}, params::Poly_parameters, filename) where {T}
    return Wilson_loop_measurement(
        U,
        filename = filename,
        verbose_level = params.verbose_level,
        printvalues = params.printvalues,
    )
end




function measure(m::M, itrj, U; additional_string = "") where {M<:Wilson_loop_measurement}
    temps = get_temporary_gaugefields(m)
    poly = calculate_Polyakov_loop(U, temps[1], temps[2])
    measurestring = ""

    if m.printvalues
        #println_verbose_level2(U[1],"-----------------")
        measurestring = "$itrj $additional_string $(real(poly)) $(imag(poly)) # poly"
        println_verbose_level2(m.verbose_print, measurestring)
        #println_verbose_level2(U[1],"-----------------")
    end

    return poly, measurestring
end

function calc_Wilson_loop(U::Array{T,1}, Lt, Ls) where {T<:GaugeFields}
    # Making a ( Ls × Lt) Wilson loop operator for potential calculations
    WL = 0.0 + 0.0im
    NV = U[1].NV
    NC = U[1].NC
    Wmat = Array{GaugeFields_1d,2}(undef, 4, 4)
    #
    calc_large_wiloson_loop!(Wmat, Lt, Ls, U) # make wilon loop operator and evaluate as a field, not traced.
    WL = calc_Wilson_loop_core(Wmat, U, NV) # tracing over color and average over spacetime and x,y,z.
    NDir = 3.0 # in 4 diemension, 3 associated staples. t-x plane, t-y plane, t-z plane
    return real(WL) / NV / NDir / NC
end
function calc_Wilson_loop_core(Wmat, U::Array{GaugeFields{S},1}, NV) where {S<:SUn}
    if S == SU3
        NC = 3
    elseif S == SU2
        NC = 2
    else
        NC = U[1].NC
        #error("NC != 2,3 is not supported")
    end
    W = 0.0 + 0.0im
    for n = 1:NV
        for μ = 1:3 # spatial directions
            ν = 4  # T-direction is not summed over
            W += tr(Wmat[μ, ν][:, :, n])
        end
    end
    return W
end
function calc_large_wiloson_loop!(Wmat, Lt, Ls, U)
    W_operator = make_Wilson_loop(Lt, Ls)
    calc_large_wiloson_loop!(Wmat, W_operator, U)
    return
end
function make_Wilson_loop(Lt, Ls, Dim)
    #= Making a Wilson loop operator for potential calculations
        Ls × Lt
        ν=4
        ↑
        +--+ 
        |  |
        |  |
        |  |
        +--+ → μ=1,2,3 (averaged)
    =#
    Wmatset = Array{Wilsonline{Dim},2}(undef, 4, 4)
    for μ = 1:3 # spatial directions
        ν = 4 # T-direction is not summed over
        loops = Wilsonline{Dim}[]
        loop = Wilsonline([(μ, Ls), (ν, Lt), (μ, -Ls), (ν, -Lt)])
        push!(loops, loop)
        Wmatset[μ, ν] = loops
    end
    return Wmatset
end
function calc_large_wiloson_loop!(temp_Wmat, loops_μν, U)
    W = temp_Wmat
    for μ = 1:3 # spatial directions
        ν = 4 # T-direction is not summed over
        loopset = Loops(U, loops_μν[μ, ν])
        W[μ, ν] = evaluate_loops(loopset, U)
    end
    return
end
