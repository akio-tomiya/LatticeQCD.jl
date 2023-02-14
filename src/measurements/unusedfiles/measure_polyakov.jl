mutable struct Polyakov_measurement{Dim,TG} <: AbstractMeasurement
    filename::Union{Nothing,String}
    _temporary_gaugefields::Vector{TG}
    Dim::Int8
    #factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool


    function Polyakov_measurement(
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

function Polyakov_measurement(U::Vector{T}, params::Poly_parameters, filename) where {T}
    return Polyakov_measurement(
        U,
        filename = filename,
        verbose_level = params.verbose_level,
        printvalues = params.printvalues,
    )
end




function measure(m::M, itrj, U; additional_string = "") where {M<:Polyakov_measurement}
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
