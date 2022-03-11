mutable struct Plaquette_measurement{Dim,TG} <: AbstractMeasurement
    filename::String
    _temporary_gaugefields::Vector{TG}
    Dim::Int8
    factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool

    function Plaquette_measurement(U::Vector{T};
            filename = nothing,
            verbose_level = 2,
            printvalues = true) where T
        myrank = get_myrank(U)
        #=
        if U[1].mpi == false
            myrank = 0
        else
            myrank = U[1].myrank
        end
        =#
        if printvalues
            verbose_print = Verbose_print(verbose_level,myid = myrank,filename=filename)
        else
            verbose_print = nothing
        end
        Dim = length(U)

        if Dim == 4
            comb = 6
        elseif Dim == 2
            comb = 1
        else
            error("Dim = $Dim is not supported in Plaquette_measurement")
        end
        factor = 1/(comb*U[1].NV*U[1].NC)

        numg = 2
        _temporary_gaugefields = Vector{T}(undef,numg)
        _temporary_gaugefields[1] = similar(U[1])
        for i=2:numg
            _temporary_gaugefields[i] = similar(U[1])
        end

        return new{Dim,T}(filename,_temporary_gaugefields,Dim,factor,verbose_print,printvalues)

    end



end

function measure(m::M,itrj,U) where M <: Plaquette_measurement
    temps = get_temporary_gaugefields(m)
    plaq = real(calculate_Plaquette(U,temps[1],temps[2])*m.factor)

    if m.printvalues
        println_verbose_level2(U[1],"-----------------")
        println_verbose_level2(m.verbose_print,"$itrj $plaq # plaq")
        println_verbose_level2(U[1],"-----------------")
    end

    return plaq
end