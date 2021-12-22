module Measure_plaquet_module
    import ..AbstractMeasurement_module:AbstractMeasurement,measure
    import ..Gaugefield:AbstractGaugefields
    import ..Gaugefield:calculate_Plaquette
    import ..Gaugefield:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3,println_verbose2,println_verbose1,
            print_verbose1,print_verbose2,print_verbose3

    mutable struct Measure_plaquet{T} <: AbstractMeasurement
        filename::String
        fp::IOStream
        tempU::Array{T,1}
        factor::Float64
        printvalues::Bool

        function Measure_plaquet(filename,U::Array{T,1},printvalues = true) where T
            fp = open(filename,"w")
            
            tempU = Array{T,1}(undef,2)
            for i=1:2
                tempU[i] = similar(U[1])
            end
            Dim = length(U)

            if Dim == 4
                comb = 6 #4*3/2
            elseif Dim == 3
                comb = 3
            elseif Dim == 2
                comb = 1
            else
                error("dimension $Dim is not supported")
            end
            factor = 1/(comb*U[1].NV*U[1].NC)


            m = new{T}(filename,fp,tempU,factor,printvalues)
            finalizer(m) do m
                close(m.fp)
            end
            return m
        end

    end

    function measure(measurement::M,itrj,U::Array{T,1};verbose = Verbose_2()) where {M <: Measure_plaquet,T <: AbstractGaugefields}
        temp1 = measurement.tempU[1]
        temp2 = measurement.tempU[2]
        
        printvalues = measurement.printvalues
        fp = measurement.fp

        plaq = real(calculate_Plaquette(U,temp1,temp2))*measurement.factor
        

        if printvalues
            println_verbose1(verbose,"-----------------")
            println_verbose1(verbose,"$itrj $plaq # plaq")
            println(fp,"$itrj $plaq # plaq")
            println_verbose1(verbose,"-----------------")
        end

        return plaq
    end
end