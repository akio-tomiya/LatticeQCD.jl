module Measure_polyakov_loop_module
    import ..AbstractMeasurement_module:AbstractMeasurement,measure
    import ..Gaugefield:AbstractGaugefields
    import ..Gaugefields:calculate_Plaquette,calculate_Polyakov_loop
    import ..Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3,println_verbose2,println_verbose1,
            print_verbose1,print_verbose2,print_verbose3

    mutable struct Measure_polyakov{T} <: AbstractMeasurement
        filename::String
        fp::IOStream
        tempU::Array{T,1}
        printvalues::Bool

        function Measure_polyakov(filename,U::Array{T,1},printvalues = true) where T
            fp = open(filename,"w")
            
            tempU = Array{T,1}(undef,2)
            for i=1:2
                tempU[i] = similar(U[1])
            end

            m = new{T}(filename,fp,tempU,printvalues)
            finalizer(m) do m
                close(m.fp)
            end
            return m
        end

    end

    function measure(measurement::M,itrj,U::Array{T,1};verbose = Verbose_2()) where {M <: Measure_polyakov,T <: AbstractGaugefields}
        #error("measure is not implemented in $(typeof(measurement))")

        
        temp1 = measurement.tempU[1]
        temp2 = measurement.tempU[2]

        printvalues = measurement.printvalues
        fp = measurement.fp

        poly = calculate_Polyakov_loop(U,temp1,temp2) 
        

        if printvalues
            println_verbose1(verbose,"-----------------")
            println_verbose1(verbose,"$itrj $(real(poly)) $(imag(poly)) # poly")
            println(fp,"$itrj $(real(poly)) $(imag(poly)) # poly")
            println_verbose1(verbose,"-----------------")
        end

        return poly
    end
end