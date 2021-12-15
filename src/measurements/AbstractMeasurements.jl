module AbstractMeasurement_module
    import ..Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3,println_verbose2,println_verbose1,
            print_verbose1,print_verbose2,print_verbose3
    import ..AbstractGaugefields_module:AbstractGaugefields

    abstract type AbstractMeasurement end

    function measure(measurement::M,itrj,U::Array{T,1};verbose = Verbose_2()) where {M <: AbstractMeasurement,T <: AbstractGaugefields}
        error("measure is not implemented in $(typeof(measurement))")
    end
end