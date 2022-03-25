mutable struct Pion_correlator_measurement{Dim,TG} <: AbstractMeasurement
    filename::String
    _temporary_gaugefields::Vector{TG}
    Dim::Int8
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool

    function Pion_correlator_measurement(U::Vector{T};
        filename = nothing,
        verbose_level = 2,
        printvalues = true) where T

        return 
    end

end