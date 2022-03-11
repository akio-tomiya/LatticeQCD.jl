module AbstractMeasurement_module
using Gaugefields
import Gaugefields.Verbose_print
import Gaugefields.println_verbose_level2
import Gaugefields.AbstractGaugefields_module.get_myrank

abstract type AbstractMeasurement end

function get_temporary_gaugefields(m::AbstractMeasurement)
    return m._temporary_gaugefields
end

include("measure_plaquette.jl")

function measure(measurement::M,itrj,U) where M <: AbstractMeasurement
    error("measure with a type $M is not supported")
end

end