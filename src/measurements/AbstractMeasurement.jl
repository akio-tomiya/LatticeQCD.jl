module AbstractMeasurement_module
using Wilsonloop
using Gaugefields
import Gaugefields.Verbose_print
import Gaugefields.println_verbose_level2
import Gaugefields.AbstractGaugefields_module.get_myrank
import Gaugefields:AbstractGaugefields
using LinearAlgebra
import Wilsonloop.make_cloverloops

abstract type AbstractMeasurement end

function get_temporary_gaugefields(m::AbstractMeasurement)
    return m._temporary_gaugefields
end

include("measure_plaquette.jl")
include("measure_polyakov.jl")
include("measure_topological_charge.jl")
include("measure_energy_density.jl")
include("Measurement_set.jl")

function measure(measurement::M,itrj,U) where M <: AbstractMeasurement
    error("measure with a type $M is not supported")
end



end