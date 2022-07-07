module AbstractMeasurement_module
using Wilsonloop
using Gaugefields
import Gaugefields.Verbose_print
import Gaugefields.println_verbose_level2
import Gaugefields.AbstractGaugefields_module.get_myrank
import Gaugefields: AbstractGaugefields
import Gaugefields.Abstractsmearing
using LinearAlgebra
using LatticeDiracOperators
import LatticeDiracOperators.Dirac_operators:
    clear_fermion!, AbstractFermionfields_4D, Z4_distribution_fermi!
import Wilsonloop.make_cloverloops
import ..Parameter_structs:Plaq_parameters,Poly_parameters,
TopologicalCharge_parameters,
ChiralCondensate_parameters,
Pion_parameters,
Energy_density_parameters

abstract type AbstractMeasurement end

function get_temporary_gaugefields(m::AbstractMeasurement)
    return m._temporary_gaugefields
end

function get_temporary_fermionfields(m::AbstractMeasurement)
    return m._temporary_fermionfields
end

include("measure_plaquette.jl")
include("measure_polyakov.jl")
include("measure_topological_charge.jl")
include("measure_energy_density.jl")
include("measure_chiral_condensate.jl")
include("measure_Pion_correlator.jl")

include("measurement_parameters_set.jl")
include("Measurement_set.jl")


function measure(measurement::M, itrj, U) where {M<:AbstractMeasurement}
    error("measure with a type $M is not supported")
end



end
