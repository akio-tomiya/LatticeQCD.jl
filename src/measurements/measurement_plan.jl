module MeasurementPlan_module

import QCDMeasurements
import Gaugefields

import ..LQCDCommunication: is_root

import ..System_parameters: Params
import ..LQCDConfig_module: show_config
import ..Simulation_module:
    AbstractConfiguration,
    GaugeConfiguration,
    Simulation,
    copy_configuration!,
    update!

"""Common supertype for serializable observable settings."""
abstract type AbstractObservableConfig end

"""Measure the normalized plaquette."""
struct PlaquetteObservableConfig <: AbstractObservableConfig end

"""Measure the Polyakov loop in the final lattice direction."""
struct PolyakovLoopObservableConfig <: AbstractObservableConfig end

"""Measure one or more gluonic definitions of topological charge."""
struct TopologicalChargeObservableConfig{K<:Tuple} <:
       AbstractObservableConfig
    kinds::K
    improved_definition::Symbol

    function TopologicalChargeObservableConfig(
        kinds::K,
        improved_definition::Union{Symbol,AbstractString}=:alexandrou,
    ) where {K<:Tuple}
        isempty(kinds) && throw(ArgumentError(
            "at least one topological-charge definition is required",
        ))
        normalized_kinds = map(kinds) do kind
            value = Symbol(lowercase(strip(String(kind))))
            value in (:plaquette, :clover) || throw(ArgumentError(
                "topological-charge kind must be plaquette or clover; " *
                "got $(repr(kind))",
            ))
            value
        end
        definition = Symbol(lowercase(replace(
            strip(String(improved_definition)),
            "-" => "_",
        )))
        definition in (:alexandrou, :bilson_thompson, :both) ||
            throw(ArgumentError(
                "improved topological-charge definition must be " *
                "alexandrou, bilson_thompson, or both; got " *
                repr(improved_definition),
            ))
        return new{typeof(normalized_kinds)}(
            normalized_kinds,
            definition,
        )
    end
end

TopologicalChargeObservableConfig(
    kinds::AbstractVector,
    improved_definition::Union{Symbol,AbstractString}=:alexandrou,
) = TopologicalChargeObservableConfig(Tuple(kinds), improved_definition)

"""Measure rectangular Wilson loops up to `Tmax` by `Rmax`."""
struct WilsonLoopObservableConfig <: AbstractObservableConfig
    Tmax::Int
    Rmax::Int

    function WilsonLoopObservableConfig(Tmax::Integer, Rmax::Integer)
        Tmax >= 0 || throw(ArgumentError(
            "Tmax must be nonnegative; got $Tmax",
        ))
        Rmax >= 0 || throw(ArgumentError(
            "Rmax must be nonnegative; got $Rmax",
        ))
        return new(Int(Tmax), Int(Rmax))
    end
end

"""Measure the clover energy density."""
struct EnergyDensityObservableConfig <: AbstractObservableConfig end

"""Common supertype for fermion operators used only by measurements."""
abstract type AbstractMeasurementFermionConfig end

"""Wilson or Wilson-clover parameters used by a fermionic observable."""
struct WilsonMeasurementFermionConfig{
    H<:Real,
    R<:Real,
    C<:Real,
} <: AbstractMeasurementFermionConfig
    hopping::H
    r::R
    clover_coefficient::C
    clover::Bool
end

"""Staggered parameters used by a fermionic observable."""
struct StaggeredMeasurementFermionConfig{M<:Real} <:
       AbstractMeasurementFermionConfig
    mass::M
    flavors::Int

    function StaggeredMeasurementFermionConfig(
        mass::M,
        flavors::Integer,
    ) where {M<:Real}
        flavors > 0 || throw(ArgumentError(
            "the number of fermion flavors must be positive; got $flavors",
        ))
        return new{M}(mass, Int(flavors))
    end
end

"""CG settings owned by a fermionic observable."""
struct MeasurementSolverConfig{T<:Real}
    tolerance::T
    maximum_iterations::Int
    method::Symbol

    function MeasurementSolverConfig(
        tolerance::T,
        maximum_iterations::Integer,
        method::Union{Symbol,AbstractString}=:bicg,
    ) where {T<:Real}
        isfinite(tolerance) && tolerance > 0 || throw(ArgumentError(
            "measurement solver tolerance must be positive and finite; " *
            "got $tolerance",
        ))
        maximum_iterations > 0 || throw(ArgumentError(
            "measurement solver iterations must be positive; " *
            "got $maximum_iterations",
        ))
        normalized_method = Symbol(lowercase(strip(String(method))))
        normalized_method in (:bicg, :bicgstab, :preconditiond_bicgstab) ||
            throw(ArgumentError(
                "measurement solver method must be bicg, bicgstab, or " *
                "preconditiond_bicgstab; got $(repr(method))",
            ))
        return new{T}(
            tolerance,
            Int(maximum_iterations),
            normalized_method,
        )
    end
end

"""Common supertype for fermion smearing used by measurements."""
abstract type AbstractMeasurementSmearingConfig end

"""Use the gauge field without smearing."""
struct NoMeasurementSmearingConfig <:
       AbstractMeasurementSmearingConfig end

"""Apply the configured stout layers before a fermionic measurement."""
struct StoutMeasurementSmearingConfig{R<:Tuple,L<:Tuple} <:
       AbstractMeasurementSmearingConfig
    layers::Int
    coefficients::R
    loops::L

    function StoutMeasurementSmearingConfig(
        layers::Integer,
        coefficients::R,
        loops::L,
    ) where {R<:Tuple,L<:Tuple}
        layers > 0 || throw(ArgumentError(
            "the number of stout layers must be positive; got $layers",
        ))
        isempty(coefficients) && throw(ArgumentError(
            "stout smearing requires at least one coefficient",
        ))
        isempty(loops) && throw(ArgumentError(
            "stout smearing requires at least one loop",
        ))
        all(coefficient -> coefficient isa Real, coefficients) ||
            throw(ArgumentError("stout coefficients must be real"))
        all(loop -> loop isa Symbol || loop isa AbstractString, loops) ||
            throw(ArgumentError("stout loops must be strings or symbols"))
        return new{R,L}(Int(layers), coefficients, loops)
    end
end

StoutMeasurementSmearingConfig(
    layers::Integer,
    coefficients::AbstractVector,
    loops::AbstractVector,
) = StoutMeasurementSmearingConfig(
    layers,
    Tuple(coefficients),
    Tuple(loops),
)

"""Stochastic chiral-condensate measurement settings."""
struct ChiralCondensateObservableConfig{
    F<:AbstractMeasurementFermionConfig,
    S<:MeasurementSolverConfig,
    M<:AbstractMeasurementSmearingConfig,
} <: AbstractObservableConfig
    fermion::F
    solver::S
    smearing::M
    noise_vectors::Int

    function ChiralCondensateObservableConfig(
        fermion::F,
        solver::S,
        smearing::M,
        noise_vectors::Integer,
    ) where {
        F<:AbstractMeasurementFermionConfig,
        S<:MeasurementSolverConfig,
        M<:AbstractMeasurementSmearingConfig,
    }
        noise_vectors > 0 || throw(ArgumentError(
            "the number of chiral-condensate noise vectors must be " *
            "positive; got $noise_vectors",
        ))
        return new{F,S,M}(
            fermion,
            solver,
            smearing,
            Int(noise_vectors),
        )
    end
end

"""Point-source pion-correlator measurement settings."""
struct PionCorrelatorObservableConfig{
    F<:AbstractMeasurementFermionConfig,
    S<:MeasurementSolverConfig,
    M<:AbstractMeasurementSmearingConfig,
} <: AbstractObservableConfig
    fermion::F
    solver::S
    smearing::M
end

"""Legacy-compatible periodic scheduling on trajectory numbers."""
struct PeriodicSchedule
    every::Int
    start::Int

    function PeriodicSchedule(every::Integer, start::Integer=0)
        every > 0 || throw(ArgumentError(
            "measurement interval must be positive; got $every",
        ))
        start >= 0 || throw(ArgumentError(
            "measurement start must be nonnegative; got $start",
        ))
        return new(Int(every), Int(start))
    end
end

is_due(schedule::PeriodicSchedule, trajectory::Integer) =
    trajectory >= schedule.start && mod(trajectory, schedule.every) == 0

"""One observable and the independent schedule that selects it."""
struct ScheduledMeasurementConfig{
    O<:AbstractObservableConfig,
    S,
}
    observable::O
    schedule::S
end

"""A heterogeneous, fully concrete tuple of scheduled measurements."""
struct MeasurementPlan{T<:Tuple}
    measurements::T
end

MeasurementPlan(measurements::ScheduledMeasurementConfig...) =
    MeasurementPlan(measurements)

"""Marker used when a measurement program has no gradient-flow stage."""
struct NoGradientFlowMeasurementConfig end

"""
Measurements made along a Wilson gradient-flow trajectory.

`samples` is the number of reported flow points. Each reported point advances
the integrator by `integration_steps` steps of size `step_size`. Schedules in
`measurements` are evaluated against the reported flow-point index, preserving
the meaning of `measure_every` in the legacy input format.
"""
struct GradientFlowMeasurementConfig{
    T<:Real,
    M<:MeasurementPlan,
}
    step_size::T
    samples::Int
    integration_steps::Int
    measurements::M

    function GradientFlowMeasurementConfig(
        step_size::T,
        samples::Integer,
        integration_steps::Integer,
        measurements::M,
    ) where {T<:Real,M<:MeasurementPlan}
        isfinite(step_size) && step_size > 0 || throw(ArgumentError(
            "gradient-flow step size must be positive and finite; " *
            "got $step_size",
        ))
        samples > 0 || throw(ArgumentError(
            "the number of gradient-flow samples must be positive; " *
            "got $samples",
        ))
        integration_steps > 0 || throw(ArgumentError(
            "gradient-flow integration steps must be positive; " *
            "got $integration_steps",
        ))
        isempty(measurements.measurements) && throw(ArgumentError(
            "gradient flow requires at least one measurement",
        ))
        return new{T,M}(
            step_size,
            Int(samples),
            Int(integration_steps),
            measurements,
        )
    end
end

"""Direct and optional gradient-flow measurements for one simulation."""
struct MeasurementProgram{D<:MeasurementPlan,F}
    direct::D
    gradient_flow::F

    function MeasurementProgram(
        direct::D,
        gradient_flow::F,
    ) where {
        D<:MeasurementPlan,
        F<:Union{
            NoGradientFlowMeasurementConfig,
            GradientFlowMeasurementConfig,
        },
    }
        return new{D,F}(direct, gradient_flow)
    end
end

MeasurementProgram(direct::MeasurementPlan) =
    MeasurementProgram(direct, NoGradientFlowMeasurementConfig())

function fermion_config(parameters::QCDMeasurements.ChiralCondensateParameters)
    kind = String(parameters.fermiontype)
    if kind == "Staggered"
        return StaggeredMeasurementFermionConfig(
            parameters.mass,
            parameters.Nf,
        )
    elseif kind == "Wilson"
        return WilsonMeasurementFermionConfig(
            parameters.hop,
            parameters.r,
            parameters.Clover_coefficient,
            false,
        )
    end
    throw(ArgumentError(
        "Chiral_condensate supports Wilson and Staggered fermions; " *
        "got $(repr(kind))",
    ))
end

function fermion_config(parameters::QCDMeasurements.PionCorrelatorParameters)
    kind = String(parameters.fermiontype)
    fermion = parameters.fermion_parameters
    if kind == "Staggered"
        return StaggeredMeasurementFermionConfig(
            fermion.mass,
            fermion.Nf,
        )
    elseif kind in ("Wilson", "WilsonClover")
        return WilsonMeasurementFermionConfig(
            fermion.hop,
            fermion.r,
            fermion.Clover_coefficient,
            kind == "WilsonClover",
        )
    end
    throw(ArgumentError(
        "Pion_correlator supports Wilson, WilsonClover, and Staggered " *
        "fermions; got $(repr(kind))",
    ))
end

function smearing_config(parameters)
    kind = lowercase(strip(String(parameters.smearing_for_fermion)))
    kind == "nothing" && return NoMeasurementSmearingConfig()
    kind == "stout" || throw(ArgumentError(
        "measurement smearing must be nothing or stout; got $(repr(kind))",
    ))
    return StoutMeasurementSmearingConfig(
        parameters.stout_numlayers,
        parameters.stout_ρ,
        parameters.stout_loops,
    )
end

observable_config(::QCDMeasurements.PlaquetteParameters) =
    PlaquetteObservableConfig()

observable_config(::QCDMeasurements.PolyakovParameters) =
    PolyakovLoopObservableConfig()

function observable_config(
    parameters::QCDMeasurements.TopologicalChargeParameters,
)
    return TopologicalChargeObservableConfig(
        parameters.kinds_of_topological_charge,
        parameters.improved_topological_charge_definition,
    )
end

function observable_config(
    parameters::QCDMeasurements.ChiralCondensateParameters,
)
    return ChiralCondensateObservableConfig(
        fermion_config(parameters),
        MeasurementSolverConfig(
            parameters.eps,
            parameters.MaxCGstep,
            :bicg,
        ),
        smearing_config(parameters),
        parameters.Nr,
    )
end

function observable_config(
    parameters::QCDMeasurements.PionCorrelatorParameters,
)
    return PionCorrelatorObservableConfig(
        fermion_config(parameters),
        MeasurementSolverConfig(
            parameters.eps,
            parameters.MaxCGstep,
            parameters.method_CG,
        ),
        smearing_config(parameters),
    )
end

observable_config(parameters::QCDMeasurements.WilsonLoopParameters) =
    WilsonLoopObservableConfig(parameters.Tmax, parameters.Rmax)

observable_config(::QCDMeasurements.EnergyDensityParameters) =
    EnergyDensityObservableConfig()

"""Convert one Wizard/legacy measurement dictionary to typed settings."""
function measurement_config(values::AbstractDict; start::Integer=0)
    parameters =
        QCDMeasurements.construct_Measurement_parameters_from_dict(values)
    observable = observable_config(parameters)
    schedule = PeriodicSchedule(parameters.measure_every, start)
    return ScheduledMeasurementConfig(observable, schedule)
end

"""Convert a list of Wizard/legacy dictionaries to a typed plan."""
function measurement_plan(
    values::AbstractVector{<:AbstractDict};
    start::Integer=0,
)
    measurements = Tuple(
        measurement_config(value; start) for value in values
    )
    return MeasurementPlan(measurements)
end

"""
Convert the direct (unflowed) measurements in legacy `Params` to a typed plan.
Use [`measurement_program`](@ref) when gradient-flow settings must be included.
"""
function legacy_measurement_plan(parameters)
    start = parameters.update_method == "Fileloading" ?
            0 : max(0, parameters.Nthermalization)
    return measurement_plan(parameters.measurement_methods; start)
end


measurement_plan(parameters::Params) = legacy_measurement_plan(parameters)

"""Convert the gradient-flow part of legacy `Params` to typed settings."""
function legacy_gradient_flow_measurement_config(parameters)
    parameters.hasgradientflow || return NoGradientFlowMeasurementConfig()
    measurements = measurement_plan(parameters.measurements_for_flow)
    return GradientFlowMeasurementConfig(
        parameters.eps_flow,
        parameters.numflow,
        parameters.Nflow,
        measurements,
    )
end


gradient_flow_measurement_config(parameters::Params) =
    legacy_gradient_flow_measurement_config(parameters)

"""Convert all direct and gradient-flow measurements in legacy `Params`."""
function legacy_measurement_program(parameters)
    return MeasurementProgram(
        legacy_measurement_plan(parameters),
        legacy_gradient_flow_measurement_config(parameters),
    )
end


measurement_program(parameters::Params) = legacy_measurement_program(parameters)

observable_name(::PlaquetteObservableConfig) = :plaquette
observable_name(::PolyakovLoopObservableConfig) = :polyakov_loop
observable_name(::TopologicalChargeObservableConfig) = :topological_charge
observable_name(::ChiralCondensateObservableConfig) = :chiral_condensate
observable_name(::PionCorrelatorObservableConfig) = :pion_correlator
observable_name(::WilsonLoopObservableConfig) = :wilson_loop
observable_name(::EnergyDensityObservableConfig) = :energy_density

function measurement_dictionary(::PlaquetteObservableConfig)
    return Dict{String,Any}("methodname" => "Plaquette")
end

function measurement_dictionary(::PolyakovLoopObservableConfig)
    return Dict{String,Any}("methodname" => "Polyakov_loop")
end

function measurement_dictionary(config::TopologicalChargeObservableConfig)
    return Dict{String,Any}(
        "methodname" => "Topological_charge",
        "kinds_of_topological_charge" => collect(String.(config.kinds)),
        "improved_topological_charge_definition" =>
            String(config.improved_definition),
    )
end

function measurement_dictionary(config::WilsonLoopObservableConfig)
    return Dict{String,Any}(
        "methodname" => "Wilson_loop",
        "Tmax" => config.Tmax,
        "Rmax" => config.Rmax,
    )
end

function measurement_dictionary(::EnergyDensityObservableConfig)
    return Dict{String,Any}("methodname" => "Energy_density")
end

function add_fermion_parameters!(
    values::Dict{String,Any},
    config::WilsonMeasurementFermionConfig,
)
    values["fermiontype"] = config.clover ? "WilsonClover" : "Wilson"
    values["hop"] = config.hopping
    values["r"] = config.r
    values["Clover_coefficient"] = config.clover_coefficient
    return values
end

function add_fermion_parameters!(
    values::Dict{String,Any},
    config::StaggeredMeasurementFermionConfig,
)
    values["fermiontype"] = "Staggered"
    values["mass"] = config.mass
    values["Nf"] = config.flavors
    return values
end

function add_smearing_parameters!(
    values::Dict{String,Any},
    ::NoMeasurementSmearingConfig,
)
    values["smearing_for_fermion"] = "nothing"
    return values
end

function add_smearing_parameters!(
    values::Dict{String,Any},
    config::StoutMeasurementSmearingConfig,
)
    values["smearing_for_fermion"] = "stout"
    values["stout_numlayers"] = config.layers
    values["stout_ρ"] = collect(config.coefficients)
    values["stout_loops"] = collect(String.(config.loops))
    return values
end

function measurement_dictionary(config::ChiralCondensateObservableConfig)
    values = Dict{String,Any}(
        "methodname" => "Chiral_condensate",
        "eps" => config.solver.tolerance,
        "MaxCGstep" => config.solver.maximum_iterations,
        "Nr" => config.noise_vectors,
    )
    add_fermion_parameters!(values, config.fermion)
    add_smearing_parameters!(values, config.smearing)
    return values
end

function measurement_dictionary(config::PionCorrelatorObservableConfig)
    values = Dict{String,Any}(
        "methodname" => "Pion_correlator",
        "eps" => config.solver.tolerance,
        "MaxCGstep" => config.solver.maximum_iterations,
        "method_CG" => String(config.solver.method),
    )
    add_fermion_parameters!(values, config.fermion)
    add_smearing_parameters!(values, config.smearing)
    return values
end

function build_qcd_measurement(config::AbstractObservableConfig, gauge)
    values = measurement_dictionary(config)
    values["printvalues"] = false
    return QCDMeasurements.prepare_measurement(gauge, values)
end

"""One concrete QCDMeasurements object and its typed settings."""
struct ScheduledMeasurementRuntime{C,M,S}
    config::C
    measurement::M
    schedule::S
end

"""All measurement workspaces belonging to one simulation."""
struct MeasurementRuntime{T<:Tuple}
    measurements::T
end

"""Runtime marker for a program without gradient-flow measurements."""
struct NoGradientFlowMeasurementRuntime end

"""Preallocated gradient-flow field, integrator, and measurements."""
struct GradientFlowMeasurementRuntime{
    C<:GradientFlowMeasurementConfig,
    W<:GaugeConfiguration,
    F,
    M<:MeasurementRuntime,
}
    config::C
    configuration::W
    flow::F
    measurements::M
end

"""Runtime resources for direct and optional gradient-flow measurements."""
struct MeasurementProgramRuntime{D<:MeasurementRuntime,F}
    direct::D
    gradient_flow::F

    function MeasurementProgramRuntime(
        direct::D,
        gradient_flow::F,
    ) where {
        D<:MeasurementRuntime,
        F<:Union{
            NoGradientFlowMeasurementRuntime,
            GradientFlowMeasurementRuntime,
        },
    }
        return new{D,F}(direct, gradient_flow)
    end
end

function build_measurement_runtime(
    plan::MeasurementPlan,
    configuration::AbstractConfiguration,
)
    runtimes = map(plan.measurements) do scheduled
        measurement = build_qcd_measurement(
            scheduled.observable,
            configuration.gauge,
        )
        ScheduledMeasurementRuntime(
            scheduled.observable,
            measurement,
            scheduled.schedule,
        )
    end
    return MeasurementRuntime(runtimes)
end

build_measurement_runtime(plan::MeasurementPlan, simulation::Simulation) =
    build_measurement_runtime(plan, simulation.configuration)

build_gradient_flow_runtime(
    ::NoGradientFlowMeasurementConfig,
    ::AbstractConfiguration,
) = NoGradientFlowMeasurementRuntime()

function build_gradient_flow_runtime(
    config::GradientFlowMeasurementConfig,
    configuration::AbstractConfiguration,
)
    flowed_configuration = GaugeConfiguration([
        similar(link) for link in configuration.gauge
    ])
    copy_configuration!(flowed_configuration, configuration)
    flow = Gaugefields.gradient_flow(
        flowed_configuration.gauge;
        steps=config.integration_steps,
        step_size=config.step_size,
    )
    measurements = build_measurement_runtime(
        config.measurements,
        flowed_configuration,
    )
    return GradientFlowMeasurementRuntime(
        config,
        flowed_configuration,
        flow,
        measurements,
    )
end

function build_measurement_runtime(
    program::MeasurementProgram,
    configuration::AbstractConfiguration,
)
    return MeasurementProgramRuntime(
        build_measurement_runtime(program.direct, configuration),
        build_gradient_flow_runtime(program.gradient_flow, configuration),
    )
end

build_measurement_runtime(program::MeasurementProgram, simulation::Simulation) =
    build_measurement_runtime(program, simulation.configuration)

"""Coordinates identifying the configuration and optional flow time."""
struct MeasurementPoint{F}
    trajectory::Int
    flow_time::F

    function MeasurementPoint(trajectory::Integer, flow_time::F) where {F}
        trajectory >= 0 || throw(ArgumentError(
            "measurement trajectory must be nonnegative; got $trajectory",
        ))
        return new{F}(Int(trajectory), flow_time)
    end
end

MeasurementPoint(trajectory::Integer) = MeasurementPoint(trajectory, nothing)

"""A numerical measurement result with stable observable metadata."""
struct MeasurementRecord{N,V,P<:MeasurementPoint}
    name::N
    value::V
    point::P
end

snapshot_measurement_value(value::Number) = value
snapshot_measurement_value(value::AbstractString) = value
snapshot_measurement_value(value) = deepcopy(value)

function measure_one!(
    runtime::ScheduledMeasurementRuntime,
    gauge,
    point::MeasurementPoint,
)
    output = QCDMeasurements.measure(runtime.measurement, gauge)
    value = snapshot_measurement_value(QCDMeasurements.get_value(output))
    return MeasurementRecord(
        observable_name(runtime.config),
        value,
        point,
    )
end

"""Measure every configured observable, independently of its schedule."""
function measure_now!(
    runtime::MeasurementRuntime,
    simulation::Simulation;
    trajectory::Integer=simulation.state.trajectory,
    flow_time=nothing,
)
    point = MeasurementPoint(trajectory, flow_time)
    gauge = simulation.configuration.gauge
    return map(
        measurement -> measure_one!(measurement, gauge, point),
        runtime.measurements,
    )
end

function measure_due_measurements!(
    measurements::Tuple,
    gauge,
    point::MeasurementPoint,
    schedule_index::Integer,
)
    isempty(measurements) && return ()
    current = first(measurements)
    remaining = measure_due_measurements!(
        Base.tail(measurements),
        gauge,
        point,
        schedule_index,
    )
    if is_due(current.schedule, schedule_index)
        return (measure_one!(current, gauge, point), remaining...)
    end
    return remaining
end

"""Measure only observables due at the selected trajectory."""
function measure_due!(
    runtime::MeasurementRuntime,
    simulation::Simulation;
    trajectory::Integer=simulation.state.trajectory,
    flow_time=nothing,
)
    point = MeasurementPoint(trajectory, flow_time)
    return measure_due_measurements!(
        runtime.measurements,
        simulation.configuration.gauge,
        point,
        trajectory,
    )
end

function append_measurement_records!(
    destination::Vector{MeasurementRecord},
    records,
)
    append!(destination, records)
    return destination
end

measure_gradient_flow_due!(
    ::NoGradientFlowMeasurementRuntime,
    ::Simulation;
    trajectory::Integer,
) = MeasurementRecord[]

function measure_gradient_flow_due!(
    runtime::GradientFlowMeasurementRuntime,
    simulation::Simulation;
    trajectory::Integer,
)
    copy_configuration!(
        runtime.configuration,
        simulation.configuration,
    )
    records = MeasurementRecord[]
    for flow_index in 1:runtime.config.samples
        Gaugefields.flow!(runtime.configuration.gauge, runtime.flow)
        flow_time = flow_index *
                    runtime.config.integration_steps *
                    runtime.config.step_size
        point = MeasurementPoint(trajectory, flow_time)
        current = measure_due_measurements!(
            runtime.measurements.measurements,
            runtime.configuration.gauge,
            point,
            flow_index,
        )
        append_measurement_records!(records, current)
    end
    return records
end

measure_gradient_flow_now!(
    ::NoGradientFlowMeasurementRuntime,
    ::Simulation;
    trajectory::Integer,
) = MeasurementRecord[]

function measure_gradient_flow_now!(
    runtime::GradientFlowMeasurementRuntime,
    simulation::Simulation;
    trajectory::Integer,
)
    copy_configuration!(
        runtime.configuration,
        simulation.configuration,
    )
    records = MeasurementRecord[]
    for flow_index in 1:runtime.config.samples
        Gaugefields.flow!(runtime.configuration.gauge, runtime.flow)
        flow_time = flow_index *
                    runtime.config.integration_steps *
                    runtime.config.step_size
        point = MeasurementPoint(trajectory, flow_time)
        current = map(
            measurement -> measure_one!(
                measurement,
                runtime.configuration.gauge,
                point,
            ),
            runtime.measurements.measurements,
        )
        append_measurement_records!(records, current)
    end
    return records
end

function combined_measurement_records(
    direct::Tuple,
    flowed::Vector{MeasurementRecord},
)
    records = MeasurementRecord[]
    append_measurement_records!(records, direct)
    append_measurement_records!(records, flowed)
    return records
end

function measure_due!(
    runtime::MeasurementProgramRuntime{D,NoGradientFlowMeasurementRuntime},
    simulation::Simulation;
    trajectory::Integer=simulation.state.trajectory,
    flow_time=nothing,
) where {D<:MeasurementRuntime}
    flow_time === nothing || throw(ArgumentError(
        "flow_time is determined by the measurement program",
    ))
    return measure_due!(runtime.direct, simulation; trajectory)
end

function measure_due!(
    runtime::MeasurementProgramRuntime,
    simulation::Simulation;
    trajectory::Integer=simulation.state.trajectory,
    flow_time=nothing,
)
    flow_time === nothing || throw(ArgumentError(
        "flow_time is determined by the measurement program",
    ))
    direct = measure_due!(runtime.direct, simulation; trajectory)
    flowed = measure_gradient_flow_due!(
        runtime.gradient_flow,
        simulation;
        trajectory,
    )
    return combined_measurement_records(direct, flowed)
end

function measure_now!(
    runtime::MeasurementProgramRuntime{D,NoGradientFlowMeasurementRuntime},
    simulation::Simulation;
    trajectory::Integer=simulation.state.trajectory,
    flow_time=nothing,
) where {D<:MeasurementRuntime}
    flow_time === nothing || throw(ArgumentError(
        "flow_time is determined by the measurement program",
    ))
    return measure_now!(runtime.direct, simulation; trajectory)
end

function measure_now!(
    runtime::MeasurementProgramRuntime,
    simulation::Simulation;
    trajectory::Integer=simulation.state.trajectory,
    flow_time=nothing,
)
    flow_time === nothing || throw(ArgumentError(
        "flow_time is determined by the measurement program",
    ))
    direct = measure_now!(runtime.direct, simulation; trajectory)
    flowed = measure_gradient_flow_now!(
        runtime.gradient_flow,
        simulation;
        trajectory,
    )
    return combined_measurement_records(direct, flowed)
end

"""Default sink: return records without I/O or printing."""
struct NoMeasurementSink end

"""Forward every produced record to a concrete callback."""
struct FunctionMeasurementSink{F}
    callback::F
end

emit_measurements!(
    ::NoMeasurementSink,
    records::Tuple,
    ::Simulation,
) = records

function is_measurement_root(simulation::Simulation)
    communicator = Gaugefields.gauge_communicator(
        simulation.configuration.gauge,
    )
    return is_root(communicator)
end

function emit_measurements!(
    sink::FunctionMeasurementSink,
    records::Tuple,
    simulation::Simulation,
)
    is_measurement_root(simulation) && foreach(sink.callback, records)
    return records
end

function emit_measurements!(
    sink::FunctionMeasurementSink,
    records::AbstractVector,
    simulation::Simulation,
)
    is_measurement_root(simulation) && foreach(sink.callback, records)
    return records
end

emit_measurements!(
    ::NoMeasurementSink,
    records::AbstractVector,
    ::Simulation,
) = records

"""A simulation composed with measurement workspaces and an output sink."""
struct SimulationRunner{S,M,K}
    simulation::S
    measurements::M
    sink::K
end

function build_simulation_runner(
    simulation::Simulation,
    plan::MeasurementPlan;
    sink=NoMeasurementSink(),
)
    measurements = build_measurement_runtime(plan, simulation)
    return SimulationRunner(simulation, measurements, sink)
end

function build_simulation_runner(
    simulation::Simulation,
    program::MeasurementProgram;
    sink=NoMeasurementSink(),
)
    measurements = build_measurement_runtime(program, simulation)
    return SimulationRunner(simulation, measurements, sink)
end

function build_simulation_runner(
    simulation::Simulation,
    parameters::Params;
    sink=NoMeasurementSink(),
)
    return build_simulation_runner(
        simulation,
        measurement_program(parameters);
        sink,
    )
end

"""Update diagnostics and the measurements produced after that update."""
struct SimulationStepResult{U,M}
    update::U
    measurements::M
end

function step!(runner::SimulationRunner)
    update_result = update!(runner.simulation)
    records = measure_due!(runner.measurements, runner.simulation)
    emit_measurements!(runner.sink, records, runner.simulation)
    return SimulationStepResult(update_result, records)
end

function measure_now!(
    runner::SimulationRunner;
    trajectory::Integer=runner.simulation.state.trajectory,
    flow_time=nothing,
)
    records = measure_now!(
        runner.measurements,
        runner.simulation;
        trajectory,
        flow_time,
    )
    emit_measurements!(runner.sink, records, runner.simulation)
    return records
end

function measure_due!(
    runner::SimulationRunner;
    trajectory::Integer=runner.simulation.state.trajectory,
    flow_time=nothing,
)
    records = measure_due!(
        runner.measurements,
        runner.simulation;
        trajectory,
        flow_time,
    )
    emit_measurements!(runner.sink, records, runner.simulation)
    return records
end

"""Run a fixed number of update-and-measure steps without retaining history."""
function run!(runner::SimulationRunner; steps::Integer)
    steps >= 0 || throw(ArgumentError("steps must be nonnegative; got $steps"))
    for _ in 1:steps
        step!(runner)
    end
    return runner
end

function print_observable(io::IO, observable::AbstractObservableConfig)
    print(io, observable_name(observable))
    return nothing
end

function print_observable(io::IO, observable::TopologicalChargeObservableConfig)
    print(io, observable_name(observable), " kinds=")
    show(io, observable.kinds)
    print(io, " improved=", observable.improved_definition)
    return nothing
end

function print_observable(io::IO, observable::WilsonLoopObservableConfig)
    print(
        io,
        observable_name(observable),
        " Tmax=",
        observable.Tmax,
        " Rmax=",
        observable.Rmax,
    )
    return nothing
end

function print_observable(
    io::IO,
    observable::Union{
        ChiralCondensateObservableConfig,
        PionCorrelatorObservableConfig,
    },
)
    print(io, observable_name(observable), " fermion=")
    if observable.fermion isa WilsonMeasurementFermionConfig
        print(io, observable.fermion.clover ? "WilsonClover" : "Wilson")
    else
        print(io, "Staggered")
    end
    return nothing
end

"""Print a human-readable typed measurement plan."""
function show_config(io::IO, plan::MeasurementPlan)
    println(io, "MeasurementPlan")
    if isempty(plan.measurements)
        println(io, "  measurements: none")
        return nothing
    end
    for scheduled in plan.measurements
        print(io, "  ")
        print_observable(io, scheduled.observable)
        println(
            io,
            " every=",
            scheduled.schedule.every,
            " start=",
            scheduled.schedule.start,
        )
    end
    return nothing
end

function Base.show(io::IO, plan::MeasurementPlan)
    print(io, "MeasurementPlan(measurements=", length(plan.measurements), ")")
end

function Base.show(io::IO, ::MIME"text/plain", plan::MeasurementPlan)
    return show_config(io, plan)
end

function show_config(io::IO, ::NoGradientFlowMeasurementConfig)
    println(io, "GradientFlowMeasurements: disabled")
    return nothing
end

function show_config(io::IO, config::GradientFlowMeasurementConfig)
    println(io, "GradientFlowMeasurements")
    println(io, "  step_size: ", config.step_size)
    println(io, "  samples: ", config.samples)
    println(io, "  integration_steps: ", config.integration_steps)
    for scheduled in config.measurements.measurements
        print(io, "  ")
        print_observable(io, scheduled.observable)
        println(
            io,
            " every_flow_sample=",
            scheduled.schedule.every,
        )
    end
    return nothing
end

function show_config(io::IO, program::MeasurementProgram)
    println(io, "MeasurementProgram")
    println(io, "Direct measurements:")
    if isempty(program.direct.measurements)
        println(io, "  none")
    else
        for scheduled in program.direct.measurements
            print(io, "  ")
            print_observable(io, scheduled.observable)
            println(
                io,
                " every=",
                scheduled.schedule.every,
                " start=",
                scheduled.schedule.start,
            )
        end
    end
    show_config(io, program.gradient_flow)
    return nothing
end

function Base.show(io::IO, config::GradientFlowMeasurementConfig)
    print(
        io,
        "GradientFlowMeasurementConfig(samples=",
        config.samples,
        ", measurements=",
        length(config.measurements.measurements),
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", config::GradientFlowMeasurementConfig)
    return show_config(io, config)
end

function Base.show(io::IO, program::MeasurementProgram)
    print(
        io,
        "MeasurementProgram(direct=",
        length(program.direct.measurements),
        ", gradient_flow=",
        program.gradient_flow isa NoGradientFlowMeasurementConfig ?
            0 : length(program.gradient_flow.measurements.measurements),
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", program::MeasurementProgram)
    return show_config(io, program)
end

end
