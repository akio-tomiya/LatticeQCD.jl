module SimulationSession_module

import Gaugefields
import JLD2
import Serialization

import ..LQCDCommunication: broadcast!, comm_rank, is_root

import ..System_parameters: Params
import ..LQCDConfig_module:
    ConfigurationSequenceConfig,
    HMCConfig,
    HISQDiracConfig,
    LQCDConfig,
    SLHMCConfig,
    show_config
import ..Simulation_module:
    ConfigurationSequenceUpdateResult,
    ConfigurationSequenceState,
    GaugefieldsEnvironment,
    HeatbathUpdateResult,
    HMCUpdateResult,
    HMCState,
    SLHMCUpdateResult,
    Simulation,
    build_simulation,
    current_configuration_path,
    load_configuration!,
    save_configuration,
    update!
import ..MeasurementPlan_module:
    MeasurementPlan,
    MeasurementProgram,
    MeasurementRecord,
    NoGradientFlowMeasurementConfig,
    build_measurement_runtime,
    legacy_measurement_program,
    measure_due!,
    measurement_program,
    step!,
    run!

"""A run without direct or gradient-flow measurements."""
empty_measurement_program() = MeasurementProgram(
    MeasurementPlan(),
    NoGradientFlowMeasurementConfig(),
)

"""
Typed execution counts and measurement settings.

`initial_trajectory` is the trajectory number of the input configuration.
The first completed update is therefore `initial_trajectory + 1`.
"""
struct SimulationSchedule{M<:MeasurementProgram}
    thermalization_steps::Int
    production_steps::Int
    initial_trajectory::Int
    measurements::M

    function SimulationSchedule(
        thermalization_steps::Integer,
        production_steps::Integer,
        initial_trajectory::Integer,
        measurements::M,
    ) where {M<:MeasurementProgram}
        thermalization_steps >= 0 || throw(ArgumentError(
            "thermalization_steps must be nonnegative; " *
            "got $thermalization_steps",
        ))
        production_steps >= 0 || throw(ArgumentError(
            "production_steps must be nonnegative; got $production_steps",
        ))
        initial_trajectory >= 0 || throw(ArgumentError(
            "initial_trajectory must be nonnegative; " *
            "got $initial_trajectory",
        ))
        return new{M}(
            Int(thermalization_steps),
            Int(production_steps),
            Int(initial_trajectory),
            measurements,
        )
    end
end

SimulationSchedule(
    thermalization_steps::Integer,
    production_steps::Integer;
    initial_trajectory::Integer=0,
    measurements=empty_measurement_program(),
) = SimulationSchedule(
    thermalization_steps,
    production_steps,
    initial_trajectory,
    measurements,
)

"""Do not write gauge configurations during a run."""
struct NoConfigurationOutput end

"""Portable, backend-independent JLD2 configuration checkpoints."""
struct JLD2ConfigurationOutput{
    D<:AbstractString,
    P<:AbstractString,
}
    directory::D
    prefix::P
    every::Int
    width::Int

    function JLD2ConfigurationOutput(
        directory::D,
        prefix::P,
        every::Integer,
        width::Integer,
    ) where {D<:AbstractString,P<:AbstractString}
        isempty(strip(directory)) && throw(ArgumentError(
            "the JLD2 output directory must not be empty",
        ))
        isempty(strip(prefix)) && throw(ArgumentError(
            "the JLD2 filename prefix must not be empty",
        ))
        every > 0 || throw(ArgumentError(
            "the JLD2 save interval must be positive; got $every",
        ))
        width > 0 || throw(ArgumentError(
            "the JLD2 trajectory width must be positive; got $width",
        ))
        return new{D,P}(
            directory,
            prefix,
            Int(every),
            Int(width),
        )
    end
end

JLD2ConfigurationOutput(
    directory::AbstractString;
    prefix::AbstractString="conf_",
    every::Integer=1,
    width::Integer=8,
) = JLD2ConfigurationOutput(directory, prefix, every, width)

"""Do not write restart checkpoints."""
struct NoCheckpointOutput end

"""
Portable HMC restart checkpoints written at completed trajectory boundaries.

Each checkpoint is first assembled in a `.pending` file.  Rank zero appends
the HMC/session state and atomically renames the completed file, so a failed
write never replaces the last complete checkpoint.
"""
struct JLD2CheckpointOutput{
    D<:AbstractString,
    P<:AbstractString,
}
    directory::D
    prefix::P
    every::Int
    width::Int

    function JLD2CheckpointOutput(
        directory::D,
        prefix::P,
        every::Integer,
        width::Integer,
    ) where {D<:AbstractString,P<:AbstractString}
        isempty(strip(directory)) && throw(ArgumentError(
            "the checkpoint directory must not be empty",
        ))
        isempty(strip(prefix)) && throw(ArgumentError(
            "the checkpoint filename prefix must not be empty",
        ))
        every > 0 || throw(ArgumentError(
            "the checkpoint interval must be positive; got $every",
        ))
        width > 0 || throw(ArgumentError(
            "the checkpoint trajectory width must be positive; got $width",
        ))
        return new{D,P}(directory, prefix, Int(every), Int(width))
    end
end

JLD2CheckpointOutput(
    directory::AbstractString;
    prefix::AbstractString="restart_",
    every::Integer=10,
    width::Integer=8,
) = JLD2CheckpointOutput(directory, prefix, every, width)

"""Legacy Bridge-text configuration checkpoints."""
struct BridgeTextConfigurationOutput{
    D<:AbstractString,
    P<:AbstractString,
}
    directory::D
    prefix::P
    every::Int
    width::Int

    function BridgeTextConfigurationOutput(
        directory::D,
        prefix::P,
        every::Integer,
        width::Integer,
    ) where {D<:AbstractString,P<:AbstractString}
        isempty(strip(directory)) && throw(ArgumentError(
            "the Bridge-text output directory must not be empty",
        ))
        isempty(strip(prefix)) && throw(ArgumentError(
            "the Bridge-text filename prefix must not be empty",
        ))
        every > 0 || throw(ArgumentError(
            "the Bridge-text save interval must be positive; got $every",
        ))
        width > 0 || throw(ArgumentError(
            "the Bridge-text trajectory width must be positive; got $width",
        ))
        return new{D,P}(directory, prefix, Int(every), Int(width))
    end
end

BridgeTextConfigurationOutput(
    directory::AbstractString;
    prefix::AbstractString="conf_",
    every::Integer=1,
    width::Integer=8,
) = BridgeTextConfigurationOutput(directory, prefix, every, width)

"""Legacy ILDG configuration checkpoints."""
struct ILDGConfigurationOutput{
    D<:AbstractString,
    P<:AbstractString,
}
    directory::D
    prefix::P
    every::Int
    width::Int

    function ILDGConfigurationOutput(
        directory::D,
        prefix::P,
        every::Integer,
        width::Integer,
    ) where {D<:AbstractString,P<:AbstractString}
        isempty(strip(directory)) && throw(ArgumentError(
            "the ILDG output directory must not be empty",
        ))
        isempty(strip(prefix)) && throw(ArgumentError(
            "the ILDG filename prefix must not be empty",
        ))
        every > 0 || throw(ArgumentError(
            "the ILDG save interval must be positive; got $every",
        ))
        width > 0 || throw(ArgumentError(
            "the ILDG trajectory width must be positive; got $width",
        ))
        return new{D,P}(directory, prefix, Int(every), Int(width))
    end
end

ILDGConfigurationOutput(
    directory::AbstractString;
    prefix::AbstractString="conf_",
    every::Integer=1,
    width::Integer=8,
) = ILDGConfigurationOutput(directory, prefix, every, width)

"""All side-effect settings, kept separate from the physics config."""
struct OutputConfig{C,K}
    configurations::C
    checkpoints::K

    function OutputConfig(configurations::C, checkpoints::K) where {
        C<:Union{
            NoConfigurationOutput,
            JLD2ConfigurationOutput,
            BridgeTextConfigurationOutput,
            ILDGConfigurationOutput,
        },
        K<:Union{
            NoCheckpointOutput,
            JLD2CheckpointOutput,
        },
    }
        return new{C,K}(configurations, checkpoints)
    end
end

OutputConfig(configurations) = OutputConfig(
    configurations,
    NoCheckpointOutput(),
)

OutputConfig(;
    configurations=NoConfigurationOutput(),
    checkpoints=NoCheckpointOutput(),
) = OutputConfig(configurations, checkpoints)

"""Serializable, GUI-neutral specification for one simulation run."""
struct SimulationSpec{
    C<:LQCDConfig,
    S<:SimulationSchedule,
    O<:OutputConfig,
}
    config::C
    schedule::S
    output::O
end

SimulationSpec(config::LQCDConfig, schedule::SimulationSchedule) =
    SimulationSpec(config, schedule, OutputConfig())

"""One field-addressable validation problem suitable for a GUI form."""
struct ValidationIssue{
    P<:AbstractVector{String},
    C<:Symbol,
    M<:AbstractString,
}
    path::P
    code::C
    message::M
end

function validation_issue(path, code::Symbol, message::AbstractString)
    return ValidationIssue(String.(collect(path)), code, String(message))
end

"""Return cross-component validation issues without throwing."""
function validate(spec::SimulationSpec)
    issues = ValidationIssue{Vector{String},Symbol,String}[]
    for (index, fermion) in enumerate(spec.config.fermions)
        if fermion.operator isa HISQDiracConfig
            if spec.config.gauge.halo < 3
                push!(issues, validation_issue(
                    ("config", "gauge", "halo"),
                    :minimum,
                    "HISQ fermion $index requires a halo width of at least 3",
                ))
            end
        end
    end
    if spec.config.update isa ConfigurationSequenceConfig
        if spec.schedule.thermalization_steps != 0
            push!(issues, validation_issue(
                ("schedule", "thermalization_steps"),
                :not_supported,
                "configuration-file sequences cannot be thermalized",
            ))
        end
        if spec.schedule.initial_trajectory != 0
            push!(issues, validation_issue(
                ("schedule", "initial_trajectory"),
                :must_be_zero,
                "configuration-file sequences determine their own position",
            ))
        end
        if !(spec.output.configurations isa NoConfigurationOutput)
            push!(issues, validation_issue(
                ("output", "configurations"),
                :not_supported,
                "configuration-file sequences are input-only and are not " *
                "saved again by the runner",
            ))
        end
    end
    if !(spec.output.checkpoints isa NoCheckpointOutput) &&
       !(spec.config.update isa Union{HMCConfig,SLHMCConfig})
        push!(issues, validation_issue(
            ("output", "checkpoints"),
            :not_supported,
            "restart checkpoints currently require HMC or SLHMC",
        ))
    end
    return issues
end

function throw_validation_issues(spec::SimulationSpec)
    issues = validate(spec)
    isempty(issues) && return nothing
    messages = map(issues) do issue
        "$(join(issue.path, '.')) [$(issue.code)]: $(issue.message)"
    end
    throw(ArgumentError(join(messages, "\n")))
end

"""
Translate legacy trajectory thresholds without changing their behavior.

Legacy `Params` updates trajectories `initialtrj:Nsteps` and begins measuring
when the trajectory number reaches `Nthermalization`.  The typed schedule uses
explicit counts, so the threshold is converted into a thermalization count.
"""
function legacy_simulation_schedule(parameters)
    first_trajectory = parameters.initialtrj
    last_trajectory = parameters.Nsteps
    first_trajectory >= 1 || throw(ArgumentError(
        "typed update schedules require initialtrj >= 1; " *
        "got $first_trajectory",
    ))
    total_steps = max(0, last_trajectory - first_trajectory + 1)
    first_production = max(first_trajectory, parameters.Nthermalization)
    thermalization_steps = clamp(
        first_production - first_trajectory,
        0,
        total_steps,
    )
    production_steps = total_steps - thermalization_steps
    return SimulationSchedule(
        thermalization_steps,
        production_steps,
        first_trajectory - 1,
        legacy_measurement_program(parameters),
    )
end


simulation_schedule(parameters::Params) = legacy_simulation_schedule(parameters)

function legacy_output_config(parameters)
    # The legacy Fileloading runner deliberately ignores saveU_* settings.
    parameters.update_method == "Fileloading" && return OutputConfig()
    configurations = if isnothing(parameters.saveU_format)
        NoConfigurationOutput()
    else
        format = lowercase(replace(
            strip(String(parameters.saveU_format)),
            "_" => "",
            "-" => "",
        ))
        if format in ("jld", "jld2")
            JLD2ConfigurationOutput(
                parameters.saveU_dir;
                every=parameters.saveU_every,
            )
        elseif format in ("bridgetext", "bridge", "text")
            BridgeTextConfigurationOutput(
                parameters.saveU_dir;
                every=parameters.saveU_every,
            )
        elseif format == "ildg"
            ILDGConfigurationOutput(
                parameters.saveU_dir;
                every=parameters.saveU_every,
            )
        else
            throw(ArgumentError(
                "unsupported saveU_format=$(repr(parameters.saveU_format)); " *
                "expected JLD, ILDG, or BridgeText",
            ))
        end
    end

    checkpoint_every = hasproperty(parameters, :checkpoint_every) ?
        Int(parameters.checkpoint_every) : 0
    checkpoint_every >= 0 || throw(ArgumentError(
        "checkpoint_every must be nonnegative; got $checkpoint_every",
    ))
    checkpoints = if iszero(checkpoint_every)
        NoCheckpointOutput()
    else
        checkpoint_dir = String(parameters.checkpoint_dir)
        JLD2CheckpointOutput(checkpoint_dir; every=checkpoint_every)
    end
    return OutputConfig(configurations, checkpoints)
end


output_config(parameters::Params) = legacy_output_config(parameters)

SimulationSpec(parameters::Params) = SimulationSpec(
    LQCDConfig(parameters),
    simulation_schedule(parameters),
    output_config(parameters),
)

abstract type AbstractSimulationEvent end

struct RunStarted{S<:SimulationSchedule} <: AbstractSimulationEvent
    schedule::S
    trajectory::Int
end

struct ThermalizationStepFinished{U} <: AbstractSimulationEvent
    trajectory::Int
    step::Int
    update::U
end

struct TrajectoryFinished{U} <: AbstractSimulationEvent
    trajectory::Int
    step::Int
    update::U
end

struct MeasurementFinished{R<:MeasurementRecord} <:
       AbstractSimulationEvent
    record::R
end

struct ConfigurationSaved{P<:AbstractString} <: AbstractSimulationEvent
    trajectory::Int
    path::P
    format::Symbol
end

"""A complete restart checkpoint became visible at a trajectory boundary."""
struct CheckpointSaved{P<:AbstractString} <: AbstractSimulationEvent
    trajectory::Int
    path::P
end

"""A session was restored from a restart checkpoint."""
struct CheckpointLoaded{P<:AbstractString} <: AbstractSimulationEvent
    trajectory::Int
    path::P
end

"""A configuration from a finite Fileloading input sequence became current."""
struct ConfigurationLoaded{P<:AbstractString} <: AbstractSimulationEvent
    trajectory::Int
    current_index::Int
    path::P
end

struct SimulationRunSummary
    initial_trajectory::Int
    final_trajectory::Int
    thermalization_completed::Int
    production_completed::Int
    saved_configurations::Int
    saved_checkpoints::Int
    stopped::Bool
end

struct RunStopped <: AbstractSimulationEvent
    summary::SimulationRunSummary
end

struct RunFinished <: AbstractSimulationEvent
    summary::SimulationRunSummary
end

struct RunFailed{M<:AbstractString} <: AbstractSimulationEvent
    trajectory::Int
    message::M
end

struct NoSimulationEventSink end

"""Print rank-zero simulation progress in a terminal or notebook."""
struct ConsoleSimulationEventSink{I<:IO}
    io::I
end

ConsoleSimulationEventSink() = ConsoleSimulationEventSink(stdout)

struct FunctionSimulationEventSink{F}
    callback::F
end

mutable struct RecordingSimulationEventSink{E<:AbstractVector}
    events::E
end

RecordingSimulationEventSink() = RecordingSimulationEventSink(Any[])

struct CompositeSimulationEventSink{S<:Tuple}
    sinks::S
end

CompositeSimulationEventSink(sinks...) = CompositeSimulationEventSink(sinks)

deliver_event!(::NoSimulationEventSink, event) = event

function print_update_summary(io::IO, update::HMCUpdateResult)
    print(
        io,
        "accepted=", update.accepted,
        " H_initial=", update.initial_hamiltonian,
        " H_final=", update.final_hamiltonian,
        " delta_H=", update.delta_hamiltonian,
    )
end

function print_update_summary(io::IO, update::SLHMCUpdateResult)
    print(
        io,
        "accepted=", update.accepted,
        " target_H_initial=", update.initial_hamiltonian,
        " target_H_final=", update.final_hamiltonian,
        " target_delta_H=", update.delta_hamiltonian,
        " md_delta_H=", update.md_delta_hamiltonian,
    )
end

function print_update_summary(io::IO, update::HeatbathUpdateResult)
    print(
        io,
        "heatbath_sweep=", update.heatbath_sweep,
        " overrelaxation_sweep=", update.overrelaxation_sweep,
    )
end

function print_update_summary(
    io::IO,
    update::ConfigurationSequenceUpdateResult,
)
    print(
        io,
        "configuration=", update.current_index,
        " path=", repr(update.path),
    )
end

print_update_summary(io::IO, update) = print(io, nameof(typeof(update)))

function print_simulation_event(io::IO, event::RunStarted)
    println(
        io,
        "# run started: trajectory=", event.trajectory,
        " thermalization=", event.schedule.thermalization_steps,
        " production=", event.schedule.production_steps,
    )
end

function print_simulation_event(io::IO, event::ThermalizationStepFinished)
    print(
        io,
        "# thermalization step=", event.step,
        " trajectory=", event.trajectory,
        " ",
    )
    print_update_summary(io, event.update)
    println(io)
end

function print_simulation_event(io::IO, event::TrajectoryFinished)
    print(io, "# trajectory=", event.trajectory, " step=", event.step, " ")
    print_update_summary(io, event.update)
    println(io)
end

function print_simulation_event(io::IO, event::MeasurementFinished)
    record = event.record
    print(io, "# measurement: trajectory=", record.point.trajectory)
    record.point.flow_time === nothing ||
        print(io, " flow_time=", record.point.flow_time)
    print(io, " ", record.name, " = ")
    show(IOContext(io, :compact => true, :limit => true), record.value)
    println(io)
end

function print_simulation_event(io::IO, event::ConfigurationSaved)
    println(
        io,
        "# configuration saved: trajectory=", event.trajectory,
        " format=", event.format,
        " path=", repr(event.path),
    )
end

function print_simulation_event(io::IO, event::CheckpointSaved)
    println(
        io,
        "# restart checkpoint saved: trajectory=", event.trajectory,
        " path=", repr(event.path),
    )
end

function print_simulation_event(io::IO, event::CheckpointLoaded)
    println(
        io,
        "# restart checkpoint loaded: trajectory=", event.trajectory,
        " path=", repr(event.path),
    )
end

function print_simulation_event(io::IO, event::ConfigurationLoaded)
    println(
        io,
        "# configuration loaded: trajectory=", event.trajectory,
        " index=", event.current_index,
        " path=", repr(event.path),
    )
end

function print_simulation_event(io::IO, event::RunStopped)
    println(io, "# run stopped: ", event.summary)
end

function print_simulation_event(io::IO, event::RunFinished)
    println(io, "# run finished: ", event.summary)
end

function print_simulation_event(io::IO, event::RunFailed)
    println(
        io,
        "# run failed: trajectory=", event.trajectory,
        " error=", event.message,
    )
end

function deliver_event!(sink::ConsoleSimulationEventSink, event)
    print_simulation_event(sink.io, event)
    flush(sink.io)
    return event
end

function deliver_event!(sink::FunctionSimulationEventSink, event)
    sink.callback(event)
    return event
end

function deliver_event!(sink::RecordingSimulationEventSink, event)
    push!(sink.events, event)
    return event
end

function deliver_event!(sink::CompositeSimulationEventSink, event)
    foreach(child -> deliver_event!(child, event), sink.sinks)
    return event
end

"""Mutable progress and cooperative stop state; no physics objects live here."""
mutable struct SimulationSessionState
    thermalization_completed::Int
    production_completed::Int
    saved_configurations::Int
    saved_checkpoints::Int
    sequence_initial_pending::Bool
    stop_requested::Base.Threads.Atomic{Bool}
    running::Base.Threads.Atomic{Bool}
end

SimulationSessionState(sequence_initial_pending::Bool=false) = SimulationSessionState(
    0,
    0,
    0,
    0,
    sequence_initial_pending,
    Base.Threads.Atomic{Bool}(false),
    Base.Threads.Atomic{Bool}(false),
)

"""Concrete runtime used by CLI, notebooks, GUIs, and remote frontends."""
struct SimulationSession{S,C,M,O,K,R}
    simulation::S
    schedule::C
    measurements::M
    output::O
    sink::K
    state::R
end


function resolve_session_schedule(
    schedule::SimulationSchedule,
    simulation::Simulation,
)
    simulation.state isa ConfigurationSequenceState || return schedule
    return SimulationSchedule(
        0,
        length(simulation.updater.paths),
        0,
        schedule.measurements,
    )
end

function initialize_trajectory!(
    simulation::Simulation,
    schedule::SimulationSchedule,
)
    if simulation.state isa ConfigurationSequenceState
        schedule.initial_trajectory == 0 || throw(ArgumentError(
            "configuration-file sequences require initial_trajectory=0",
        ))
        return simulation
    end
    simulation.state.trajectory = schedule.initial_trajectory
    return simulation
end

function build_simulation_session(
    spec::SimulationSpec,
    environment::GaugefieldsEnvironment=GaugefieldsEnvironment();
    sink=NoSimulationEventSink(),
)
    throw_validation_issues(spec)
    simulation = build_simulation(spec.config, environment)
    initialize_trajectory!(simulation, spec.schedule)
    schedule = resolve_session_schedule(spec.schedule, simulation)
    measurements = build_measurement_runtime(
        schedule.measurements,
        simulation,
    )
    return SimulationSession(
        simulation,
        schedule,
        measurements,
        spec.output,
        sink,
        SimulationSessionState(
            simulation.state isa ConfigurationSequenceState,
        ),
    )
end

build_simulation(
    spec::SimulationSpec,
    environment::GaugefieldsEnvironment=GaugefieldsEnvironment();
    sink=NoSimulationEventSink(),
) = build_simulation_session(spec, environment; sink)

function session_status(session::SimulationSession)
    is_running(session) && return :running
    is_finished(session) && return :finished
    session.state.stop_requested[] && return :stop_requested
    return :ready
end

function Base.show(io::IO, session::SimulationSession)
    print(io, "SimulationSession(")
    show(io, session.simulation)
    print(
        io,
        ", thermalization=",
        session.state.thermalization_completed,
        "/",
        session.schedule.thermalization_steps,
        ", production=",
        session.state.production_completed,
        "/",
        session.schedule.production_steps,
        ", status=",
        session_status(session),
        ")",
    )
end

function show_config(io::IO, session::SimulationSession)
    println(io, "SimulationSession")
    print(io, "  simulation: ")
    show(io, session.simulation)
    println(io)
    println(
        io,
        "  thermalization: ",
        session.state.thermalization_completed,
        " / ",
        session.schedule.thermalization_steps,
    )
    println(
        io,
        "  production: ",
        session.state.production_completed,
        " / ",
        session.schedule.production_steps,
    )
    println(
        io,
        "  saved configurations: ",
        session.state.saved_configurations,
    )
    println(
        io,
        "  saved restart checkpoints: ",
        session.state.saved_checkpoints,
    )
    println(io, "  status: ", session_status(session))
    print(io, "  measurements: ")
    show(io, session.schedule.measurements)
    println(io)
    print(io, "  output: ")
    show(io, session.output)
    println(io)
    return nothing
end

Base.show(io::IO, ::MIME"text/plain", session::SimulationSession) =
    show_config(io, session)

function is_session_root(session::SimulationSession)
    gauge = session.simulation.configuration.gauge
    communicator = Gaugefields.gauge_communicator(gauge)
    return is_root(communicator)
end

function emit_event!(session::SimulationSession, event)
    is_session_root(session) && deliver_event!(session.sink, event)
    return event
end

function configuration_output_path(
    output::Union{
        JLD2ConfigurationOutput,
        BridgeTextConfigurationOutput,
        ILDGConfigurationOutput,
    },
    trajectory::Integer,
)
    number = lpad(string(trajectory), output.width, '0')
    return joinpath(
        output.directory,
        "$(output.prefix)$number$(configuration_output_extension(output))",
    )
end


configuration_output_extension(::JLD2ConfigurationOutput) = ".jld2"
configuration_output_extension(::BridgeTextConfigurationOutput) = ".txt"
configuration_output_extension(::ILDGConfigurationOutput) = ".ildg"

configuration_output_format(::JLD2ConfigurationOutput) = :jld2
configuration_output_format(::BridgeTextConfigurationOutput) = :bridge
configuration_output_format(::ILDGConfigurationOutput) = :ildg

save_configuration_if_due!(
    ::NoConfigurationOutput,
    ::SimulationSession,
    ::Integer,
) = nothing

function save_configuration_if_due!(
    output::Union{
        JLD2ConfigurationOutput,
        BridgeTextConfigurationOutput,
        ILDGConfigurationOutput,
    },
    session::SimulationSession,
    trajectory::Integer,
)
    mod(trajectory, output.every) == 0 || return nothing
    if is_session_root(session)
        mkpath(output.directory)
    end
    path = configuration_output_path(output, trajectory)
    format = configuration_output_format(output)
    save_configuration(
        path,
        session.simulation.configuration;
        format,
    )
    session.state.saved_configurations += 1
    emit_event!(
        session,
        ConfigurationSaved(Int(trajectory), path, format),
    )
    return path
end

function save_configuration_if_due!(
    output::OutputConfig,
    session::SimulationSession,
    trajectory::Integer,
)
    return save_configuration_if_due!(
        output.configurations,
        session,
        trajectory,
    )
end

const RESTART_CHECKPOINT_FORMAT =
    "LatticeQCD.jl trajectory-boundary restart checkpoint"
const RESTART_CHECKPOINT_VERSION = 1

checkpoint_output_path(
    output::JLD2CheckpointOutput,
    trajectory::Integer,
) = joinpath(
    output.directory,
    "$(output.prefix)$(lpad(string(trajectory), output.width, '0')).jld2",
)

function checkpoint_communicator(session::SimulationSession)
    return Gaugefields.gauge_communicator(
        session.simulation.configuration.gauge,
    )
end

function root_operation(operation, communicator, description::AbstractString)
    success = Ref(true)
    message = ""
    if comm_rank(communicator) == 0
        try
            operation()
        catch exception
            success[] = false
            message = sprint(showerror, exception, catch_backtrace())
        end
    end
    broadcast!(success, 0, communicator)
    success[] || error(
        comm_rank(communicator) == 0 ? message :
        "$description failed on rank zero",
    )
    return nothing
end

function checkpoint_metadata(session::SimulationSession)
    simulation_state = session.simulation.state
    simulation_state isa HMCState || throw(ArgumentError(
        "restart checkpoints currently require HMC or SLHMC state",
    ))
    return (
        trajectory=simulation_state.trajectory,
        accepted=simulation_state.accepted,
        metropolis_rng=deepcopy(simulation_state.metropolis_rng),
        initial_trajectory=session.schedule.initial_trajectory,
        thermalization_steps=session.schedule.thermalization_steps,
        production_steps=session.schedule.production_steps,
        thermalization_completed=session.state.thermalization_completed,
        production_completed=session.state.production_completed,
        saved_configurations=session.state.saved_configurations,
        saved_checkpoints=session.state.saved_checkpoints,
    )
end

function append_checkpoint_metadata!(path::AbstractString, metadata)
    JLD2.jldopen(path, "a+") do file
        file["latticeqcd_checkpoint_format"] = RESTART_CHECKPOINT_FORMAT
        file["latticeqcd_checkpoint_version"] = RESTART_CHECKPOINT_VERSION
        file["latticeqcd_trajectory"] = metadata.trajectory
        file["latticeqcd_accepted"] = metadata.accepted
        file["latticeqcd_metropolis_rng"] = metadata.metropolis_rng
        file["latticeqcd_initial_trajectory"] = metadata.initial_trajectory
        file["latticeqcd_thermalization_steps"] =
            metadata.thermalization_steps
        file["latticeqcd_production_steps"] = metadata.production_steps
        file["latticeqcd_thermalization_completed"] =
            metadata.thermalization_completed
        file["latticeqcd_production_completed"] =
            metadata.production_completed
        file["latticeqcd_saved_configurations"] =
            metadata.saved_configurations
        file["latticeqcd_saved_checkpoints"] = metadata.saved_checkpoints
    end
    return path
end

function read_checkpoint_metadata(path::AbstractString)
    return JLD2.jldopen(path, "r") do file
        get(file, "latticeqcd_checkpoint_format", nothing) ==
            RESTART_CHECKPOINT_FORMAT || throw(ArgumentError(
            "$path is not a LatticeQCD restart checkpoint",
        ))
        version = Int(file["latticeqcd_checkpoint_version"])
        version == RESTART_CHECKPOINT_VERSION || throw(ArgumentError(
            "unsupported restart checkpoint version $version in $path",
        ))
        return (
            trajectory=Int(file["latticeqcd_trajectory"]),
            accepted=Int(file["latticeqcd_accepted"]),
            metropolis_rng=file["latticeqcd_metropolis_rng"],
            initial_trajectory=Int(file["latticeqcd_initial_trajectory"]),
            thermalization_steps=Int(
                file["latticeqcd_thermalization_steps"],
            ),
            production_steps=Int(file["latticeqcd_production_steps"]),
            thermalization_completed=Int(
                file["latticeqcd_thermalization_completed"],
            ),
            production_completed=Int(
                file["latticeqcd_production_completed"],
            ),
            saved_configurations=Int(
                file["latticeqcd_saved_configurations"],
            ),
            saved_checkpoints=Int(file["latticeqcd_saved_checkpoints"]),
        )
    end
end

function serialized_root_value(operation, communicator, description)
    success = Ref(true)
    bytes = UInt8[]
    message = ""
    if comm_rank(communicator) == 0
        try
            buffer = IOBuffer()
            Serialization.serialize(buffer, operation())
            bytes = take!(buffer)
        catch exception
            success[] = false
            message = sprint(showerror, exception, catch_backtrace())
        end
    end
    broadcast!(success, 0, communicator)
    success[] || error(
        comm_rank(communicator) == 0 ? message :
        "$description failed on rank zero",
    )
    length_buffer = Ref(length(bytes))
    broadcast!(length_buffer, 0, communicator)
    comm_rank(communicator) == 0 || resize!(bytes, length_buffer[])
    broadcast!(bytes, 0, communicator)
    return Serialization.deserialize(IOBuffer(bytes))
end

function validate_checkpoint_metadata(
    session::SimulationSession,
    metadata,
)
    session.simulation.state isa HMCState || throw(ArgumentError(
        "restart checkpoints currently require HMC or SLHMC state",
    ))
    metadata.initial_trajectory == session.schedule.initial_trajectory ||
        throw(ArgumentError(
            "checkpoint initial trajectory $(metadata.initial_trajectory) " *
            "does not match session initial trajectory " *
            "$(session.schedule.initial_trajectory)",
        ))
    metadata.thermalization_steps == session.schedule.thermalization_steps ||
        throw(ArgumentError(
            "checkpoint thermalization length " *
            "$(metadata.thermalization_steps) does not match session " *
            "length $(session.schedule.thermalization_steps)",
        ))
    metadata.production_steps == session.schedule.production_steps ||
        throw(ArgumentError(
            "checkpoint production length $(metadata.production_steps) " *
            "does not match session length " *
            "$(session.schedule.production_steps)",
        ))
    0 <= metadata.thermalization_completed <=
          session.schedule.thermalization_steps || throw(ArgumentError(
        "checkpoint thermalization progress is outside the session schedule",
    ))
    metadata.production_completed >= 0 || throw(ArgumentError(
        "checkpoint production progress must be nonnegative",
    ))
    metadata.production_completed <= session.schedule.production_steps ||
        throw(ArgumentError(
            "checkpoint has completed $(metadata.production_completed) " *
            "production steps, but the session only schedules " *
            "$(session.schedule.production_steps)",
        ))
    expected_trajectory = metadata.initial_trajectory +
        metadata.thermalization_completed + metadata.production_completed
    metadata.trajectory == expected_trajectory || throw(ArgumentError(
        "checkpoint trajectory $(metadata.trajectory) is inconsistent with " *
        "its completed-step counters (expected $expected_trajectory)",
    ))
    metadata.accepted >= 0 || throw(ArgumentError(
        "checkpoint accepted count must be nonnegative",
    ))
    typeof(metadata.metropolis_rng) ==
        typeof(session.simulation.state.metropolis_rng) ||
        throw(ArgumentError(
            "checkpoint Metropolis RNG type " *
            "$(typeof(metadata.metropolis_rng)) does not match session type " *
            "$(typeof(session.simulation.state.metropolis_rng))",
        ))
    return nothing
end

function _save_checkpoint(path::AbstractString, session::SimulationSession)
    metadata = checkpoint_metadata(session)
    validate_checkpoint_metadata(session, metadata)
    communicator = checkpoint_communicator(session)
    pending_path = path * ".pending"
    root_operation(communicator, "checkpoint preparation") do
        mkpath(dirname(path))
        rm(pending_path; force=true)
    end
    save_configuration(
        pending_path,
        session.simulation.configuration;
        format=:jld2,
    )
    root_operation(communicator, "checkpoint finalization") do
        append_checkpoint_metadata!(pending_path, metadata)
        Base.Filesystem.rename(pending_path, path)
    end
    return path
end

"""
    save_checkpoint(path, session)

Safely save a complete HMC/SLHMC restart checkpoint.  Call this only between
trajectories; a running session is rejected.  Scheduled checkpoint output uses
the same implementation internally at completed trajectory boundaries.
"""
function save_checkpoint(path::AbstractString, session::SimulationSession)
    is_running(session) && throw(ArgumentError(
        "cannot save a restart checkpoint while a trajectory is running",
    ))
    return _save_checkpoint(path, session)
end

"""
    load_checkpoint!(session, path)

Restore gauge links, HMC counters, the Metropolis RNG, and session progress.
The session must be built from a compatible specification on every rank.
"""
function load_checkpoint!(
    session::SimulationSession,
    path::AbstractString,
)
    is_running(session) && throw(ArgumentError(
        "cannot restore a checkpoint while the session is running",
    ))
    communicator = checkpoint_communicator(session)
    metadata = serialized_root_value(
        () -> read_checkpoint_metadata(path),
        communicator,
        "checkpoint metadata read",
    )
    validate_checkpoint_metadata(session, metadata)
    load_configuration!(session.simulation.configuration, path; format=:jld2)

    simulation_state = session.simulation.state
    simulation_state.trajectory = metadata.trajectory
    simulation_state.accepted = metadata.accepted
    simulation_state.metropolis_rng = metadata.metropolis_rng
    progress = session.state
    progress.thermalization_completed = metadata.thermalization_completed
    progress.production_completed = metadata.production_completed
    progress.saved_configurations = metadata.saved_configurations
    progress.saved_checkpoints = metadata.saved_checkpoints
    progress.sequence_initial_pending = false
    progress.stop_requested[] = false
    emit_event!(session, CheckpointLoaded(metadata.trajectory, path))
    return session
end

save_checkpoint_if_due!(
    ::NoCheckpointOutput,
    ::SimulationSession,
    ::Integer,
) = nothing

function save_checkpoint_if_due!(
    output::JLD2CheckpointOutput,
    session::SimulationSession,
    trajectory::Integer,
)
    mod(trajectory, output.every) == 0 || return nothing
    path = checkpoint_output_path(output, trajectory)
    session.state.saved_checkpoints += 1
    try
        _save_checkpoint(path, session)
    catch
        session.state.saved_checkpoints -= 1
        rethrow()
    end
    emit_event!(session, CheckpointSaved(Int(trajectory), path))
    return path
end

function save_checkpoint_if_due!(
    output::OutputConfig,
    session::SimulationSession,
    trajectory::Integer,
)
    return save_checkpoint_if_due!(
        output.checkpoints,
        session,
        trajectory,
    )
end

struct SessionStepResult{P<:Symbol,U,M,S,C}
    phase::P
    update::U
    measurements::M
    saved_path::S
    checkpoint_path::C
end

SessionStepResult(phase, update, measurements, saved_path) =
    SessionStepResult(phase, update, measurements, saved_path, nothing)

function emit_measurement_events!(session::SimulationSession, records)
    foreach(records) do record
        emit_event!(session, MeasurementFinished(record))
    end
    return records
end


function measure_and_save_production!(
    session::SimulationSession,
    trajectory::Integer,
    update_result,
)
    records = measure_due!(
        session.measurements,
        session.simulation;
        trajectory,
    )
    emit_measurement_events!(session, records)
    path = save_configuration_if_due!(session.output, session, trajectory)
    checkpoint_path = save_checkpoint_if_due!(
        session.output,
        session,
        trajectory,
    )
    return SessionStepResult(
        :production,
        update_result,
        records,
        path,
        checkpoint_path,
    )
end

function advance_session_step!(session::SimulationSession)
    progress = session.state
    schedule = session.schedule
    if progress.sequence_initial_pending
        progress.sequence_initial_pending = false
        progress.production_completed += 1
        trajectory = session.simulation.state.trajectory
        path = current_configuration_path(session.simulation)
        emit_event!(session, ConfigurationLoaded(
            trajectory,
            session.simulation.state.current_index,
            path,
        ))
        return measure_and_save_production!(
            session,
            trajectory,
            nothing,
        )
    end
    if progress.thermalization_completed < schedule.thermalization_steps
        update_result = update!(session.simulation)
        progress.thermalization_completed += 1
        trajectory = session.simulation.state.trajectory
        emit_event!(session, ThermalizationStepFinished(
            trajectory,
            progress.thermalization_completed,
            update_result,
        ))
        checkpoint_path = save_checkpoint_if_due!(
            session.output,
            session,
            trajectory,
        )
        return SessionStepResult(
            :thermalization,
            update_result,
            (),
            nothing,
            checkpoint_path,
        )
    end

    progress.production_completed < schedule.production_steps ||
        throw(EOFError())
    update_result = update!(session.simulation)
    progress.production_completed += 1
    trajectory = session.simulation.state.trajectory
    if update_result isa ConfigurationSequenceUpdateResult
        emit_event!(session, ConfigurationLoaded(
            trajectory,
            update_result.current_index,
            update_result.path,
        ))
    else
        emit_event!(session, TrajectoryFinished(
            trajectory,
            progress.production_completed,
            update_result,
        ))
    end
    return measure_and_save_production!(
        session,
        trajectory,
        update_result,
    )
end

function step!(session::SimulationSession)
    was_running = Base.Threads.atomic_cas!(
        session.state.running,
        false,
        true,
    )
    was_running == false || throw(ArgumentError(
        "this SimulationSession is already running",
    ))
    try
        return advance_session_step!(session)
    catch exception
        if !(exception isa EOFError)
            emit_event!(session, RunFailed(
                session.simulation.state.trajectory,
                sprint(showerror, exception),
            ))
        end
        rethrow()
    finally
        session.state.running[] = false
    end
end

function is_finished(session::SimulationSession)
    return session.state.thermalization_completed >=
           session.schedule.thermalization_steps &&
           session.state.production_completed >=
           session.schedule.production_steps
end

is_running(session::SimulationSession) = session.state.running[]

function request_stop!(session::SimulationSession)
    session.state.stop_requested[] = true
    return session
end

function run_summary(session::SimulationSession; stopped::Bool=false)
    return SimulationRunSummary(
        session.schedule.initial_trajectory,
        session.simulation.state.trajectory,
        session.state.thermalization_completed,
        session.state.production_completed,
        session.state.saved_configurations,
        session.state.saved_checkpoints,
        stopped,
    )
end

"""
Run until the schedule is complete or a cooperative stop is requested.

Stop requests are checked only between complete trajectories, after HMC
accept/reject and rollback and after any due measurement/checkpoint.
"""
function session_with_sink(session::SimulationSession, sink)
    return SimulationSession(
        session.simulation,
        session.schedule,
        session.measurements,
        session.output,
        sink,
        session.state,
    )
end

"""
Run a complete simulation schedule.

Progress is printed on rank zero by default. Pass `verbose=false` for GUI,
batch, or library use; events are still delivered to the session's own sink.
"""
function run!(
    session::SimulationSession;
    verbose::Bool=true,
    io::IO=stdout,
)
    active_session = if verbose
        session_with_sink(
            session,
            CompositeSimulationEventSink(
                session.sink,
                ConsoleSimulationEventSink(io),
            ),
        )
    else
        session
    end
    was_running = Base.Threads.atomic_cas!(
        active_session.state.running,
        false,
        true,
    )
    was_running == false || throw(ArgumentError(
        "this SimulationSession is already running",
    ))
    active_session.state.stop_requested[] = false
    try
        emit_event!(active_session, RunStarted(
            active_session.schedule,
            active_session.simulation.state.trajectory,
        ))
        while !is_finished(active_session) &&
              !active_session.state.stop_requested[]
            advance_session_step!(active_session)
        end
        stopped = active_session.state.stop_requested[] &&
                  !is_finished(active_session)
        summary = run_summary(active_session; stopped)
        if stopped
            emit_event!(active_session, RunStopped(summary))
        else
            emit_event!(active_session, RunFinished(summary))
        end
        return summary
    catch exception
        message = sprint(showerror, exception)
        emit_event!(active_session, RunFailed(
            active_session.simulation.state.trajectory,
            message,
        ))
        rethrow()
    finally
        active_session.state.running[] = false
    end
end

function Base.show(io::IO, summary::SimulationRunSummary)
    print(
        io,
        "SimulationRunSummary(initial=", summary.initial_trajectory,
        ", final=", summary.final_trajectory,
        ", thermalization=", summary.thermalization_completed,
        ", production=", summary.production_completed,
        ", saved=", summary.saved_configurations,
        ", checkpoints=", summary.saved_checkpoints,
        ", stopped=", summary.stopped,
        ")",
    )
end

function Base.show(io::IO, ::MIME"text/plain", summary::SimulationRunSummary)
    println(io, "SimulationRunSummary")
    println(io, "  initial trajectory: ", summary.initial_trajectory)
    println(io, "  final trajectory: ", summary.final_trajectory)
    println(io, "  thermalization completed: ", summary.thermalization_completed)
    println(io, "  production completed: ", summary.production_completed)
    println(io, "  saved configurations: ", summary.saved_configurations)
    println(io, "  saved restart checkpoints: ", summary.saved_checkpoints)
    print(io, "  stopped: ", summary.stopped)
end

function show_config(io::IO, schedule::SimulationSchedule)
    println(io, "SimulationSchedule")
    println(io, "  initial trajectory: ", schedule.initial_trajectory)
    println(io, "  thermalization steps: ", schedule.thermalization_steps)
    println(io, "  production steps: ", schedule.production_steps)
    println(
        io,
        "  direct measurements: ",
        length(schedule.measurements.direct.measurements),
    )
    print(io, "  gradient flow: ")
    println(
        io,
        schedule.measurements.gradient_flow isa
        NoGradientFlowMeasurementConfig ? "disabled" : "enabled",
    )
    return nothing
end

function show_config(io::IO, output::OutputConfig)
    println(io, "OutputConfig")
    if output.configurations isa NoConfigurationOutput
        println(io, "  configurations: disabled")
    else
        config = output.configurations
        format = configuration_output_format(config)
        description = format === :jld2 ? "portable JLD2" : String(format)
        println(io, "  configurations: ", description)
        println(io, "  directory: ", repr(config.directory))
        println(io, "  every: ", config.every)
        println(io, "  prefix: ", repr(config.prefix))
    end
    if output.checkpoints isa NoCheckpointOutput
        println(io, "  restart checkpoints: disabled")
    else
        checkpoint = output.checkpoints
        println(io, "  restart checkpoints: portable JLD2")
        println(io, "  checkpoint directory: ", repr(checkpoint.directory))
        println(io, "  checkpoint every: ", checkpoint.every)
        println(io, "  checkpoint prefix: ", repr(checkpoint.prefix))
    end
    return nothing
end

function show_config(io::IO, spec::SimulationSpec)
    println(io, "SimulationSpec")
    show_config(io, spec.config)
    show_config(io, spec.schedule)
    show_config(io, spec.output)
    return nothing
end

function Base.show(io::IO, schedule::SimulationSchedule)
    print(
        io,
        "SimulationSchedule(thermalization=",
        schedule.thermalization_steps,
        ", production=",
        schedule.production_steps,
        ", initial=",
        schedule.initial_trajectory,
        ")",
    )
end

Base.show(io::IO, ::MIME"text/plain", schedule::SimulationSchedule) =
    show_config(io, schedule)

function Base.show(io::IO, output::OutputConfig)
    print(io, "OutputConfig(configurations=")
    print(io, nameof(typeof(output.configurations)))
    print(io, ", checkpoints=")
    print(io, nameof(typeof(output.checkpoints)))
    print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", output::OutputConfig) =
    show_config(io, output)

function Base.show(io::IO, spec::SimulationSpec)
    print(io, "SimulationSpec(config=")
    print(io, nameof(typeof(spec.config)))
    print(io, ", schedule=")
    show(io, spec.schedule)
    print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", spec::SimulationSpec) =
    show_config(io, spec)

export SimulationSchedule,
    NoConfigurationOutput,
    JLD2ConfigurationOutput,
    BridgeTextConfigurationOutput,
    ILDGConfigurationOutput,
    NoCheckpointOutput,
    JLD2CheckpointOutput,
    OutputConfig,
    SimulationSpec,
    ValidationIssue,
    validate,
    simulation_schedule,
    output_config,
    AbstractSimulationEvent,
    RunStarted,
    ThermalizationStepFinished,
    TrajectoryFinished,
    MeasurementFinished,
    ConfigurationSaved,
    CheckpointSaved,
    CheckpointLoaded,
    ConfigurationLoaded,
    SimulationRunSummary,
    RunStopped,
    RunFinished,
    RunFailed,
    NoSimulationEventSink,
    ConsoleSimulationEventSink,
    FunctionSimulationEventSink,
    RecordingSimulationEventSink,
    CompositeSimulationEventSink,
    SimulationSessionState,
    SimulationSession,
    SessionStepResult,
    build_simulation_session,
    emit_event!,
    configuration_output_path,
    checkpoint_output_path,
    save_checkpoint,
    load_checkpoint!,
    step!,
    is_finished,
    is_running,
    request_stop!,
    run_summary,
    run!

end
