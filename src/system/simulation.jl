module Simulation_module

import Random
import Gaugefields
import Gaugefields: load_configuration!, save_configuration
import LatticeDiracOperators

import ..LQCDCommunication: broadcast!, comm_rank, is_distributed

import ..LQCDConfig_module:
    AbstractDiracOperatorConfig,
    AbstractUpdateConfig,
    ConfigurationSequenceConfig,
    ColdStartConfig,
    DomainwallDiracConfig,
    DirectorySourceConfig,
    EmbeddedInstantonConfig,
    FermionActionConfig,
    FermionSolverConfig,
    FileStartConfig,
    GaugeActionConfig,
    HeatbathConfig,
    HMCConfig,
    HotStartConfig,
    HISQDiracConfig,
    InstantonConfig,
    LQCDConfig,
    LeapfrogConfig,
    ManifestSourceConfig,
    MobiusDomainwallDiracConfig,
    NoFermionSmearingConfig,
    PQPConfig,
    PseudofermionRefreshConfig,
    QPQConfig,
    SLHMCConfig,
    SextonWeingartenConfig,
    StaggeredDiracConfig,
    StoutFermionSmearingConfig,
    WilsonCloverDiracConfig,
    WilsonDiracConfig,
    canonical_configuration_format,
    fermion_action_names,
    md_trajectory_length,
    momentum_denominator

"""Common supertype for runtime field configurations."""
abstract type AbstractConfiguration end

"""A runtime configuration containing gauge fields only."""
struct GaugeConfiguration{G} <: AbstractConfiguration
    gauge::G
end

"""A runtime configuration containing gauge and named fermion fields."""
struct GaugeFermionConfiguration{G,F<:NamedTuple} <: AbstractConfiguration
    gauge::G
    fermions::F
end

gauge_links(configuration::GaugeConfiguration) = configuration.gauge
gauge_links(configuration::GaugeFermionConfiguration) = configuration.gauge

"""
    save_configuration(filename, configuration; format=:jld2, kwargs...)

Save the gauge links contained in a runtime configuration.  Dynamical
pseudofermion fields are trajectory workspaces and are deliberately not part
of a configuration checkpoint.  Portable JLD2 output is gathered and written
by rank 0 through the Gaugefields v1 API; every rank owning a distributed
configuration must call this function.
"""
function save_configuration(
    filename::AbstractString,
    configuration::AbstractConfiguration;
    format::Union{Symbol,AbstractString}=:jld2,
    kwargs...,
)
    canonical_format = canonical_configuration_format(format)
    Gaugefields.save_configuration(
        filename,
        gauge_links(configuration);
        format=canonical_format,
        kwargs...,
    )
    return filename
end

"""
    load_configuration!(configuration, filename; format=:jld2)

Replace the gauge links in an existing runtime configuration.  Fermion fields
remain allocated on the current backend and are refreshed by the next HMC
trajectory.
"""
function load_configuration!(
    configuration::AbstractConfiguration,
    filename::AbstractString;
    format::Union{Symbol,AbstractString}=:jld2,
)
    canonical_format = canonical_configuration_format(format)
    Gaugefields.load_configuration!(
        gauge_links(configuration),
        filename;
        format=canonical_format,
    )
    return configuration
end

"""
Gaugefields v1 execution choices.

Backend initialization remains under the control of the application or
notebook. Pass `Gaugefields.SerialCommunicator()` for an MPI-free single-rank
run, or an initialized MPI communicator for a distributed run. The default
`nothing` retains Gaugefields' automatic communicator selection.
"""
struct GaugefieldsEnvironment{B,P,C,E}
    backend::B
    process_grid::P
    communicator::C
    element_type::Type{E}
    verbose::Int
end

function GaugefieldsEnvironment(;
    backend=Gaugefields.LatticeMatricesBackend(),
    process_grid=nothing,
    communicator=nothing,
    element_type::Type{E}=ComplexF64,
    verbose::Integer=0,
) where {E}
    return GaugefieldsEnvironment{
        typeof(backend),
        typeof(process_grid),
        typeof(communicator),
        E,
    }(
        backend,
        process_grid,
        communicator,
        element_type,
        Int(verbose),
    )
end

function unsigned_seed(seed::Integer)
    seed >= 0 || throw(ArgumentError("random seeds must be nonnegative; got $seed"))
    return UInt64(seed)
end

function random_stream_seed(random_config)
    seed = unsigned_seed(random_config.seed)
    if random_config.name == :momentum
        return seed ⊻ UInt64(0x243f6a8885a308d3)
    elseif random_config.name == :metropolis
        return seed ⊻ UInt64(0x13198a2e03707344)
    elseif random_config.name == :heatbath
        return seed ⊻ UInt64(0xa4093822299f31d0)
    elseif random_config.name == :pseudofermion
        return seed ⊻ UInt64(0x082efa98ec4e6c89)
    end
    throw(ArgumentError(
        "unsupported random stream name $(repr(random_config.name))",
    ))
end

configuration_format(format::Union{Symbol,AbstractString}) =
    canonical_configuration_format(format)

function configuration_path_matches(path::AbstractString, format::Symbol)
    name = lowercase(basename(path))
    format === :jld2 && return endswith(name, ".jld2")
    format === :ildg && return endswith(name, ".ildg")
    format === :bridge && return endswith(name, ".txt")
    throw(ArgumentError("unsupported configuration format $(repr(format))"))
end

function configuration_paths(
    source::DirectorySourceConfig,
    format::Symbol,
)
    directory = abspath(source.directory)
    isdir(directory) || throw(ArgumentError(
        "configuration directory does not exist: $directory",
    ))
    paths = sort(filter(
        path -> isfile(path) && configuration_path_matches(path, format),
        readdir(directory; join=true),
    ))
    isempty(paths) && throw(ArgumentError(
        "configuration directory contains no $format files: $directory",
    ))
    return paths
end

function configuration_paths(
    source::ManifestSourceConfig,
    format::Symbol,
)
    directory = abspath(source.directory)
    isdir(directory) || throw(ArgumentError(
        "configuration directory does not exist: $directory",
    ))
    manifest = isabspath(source.manifest) ?
               normpath(source.manifest) :
               normpath(joinpath(directory, source.manifest))
    isfile(manifest) || throw(ArgumentError(
        "configuration manifest does not exist: $manifest",
    ))

    paths = String[]
    for (line_number, line) in enumerate(eachline(manifest))
        entry = strip(first(split(line, '#'; limit=2)))
        isempty(entry) && continue
        path = if isabspath(entry)
            normpath(entry)
        else
            normpath(joinpath(directory, entry))
        end
        configuration_path_matches(path, format) || throw(ArgumentError(
            "configuration manifest entry $line_number has the wrong " *
            "extension for $format: $entry",
        ))
        isfile(path) || throw(ArgumentError(
            "configuration manifest entry $line_number does not exist: $path",
        ))
        push!(paths, path)
    end
    isempty(paths) && throw(ArgumentError(
        "configuration manifest contains no $format files: $manifest",
    ))
    return paths
end

configuration_paths(config::ConfigurationSequenceConfig) =
    configuration_paths(config.source, config.format)

function validate_one_instanton(config::LQCDConfig)
    dimensions = length(config.lattice.L)
    dimensions in (2, 4) || throw(ArgumentError(
        "one-instanton initialization supports two- and four-dimensional " *
        "lattices; got $dimensions dimensions",
    ))
    config.gauge.NC == 2 || throw(ArgumentError(
        "legacy one-instanton initialization requires NC=2; " *
        "got NC=$(config.gauge.NC)",
    ))
    return nothing
end

function validate_embedded_instanton(
    config::LQCDConfig,
    initialization::EmbeddedInstantonConfig,
)
    dimensions = length(config.lattice.L)
    dimensions == 4 || throw(ArgumentError(
        "embedded-instanton initialization requires a four-dimensional " *
        "lattice; got $dimensions dimensions",
    ))
    config.gauge.NC >= 2 || throw(ArgumentError(
        "embedded-instanton initialization requires NC >= 2; " *
        "got NC=$(config.gauge.NC)",
    ))
    first_color, second_color = initialization.block
    first_color != second_color || throw(ArgumentError(
        "embedded-instanton color indices must be distinct",
    ))
    for color in initialization.block
        1 <= color <= config.gauge.NC || throw(ArgumentError(
            "embedded-instanton color index $color is outside " *
            "1:$(config.gauge.NC)",
        ))
    end
    return nothing
end

function build_one_instanton(
    config::LQCDConfig,
    environment::GaugefieldsEnvironment{B},
) where {B<:Gaugefields.LatticeMatricesBackend}
    validate_one_instanton(config)
    environment.element_type in (ComplexF32, ComplexF64) || throw(ArgumentError(
        "one-instanton initialization on the LatticeMatrices backend " *
        "requires ComplexF32 or ComplexF64 fields; " *
        "got $(environment.element_type)",
    ))
    return Gaugefields.Oneinstanton(
        config.gauge.NC,
        config.gauge.halo,
        config.lattice.L...;
        isMPILattice=true,
        PEs=environment.process_grid,
        elementtype=environment.element_type,
        verbose_level=environment.verbose,
    )
end

function build_one_instanton(
    config::LQCDConfig,
    environment::GaugefieldsEnvironment{B},
) where {B<:Gaugefields.LegacyBackend}
    validate_one_instanton(config)
    environment.process_grid === nothing || throw(ArgumentError(
        "LegacyBackend does not accept a process grid",
    ))
    environment.element_type === ComplexF64 || throw(ArgumentError(
        "LegacyBackend one-instanton initialization requires ComplexF64",
    ))
    return Gaugefields.Oneinstanton(
        config.gauge.NC,
        config.gauge.halo,
        config.lattice.L...;
        verbose_level=environment.verbose,
    )
end

function build_one_instanton(
    ::LQCDConfig,
    environment::GaugefieldsEnvironment,
)
    throw(ArgumentError(
        "one-instanton initialization does not support backend " *
        "$(typeof(environment.backend))",
    ))
end

function embedded_instanton_parameters(
    config::LQCDConfig,
    initialization::EmbeddedInstantonConfig,
)
    lattice = config.lattice.L
    center = initialization.center === nothing ?
             ntuple(index -> lattice[index] / 2 + 0.5, 4) :
             initialization.center
    radius = initialization.radius === nothing ?
             div(lattice[1], 2) :
             initialization.radius
    return (; center, radius)
end

function build_embedded_instanton(
    config::LQCDConfig,
    environment::GaugefieldsEnvironment{B},
    initialization::EmbeddedInstantonConfig,
) where {B<:Gaugefields.LatticeMatricesBackend}
    validate_embedded_instanton(config, initialization)
    environment.element_type in (ComplexF32, ComplexF64) || throw(ArgumentError(
        "embedded-instanton initialization on the LatticeMatrices backend " *
        "requires ComplexF32 or ComplexF64 fields; " *
        "got $(environment.element_type)",
    ))
    parameters = embedded_instanton_parameters(config, initialization)
    return Gaugefields.Oneinstanton_SUN_embedded(
        config.gauge.NC,
        config.lattice.L...;
        NDW=config.gauge.halo,
        center=parameters.center,
        radius=parameters.radius,
        sign=initialization.sign,
        block=initialization.block,
        isMPILattice=true,
        PEs=environment.process_grid,
        elementtype=environment.element_type,
        verbose_level=environment.verbose,
    )
end

function build_embedded_instanton(
    config::LQCDConfig,
    environment::GaugefieldsEnvironment{B},
    initialization::EmbeddedInstantonConfig,
) where {B<:Gaugefields.LegacyBackend}
    validate_embedded_instanton(config, initialization)
    environment.process_grid === nothing || throw(ArgumentError(
        "LegacyBackend does not accept a process grid",
    ))
    environment.element_type === ComplexF64 || throw(ArgumentError(
        "LegacyBackend embedded-instanton initialization requires ComplexF64",
    ))
    parameters = embedded_instanton_parameters(config, initialization)
    return Gaugefields.Oneinstanton_SUN_embedded(
        config.gauge.NC,
        config.lattice.L...;
        NDW=config.gauge.halo,
        center=parameters.center,
        radius=parameters.radius,
        sign=initialization.sign,
        block=initialization.block,
        verbose_level=environment.verbose,
    )
end

function build_embedded_instanton(
    ::LQCDConfig,
    environment::GaugefieldsEnvironment,
    ::EmbeddedInstantonConfig,
)
    throw(ArgumentError(
        "embedded-instanton initialization does not support backend " *
        "$(typeof(environment.backend))",
    ))
end

function build_standard_configuration(
    config::LQCDConfig,
    environment::GaugefieldsEnvironment,
    start::Symbol,
    seed,
)
    return Gaugefields.gauge_configuration(
        config.lattice.L;
        backend=environment.backend,
        colors=config.gauge.NC,
        halo=config.gauge.halo,
        start,
        seed,
        process_grid=environment.process_grid,
        comm=environment.communicator,
        boundary=:periodic,
        eltype=environment.element_type,
        verbose=environment.verbose,
    )
end

function build_initial_gauge(
    config::LQCDConfig,
    environment::GaugefieldsEnvironment,
    ::ColdStartConfig,
)
    return build_standard_configuration(config, environment, :cold, nothing)
end

function build_initial_gauge(
    config::LQCDConfig,
    environment::GaugefieldsEnvironment,
    initialization::HotStartConfig,
)
    seed = initialization.seed === nothing ?
           nothing : unsigned_seed(initialization.seed)
    return build_standard_configuration(config, environment, :hot, seed)
end

function build_initial_gauge(
    config::LQCDConfig,
    environment::GaugefieldsEnvironment,
    initialization::FileStartConfig,
)
    initialization.format === nothing && throw(ArgumentError(
        "load format is required when initialization is a configuration file",
    ))
    gauge = build_standard_configuration(config, environment, :cold, nothing)
    Gaugefields.load_configuration!(
        gauge,
        initialization.path;
        format=configuration_format(String(initialization.format)),
    )
    return gauge
end

function build_initial_gauge(
    config::LQCDConfig,
    environment::GaugefieldsEnvironment,
    ::InstantonConfig,
)
    return build_one_instanton(config, environment)
end

function build_initial_gauge(
    config::LQCDConfig,
    environment::GaugefieldsEnvironment,
    initialization::EmbeddedInstantonConfig,
)
    return build_embedded_instanton(config, environment, initialization)
end

fermion_field_family(::WilsonDiracConfig) = "Wilson"
fermion_field_family(::WilsonCloverDiracConfig) = "Wilson"
fermion_field_family(::StaggeredDiracConfig) = "staggered"
fermion_field_family(::HISQDiracConfig) = "staggered"
fermion_field_family(::DomainwallDiracConfig) = "Domainwall"
fermion_field_family(::MobiusDomainwallDiracConfig) = "MobiusDomainwall"

fermion_fifth_dimension(::AbstractDiracOperatorConfig) = 2
fermion_fifth_dimension(config::DomainwallDiracConfig) =
    config.fifth_dimension
fermion_fifth_dimension(config::MobiusDomainwallDiracConfig) =
    config.fifth_dimension

fermion_field_boundary_conditions(
    ::AbstractDiracOperatorConfig,
    boundary_conditions::NTuple{4,<:Integer},
) = boundary_conditions

function fermion_field_boundary_conditions(
    ::Union{DomainwallDiracConfig,MobiusDomainwallDiracConfig},
    boundary_conditions::NTuple{4,<:Integer},
)
    # LatticeMatrices stores domain-wall fields as a five-dimensional lattice.
    # The physical four-dimensional phases come from the solver configuration;
    # the undistributed fifth direction is periodic.
    return (boundary_conditions..., 1)
end

function build_fermion_field(config::FermionActionConfig, gauge)
    boundary_conditions = config.solver.boundary_conditions
    length(boundary_conditions) == 4 || throw(ArgumentError(
        "dynamical fermions require four boundary-condition phases; " *
        "got $boundary_conditions",
    ))
    return LatticeDiracOperators.Initialize_pseudofermion_fields(
        gauge[1],
        fermion_field_family(config.operator);
        L5=fermion_fifth_dimension(config.operator),
        nowing=true,
        boundarycondition=fermion_field_boundary_conditions(
            config.operator,
            boundary_conditions,
        ),
    )
end

function build_fermion_fields(fermions::Tuple, gauge)
    names = fermion_action_names(fermions)
    length(unique(names)) == length(names) || throw(ArgumentError(
        "fermion action names must be unique; got $names",
    ))
    fields = map(config -> build_fermion_field(config, gauge), fermions)
    return NamedTuple{names}(fields)
end

"""Construct gauge fields exclusively through the Gaugefields v1 API."""
function build_configuration(
    config::LQCDConfig,
    environment::GaugefieldsEnvironment=GaugefieldsEnvironment(),
)
    gauge = build_initial_gauge(
        config,
        environment,
        config.gauge.initialization,
    )
    configuration = if isempty(config.fermions)
        GaugeConfiguration(gauge)
    else
        GaugeFermionConfiguration(
            gauge,
            build_fermion_fields(config.fermions, gauge),
        )
    end
    initialize_configuration!(configuration, config.update)
    return configuration
end

function initialize_configuration!(
    configuration::AbstractConfiguration,
    ::AbstractUpdateConfig,
)
    return configuration
end

function initialize_configuration!(
    configuration::AbstractConfiguration,
    config::ConfigurationSequenceConfig,
)
    paths = configuration_paths(config)
    Gaugefields.load_configuration!(
        gauge_links(configuration),
        first(paths);
        format=config.format,
    )
    return configuration
end

function action_loops(loop, dimensions::Int)
    if loop isa AbstractString || loop isa Symbol
        return Gaugefields.make_loops_fromname(String(loop); Dim=dimensions)
    elseif loop isa AbstractVector
        return collect(loop)
    end
    return [loop]
end

"""Construct a Gaugefields v1 `GaugeAction` from typed action terms."""
function build_gauge_action(
    config::GaugeActionConfig,
    configuration::AbstractConfiguration,
)
    gauge = gauge_links(configuration)
    action = Gaugefields.GaugeAction(gauge)
    for term in config.terms
        loops = action_loops(term.loop, length(gauge))
        coefficient = term.coupling
        if term.include_adjoint
            append!(loops, loops')
            coefficient /= 2
        end
        push!(action, coefficient, loops)
    end
    return action
end

"""Copy gauge links without moving them to the host."""
function copy_configuration!(
    destination::GaugeConfiguration,
    source::GaugeConfiguration,
)
    length(destination.gauge) == length(source.gauge) || throw(ArgumentError(
        "source and destination must have the same number of directions",
    ))
    for direction in eachindex(destination.gauge, source.gauge)
        Gaugefields.substitute_U!(
            destination.gauge[direction],
            source.gauge[direction],
        )
    end
    return destination
end


function copy_configuration!(
    destination::GaugeFermionConfiguration,
    source::GaugeFermionConfiguration,
)
    keys(destination.fermions) == keys(source.fermions) || throw(ArgumentError(
        "source and destination must contain the same fermion fields",
    ))
    copy_gauge_links!(destination.gauge, source.gauge)
    for name in keys(destination.fermions)
        LatticeDiracOperators.substitute_fermion!(
            getproperty(destination.fermions, name),
            getproperty(source.fermions, name),
        )
    end
    return destination
end

function copy_gauge_links!(destination, source)
    length(destination) == length(source) || throw(ArgumentError(
        "source and destination must have the same number of directions",
    ))
    for direction in eachindex(destination, source)
        Gaugefields.substitute_U!(destination[direction], source[direction])
    end
    return destination
end

function copy_configuration!(
    destination::GaugeConfiguration,
    source::GaugeFermionConfiguration,
)
    copy_gauge_links!(destination.gauge, source.gauge)
    return destination
end

function copy_configuration!(
    destination::GaugeFermionConfiguration,
    source::GaugeConfiguration,
)
    copy_gauge_links!(destination.gauge, source.gauge)
    return destination
end

function common_dirac_parameters(config::FermionActionConfig)
    solver = config.solver
    return Dict{String,Any}(
        "eps_CG" => solver.tolerance,
        "MaxCGstep" => solver.max_steps,
        "verbose_level" => solver.verbose,
        "boundarycondition" => collect(solver.boundary_conditions),
    )
end

function dirac_parameters(config::FermionActionConfig{O}) where {
    O<:WilsonDiracConfig,
}
    parameters = common_dirac_parameters(config)
    parameters["Dirac_operator"] = "Wilson"
    parameters["κ"] = config.operator.hopping_parameter
    parameters["r"] = config.operator.wilson_parameter
    return parameters
end

function dirac_parameters(config::FermionActionConfig{O}) where {
    O<:WilsonCloverDiracConfig,
}
    parameters = common_dirac_parameters(config)
    parameters["Dirac_operator"] = "WilsonClover"
    parameters["κ"] = config.operator.hopping_parameter
    parameters["r"] = config.operator.wilson_parameter
    parameters["cSW"] = config.operator.clover_coefficient
    return parameters
end

function dirac_parameters(config::FermionActionConfig{O}) where {
    O<:StaggeredDiracConfig,
}
    parameters = common_dirac_parameters(config)
    parameters["Dirac_operator"] = "staggered"
    parameters["mass"] = config.operator.mass
    return parameters
end

function dirac_parameters(config::FermionActionConfig{O}) where {
    O<:HISQDiracConfig,
}
    parameters = common_dirac_parameters(config)
    parameters["Dirac_operator"] = "HISQ"
    parameters["mass"] = config.operator.mass
    parameters["naik_epsilon"] = config.operator.naik_epsilon
    return parameters
end

function dirac_parameters(config::FermionActionConfig{O}) where {
    O<:DomainwallDiracConfig,
}
    parameters = common_dirac_parameters(config)
    parameters["Dirac_operator"] = "Domainwall"
    parameters["mass"] = config.operator.mass
    parameters["M"] = config.operator.domainwall_height
    parameters["L5"] = config.operator.fifth_dimension
    return parameters
end


function dirac_parameters(config::FermionActionConfig{O}) where {
    O<:MobiusDomainwallDiracConfig,
}
    parameters = common_dirac_parameters(config)
    parameters["Dirac_operator"] = "MobiusDomainwall"
    parameters["mass"] = config.operator.mass
    parameters["M"] = config.operator.domainwall_height
    parameters["L5"] = config.operator.fifth_dimension
    parameters["b"] = config.operator.b
    parameters["c"] = config.operator.c
    return parameters
end

build_fermion_smearing(::NoFermionSmearingConfig, gauge) = nothing

function build_fermion_smearing(
    config::StoutFermionSmearingConfig,
    gauge,
)
    dimension = length(gauge)
    numtemps = max(
        4 + config.layers * dimension,
        3 * dimension,
    )
    smearing = Gaugefields.CovNeuralnet(first(gauge); numtemps)
    loops = collect(config.loops)
    coefficients = collect(config.coefficients)
    for _ in 1:config.layers
        push!(
            smearing,
            Gaugefields.STOUT_Layer(loops, coefficients, gauge),
        )
    end
    return smearing
end

function _build_fermion_action_and_smearing(
    config::FermionActionConfig,
    gauge,
    pseudofermion,
)
    operator = LatticeDiracOperators.Dirac_operator(
        gauge,
        pseudofermion,
        dirac_parameters(config),
    )
    smearing = build_fermion_smearing(config.smearing, gauge)
    action_parameters = Dict{String,Any}("Nf" => config.flavors)
    action = LatticeDiracOperators.FermiAction(
        operator,
        action_parameters;
        covneuralnet=smearing,
    )
    return action, smearing
end

function build_fermion_action(
    config::FermionActionConfig,
    gauge,
    pseudofermion,
)
    action, _ = _build_fermion_action_and_smearing(
        config,
        gauge,
        pseudofermion,
    )
    return action
end

"""Concrete resources needed to refresh one pseudofermion action term."""
struct PseudofermionRefreshRuntime{P,N,S<:Integer}
    provider::P
    noise::N
    seed::S
    subgroup::Int
end

function pseudofermion_refresh_config(
    refreshes::Tuple,
    action_name::Symbol,
)
    matches = filter(config -> config.action_name === action_name, refreshes)
    length(matches) == 1 || throw(ArgumentError(
        "expected one pseudofermion refresh for $action_name; " *
        "found $(length(matches))",
    ))
    return only(matches)
end

function build_fermion_provider(
    config::FermionActionConfig,
    configuration::GaugeFermionConfiguration,
)
    gauge = configuration.gauge
    pseudofermion = getproperty(configuration.fermions, config.name)
    action, smearing = _build_fermion_action_and_smearing(
        config,
        gauge,
        pseudofermion,
    )
    return LatticeDiracOperators.PseudofermionMDAction(
        action,
        pseudofermion,
        smearing,
    )
end

function build_fermion_runtime(
    config::FermionActionConfig,
    refreshes::Tuple,
    configuration::GaugeFermionConfiguration,
)
    pseudofermion = getproperty(configuration.fermions, config.name)
    provider = build_fermion_provider(config, configuration)
    refresh_config = pseudofermion_refresh_config(refreshes, config.name)
    runtime = PseudofermionRefreshRuntime(
        provider,
        similar(pseudofermion),
        random_stream_seed(refresh_config.random),
        refresh_config.subgroup,
    )
    return provider, runtime
end

function named_md_action_set(gauge_action, fermions::Tuple, providers::Tuple)
    fermion_names = fermion_action_names(fermions)
    names = (:gauge, fermion_names...)
    return Gaugefields.MDActionSet(
        NamedTuple{names}((gauge_action, providers...)),
    )
end

function build_md_action_set(
    gauge_action,
    fermions::Tuple,
    configuration::GaugeFermionConfiguration,
)
    providers = map(fermions) do config
        build_fermion_provider(config, configuration)
    end
    return named_md_action_set(gauge_action, fermions, providers)
end

build_md_action_set(gauge_action, ::Tuple{}, ::GaugeConfiguration) =
    gauge_action

function build_md_actions(
    gauge_action,
    fermions::Tuple,
    refreshes::Tuple,
    configuration::GaugeFermionConfiguration,
)
    fermion_names = fermion_action_names(fermions)
    refresh_names = map(config -> config.action_name, refreshes)
    length(unique(refresh_names)) == length(refresh_names) ||
        throw(ArgumentError(
            "pseudofermion refresh names must be unique; got $refresh_names",
        ))
    Set(refresh_names) == Set(fermion_names) || throw(ArgumentError(
        "pseudofermion refresh names $refresh_names do not match fermion " *
        "actions $fermion_names",
    ))
    runtimes = map(fermions) do config
        build_fermion_runtime(config, refreshes, configuration)
    end
    providers = map(first, runtimes)
    refresh_runtimes = map(last, runtimes)
    action_set = named_md_action_set(gauge_action, fermions, providers)
    return action_set, refresh_runtimes
end

function build_md_actions(
    gauge_action,
    ::Tuple{},
    refreshes::Tuple,
    ::GaugeConfiguration,
)
    isempty(refreshes) || throw(ArgumentError(
        "a gauge-only simulation must not configure pseudofermion refreshes",
    ))
    return gauge_action, ()
end

"""Preallocated runtime objects for one gauge or pseudofermion HMC updater."""
mutable struct HMCUpdater{
    D,
    P,
    B<:GaugeConfiguration,
    T<:Real,
    S<:Integer,
    R<:Tuple,
}
    md_driver::D
    momentum::P
    backup::B
    momentum_sigma::T
    momentum_seed::S
    pseudofermion_refreshes::R
end

"""Target-Hamiltonian and proposal-Hamiltonian drivers for one SLHMC chain."""
mutable struct SLHMCUpdater{
    D,
    E,
    P,
    B<:GaugeConfiguration,
    T<:Real,
    S<:Integer,
    R<:Tuple,
    Q<:Tuple,
}
    md_driver::D
    target_driver::E
    momentum::P
    backup::B
    momentum_sigma::T
    momentum_seed::S
    pseudofermion_refreshes::R
    md_trajectory_state_providers::Q
end

"""Mutable, checkpointable state of an HMC chain."""
mutable struct HMCState{R}
    trajectory::Int
    accepted::Int
    metropolis_rng::R
end

"""Diagnostics and decision from one HMC trajectory."""
struct HMCUpdateResult{T<:Real}
    trajectory::Int
    accepted::Bool
    initial_hamiltonian::T
    final_hamiltonian::T
    delta_hamiltonian::T
end


"""
Exact target diagnostics plus the MD-action diagnostics for one SLHMC proposal.
"""
struct SLHMCUpdateResult{T<:Real,M<:Real}
    trajectory::Int
    accepted::Bool
    initial_hamiltonian::T
    final_hamiltonian::T
    delta_hamiltonian::T
    md_initial_hamiltonian::M
    md_final_hamiltonian::M
    md_delta_hamiltonian::M
end

"""Gaugefields heatbath storage and the number of OR sweeps per update."""
mutable struct HeatbathUpdater{K}
    kernel::K
    overrelaxation_steps::Int
end

"""Portable counters needed to resume a heatbath Markov chain."""
mutable struct HeatbathState
    trajectory::Int
    heatbath_sweep::UInt64
    overrelaxation_sweep::UInt64
end

"""Counters after one completed heatbath plus optional OR update."""
struct HeatbathUpdateResult
    trajectory::Int
    heatbath_sweep::UInt64
    overrelaxation_sweep::UInt64
end

"""Resolved configuration paths and their Gaugefields v1 input format."""
struct ConfigurationSequenceUpdater{P<:AbstractVector,F<:Symbol}
    paths::P
    format::F
end

"""Current position in a finite sequence of configurations."""
mutable struct ConfigurationSequenceState
    trajectory::Int
    current_index::Int
end

"""The configuration selected by one successful sequence update."""
struct ConfigurationSequenceUpdateResult{P<:AbstractString}
    trajectory::Int
    current_index::Int
    path::P
end

"""Runtime simulation assembled from typed input and concrete resources."""
struct Simulation{C<:AbstractConfiguration,A,U,S,E,I}
    configuration::C
    action::A
    updater::U
    state::S
    environment::E
    input::I
end

force_group_names(group) = Tuple(Symbol.(group.names))

function validate_force_partition(configured::Tuple, available::Tuple)
    length(unique(configured)) == length(configured) || throw(ArgumentError(
        "force names must not be duplicated; got $configured",
    ))
    Set(configured) == Set(available) || throw(ArgumentError(
        "configured forces $configured do not match available actions " *
        "$available",
    ))
    return nothing
end

function gaugefields_integrator(config::LeapfrogConfig, available::Tuple)
    validate_force_partition(force_group_names(config.forces), available)
    return config.ordering isa QPQConfig ? Gaugefields.QPQ() : Gaugefields.PQP()
end

function gaugefields_integrator(
    config::SextonWeingartenConfig,
    available::Tuple,
)
    config.fast_integrator isa LeapfrogConfig || throw(ArgumentError(
        "the Gaugefields Sexton-Weingarten driver currently supports one " *
        "nested leapfrog fast integrator",
    ))
    config.fast_integrator.ordering isa QPQConfig || throw(ArgumentError(
        "the Sexton-Weingarten fast integrator must use QPQ ordering",
    ))
    slow = force_group_names(config.slow_forces)
    fast = force_group_names(config.fast_integrator.forces)
    validate_force_partition((slow..., fast...), available)
    ordering = config.ordering isa QPQConfig ?
               Gaugefields.QPQ() : Gaugefields.PQP()
    return Gaugefields.SextonWeingarten(
        slow=slow,
        fast=fast,
        n_fast=config.fast_steps,
        ordering=ordering,
    )
end

function build_md_driver(md_config, momentum_config, gauge, action)
    available = action isa Gaugefields.MDActionSet ?
                Tuple(keys(action.terms)) : (:gauge,)
    return Gaugefields.md_driver(
        gauge,
        action;
        steps=md_config.steps,
        trajectory_length=md_trajectory_length(md_config),
        integrator=gaugefields_integrator(md_config.integrator, available),
        momentum_denominator=momentum_denominator(momentum_config),
    )
end

function build_hmc_updater(
    config::HMCConfig,
    configuration::AbstractConfiguration,
    action,
    refreshes::Tuple,
)
    gauge = gauge_links(configuration)
    driver = build_md_driver(config.md, config.momentum, gauge, action)
    momentum = Gaugefields.gauge_momenta(gauge)
    backup = GaugeConfiguration([similar(link) for link in gauge])
    return HMCUpdater(
        driver,
        momentum,
        backup,
        config.momentum.sigma,
        random_stream_seed(config.momentum.random),
        refreshes,
    )
end

function validate_slhmc_fermions(target::Tuple, md::Tuple)
    target_names = fermion_action_names(target)
    md_names = fermion_action_names(md)
    target_names == md_names || throw(ArgumentError(
        "SLHMC target fermion names $target_names must match MD-action " *
        "fermion names $md_names in the same order",
    ))
    for (target_term, md_term) in zip(target, md)
        target_family = fermion_field_family(target_term.operator)
        md_family = fermion_field_family(md_term.operator)
        target_family == md_family || throw(ArgumentError(
            "SLHMC action $(target_term.name) uses incompatible target " *
            "and MD fermion-field families $target_family and $md_family",
        ))
        target_l5 = fermion_fifth_dimension(target_term.operator)
        md_l5 = fermion_fifth_dimension(md_term.operator)
        target_l5 == md_l5 || throw(ArgumentError(
            "SLHMC action $(target_term.name) uses target L5=$target_l5 " *
            "but MD-action L5=$md_l5",
        ))
    end
    return nothing
end

function build_slhmc_updater(
    config::SLHMCConfig,
    configuration::AbstractConfiguration,
    target_action,
    md_action,
    refreshes::Tuple,
)
    gauge = gauge_links(configuration)
    md_driver = build_md_driver(config.md, config.momentum, gauge, md_action)
    target_driver = build_md_driver(
        config.md,
        config.momentum,
        gauge,
        target_action,
    )
    momentum = Gaugefields.gauge_momenta(gauge)
    backup = GaugeConfiguration([similar(link) for link in gauge])
    md_trajectory_state_providers = md_action isa Gaugefields.MDActionSet ?
        Tuple(
            provider for (name, provider) in pairs(md_action.terms)
            if name !== :gauge
        ) : ()
    return SLHMCUpdater(
        md_driver,
        target_driver,
        momentum,
        backup,
        config.momentum.sigma,
        random_stream_seed(config.momentum.random),
        refreshes,
        md_trajectory_state_providers,
    )
end

function plaquette_heatbath_beta(config::GaugeActionConfig)
    length(config.terms) == 1 || throw(ArgumentError(
        "checkerboard heatbath requires exactly one plaquette action term; " *
        "got $(length(config.terms)) terms",
    ))
    term = only(config.terms)
    loop_name = if term.loop isa Union{Symbol,AbstractString}
        lowercase(strip(String(term.loop)))
    else
        nothing
    end
    loop_name == "plaquette" || throw(ArgumentError(
        "checkerboard heatbath requires a plaquette loop; got " *
        repr(term.loop),
    ))
    term.include_adjoint || throw(ArgumentError(
        "checkerboard heatbath requires both plaquette orientations",
    ))
    term.coupling isa Real || throw(ArgumentError(
        "checkerboard heatbath requires a real beta; got " *
        repr(term.coupling),
    ))
    return term.coupling
end

function build_heatbath_updater(
    config::HeatbathConfig,
    input::LQCDConfig,
    configuration::GaugeConfiguration,
    action,
)
    gauge = configuration.gauge
    seed = random_stream_seed(config.random)
    kernel = if config.even_odd
        beta = plaquette_heatbath_beta(input.gauge_action)
        Gaugefields.heatbath_updater(
            gauge;
            beta,
            ITERATION_MAX=config.max_iterations,
            seed,
            sweep=0,
            overrelaxation_sweep=0,
        )
    else
        Gaugefields.heatbath_updater(
            gauge,
            action;
            ITERATION_MAX=config.max_iterations,
            seed,
            sweep=0,
            overrelaxation_sweep=0,
        )
    end
    return HeatbathUpdater(kernel, config.overrelaxation_steps)
end

function build_update_runtime(
    config::HMCConfig,
    ::LQCDConfig,
    configuration::AbstractConfiguration,
    action,
    refreshes::Tuple,
)
    updater = build_hmc_updater(config, configuration, action, refreshes)
    state = HMCState(
        0,
        0,
        Random.MersenneTwister(
            random_stream_seed(config.acceptance.random),
        ),
    )
    return updater, state
end

function build_update_runtime(
    config::SLHMCConfig,
    input::LQCDConfig,
    configuration::AbstractConfiguration,
    target_action,
    refreshes::Tuple,
)
    validate_slhmc_fermions(input.fermions, config.md_fermions)
    md_gauge_action = build_gauge_action(
        config.md_gauge_action,
        configuration,
    )
    md_action = build_md_action_set(
        md_gauge_action,
        config.md_fermions,
        configuration,
    )
    updater = build_slhmc_updater(
        config,
        configuration,
        target_action,
        md_action,
        refreshes,
    )
    state = HMCState(
        0,
        0,
        Random.MersenneTwister(
            random_stream_seed(config.acceptance.random),
        ),
    )
    return updater, state
end

function build_update_runtime(
    config::HeatbathConfig,
    input::LQCDConfig,
    configuration::GaugeConfiguration,
    action,
    ::Tuple{},
)
    updater = build_heatbath_updater(
        config,
        input,
        configuration,
        action,
    )
    state = HeatbathState(0, 0, 0)
    return updater, state
end

function build_update_runtime(
    config::ConfigurationSequenceConfig,
    ::LQCDConfig,
    ::GaugeConfiguration,
    ::Any,
    ::Tuple{},
)
    paths = configuration_paths(config)
    updater = ConfigurationSequenceUpdater(paths, config.format)
    state = ConfigurationSequenceState(0, 1)
    return updater, state
end

function validate_runtime_input(config::LQCDConfig)
    if !isempty(config.fermions) &&
       !(config.update isa Union{HMCConfig,SLHMCConfig})
        throw(ArgumentError(
            "dynamical fermion actions require an HMC or SLHMC update; got " *
            String(nameof(typeof(config.update))),
        ))
    end
    for fermion in config.fermions
        if fermion.operator isa HISQDiracConfig
            config.gauge.halo >= 3 || throw(ArgumentError(
                "HISQ requires a gauge-field halo width of at least 3; " *
                "got $(config.gauge.halo)",
            ))
        end
    end
    return nothing
end


pseudofermion_configs(config::HMCConfig) = config.pseudofermions
pseudofermion_configs(config::SLHMCConfig) = config.pseudofermions
pseudofermion_configs(::AbstractUpdateConfig) = ()

"""Build a gauge or dynamical-fermion runtime using the v1 package APIs."""
function build_simulation(
    config::LQCDConfig,
    environment::GaugefieldsEnvironment=GaugefieldsEnvironment(),
)
    validate_runtime_input(config)
    configuration = build_configuration(config, environment)
    gauge_action = build_gauge_action(config.gauge_action, configuration)
    action, refreshes = build_md_actions(
        gauge_action,
        config.fermions,
        pseudofermion_configs(config.update),
        configuration,
    )
    updater, state = build_update_runtime(
        config.update,
        config,
        configuration,
        action,
        refreshes,
    )
    return Simulation(
        configuration,
        action,
        updater,
        state,
        environment,
        config,
    )
end

function check_heatbath_counters(
    updater::HeatbathUpdater,
    state::HeatbathState,
)
    updater.kernel.sweep == state.heatbath_sweep || throw(ArgumentError(
        "heatbath sweep counter mismatch: runtime=$(updater.kernel.sweep), " *
        "state=$(state.heatbath_sweep)",
    ))
    updater.kernel.overrelaxation_sweep == state.overrelaxation_sweep ||
        throw(ArgumentError(
            "overrelaxation sweep counter mismatch: runtime=" *
            "$(updater.kernel.overrelaxation_sweep), " *
            "state=$(state.overrelaxation_sweep)",
        ))
    return nothing
end

function check_configuration_sequence_state(
    updater::ConfigurationSequenceUpdater,
    state::ConfigurationSequenceState,
)
    1 <= state.current_index <= length(updater.paths) || throw(ArgumentError(
        "configuration sequence index $(state.current_index) is outside " *
        "1:$(length(updater.paths))",
    ))
    state.trajectory == state.current_index - 1 || throw(ArgumentError(
        "configuration sequence counter mismatch: trajectory=" *
        "$(state.trajectory), current_index=$(state.current_index)",
    ))
    return nothing
end

function has_next_configuration(
    updater::ConfigurationSequenceUpdater,
    state::ConfigurationSequenceState,
)
    check_configuration_sequence_state(updater, state)
    return state.current_index < length(updater.paths)
end

function current_configuration_path(
    updater::ConfigurationSequenceUpdater,
    state::ConfigurationSequenceState,
)
    check_configuration_sequence_state(updater, state)
    return updater.paths[state.current_index]
end

"""Run one heatbath sweep followed by the configured number of OR sweeps."""
function update!(
    updater::HeatbathUpdater,
    configuration::GaugeConfiguration,
    state::HeatbathState,
    ::GaugefieldsEnvironment,
)
    check_heatbath_counters(updater, state)
    trajectory = state.trajectory
    gauge = configuration.gauge

    Gaugefields.heatbath!(gauge, updater.kernel)
    state.heatbath_sweep = updater.kernel.sweep
    for _ in 1:updater.overrelaxation_steps
        Gaugefields.overrelaxation!(gauge, updater.kernel)
        state.overrelaxation_sweep = updater.kernel.overrelaxation_sweep
    end

    check_heatbath_counters(updater, state)
    state.trajectory += 1
    return HeatbathUpdateResult(
        trajectory,
        state.heatbath_sweep,
        state.overrelaxation_sweep,
    )
end

"""Load the next configuration in a finite sequence."""
function update!(
    updater::ConfigurationSequenceUpdater,
    configuration::GaugeConfiguration,
    state::ConfigurationSequenceState,
    ::GaugefieldsEnvironment,
)
    check_configuration_sequence_state(updater, state)
    has_next_configuration(updater, state) || throw(EOFError())
    next_index = state.current_index + 1
    path = updater.paths[next_index]

    Gaugefields.load_configuration!(
        configuration.gauge,
        path;
        format=updater.format,
    )
    trajectory = state.trajectory
    state.current_index = next_index
    state.trajectory += 1
    check_configuration_sequence_state(updater, state)
    return ConfigurationSequenceUpdateResult(
        trajectory,
        next_index,
        path,
    )
end

function metropolis_rule(delta_hamiltonian::Real, uniform_random::Real)
    isfinite(delta_hamiltonian) || return false
    0 < uniform_random <= 1 || throw(ArgumentError(
        "the Metropolis random number must be in (0, 1]; got $uniform_random",
    ))
    return log(uniform_random) < min(zero(delta_hamiltonian), -delta_hamiltonian)
end

function metropolis_decision(delta_hamiltonian, state::HMCState, gauge)
    communicator = Gaugefields.gauge_communicator(gauge)
    if !is_distributed(communicator)
        return metropolis_rule(
            delta_hamiltonian,
            Random.rand(state.metropolis_rng),
        )
    end

    decision = Ref(false)
    if comm_rank(communicator) == 0
        decision[] = metropolis_rule(
            delta_hamiltonian,
            Random.rand(state.metropolis_rng),
        )
    end
    broadcast!(decision, 0, communicator)
    return decision[]
end

function apply_metropolis_decision!(
    configuration::AbstractConfiguration,
    backup::GaugeConfiguration,
    accepted::Bool,
)
    accepted || copy_configuration!(configuration, backup)
    return accepted
end

function refresh_pseudofermions!(refreshes::Tuple, gauge, trajectory::Integer)
    for refresh in refreshes
        LatticeDiracOperators.refresh_pseudofermion!(
            refresh.provider,
            gauge,
            refresh.noise;
            seed=refresh.seed,
            sweep=trajectory,
            subgroup=refresh.subgroup,
        )
    end
    return nothing
end

"""
Discard chronological solver guesses before each pseudofermion trajectory.

LDO action workspaces are intentionally not part of a portable checkpoint.
Starting every trajectory from cleared Krylov solution buffers makes a
continuous run and a newly constructed session follow the same solver path.
"""
function reset_pseudofermion_solver_state!(refreshes::Tuple)
    for refresh in refreshes
        LatticeDiracOperators.reset_trajectory_state!(refresh.provider)
    end
    return nothing
end

function reset_md_trajectory_state!(providers::Tuple)
    for provider in providers
        LatticeDiracOperators.reset_trajectory_state!(provider)
    end
    return nothing
end

"""Run one gauge or dynamical-fermion HMC trajectory."""
function update!(
    updater::HMCUpdater,
    configuration::AbstractConfiguration,
    state::HMCState,
    ::GaugefieldsEnvironment,
)
    gauge = gauge_links(configuration)
    trajectory = state.trajectory
    reset_pseudofermion_solver_state!(updater.pseudofermion_refreshes)
    refresh_pseudofermions!(
        updater.pseudofermion_refreshes,
        gauge,
        trajectory,
    )
    updater.momentum = Gaugefields.gaussian_momenta(
        gauge;
        sigma=updater.momentum_sigma,
        seed=updater.momentum_seed,
        sweep=trajectory,
    )
    copy_configuration!(updater.backup, configuration)

    diagnostics = try
        Gaugefields.md_trajectory!(gauge, updater.momentum, updater.md_driver)
    catch
        copy_configuration!(configuration, updater.backup)
        rethrow()
    end
    accepted = metropolis_decision(
        diagnostics.delta_hamiltonian,
        state,
        gauge,
    )
    apply_metropolis_decision!(configuration, updater.backup, accepted)

    state.trajectory += 1
    state.accepted += accepted
    return HMCUpdateResult(
        trajectory,
        accepted,
        diagnostics.initial_hamiltonian,
        diagnostics.final_hamiltonian,
        diagnostics.delta_hamiltonian,
    )
end

"""
Run one self-learning HMC trajectory.

The pseudofermion and momentum refreshes, the initial/final Hamiltonians, and
the Metropolis decision use the target action.  Only the reversible MD
trajectory is generated with the approximate action.
"""
function update!(
    updater::SLHMCUpdater,
    configuration::AbstractConfiguration,
    state::HMCState,
    ::GaugefieldsEnvironment,
)
    gauge = gauge_links(configuration)
    trajectory = state.trajectory
    reset_pseudofermion_solver_state!(updater.pseudofermion_refreshes)
    reset_md_trajectory_state!(updater.md_trajectory_state_providers)
    refresh_pseudofermions!(
        updater.pseudofermion_refreshes,
        gauge,
        trajectory,
    )
    updater.momentum = Gaugefields.gaussian_momenta(
        gauge;
        sigma=updater.momentum_sigma,
        seed=updater.momentum_seed,
        sweep=trajectory,
    )
    copy_configuration!(updater.backup, configuration)

    target_initial = Gaugefields.md_hamiltonian(
        gauge,
        updater.momentum,
        updater.target_driver,
    )
    md_diagnostics = try
        diagnostics = Gaugefields.md_trajectory!(
            gauge,
            updater.momentum,
            updater.md_driver,
        )
        target_final = Gaugefields.md_hamiltonian(
            gauge,
            updater.momentum,
            updater.target_driver,
        )
        (diagnostics, target_final)
    catch
        copy_configuration!(configuration, updater.backup)
        rethrow()
    end
    proposal_diagnostics, target_final = md_diagnostics
    target_delta = target_final - target_initial
    accepted = metropolis_decision(target_delta, state, gauge)
    apply_metropolis_decision!(configuration, updater.backup, accepted)

    state.trajectory += 1
    state.accepted += accepted
    return SLHMCUpdateResult(
        trajectory,
        accepted,
        target_initial,
        target_final,
        target_delta,
        proposal_diagnostics.initial_hamiltonian,
        proposal_diagnostics.final_hamiltonian,
        proposal_diagnostics.delta_hamiltonian,
    )
end

function update!(simulation::Simulation)
    return update!(
        simulation.updater,
        simulation.configuration,
        simulation.state,
        simulation.environment,
    )
end

has_next_configuration(simulation::Simulation) = has_next_configuration(
    simulation.updater,
    simulation.state,
)

current_configuration_path(simulation::Simulation) =
    current_configuration_path(simulation.updater, simulation.state)

function Base.show(io::IO, simulation::Simulation)
    print(io, "Simulation(configuration=")
    print(io, nameof(typeof(simulation.configuration)))
    print(io, ", updater=", nameof(typeof(simulation.updater)))
    print(io, ", trajectory=", simulation.state.trajectory, ")")
end

end
