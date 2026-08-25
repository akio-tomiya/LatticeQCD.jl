module SimulationInput_module

import TOML
import QCDMeasurements

import ..Parameter_structs:
    Print_Fermions_parameters,
    Print_Gradientflow_parameters,
    Print_HMCrelated_parameters,
    Print_Physical_parameters,
    Print_System_control_parameters
import ..LQCDConfig_module: lqcd_config
import ..LQCDConfig_module
import ..MeasurementPlan_module
import ..SimulationSession_module:
    SimulationSpec,
    legacy_output_config,
    legacy_simulation_schedule
import ..SimulationSession_module

const LEGACY_PHYSICAL = "Physical setting"
const LEGACY_FERMIONS = "Physical setting(fermions)"
const LEGACY_CONTROL = "System Control"
const LEGACY_HMC = "HMC related"
const LEGACY_MEASUREMENTS = "Measurement set"
const LEGACY_FLOW = "gradientflow_measurements"

function legacy_section(document::AbstractDict, name::AbstractString)
    section = get(document, name, Dict{String,Any}())
    section isa AbstractDict || throw(ArgumentError(
        "legacy TOML section $(repr(name)) must be a table",
    ))
    return section
end

function legacy_value(
    section::AbstractDict,
    name::AbstractString,
    default;
    preserve_nothing_string::Bool=false,
)
    value = get(section, name, default)
    if !preserve_nothing_string && value == "nothing"
        return nothing
    end
    return value
end

function legacy_value(
    sections::Tuple,
    name::AbstractString,
    default;
    preserve_nothing_string::Bool=false,
)
    for section in sections
        if haskey(section, name)
            return legacy_value(
                section,
                name,
                default;
                preserve_nothing_string,
            )
        end
    end
    return default
end

function shallow_property_dictionary(value)
    return Dict{String,Any}(
        String(name) => getproperty(value, name) for name in propertynames(value)
    )
end

function normalize_legacy_measurement(method::AbstractDict)
    normalized = Dict{String,Any}()
    for (raw_name, value) in method
        name = String(raw_name)
        if name == "fermion_parameters"
            value isa AbstractDict || throw(ArgumentError(
                "measurement fermion_parameters must be a table",
            ))
            fermion_type = get(value, "Dirac_operator", nothing)
            isnothing(fermion_type) && throw(ArgumentError(
                "measurement fermion_parameters requires Dirac_operator",
            ))
            defaults = shallow_property_dictionary(
                QCDMeasurements.initialize_fermion_parameters(fermion_type),
            )
            for (parameter_name, parameter_default) in defaults
                normalized[parameter_name] = get(
                    value,
                    parameter_name,
                    parameter_default,
                )
            end
        else
            normalized[name] = value
        end
    end
    return normalized
end

function normalize_legacy_measurements(methods::AbstractDict)
    return Dict[
        normalize_legacy_measurement(method) for (_, method) in methods
    ]
end

function normalize_legacy_measurements(methods::AbstractVector)
    return Dict[
        normalize_legacy_measurement(method) for method in methods
    ]
end

function normalize_legacy_measurements(methods)
    throw(ArgumentError(
        "legacy measurement_methods must be a table or array of tables; " *
        "got $(typeof(methods))",
    ))
end

function legacy_measurements(section::AbstractDict, name::AbstractString)
    methods = get(section, name, Dict{String,Any}())
    return normalize_legacy_measurements(methods)
end

function legacy_domainwall_value(
    fermions::AbstractDict,
    current_name::AbstractString,
    old_name::AbstractString,
    default,
)
    if haskey(fermions, current_name)
        return legacy_value(fermions, current_name, default)
    end
    return legacy_value(fermions, old_name, default)
end

"""
Return the side-effect-free legacy values consumed by the typed decoders.

Unlike `construct_Params_from_TOML`, this function does not open a logfile,
create output directories, or mutate the parsed TOML dictionary. The returned
`NamedTuple` is only an input adapter; dictionaries are not retained in the
resulting `SimulationSpec`.
"""
function legacy_simulation_values(document::AbstractDict)
    physical = legacy_section(document, LEGACY_PHYSICAL)
    fermions = legacy_section(document, LEGACY_FERMIONS)
    control = legacy_section(document, LEGACY_CONTROL)
    hmc = legacy_section(document, LEGACY_HMC)
    measurements = legacy_section(document, LEGACY_MEASUREMENTS)
    flow = legacy_section(document, LEGACY_FLOW)
    slhmc = legacy_section(document, "SLHMC related")
    slmc = legacy_section(document, "SLMC related")

    physical_defaults = Print_Physical_parameters()
    fermion_defaults = Print_Fermions_parameters()
    control_defaults = Print_System_control_parameters()
    hmc_defaults = Print_HMCrelated_parameters()
    flow_defaults = Print_Gradientflow_parameters()

    return (
        L=Tuple(Int.(legacy_value(physical, "L", physical_defaults.L))),
        β=legacy_value(physical, "β", physical_defaults.β),
        NC=Int(legacy_value(physical, "NC", physical_defaults.NC)),
        Nthermalization=Int(legacy_value(
            physical,
            "Nthermalization",
            physical_defaults.Nthermalization,
        )),
        Nsteps=Int(legacy_value(physical, "Nsteps", physical_defaults.Nsteps)),
        initial=String(legacy_value(
            physical,
            "initial",
            physical_defaults.initial,
            preserve_nothing_string=true,
        )),
        initialtrj=Int(legacy_value(
            physical,
            "initialtrj",
            physical_defaults.initialtrj,
        )),
        update_method=String(legacy_value(
            physical,
            "update_method",
            physical_defaults.update_method,
            preserve_nothing_string=true,
        )),
        useOR=Bool(legacy_value(physical, "useOR", physical_defaults.useOR)),
        numOR=Int(legacy_value(physical, "numOR", physical_defaults.numOR)),
        Nwing=Int(legacy_value(physical, "Nwing", physical_defaults.Nwing)),
        quench=Bool(legacy_value(fermions, "quench", fermion_defaults.quench)),
        Dirac_operator=legacy_value(
            fermions,
            "Dirac_operator",
            fermion_defaults.Dirac_operator,
        ),
        Clover_coefficient=legacy_value(
            fermions,
            "Clover_coefficient",
            fermion_defaults.Clover_coefficient,
        ),
        r=legacy_value(fermions, "r", fermion_defaults.r),
        hop=legacy_value(fermions, "hop", fermion_defaults.hop),
        Nf=Int(legacy_value(fermions, "Nf", fermion_defaults.Nf)),
        mass=legacy_value(fermions, "mass", fermion_defaults.mass),
        naik_epsilon=legacy_value(
            fermions,
            "naik_epsilon",
            fermion_defaults.naik_epsilon,
        ),
        b=legacy_value(fermions, "b", fermion_defaults.b),
        c=legacy_value(fermions, "c", fermion_defaults.c),
        Domainwall_M=legacy_domainwall_value(
            fermions,
            "Domainwall_M",
            "M",
            fermion_defaults.Domainwall_M,
        ),
        Domainwall_m=legacy_domainwall_value(
            fermions,
            "Domainwall_m",
            "m",
            fermion_defaults.Domainwall_m,
        ),
        Domainwall_L5=legacy_domainwall_value(
            fermions,
            "Domainwall_L5",
            "N5",
            fermion_defaults.Domainwall_L5,
        ),
        BoundaryCondition=legacy_value(
            fermions,
            "BoundaryCondition",
            fermion_defaults.BoundaryCondition,
        ),
        smearing_for_fermion=String(legacy_value(
            fermions,
            "smearing_for_fermion",
            fermion_defaults.smearing_for_fermion,
            preserve_nothing_string=true,
        )),
        stout_numlayers=legacy_value(
            fermions,
            "stout_numlayers",
            fermion_defaults.stout_numlayers,
        ),
        stout_ρ=legacy_value(
            fermions,
            "stout_ρ",
            fermion_defaults.stout_ρ,
        ),
        stout_loops=legacy_value(
            fermions,
            "stout_loops",
            fermion_defaults.stout_loops,
        ),
        loadU_format=legacy_value(
            control,
            "loadU_format",
            control_defaults.loadU_format,
        ),
        loadU_dir=String(legacy_value(
            control,
            "loadU_dir",
            control_defaults.loadU_dir,
            preserve_nothing_string=true,
        )),
        loadU_fromfile=Bool(legacy_value(
            control,
            "loadU_fromfile",
            control_defaults.loadU_fromfile,
        )),
        loadU_filename=String(legacy_value(
            control,
            "loadU_filename",
            control_defaults.loadU_filename,
            preserve_nothing_string=true,
        )),
        saveU_dir=String(legacy_value(
            control,
            "saveU_dir",
            control_defaults.saveU_dir,
            preserve_nothing_string=true,
        )),
        saveU_format=legacy_value(
            control,
            "saveU_format",
            control_defaults.saveU_format,
        ),
        saveU_every=Int(legacy_value(
            control,
            "saveU_every",
            control_defaults.saveU_every,
        )),
        verboselevel=Int(legacy_value(
            control,
            "verboselevel",
            control_defaults.verboselevel,
        )),
        randomseed=Int(legacy_value(
            control,
            "randomseed",
            control_defaults.randomseed,
        )),
        isevenodd=Bool(legacy_value(
            control,
            "isevenodd",
            control_defaults.isevenodd,
        )),
        Δτ=legacy_value(hmc, "Δτ", hmc_defaults.Δτ),
        SextonWeingargten=Bool(legacy_value(
            hmc,
            "SextonWeingargten",
            hmc_defaults.SextonWeingargten,
        )),
        N_SextonWeingargten=Int(legacy_value(
            hmc,
            "N_SextonWeingargten",
            hmc_defaults.N_SextonWeingargten,
        )),
        MDsteps=Int(legacy_value(hmc, "MDsteps", hmc_defaults.MDsteps)),
        eps=legacy_value(hmc, "eps", hmc_defaults.eps),
        MaxCGstep=Int(legacy_value(
            hmc,
            "MaxCGstep",
            hmc_defaults.MaxCGstep,
        )),
        QPQ=Bool(legacy_value(hmc, "QPQ", hmc_defaults.QPQ)),
        βeff=legacy_value((slhmc, slmc), "βeff", 0.0),
        ITERATION_MAX=Int(legacy_value(
            (physical, control),
            "ITERATION_MAX",
            100_000,
        )),
        measurement_methods=legacy_measurements(
            measurements,
            "measurement_methods",
        ),
        hasgradientflow=Bool(legacy_value(
            (flow, control),
            "hasgradientflow",
            flow_defaults.hasgradientflow,
        )),
        eps_flow=legacy_value(
            (flow, control),
            "eps_flow",
            flow_defaults.eps_flow,
        ),
        numflow=Int(legacy_value(
            (flow, control),
            "numflow",
            flow_defaults.numflow,
        )),
        Nflow=Int(legacy_value(
            (flow, control),
            "Nflow",
            flow_defaults.Nflow,
        )),
        measurements_for_flow=legacy_measurements(
            flow,
            "measurements_for_flow",
        ),
    )
end

"""Convert an already-parsed legacy TOML document directly to typed input."""
function simulation_spec_from_legacy_toml(document::AbstractDict)
    values = legacy_simulation_values(document)
    return SimulationSpec(
        lqcd_config(values),
        legacy_simulation_schedule(values),
        legacy_output_config(values),
    )
end

"""
Load a Wizard/legacy TOML file directly as a typed `SimulationSpec`.

This read is side-effect free apart from reading `filename`: no `Params`,
open streams, log files, or output directories are created.
"""
const SIMULATION_SPEC_FORMAT = "LatticeQCD.SimulationSpec"
const SIMULATION_SPEC_SCHEMA_VERSION = 1

required(table::AbstractDict, name::AbstractString) = haskey(table, name) ?
    table[name] : throw(ArgumentError("required TOML key $(repr(name)) is missing"))

string_name(value) = String(value)

function initialization_dictionary(::LQCDConfig_module.ColdStartConfig)
    return Dict{String,Any}("kind" => "cold")
end

function initialization_dictionary(config::LQCDConfig_module.HotStartConfig)
    values = Dict{String,Any}("kind" => "hot")
    isnothing(config.seed) || (values["seed"] = config.seed)
    return values
end

function initialization_dictionary(config::LQCDConfig_module.FileStartConfig)
    values = Dict{String,Any}(
        "kind" => "file",
        "path" => config.path,
    )
    isnothing(config.format) || (values["format"] = String(config.format))
    return values
end

function initialization_dictionary(::LQCDConfig_module.InstantonConfig)
    return Dict{String,Any}("kind" => "instanton")
end

function initialization_dictionary(
    config::LQCDConfig_module.EmbeddedInstantonConfig,
)
    values = Dict{String,Any}(
        "kind" => "embedded_instanton",
        "sign" => config.sign,
        "block" => collect(config.block),
    )
    isnothing(config.center) || (values["center"] = collect(config.center))
    isnothing(config.radius) || (values["radius"] = config.radius)
    return values
end

function parse_initialization(values::AbstractDict)
    kind = lowercase(String(required(values, "kind")))
    kind == "cold" && return LQCDConfig_module.ColdStartConfig()
    kind == "hot" && return LQCDConfig_module.HotStartConfig(
        get(values, "seed", nothing),
    )
    kind == "file" && return LQCDConfig_module.FileStartConfig(
        String(required(values, "path")),
        get(values, "format", nothing),
    )
    kind == "instanton" && return LQCDConfig_module.InstantonConfig()
    kind == "embedded_instanton" && return (
        LQCDConfig_module.EmbeddedInstantonConfig(
            center=haskey(values, "center") ? Tuple(values["center"]) : nothing,
            radius=get(values, "radius", nothing),
            sign=Int(get(values, "sign", 1)),
            block=Tuple(get(values, "block", (1, 2))),
        )
    )
    throw(ArgumentError("unsupported gauge initialization kind=$(repr(kind))"))
end

function gauge_action_dictionary(action::LQCDConfig_module.GaugeActionConfig)
    terms = Dict{String,Any}[]
    for term in action.terms
        term.loop isa Union{Symbol,AbstractString} || throw(ArgumentError(
            "canonical TOML currently requires a named gauge loop; got " *
            "$(typeof(term.loop)) for term $(term.name)",
        ))
        push!(terms, Dict{String,Any}(
            "name" => String(term.name),
            "loop" => String(term.loop),
            "coupling" => term.coupling,
            "include_adjoint" => term.include_adjoint,
        ))
    end
    return Dict{String,Any}("terms" => terms)
end

function parse_gauge_action(values::AbstractDict)
    terms = Tuple(
        LQCDConfig_module.GaugeActionTermConfig(
            Symbol(String(required(term, "name"))),
            String(required(term, "loop")),
            required(term, "coupling"),
            Bool(get(term, "include_adjoint", true)),
        ) for term in required(values, "terms")
    )
    return LQCDConfig_module.GaugeActionConfig(terms...)
end

function operator_dictionary(config::LQCDConfig_module.WilsonDiracConfig)
    return Dict{String,Any}(
        "kind" => "wilson",
        "hopping_parameter" => config.hopping_parameter,
        "wilson_parameter" => config.wilson_parameter,
    )
end

function operator_dictionary(config::LQCDConfig_module.WilsonCloverDiracConfig)
    return Dict{String,Any}(
        "kind" => "wilson_clover",
        "hopping_parameter" => config.hopping_parameter,
        "wilson_parameter" => config.wilson_parameter,
        "clover_coefficient" => config.clover_coefficient,
    )
end

function operator_dictionary(config::LQCDConfig_module.StaggeredDiracConfig)
    return Dict{String,Any}("kind" => "staggered", "mass" => config.mass)
end

function operator_dictionary(config::LQCDConfig_module.HISQDiracConfig)
    return Dict{String,Any}(
        "kind" => "hisq",
        "mass" => config.mass,
        "naik_epsilon" => config.naik_epsilon,
    )
end

function operator_dictionary(config::LQCDConfig_module.DomainwallDiracConfig)
    return Dict{String,Any}(
        "kind" => "domainwall",
        "mass" => config.mass,
        "domainwall_height" => config.domainwall_height,
        "fifth_dimension" => config.fifth_dimension,
    )
end

function operator_dictionary(
    config::LQCDConfig_module.MobiusDomainwallDiracConfig,
)
    return Dict{String,Any}(
        "kind" => "mobius_domainwall",
        "mass" => config.mass,
        "domainwall_height" => config.domainwall_height,
        "fifth_dimension" => config.fifth_dimension,
        "b" => config.b,
        "c" => config.c,
    )
end

function parse_operator(values::AbstractDict)
    kind = lowercase(String(required(values, "kind")))
    kind == "wilson" && return LQCDConfig_module.WilsonDiracConfig(
        required(values, "hopping_parameter"),
        get(values, "wilson_parameter", 1.0),
    )
    kind == "wilson_clover" && return (
        LQCDConfig_module.WilsonCloverDiracConfig(
            required(values, "hopping_parameter"),
            get(values, "wilson_parameter", 1.0),
            required(values, "clover_coefficient"),
        )
    )
    kind == "staggered" && return LQCDConfig_module.StaggeredDiracConfig(
        required(values, "mass"),
    )
    kind == "hisq" && return LQCDConfig_module.HISQDiracConfig(
        required(values, "mass"),
        get(values, "naik_epsilon", 0.0),
    )
    kind == "domainwall" && return LQCDConfig_module.DomainwallDiracConfig(
        required(values, "mass"),
        required(values, "domainwall_height"),
        Int(required(values, "fifth_dimension")),
    )
    kind in ("mobius_domainwall", "mobiusdomainwall") && return (
        LQCDConfig_module.MobiusDomainwallDiracConfig(
            required(values, "mass"),
            required(values, "domainwall_height"),
            Int(required(values, "fifth_dimension")),
            get(values, "b", 2.0),
            get(values, "c", 1.0),
        )
    )
    throw(ArgumentError("unsupported fermion operator kind=$(repr(kind))"))
end

smearing_dictionary(::LQCDConfig_module.NoFermionSmearingConfig) =
    Dict{String,Any}("kind" => "none")

function smearing_dictionary(config::LQCDConfig_module.StoutFermionSmearingConfig)
    return Dict{String,Any}(
        "kind" => "stout",
        "layers" => config.layers,
        "coefficients" => collect(config.coefficients),
        "loops" => collect(String.(config.loops)),
    )
end

function parse_smearing(values::AbstractDict)
    kind = lowercase(String(required(values, "kind")))
    kind in ("none", "nothing") && return (
        LQCDConfig_module.NoFermionSmearingConfig()
    )
    kind == "stout" && return LQCDConfig_module.StoutFermionSmearingConfig(
        Int(required(values, "layers")),
        Tuple(required(values, "coefficients")),
        Tuple(String.(required(values, "loops"))),
    )
    throw(ArgumentError("unsupported fermion smearing kind=$(repr(kind))"))
end

function fermion_dictionary(config::LQCDConfig_module.FermionActionConfig)
    return Dict{String,Any}(
        "name" => String(config.name),
        "flavors" => config.flavors,
        "operator" => operator_dictionary(config.operator),
        "solver" => Dict{String,Any}(
            "tolerance" => config.solver.tolerance,
            "max_steps" => config.solver.max_steps,
            "verbose" => config.solver.verbose,
            "boundary_conditions" => collect(config.solver.boundary_conditions),
        ),
        "smearing" => smearing_dictionary(config.smearing),
    )
end

function parse_fermion(values::AbstractDict)
    solver = required(values, "solver")
    return LQCDConfig_module.FermionActionConfig(
        Symbol(String(required(values, "name"))),
        parse_operator(required(values, "operator")),
        Int(required(values, "flavors")),
        LQCDConfig_module.FermionSolverConfig(
            required(solver, "tolerance"),
            Int(required(solver, "max_steps")),
            Int(get(solver, "verbose", 0)),
            Tuple(Int.(required(solver, "boundary_conditions"))),
        ),
        parse_smearing(get(values, "smearing", Dict("kind" => "none"))),
    )
end

ordering_name(::LQCDConfig_module.QPQConfig) = "QPQ"
ordering_name(::LQCDConfig_module.PQPConfig) = "PQP"

function parse_ordering(value)
    name = uppercase(String(value))
    name == "QPQ" && return LQCDConfig_module.QPQConfig()
    name == "PQP" && return LQCDConfig_module.PQPConfig()
    throw(ArgumentError("integrator ordering must be QPQ or PQP; got $value"))
end

function integrator_dictionary(config::LQCDConfig_module.LeapfrogConfig)
    return Dict{String,Any}(
        "kind" => "leapfrog",
        "ordering" => ordering_name(config.ordering),
        "forces" => collect(String.(config.forces.names)),
    )
end

function integrator_dictionary(config::LQCDConfig_module.SextonWeingartenConfig)
    return Dict{String,Any}(
        "kind" => "sexton_weingarten",
        "ordering" => ordering_name(config.ordering),
        "slow_forces" => collect(String.(config.slow_forces.names)),
        "fast_steps" => config.fast_steps,
        "fast_integrator" => integrator_dictionary(config.fast_integrator),
    )
end

function parse_integrator(values::AbstractDict)
    kind = lowercase(String(required(values, "kind")))
    ordering = parse_ordering(required(values, "ordering"))
    kind == "leapfrog" && return LQCDConfig_module.LeapfrogConfig(
        ordering,
        LQCDConfig_module.ForceGroupConfig(
            Tuple(Symbol.(String.(required(values, "forces")))),
        ),
    )
    kind == "sexton_weingarten" && return (
        LQCDConfig_module.SextonWeingartenConfig(
            ordering,
            LQCDConfig_module.ForceGroupConfig(
                Tuple(Symbol.(String.(required(values, "slow_forces")))),
            ),
            parse_integrator(required(values, "fast_integrator")),
            Int(required(values, "fast_steps")),
        )
    )
    throw(ArgumentError("unsupported MD integrator kind=$(repr(kind))"))
end

function md_dictionary(config::LQCDConfig_module.MDConfig)
    return Dict{String,Any}(
        "step_size" => config.step_size,
        "steps" => config.steps,
        "integrator" => integrator_dictionary(config.integrator),
    )
end

function parse_md(values::AbstractDict)
    return LQCDConfig_module.MDConfig(
        required(values, "step_size"),
        Int(required(values, "steps")),
        parse_integrator(required(values, "integrator")),
    )
end

function random_dictionary(config::LQCDConfig_module.RandomStreamConfig)
    return Dict{String,Any}(
        "seed" => config.seed,
        "stream" => String(config.name),
    )
end

function parse_random(values::AbstractDict)
    return LQCDConfig_module.RandomStreamConfig(
        Int(required(values, "seed")),
        Symbol(String(required(values, "stream"))),
    )
end

function hmc_common_dictionary(config)
    return Dict{String,Any}(
        "md" => md_dictionary(config.md),
        "momentum" => Dict{String,Any}(
            "sigma" => config.momentum.sigma,
            "random" => random_dictionary(config.momentum.random),
        ),
        "acceptance" => random_dictionary(config.acceptance.random),
        "pseudofermions" => [
            Dict{String,Any}(
                "action_name" => String(refresh.action_name),
                "subgroup" => refresh.subgroup,
                "random" => random_dictionary(refresh.random),
            ) for refresh in config.pseudofermions
        ],
    )
end

function update_dictionary(config::LQCDConfig_module.HMCConfig)
    values = hmc_common_dictionary(config)
    values["kind"] = "HMC"
    return values
end

function update_dictionary(config::LQCDConfig_module.SLHMCConfig)
    values = hmc_common_dictionary(config)
    values["kind"] = "SLHMC"
    values["md_gauge_action"] = gauge_action_dictionary(config.md_gauge_action)
    values["md_fermions"] = [fermion_dictionary(f) for f in config.md_fermions]
    return values
end

function update_dictionary(config::LQCDConfig_module.HeatbathConfig)
    return Dict{String,Any}(
        "kind" => "Heatbath",
        "even_odd" => config.even_odd,
        "max_iterations" => config.max_iterations,
        "overrelaxation_steps" => config.overrelaxation_steps,
        "random" => random_dictionary(config.random),
    )
end

function source_dictionary(source::LQCDConfig_module.DirectorySourceConfig)
    return Dict{String,Any}(
        "kind" => "directory",
        "directory" => source.directory,
    )
end

function source_dictionary(source::LQCDConfig_module.ManifestSourceConfig)
    return Dict{String,Any}(
        "kind" => "manifest",
        "directory" => source.directory,
        "manifest" => source.manifest,
    )
end

function update_dictionary(config::LQCDConfig_module.ConfigurationSequenceConfig)
    return Dict{String,Any}(
        "kind" => "configuration_sequence",
        "format" => String(config.format),
        "source" => source_dictionary(config.source),
    )
end

function parse_hmc_common(values::AbstractDict)
    momentum = required(values, "momentum")
    pseudofermions = Tuple(
        LQCDConfig_module.PseudofermionRefreshConfig(
            Symbol(String(required(refresh, "action_name"))),
            parse_random(required(refresh, "random")),
            Int(required(refresh, "subgroup")),
        ) for refresh in get(values, "pseudofermions", Any[])
    )
    return (
        md=parse_md(required(values, "md")),
        momentum=LQCDConfig_module.GaussianMomentumConfig(
            required(momentum, "sigma"),
            parse_random(required(momentum, "random")),
        ),
        acceptance=LQCDConfig_module.RankZeroMetropolisConfig(
            parse_random(required(values, "acceptance")),
        ),
        pseudofermions=pseudofermions,
    )
end

function parse_source(values::AbstractDict)
    kind = lowercase(String(required(values, "kind")))
    kind == "directory" && return LQCDConfig_module.DirectorySourceConfig(
        String(required(values, "directory")),
    )
    kind == "manifest" && return LQCDConfig_module.ManifestSourceConfig(
        String(required(values, "directory")),
        String(required(values, "manifest")),
    )
    throw(ArgumentError("unsupported configuration source kind=$(repr(kind))"))
end

function parse_update(values::AbstractDict)
    kind = lowercase(String(required(values, "kind")))
    if kind == "hmc"
        common = parse_hmc_common(values)
        return LQCDConfig_module.HMCConfig(
            common.md,
            common.momentum,
            common.acceptance,
            common.pseudofermions,
        )
    elseif kind == "slhmc"
        common = parse_hmc_common(values)
        md_fermions = Tuple(parse_fermion(f) for f in get(
            values,
            "md_fermions",
            Any[],
        ))
        return LQCDConfig_module.SLHMCConfig(
            common.md,
            common.momentum,
            common.acceptance,
            common.pseudofermions,
            parse_gauge_action(required(values, "md_gauge_action")),
            md_fermions,
        )
    elseif kind == "heatbath"
        return LQCDConfig_module.HeatbathConfig(
            Bool(get(values, "even_odd", true)),
            Int(required(values, "max_iterations")),
            Int(get(values, "overrelaxation_steps", 0)),
            parse_random(required(values, "random")),
        )
    elseif kind == "configuration_sequence"
        return LQCDConfig_module.ConfigurationSequenceConfig(
            parse_source(required(values, "source")),
            String(required(values, "format")),
        )
    end
    throw(ArgumentError("unsupported update kind=$(repr(kind))"))
end

function scheduled_measurement_dictionary(scheduled)
    values = MeasurementPlan_module.measurement_dictionary(
        scheduled.observable,
    )
    values["every"] = scheduled.schedule.every
    values["start"] = scheduled.schedule.start
    return values
end

function parse_measurement_plan(values::AbstractVector)
    scheduled = map(values) do raw
        measurement_values = Dict{String,Any}(String(k) => v for (k, v) in raw)
        every = Int(required(measurement_values, "every"))
        start = Int(get(measurement_values, "start", 0))
        delete!(measurement_values, "every")
        delete!(measurement_values, "start")
        measurement_values["measure_every"] = every
        return MeasurementPlan_module.measurement_config(
            measurement_values;
            start,
        )
    end
    return MeasurementPlan_module.MeasurementPlan(Tuple(scheduled))
end

function measurement_program_dictionary(program)
    flow = if program.gradient_flow isa
              MeasurementPlan_module.NoGradientFlowMeasurementConfig
        Dict{String,Any}("enabled" => false)
    else
        config = program.gradient_flow
        Dict{String,Any}(
            "enabled" => true,
            "step_size" => config.step_size,
            "samples" => config.samples,
            "integration_steps" => config.integration_steps,
            "measurements" => [
                scheduled_measurement_dictionary(m) for
                m in config.measurements.measurements
            ],
        )
    end
    return Dict{String,Any}(
        "direct" => [
            scheduled_measurement_dictionary(m) for
            m in program.direct.measurements
        ],
        "gradient_flow" => flow,
    )
end

function parse_measurement_program(values::AbstractDict)
    direct = parse_measurement_plan(get(values, "direct", Any[]))
    flow_values = get(values, "gradient_flow", Dict("enabled" => false))
    flow = if Bool(get(flow_values, "enabled", false))
        MeasurementPlan_module.GradientFlowMeasurementConfig(
            required(flow_values, "step_size"),
            Int(required(flow_values, "samples")),
            Int(required(flow_values, "integration_steps")),
            parse_measurement_plan(required(flow_values, "measurements")),
        )
    else
        MeasurementPlan_module.NoGradientFlowMeasurementConfig()
    end
    return MeasurementPlan_module.MeasurementProgram(direct, flow)
end

configuration_output_dictionary(
    ::SimulationSession_module.NoConfigurationOutput,
) = Dict{String,Any}("format" => "none")

function configuration_output_dictionary(config)
    format = config isa SimulationSession_module.JLD2ConfigurationOutput ?
        "jld2" :
        config isa SimulationSession_module.BridgeTextConfigurationOutput ?
        "bridge" : "ildg"
    return Dict{String,Any}(
        "format" => format,
        "directory" => config.directory,
        "prefix" => config.prefix,
        "every" => config.every,
        "width" => config.width,
    )
end

function parse_output(values::AbstractDict)
    configurations = get(values, "configurations", Dict("format" => "none"))
    format = lowercase(String(get(configurations, "format", "none")))
    format in ("none", "nothing") && return (
        SimulationSession_module.OutputConfig()
    )
    constructor = if format == "jld2"
        SimulationSession_module.JLD2ConfigurationOutput
    elseif format in ("bridge", "bridge_text")
        SimulationSession_module.BridgeTextConfigurationOutput
    elseif format == "ildg"
        SimulationSession_module.ILDGConfigurationOutput
    else
        throw(ArgumentError("unsupported output format=$(repr(format))"))
    end
    config = constructor(
        String(required(configurations, "directory")),
        String(get(configurations, "prefix", "conf_")),
        Int(get(configurations, "every", 1)),
        Int(get(configurations, "width", 8)),
    )
    return SimulationSession_module.OutputConfig(config)
end

"""Convert a typed specification to the canonical, versioned TOML tree."""
function simulation_spec_dictionary(spec::SimulationSpec)
    config = spec.config
    return Dict{String,Any}(
        "format" => SIMULATION_SPEC_FORMAT,
        "schema_version" => SIMULATION_SPEC_SCHEMA_VERSION,
        "lattice" => Dict{String,Any}("size" => collect(config.lattice.L)),
        "gauge" => Dict{String,Any}(
            "colors" => config.gauge.NC,
            "halo" => config.gauge.halo,
            "initialization" => initialization_dictionary(
                config.gauge.initialization,
            ),
        ),
        "gauge_action" => gauge_action_dictionary(config.gauge_action),
        "fermions" => [fermion_dictionary(f) for f in config.fermions],
        "update" => update_dictionary(config.update),
        "schedule" => Dict{String,Any}(
            "thermalization_steps" => spec.schedule.thermalization_steps,
            "production_steps" => spec.schedule.production_steps,
            "initial_trajectory" => spec.schedule.initial_trajectory,
            "measurements" => measurement_program_dictionary(
                spec.schedule.measurements,
            ),
        ),
        "output" => Dict{String,Any}(
            "configurations" => configuration_output_dictionary(
                spec.output.configurations,
            ),
        ),
    )
end

"""Parse the canonical, versioned `SimulationSpec` TOML tree."""
function simulation_spec_from_toml(document::AbstractDict)
    get(document, "format", nothing) == SIMULATION_SPEC_FORMAT ||
        throw(ArgumentError("not a canonical $SIMULATION_SPEC_FORMAT document"))
    version = Int(required(document, "schema_version"))
    version == SIMULATION_SPEC_SCHEMA_VERSION || throw(ArgumentError(
        "unsupported SimulationSpec schema_version=$version; expected " *
        "$SIMULATION_SPEC_SCHEMA_VERSION",
    ))
    lattice_values = required(document, "lattice")
    gauge_values = required(document, "gauge")
    fermions = Tuple(parse_fermion(f) for f in get(document, "fermions", Any[]))
    config = LQCDConfig_module.LQCDConfig(
        LQCDConfig_module.LatticeConfig(
            Tuple(Int.(required(lattice_values, "size"))),
        ),
        LQCDConfig_module.GaugeConfig(
            Int(required(gauge_values, "colors")),
            Int(get(gauge_values, "halo", 0)),
            parse_initialization(required(gauge_values, "initialization")),
        ),
        parse_gauge_action(required(document, "gauge_action")),
        fermions,
        parse_update(required(document, "update")),
    )
    schedule_values = required(document, "schedule")
    schedule = SimulationSession_module.SimulationSchedule(
        Int(required(schedule_values, "thermalization_steps")),
        Int(required(schedule_values, "production_steps")),
        Int(get(schedule_values, "initial_trajectory", 0)),
        parse_measurement_program(get(
            schedule_values,
            "measurements",
            Dict{String,Any}(),
        )),
    )
    return SimulationSpec(
        config,
        schedule,
        parse_output(get(document, "output", Dict{String,Any}())),
    )
end

"""Parse canonical TOML, or fall back to a Wizard/legacy document."""
function parse_simulation_spec(document::AbstractDict)
    if get(document, "format", nothing) == SIMULATION_SPEC_FORMAT
        return simulation_spec_from_toml(document)
    end
    return simulation_spec_from_legacy_toml(document)
end

function load_simulation_spec(filename::AbstractString)
    return parse_simulation_spec(TOML.parsefile(filename))
end

"""Write a typed specification in the canonical, versioned TOML format."""
function write_simulation_spec(filename::AbstractString, spec::SimulationSpec)
    open(filename, "w") do io
        TOML.print(io, simulation_spec_dictionary(spec))
    end
    return filename
end

export legacy_simulation_values,
    simulation_spec_from_legacy_toml,
    SIMULATION_SPEC_FORMAT,
    SIMULATION_SPEC_SCHEMA_VERSION,
    simulation_spec_dictionary,
    simulation_spec_from_toml,
    parse_simulation_spec,
    load_simulation_spec,
    write_simulation_spec

end
