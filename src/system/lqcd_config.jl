module LQCDConfig_module

import ..System_parameters: Params

"""Geometry of a lattice, independent of fields and halo storage."""
struct LatticeConfig{Dim,T<:Integer}
    L::NTuple{Dim,T}
    function LatticeConfig{Dim,T}(
        L::NTuple{Dim,T},
    ) where {Dim,T<:Integer}
        return new{Dim,T}(L)
    end
end

function LatticeConfig(L::Tuple)
    isempty(L) && throw(ArgumentError("a lattice needs at least one extent"))
    all(extent -> extent isa Integer, L) || throw(ArgumentError(
        "lattice extents must be integers; got $L",
    ))
    extent_type = promote_type(map(typeof, L)...)
    extents = ntuple(index -> convert(extent_type, L[index]), length(L))
    return LatticeConfig{length(extents),extent_type}(extents)
end

"""Common supertype for typed gauge-field initialization settings."""
abstract type AbstractGaugeInitializationConfig end

"""Initialize every link to the identity."""
struct ColdStartConfig <: AbstractGaugeInitializationConfig end

"""Initialize links from the Gaugefields reproducible hot-start stream."""
struct HotStartConfig{S<:Union{Nothing,Integer}} <:
       AbstractGaugeInitializationConfig
    seed::S
end

HotStartConfig() = HotStartConfig(nothing)

"""Load a configuration from `path` using a Gaugefields file format."""
struct FileStartConfig{
    P<:AbstractString,
    F<:Union{Nothing,Symbol,AbstractString},
} <: AbstractGaugeInitializationConfig
    path::P
    format::F
end

"""The historical SU(2) instanton with radius equal to half the x extent."""
struct InstantonConfig <: AbstractGaugeInitializationConfig end

"""
An SU(2) instanton embedded in an SU(N) color block.

`nothing` selects the Gaugefields defaults for `center` and `radius`. `sign`
is `+1` for an instanton and `-1` for an anti-instanton.
"""
struct EmbeddedInstantonConfig{
    C<:Union{Nothing,Tuple},
    R<:Union{Nothing,Real},
    B<:Tuple,
} <: AbstractGaugeInitializationConfig
    center::C
    radius::R
    sign::Int
    block::B

    function EmbeddedInstantonConfig(
        center::C,
        radius::R,
        sign::Integer,
        block::B,
    ) where {
        C<:Union{Nothing,Tuple},
        R<:Union{Nothing,Real},
        B<:Tuple,
    }
        sign in (-1, 1) || throw(ArgumentError(
            "instanton sign must be +1 or -1; got $sign",
        ))
        center === nothing || length(center) == 4 || throw(ArgumentError(
            "an embedded-instanton center must contain four coordinates",
        ))
        radius === nothing || radius > 0 || throw(ArgumentError(
            "instanton radius must be positive; got $radius",
        ))
        length(block) == 2 || throw(ArgumentError(
            "an embedded-instanton block must contain two color indices",
        ))
        all(index -> index isa Integer, block) || throw(ArgumentError(
            "embedded-instanton block indices must be integers",
        ))
        return new{C,R,B}(center, radius, Int(sign), block)
    end
end

function EmbeddedInstantonConfig(;
    center=nothing,
    radius=nothing,
    sign::Integer=1,
    block=(1, 2),
)
    typed_center = center === nothing ? nothing : Tuple(center)
    typed_block = Tuple(block)
    return EmbeddedInstantonConfig(
        typed_center,
        radius,
        sign,
        typed_block,
    )
end

function gauge_initialization_config(
    initial::AbstractString,
    load_format,
    initial_seed,
)
    normalized = lowercase(strip(initial))
    normalized == "cold" && return ColdStartConfig()
    normalized == "hot" && return HotStartConfig(initial_seed)
    normalized in ("one instanton", "instanton") && return InstantonConfig()
    normalized in ("embedded instanton", "sun embedded instanton") &&
        return EmbeddedInstantonConfig()
    return FileStartConfig(initial, load_format)
end

"""
Settings needed to construct a gauge field.

`halo` is a property of the field representation rather than the lattice
geometry, so it is deliberately kept outside [`LatticeConfig`](@ref). The
concrete initialization setting is retained as a type parameter.
"""
struct GaugeConfig{I<:AbstractGaugeInitializationConfig}
    NC::Int
    halo::Int
    initialization::I
end

function GaugeConfig(
    NC::Integer,
    halo::Integer,
    initialization::I,
) where {I<:AbstractGaugeInitializationConfig}
    NC > 0 || throw(ArgumentError("NC must be positive; got $NC"))
    halo >= 0 || throw(ArgumentError("halo must be nonnegative; got $halo"))
    return GaugeConfig{I}(Int(NC), Int(halo), initialization)
end

function GaugeConfig(
    NC::Integer,
    halo::Integer,
    initial::AbstractString,
    load_format,
    initial_seed=nothing,
)
    initialization = gauge_initialization_config(
        initial,
        load_format,
        initial_seed,
    )
    return GaugeConfig(NC, halo, initialization)
end

"""
One term in a gauge action.

`name` is a stable identifier used by force-group settings. `loop` may be a
built-in loop name or a typed explicit loop description. This keeps the
configuration extensible to rectangle and other improved-action terms without
introducing an action-specific configuration type. When `include_adjoint` is
true, `coupling` multiplies the real part of the loop sum; the runtime builder
adds both orientations with coefficient `coupling / 2`.
"""
struct GaugeActionTermConfig{
    N<:Union{Symbol,AbstractString},
    L,
    C<:Number,
}
    name::N
    loop::L
    coupling::C
    include_adjoint::Bool
end

GaugeActionTermConfig(name, loop, coupling) =
    GaugeActionTermConfig(name, loop, coupling, true)

"""A gauge action represented by a typed, heterogeneous tuple of terms."""
struct GaugeActionConfig{T<:Tuple}
    terms::T
end

function GaugeActionConfig(terms::GaugeActionTermConfig...)
    return GaugeActionConfig(terms)
end

"""Common supertype for serializable Dirac-operator settings."""
abstract type AbstractDiracOperatorConfig end

"""Wilson operator parameters."""
struct WilsonDiracConfig{K<:Real,R<:Real} <: AbstractDiracOperatorConfig
    hopping_parameter::K
    wilson_parameter::R

    function WilsonDiracConfig(
        hopping_parameter::K,
        wilson_parameter::R=one(K),
    ) where {K<:Real,R<:Real}
        hopping_parameter > 0 || throw(ArgumentError(
            "the Wilson hopping parameter must be positive; " *
            "got $hopping_parameter",
        ))
        return new{K,R}(hopping_parameter, wilson_parameter)
    end
end

"""Wilson--clover operator parameters."""
struct WilsonCloverDiracConfig{K<:Real,R<:Real,C<:Real} <:
       AbstractDiracOperatorConfig
    hopping_parameter::K
    wilson_parameter::R
    clover_coefficient::C

    function WilsonCloverDiracConfig(
        hopping_parameter::K,
        wilson_parameter::R,
        clover_coefficient::C,
    ) where {K<:Real,R<:Real,C<:Real}
        hopping_parameter > 0 || throw(ArgumentError(
            "the Wilson hopping parameter must be positive; " *
            "got $hopping_parameter",
        ))
        return new{K,R,C}(
            hopping_parameter,
            wilson_parameter,
            clover_coefficient,
        )
    end
end

"""Staggered operator parameters."""
struct StaggeredDiracConfig{M<:Real} <: AbstractDiracOperatorConfig
    mass::M

    function StaggeredDiracConfig(mass::M) where {M<:Real}
        mass > 0 || throw(ArgumentError(
            "the staggered mass must be positive; got $mass",
        ))
        return new{M}(mass)
    end
end

"""HISQ operator parameters, including the species-dependent Naik correction."""
struct HISQDiracConfig{M<:Real,E<:Real} <: AbstractDiracOperatorConfig
    mass::M
    naik_epsilon::E

    function HISQDiracConfig(
        mass::M,
        naik_epsilon::E=zero(M),
    ) where {M<:Real,E<:Real}
        mass > 0 || throw(ArgumentError(
            "the HISQ mass must be positive; got $mass",
        ))
        isfinite(naik_epsilon) || throw(ArgumentError(
            "the HISQ Naik correction must be finite; got $naik_epsilon",
        ))
        return new{M,E}(mass, naik_epsilon)
    end
end

"""Shamir domain-wall operator parameters."""
struct DomainwallDiracConfig{M<:Real,H<:Real} <:
       AbstractDiracOperatorConfig
    mass::M
    domainwall_height::H
    fifth_dimension::Int

    function DomainwallDiracConfig(
        mass::M,
        domainwall_height::H,
        fifth_dimension::Integer,
    ) where {M<:Real,H<:Real}
        fifth_dimension > 0 || throw(ArgumentError(
            "the domain-wall fifth dimension must be positive; " *
            "got $fifth_dimension",
        ))
        return new{M,H}(
            mass,
            domainwall_height,
            Int(fifth_dimension),
        )
    end
end

"""Scalar-coefficient Möbius domain-wall operator parameters."""
struct MobiusDomainwallDiracConfig{
    M<:Real,
    H<:Real,
    B<:Real,
    C<:Real,
} <: AbstractDiracOperatorConfig
    mass::M
    domainwall_height::H
    fifth_dimension::Int
    b::B
    c::C

    function MobiusDomainwallDiracConfig(
        mass::M,
        domainwall_height::H,
        fifth_dimension::Integer,
        b::B,
        c::C,
    ) where {M<:Real,H<:Real,B<:Real,C<:Real}
        fifth_dimension > 0 || throw(ArgumentError(
            "the Möbius domain-wall fifth dimension must be positive; " *
            "got $fifth_dimension",
        ))
        all(isfinite, (mass, domainwall_height, b, c)) ||
            throw(ArgumentError(
                "the Möbius domain-wall mass, M, b, and c must be finite",
            ))
        return new{M,H,B,C}(
            mass,
            domainwall_height,
            Int(fifth_dimension),
            b,
            c,
        )
    end
end

"""Linear-solver and fermion boundary-condition settings."""
struct FermionSolverConfig{T<:Real,B<:Tuple}
    tolerance::T
    max_steps::Int
    verbose::Int
    boundary_conditions::B

    function FermionSolverConfig(
        tolerance::T,
        max_steps::Integer,
        verbose::Integer,
        boundary_conditions::B,
    ) where {T<:Real,B<:Tuple}
        tolerance > 0 || throw(ArgumentError(
            "the fermion solver tolerance must be positive; got $tolerance",
        ))
        max_steps > 0 || throw(ArgumentError(
            "the fermion solver step limit must be positive; got $max_steps",
        ))
        isempty(boundary_conditions) && throw(ArgumentError(
            "fermion boundary conditions must not be empty",
        ))
        return new{T,B}(
            tolerance,
            Int(max_steps),
            Int(verbose),
            boundary_conditions,
        )
    end
end

abstract type AbstractFermionSmearingConfig end

"""Use the original gauge links in the fermion action."""
struct NoFermionSmearingConfig <: AbstractFermionSmearingConfig end

"""Apply the same typed stout layer specification repeatedly."""
struct StoutFermionSmearingConfig{R<:Tuple,L<:Tuple} <:
       AbstractFermionSmearingConfig
    layers::Int
    coefficients::R
    loops::L

    function StoutFermionSmearingConfig(
        layers::Integer,
        coefficients::R,
        loops::L,
    ) where {R<:Tuple,L<:Tuple}
        layers > 0 || throw(ArgumentError(
            "the number of stout layers must be positive; got $layers",
        ))
        isempty(coefficients) && throw(ArgumentError(
            "stout coefficients must not be empty",
        ))
        length(coefficients) == length(loops) || throw(ArgumentError(
            "stout coefficients and loops must have the same length",
        ))
        return new{R,L}(Int(layers), coefficients, loops)
    end
end

"""
One dynamical pseudofermion action term.

The operator, solver, and smearing fields are concrete type parameters. A
tuple of these terms therefore supports several fermion species without
placing an abstractly typed field inside `LQCDConfig`.
"""
struct FermionActionConfig{
    O<:AbstractDiracOperatorConfig,
    S<:FermionSolverConfig,
    M<:AbstractFermionSmearingConfig,
}
    name::Symbol
    operator::O
    flavors::Int
    solver::S
    smearing::M

    function FermionActionConfig(
        name::Symbol,
        operator::O,
        flavors::Integer,
        solver::S,
        smearing::M=NoFermionSmearingConfig(),
    ) where {
        O<:AbstractDiracOperatorConfig,
        S<:FermionSolverConfig,
        M<:AbstractFermionSmearingConfig,
    }
        name === :gauge && throw(ArgumentError(
            "the fermion action name :gauge is reserved",
        ))
        flavors > 0 || throw(ArgumentError(
            "the number of fermion flavors must be positive; got $flavors",
        ))
        if !(operator isa Union{StaggeredDiracConfig,HISQDiracConfig}) &&
           flavors != 2
            throw(ArgumentError(
                "Wilson-family and domain-wall pseudofermion actions " *
                "currently represent two degenerate flavors; got $flavors",
            ))
        end
        return new{O,S,M}(
            name,
            operator,
            Int(flavors),
            solver,
            smearing,
        )
    end
end

"""Common supertype for serializable update settings."""
abstract type AbstractUpdateConfig end

"""Common supertype for a serializable source of gauge configurations."""
abstract type AbstractConfigurationSourceConfig end

"""Read every matching configuration in a directory, in sorted order."""
struct DirectorySourceConfig{P<:AbstractString} <:
       AbstractConfigurationSourceConfig
    directory::P

    function DirectorySourceConfig(directory::P) where {P<:AbstractString}
        isempty(strip(directory)) && throw(ArgumentError(
            "a configuration directory must not be empty",
        ))
        return new{P}(directory)
    end
end

"""Read configurations in the order specified by a manifest file."""
struct ManifestSourceConfig{
    P<:AbstractString,
    M<:AbstractString,
} <: AbstractConfigurationSourceConfig
    directory::P
    manifest::M

    function ManifestSourceConfig(
        directory::P,
        manifest::M,
    ) where {P<:AbstractString,M<:AbstractString}
        isempty(strip(directory)) && throw(ArgumentError(
            "a configuration directory must not be empty",
        ))
        isempty(strip(manifest)) && throw(ArgumentError(
            "a configuration manifest must not be empty",
        ))
        return new{P,M}(directory, manifest)
    end
end

function canonical_configuration_format(format::Union{Symbol,AbstractString})
    normalized = lowercase(replace(
        strip(String(format)),
        "_" => "",
        "-" => "",
    ))
    normalized in ("jld", "jld2") && return :jld2
    normalized in ("bridge", "bridgetext") && return :bridge
    normalized == "ildg" && return :ildg
    throw(ArgumentError(
        "configuration format must be JLD2, BridgeText, or ILDG; got " *
        repr(format),
    ))
end

"""A finite, ordered sequence of gauge configurations to load and measure."""
struct ConfigurationSequenceConfig{
    S<:AbstractConfigurationSourceConfig,
    F<:Symbol,
} <: AbstractUpdateConfig
    source::S
    format::F

    function ConfigurationSequenceConfig(
        source::S,
        format::Union{Symbol,AbstractString},
    ) where {S<:AbstractConfigurationSourceConfig}
        canonical_format = canonical_configuration_format(format)
        return new{S,typeof(canonical_format)}(source, canonical_format)
    end
end

"""Common supertype for serializable MD-integrator settings."""
abstract type AbstractMDIntegratorConfig end

"""Select the `Q(1/2) P(1) Q(1/2)` leapfrog ordering."""
struct QPQConfig end

"""Select the `P(1/2) Q(1) P(1/2)` leapfrog ordering."""
struct PQPConfig end

"""Names of action components whose forces are evaluated together."""
struct ForceGroupConfig{N<:Tuple}
    names::N

    function ForceGroupConfig(names::N) where {N<:Tuple}
        isempty(names) && throw(ArgumentError("a force group must not be empty"))
        all(name -> name isa Symbol || name isa AbstractString, names) ||
            throw(ArgumentError("force-group names must be symbols or strings"))
        return new{N}(names)
    end
end

ForceGroupConfig(names...) = ForceGroupConfig(names)

"""Single-time-scale leapfrog settings for one force group."""
struct LeapfrogConfig{
    O<:Union{QPQConfig,PQPConfig},
    G<:ForceGroupConfig,
} <: AbstractMDIntegratorConfig
    ordering::O
    forces::G
end

"""
Nested multi-time-scale settings.

`slow_forces` are updated at the outer scale. `fast_integrator` describes the
inner evolution and `fast_steps` is its subdivision count per outer step.
"""
struct SextonWeingartenConfig{
    O<:Union{QPQConfig,PQPConfig},
    S<:ForceGroupConfig,
    F<:AbstractMDIntegratorConfig,
} <: AbstractMDIntegratorConfig
    ordering::O
    slow_forces::S
    fast_integrator::F
    fast_steps::Int

    function SextonWeingartenConfig(
        ordering::O,
        slow_forces::S,
        fast_integrator::F,
        fast_steps::Integer,
    ) where {
        O<:Union{QPQConfig,PQPConfig},
        S<:ForceGroupConfig,
        F<:AbstractMDIntegratorConfig,
    }
        fast_steps > 0 || throw(ArgumentError(
            "Sexton-Weingarten fast_steps must be positive; got $fast_steps",
        ))
        return new{O,S,F}(
            ordering,
            slow_forces,
            fast_integrator,
            Int(fast_steps),
        )
    end
end

"""Molecular-dynamics step size, step count, and integration scheme."""
struct MDConfig{T<:Real,I<:AbstractMDIntegratorConfig}
    step_size::T
    steps::Int
    integrator::I

    function MDConfig(
        step_size::T,
        steps::Integer,
        integrator::I,
    ) where {T<:Real,I<:AbstractMDIntegratorConfig}
        isfinite(step_size) || throw(ArgumentError(
            "MD step_size must be finite; got $step_size",
        ))
        step_size > 0 || throw(ArgumentError(
            "MD step_size must be positive; got $step_size",
        ))
        steps > 0 || throw(ArgumentError(
            "MD steps must be positive; got $steps",
        ))
        return new{T,I}(step_size, Int(steps), integrator)
    end
end

"""A reproducible logical random-number stream."""
struct RandomStreamConfig{S<:Integer,N}
    seed::S
    name::N
end

"""
Full Gaussian momentum refresh settings.

The MD momentum denominator is derived as `sigma^2`, so the Gaussian refresh,
kinetic energy, and momentum kicks always use one consistent normalization.
"""
struct GaussianMomentumConfig{T<:Real,R<:RandomStreamConfig}
    sigma::T
    random::R

    function GaussianMomentumConfig(
        sigma::T,
        random::R,
    ) where {T<:Real,R<:RandomStreamConfig}
        isfinite(sigma) || throw(ArgumentError(
            "momentum sigma must be finite; got $sigma",
        ))
        sigma > 0 || throw(ArgumentError(
            "momentum sigma must be positive; got $sigma",
        ))
        return new{T,R}(sigma, random)
    end
end

function momentum_denominator(config::GaussianMomentumConfig)
    sigma = config.sigma
    two = one(sigma) + one(sigma)
    return sigma == sqrt(two) ? two : abs2(sigma)
end

"""Metropolis decisions drawn on rank zero and broadcast to other ranks."""
struct RankZeroMetropolisConfig{R<:RandomStreamConfig}
    random::R
end

"""Pseudofermion refresh stream associated with one named action term."""
struct PseudofermionRefreshConfig{R<:RandomStreamConfig}
    action_name::Symbol
    random::R
    subgroup::Int

    function PseudofermionRefreshConfig(
        action_name::Symbol,
        random::R,
        subgroup::Integer,
    ) where {R<:RandomStreamConfig}
        action_name === :gauge && throw(ArgumentError(
            "the pseudofermion action name :gauge is reserved",
        ))
        subgroup >= 0 || throw(ArgumentError(
            "the pseudofermion RNG subgroup must be nonnegative; " *
            "got $subgroup",
        ))
        return new{R}(action_name, random, Int(subgroup))
    end
end

"""Hybrid Monte Carlo update settings."""
struct HMCConfig{
    M<:MDConfig,
    P<:GaussianMomentumConfig,
    A<:RankZeroMetropolisConfig,
    F<:Tuple,
} <: AbstractUpdateConfig
    md::M
    momentum::P
    acceptance::A
    pseudofermions::F
end

HMCConfig(md, momentum, acceptance) =
    HMCConfig(md, momentum, acceptance, ())

"""
Self-learning HMC settings with a distinct molecular-dynamics action.

The target gauge and fermion actions remain in [`LQCDConfig`](@ref).  The
`md_gauge_action` and `md_fermions` fields describe only the approximate
action used to generate the reversible MD proposal.  Pseudofermions are
refreshed from the target action and the Metropolis decision is evaluated
with the target Hamiltonian.
"""
struct SLHMCConfig{
    M<:MDConfig,
    P<:GaussianMomentumConfig,
    A<:RankZeroMetropolisConfig,
    F<:Tuple,
    G<:GaugeActionConfig,
    D<:Tuple,
} <: AbstractUpdateConfig
    md::M
    momentum::P
    acceptance::A
    pseudofermions::F
    md_gauge_action::G
    md_fermions::D
end

SLHMCConfig(md, momentum, acceptance, md_gauge_action, md_fermions=()) =
    SLHMCConfig(
        md,
        momentum,
        acceptance,
        (),
        md_gauge_action,
        md_fermions,
    )

"""Heatbath and optional overrelaxation update settings."""
struct HeatbathConfig{R<:RandomStreamConfig} <: AbstractUpdateConfig
    even_odd::Bool
    max_iterations::Int
    overrelaxation_steps::Int
    random::R

    function HeatbathConfig(
        even_odd::Bool,
        max_iterations::Integer,
        overrelaxation_steps::Integer,
        random::R,
    ) where {R<:RandomStreamConfig}
        max_iterations > 0 || throw(ArgumentError(
            "Heatbath max_iterations must be positive; got $max_iterations",
        ))
        overrelaxation_steps >= 0 || throw(ArgumentError(
            "overrelaxation_steps must be nonnegative; " *
            "got $overrelaxation_steps",
        ))
        return new{R}(
            even_odd,
            Int(max_iterations),
            Int(overrelaxation_steps),
            random,
        )
    end
end

function fermion_smearing_config(parameters)
    normalized = lowercase(strip(parameters.smearing_for_fermion))
    normalized in ("nothing", "none", "no smearing") &&
        return NoFermionSmearingConfig()
    normalized == "stout" || throw(ArgumentError(
        "unsupported fermion smearing " * repr(parameters.smearing_for_fermion),
    ))
    parameters.stout_numlayers === nothing && throw(ArgumentError(
        "stout_numlayers is required for stout fermion smearing",
    ))
    parameters.stout_ρ === nothing && throw(ArgumentError(
        "stout_ρ is required for stout fermion smearing",
    ))
    parameters.stout_loops === nothing && throw(ArgumentError(
        "stout_loops is required for stout fermion smearing",
    ))
    return StoutFermionSmearingConfig(
        parameters.stout_numlayers,
        Tuple(parameters.stout_ρ),
        Tuple(parameters.stout_loops),
    )
end

function fermion_operator_config(parameters)
    operator_name = lowercase(strip(String(parameters.Dirac_operator)))
    operator_name == "wilson" &&
        return WilsonDiracConfig(parameters.hop, parameters.r)
    operator_name == "wilsonclover" && return WilsonCloverDiracConfig(
        parameters.hop,
        parameters.r,
        parameters.Clover_coefficient,
    )
    operator_name == "staggered" &&
        return StaggeredDiracConfig(parameters.mass)
    operator_name == "hisq" && return HISQDiracConfig(
        parameters.mass,
        parameters.naik_epsilon,
    )
    if operator_name == "domainwall"
        parameters.Domainwall_m === nothing && throw(ArgumentError(
            "Domainwall_m is required for a domain-wall action",
        ))
        parameters.Domainwall_M === nothing && throw(ArgumentError(
            "Domainwall_M is required for a domain-wall action",
        ))
        parameters.Domainwall_L5 === nothing && throw(ArgumentError(
            "Domainwall_L5 is required for a domain-wall action",
        ))
        return DomainwallDiracConfig(
            parameters.Domainwall_m,
            parameters.Domainwall_M,
            parameters.Domainwall_L5,
        )
    end
    if operator_name in ("mobiusdomainwall", "mobius_domainwall")
        parameters.Domainwall_m === nothing && throw(ArgumentError(
            "Domainwall_m is required for a Möbius domain-wall action",
        ))
        parameters.Domainwall_M === nothing && throw(ArgumentError(
            "Domainwall_M is required for a Möbius domain-wall action",
        ))
        parameters.Domainwall_L5 === nothing && throw(ArgumentError(
            "Domainwall_L5 is required for a Möbius domain-wall action",
        ))
        return MobiusDomainwallDiracConfig(
            parameters.Domainwall_m,
            parameters.Domainwall_M,
            parameters.Domainwall_L5,
            parameters.b,
            parameters.c,
        )
    end
    throw(ArgumentError(
        "unsupported Dirac_operator=" * repr(parameters.Dirac_operator),
    ))
end

function fermion_action_configs(parameters)
    (parameters.quench || isnothing(parameters.Dirac_operator)) && return ()
    solver = FermionSolverConfig(
        parameters.eps,
        parameters.MaxCGstep,
        parameters.verboselevel,
        Tuple(parameters.BoundaryCondition),
    )
    operator = fermion_operator_config(parameters)
    flavors = operator isa Union{StaggeredDiracConfig,HISQDiracConfig} ?
        parameters.Nf : 2
    smearing = fermion_smearing_config(parameters)
    return (
        FermionActionConfig(
            :fermion_1,
            operator,
            flavors,
            solver,
            smearing,
        ),
    )
end

fermion_action_names(fermions::Tuple) = map(config -> config.name, fermions)

function validate_fermion_action_names(fermions::Tuple)
    names = fermion_action_names(fermions)
    length(unique(names)) == length(names) || throw(ArgumentError(
        "fermion action names must be unique; got $names",
    ))
    return names
end

"""Convert legacy update fields in `Params` to a typed update setting."""
function update_config(parameters, fermions::Tuple=())
    seed = parameters.randomseed
    if parameters.update_method in ("HMC", "SLHMC", "SLMC")
        fermion_names = validate_fermion_action_names(fermions)
        if parameters.SextonWeingargten && isempty(fermions)
            @warn "Sexton-Weingarten is ignored for a gauge-only " *
                  "simulation because it has only one force group; " *
                  "using the configured leapfrog ordering"
        end
        ordering = parameters.QPQ ? QPQConfig() : PQPConfig()
        integrator = if parameters.SextonWeingargten && !isempty(fermions)
            iseven(parameters.N_SextonWeingargten) || throw(ArgumentError(
                "legacy N_SextonWeingargten must be even; " *
                "got $(parameters.N_SextonWeingargten)",
            ))
            SextonWeingartenConfig(
                ordering,
                ForceGroupConfig(fermion_names),
                LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge)),
                parameters.N_SextonWeingargten ÷ 2,
            )
        else
            LeapfrogConfig(
                ordering,
                ForceGroupConfig((:gauge, fermion_names...)),
            )
        end
        md = MDConfig(parameters.Δτ, parameters.MDsteps, integrator)
        momentum = GaussianMomentumConfig(
            sqrt(2.0),
            RandomStreamConfig(seed, :momentum),
        )
        acceptance = RankZeroMetropolisConfig(
            RandomStreamConfig(seed, :metropolis),
        )
        pseudofermions = ntuple(length(fermion_names)) do index
            name = fermion_names[index]
            PseudofermionRefreshConfig(
                name,
                RandomStreamConfig(seed, :pseudofermion),
                index,
            )
        end
        if parameters.update_method == "HMC"
            return HMCConfig(md, momentum, acceptance, pseudofermions)
        end
        parameters.βeff isa Real || throw(ArgumentError(
            "typed SLHMC currently requires a scalar βeff; got " *
            repr(parameters.βeff),
        ))
        isfinite(parameters.βeff) || throw(ArgumentError(
            "SLHMC βeff must be finite; got $(parameters.βeff)",
        ))
        md_gauge_action = GaugeActionConfig(
            GaugeActionTermConfig(
                :gauge_plaquette,
                "plaquette",
                parameters.βeff,
            ),
        )
        return SLHMCConfig(
            md,
            momentum,
            acceptance,
            pseudofermions,
            md_gauge_action,
            fermions,
        )
    elseif parameters.update_method == "Heatbath"
        isempty(fermions) || throw(ArgumentError(
            "heatbath updates support gauge-only simulations",
        ))
        overrelaxation_steps = parameters.useOR ? parameters.numOR : 0
        return HeatbathConfig(
            parameters.isevenodd,
            parameters.ITERATION_MAX,
            overrelaxation_steps,
            RandomStreamConfig(seed, :heatbath),
        )
    elseif parameters.update_method == "Fileloading"
        isempty(fermions) || throw(ArgumentError(
            "Fileloading is a gauge-configuration sequence and does not " *
            "construct dynamical pseudofermions",
        ))
        isnothing(parameters.loadU_format) && throw(ArgumentError(
            "loadU_format is required for update_method=\"Fileloading\"",
        ))
        source = if parameters.loadU_fromfile
            ManifestSourceConfig(
                parameters.loadU_dir,
                parameters.loadU_filename,
            )
        else
            DirectorySourceConfig(parameters.loadU_dir)
        end
        return ConfigurationSequenceConfig(source, parameters.loadU_format)
    end

    throw(ArgumentError(
        "LQCDConfig does not support update_method=" *
        repr(parameters.update_method),
    ))
end

"""
Serializable input needed to construct a simulation.

The fields contain settings only. Runtime gauge and fermion fields, action
objects, communicators, open streams, and workspaces belong to `Simulation`
and are not stored here.
"""
struct LQCDConfig{
    L<:LatticeConfig,
    G<:GaugeConfig,
    A<:GaugeActionConfig,
    F<:Tuple,
    U<:AbstractUpdateConfig,
}
    lattice::L
    gauge::G
    gauge_action::A
    fermions::F
    update::U
end

LQCDConfig(lattice, gauge, gauge_action, update) =
    LQCDConfig(lattice, gauge, gauge_action, (), update)

"""
    LQCDConfig(parameters::Params)

Extract a typed gauge and dynamical-fermion configuration from legacy
`Params`. Quenched inputs use an empty fermion tuple.
"""
function lqcd_config(parameters)
    lattice = LatticeConfig(Tuple(Int.(parameters.L)))
    fermions = fermion_action_configs(parameters)
    load_format = if isnothing(parameters.loadU_format)
        nothing
    else
        String(parameters.loadU_format)
    end
    halo = parameters.Nwing
    minimum_halo = any(
        fermion -> fermion.operator isa HISQDiracConfig,
        fermions,
    ) ? 3 : (isempty(fermions) ? 0 : 1)
    if halo < minimum_halo
        replacement = minimum_halo == 1 ?
                      "one gauge-field halo layer" :
                      "$minimum_halo gauge-field halo layers"
        @warn "Nwing=$halo in a legacy dynamical-fermion input is not " *
              "supported by LDO v1 for this operator; using $replacement"
        halo = minimum_halo
    end
    gauge = GaugeConfig(
        parameters.NC,
        halo,
        String(parameters.initial),
        load_format,
        parameters.randomseed,
    )
    gauge_action = GaugeActionConfig(
        GaugeActionTermConfig(:gauge_plaquette, "plaquette", parameters.β),
    )
    update = update_config(parameters, fermions)

    return LQCDConfig(lattice, gauge, gauge_action, fermions, update)
end

LQCDConfig(parameters::Params) = lqcd_config(parameters)

md_trajectory_length(config::MDConfig) = config.step_size * config.steps

function print_integrator_config(io::IO, config::LeapfrogConfig, indent::Int)
    padding = repeat(" ", indent)
    ordering = config.ordering isa QPQConfig ? "QPQ" : "PQP"
    println(io, padding, "scheme: ", ordering)
    print(io, padding, "forces: ")
    show(io, config.forces.names)
    println(io)
    return nothing
end

function print_integrator_config(
    io::IO,
    config::SextonWeingartenConfig,
    indent::Int,
)
    padding = repeat(" ", indent)
    ordering = config.ordering isa QPQConfig ? "QPQ" : "PQP"
    println(io, padding, "scheme: Sexton-Weingarten ", ordering)
    print(io, padding, "slow forces: ")
    show(io, config.slow_forces.names)
    println(io)
    println(io, padding, "fast steps: ", config.fast_steps)
    println(io, padding, "fast integrator:")
    print_integrator_config(io, config.fast_integrator, indent + 2)
    return nothing
end

function show_config(io::IO, config::HMCConfig)
    println(io, "HMC")
    println(io, "  MD step size: ", config.md.step_size)
    println(io, "  MD steps: ", config.md.steps)
    println(io, "  trajectory length: ", md_trajectory_length(config.md))
    println(io, "  integrator:")
    print_integrator_config(io, config.md.integrator, 4)
    println(io, "  momentum sigma: ", config.momentum.sigma)
    println(
        io,
        "  momentum denominator: ",
        momentum_denominator(config.momentum),
    )
    println(io, "  momentum seed: ", config.momentum.random.seed)
    println(io, "  momentum stream: ", config.momentum.random.name)
    if !isempty(config.pseudofermions)
        println(io, "  pseudofermion refreshes:")
        for refresh in config.pseudofermions
            println(
                io,
                "    ",
                refresh.action_name,
                ": seed=",
                refresh.random.seed,
                ", subgroup=",
                refresh.subgroup,
            )
        end
    end
    println(io, "  Metropolis policy: rank-zero broadcast")
    println(io, "  Metropolis seed: ", config.acceptance.random.seed)
    println(io, "  Metropolis stream: ", config.acceptance.random.name)
    return nothing
end

function show_config(io::IO, config::SLHMCConfig)
    println(io, "SLHMC")
    println(io, "  target action: LQCDConfig gauge_action and fermions")
    println(io, "  MD gauge-action terms: ", length(config.md_gauge_action.terms))
    print(io, "  MD fermion actions: ")
    show(io, fermion_action_names(config.md_fermions))
    println(io)
    println(io, "  MD step size: ", config.md.step_size)
    println(io, "  MD steps: ", config.md.steps)
    println(io, "  trajectory length: ", md_trajectory_length(config.md))
    println(io, "  integrator:")
    print_integrator_config(io, config.md.integrator, 4)
    println(io, "  momentum sigma: ", config.momentum.sigma)
    println(
        io,
        "  momentum denominator: ",
        momentum_denominator(config.momentum),
    )
    println(io, "  momentum seed: ", config.momentum.random.seed)
    if !isempty(config.pseudofermions)
        println(io, "  target pseudofermion refreshes:")
        for refresh in config.pseudofermions
            println(
                io,
                "    ",
                refresh.action_name,
                ": seed=",
                refresh.random.seed,
                ", subgroup=",
                refresh.subgroup,
            )
        end
    end
    println(io, "  Metropolis action: target")
    println(io, "  Metropolis policy: rank-zero broadcast")
    println(io, "  Metropolis seed: ", config.acceptance.random.seed)
    return nothing
end

function show_config(io::IO, config::HeatbathConfig)
    println(io, "Heatbath")
    println(io, "  even-odd: ", config.even_odd)
    println(io, "  maximum iterations: ", config.max_iterations)
    println(io, "  overrelaxation steps: ", config.overrelaxation_steps)
    println(io, "  random seed: ", config.random.seed)
    println(io, "  random stream: ", config.random.name)
    return nothing
end

function print_configuration_source(
    io::IO,
    source::DirectorySourceConfig,
)
    println(io, "directory scan")
    print(io, "  directory: ")
    show(io, source.directory)
    println(io)
    println(io, "  ordering: sorted")
    return nothing
end

function print_configuration_source(
    io::IO,
    source::ManifestSourceConfig,
)
    println(io, "manifest")
    print(io, "  directory: ")
    show(io, source.directory)
    println(io)
    print(io, "  manifest: ")
    show(io, source.manifest)
    println(io)
    println(io, "  ordering: manifest order")
    return nothing
end

function show_config(io::IO, config::ConfigurationSequenceConfig)
    println(io, "Configuration sequence")
    println(io, "  format: ", config.format)
    println(io, "  source:")
    source_output = sprint(print_configuration_source, config.source)
    for line in eachline(IOBuffer(source_output))
        println(io, "    ", line)
    end
    return nothing
end

function print_initialization_config(io::IO, ::ColdStartConfig)
    println(io, "cold start")
    return nothing
end

function print_initialization_config(io::IO, config::HotStartConfig)
    println(io, "hot start")
    print(io, "seed: ")
    show(io, config.seed)
    println(io)
    return nothing
end

function print_initialization_config(io::IO, config::FileStartConfig)
    println(io, "file start")
    print(io, "path: ")
    show(io, config.path)
    println(io)
    print(io, "format: ")
    show(io, config.format)
    println(io)
    return nothing
end

function print_initialization_config(io::IO, ::InstantonConfig)
    println(io, "SU(2) one-instanton start")
    println(io, "radius: half the x extent")
    return nothing
end

function print_initialization_config(
    io::IO,
    config::EmbeddedInstantonConfig,
)
    println(io, "embedded SU(2) instanton start")
    print(io, "center: ")
    show(io, config.center)
    println(io)
    print(io, "radius: ")
    show(io, config.radius)
    println(io)
    println(io, "sign: ", config.sign)
    print(io, "color block: ")
    show(io, config.block)
    println(io)
    return nothing
end

"""Print a human-readable representation of a typed LQCD configuration."""
function show_config(io::IO, config::LQCDConfig)
    println(io, "LQCDConfig")
    print(io, "  lattice: ")
    show(io, config.lattice.L)
    println(io)
    println(io, "  colors: ", config.gauge.NC)
    println(io, "  halo: ", config.gauge.halo)
    println(io, "  initialization:")
    initialization_output = sprint(
        print_initialization_config,
        config.gauge.initialization,
    )
    for line in eachline(IOBuffer(initialization_output))
        println(io, "    ", line)
    end
    println(io, "  gauge action:")
    for term in config.gauge_action.terms
        print(io, "    ", term.name, ": loop=")
        show(io, term.loop)
        println(
            io,
            ", coupling=",
            term.coupling,
            ", include adjoint=",
            term.include_adjoint,
        )
    end
    println(io, "  fermion actions:")
    if isempty(config.fermions)
        println(io, "    none")
    else
        for fermion in config.fermions
            println(
                io,
                "    ",
                fermion.name,
                ": operator=",
                nameof(typeof(fermion.operator)),
                ", flavors=",
                fermion.flavors,
            )
            println(
                io,
                "      solver tolerance=",
                fermion.solver.tolerance,
                ", max steps=",
                fermion.solver.max_steps,
            )
            println(
                io,
                "      smearing=",
                nameof(typeof(fermion.smearing)),
            )
        end
    end
    println(io, "  update:")
    update_output = sprint(show_config, config.update)
    for line in eachline(IOBuffer(update_output))
        println(io, "    ", line)
    end
    return nothing
end

show_config(config::Union{LQCDConfig,AbstractUpdateConfig}) =
    show_config(stdout, config)

function Base.show(io::IO, config::LQCDConfig)
    print(io, "LQCDConfig(lattice=")
    show(io, config.lattice.L)
    print(io, ", colors=", config.gauge.NC, ", update=")
    print(io, nameof(typeof(config.update)), ")")
end

function Base.show(io::IO, ::MIME"text/plain", config::LQCDConfig)
    return show_config(io, config)
end

end
