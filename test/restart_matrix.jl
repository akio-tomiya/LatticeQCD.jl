using LatticeQCD
using Test
import Gaugefields

const RESTART_TEST_ENVIRONMENT = GaugefieldsEnvironment(
    process_grid=(1, 1, 1, 1),
    communicator=Gaugefields.SerialCommunicator(),
    element_type=ComplexF64,
    verbose=0,
)

function restart_test_input(
    operator;
    colors=2,
    flavors=2,
    smearing=NoFermionSmearingConfig(),
    slhmc=false,
)
    uses_hisq = operator isa HISQDiracConfig
    lattice = LatticeConfig(uses_hisq ? (4, 4, 4, 4) : (2, 2, 2, 2))
    gauge = GaugeConfig(colors, uses_hisq ? 3 : 1, "cold", nothing, 0x1234)
    target_gauge_action = GaugeActionConfig(
        GaugeActionTermConfig(:gauge_plaquette, "plaquette", 1.9),
    )
    solver = FermionSolverConfig(1e-10, 2_000, 0, (1, 1, 1, -1))
    fermion = FermionActionConfig(
        :fermion_1,
        operator,
        flavors,
        solver,
        smearing,
    )
    fermions = (fermion,)
    md = MDConfig(
        0.001,
        1,
        LeapfrogConfig(
            QPQConfig(),
            ForceGroupConfig(:gauge, :fermion_1),
        ),
    )
    momentum = GaussianMomentumConfig(
        1.0,
        RandomStreamConfig(0x5678, :momentum),
    )
    acceptance = RankZeroMetropolisConfig(
        RandomStreamConfig(0x9abc, :metropolis),
    )
    refreshes = (
        PseudofermionRefreshConfig(
            :fermion_1,
            RandomStreamConfig(0xdef0, :pseudofermion),
            1,
        ),
    )
    update = if slhmc
        md_fermion = FermionActionConfig(
            fermion.name,
            WilsonDiracConfig(0.04, 1.0),
            fermion.flavors,
            fermion.solver,
            fermion.smearing,
        )
        md_gauge_action = GaugeActionConfig(
            GaugeActionTermConfig(:gauge_plaquette, "plaquette", 0.7),
        )
        SLHMCConfig(
            md,
            momentum,
            acceptance,
            refreshes,
            md_gauge_action,
            (md_fermion,),
        )
    else
        HMCConfig(md, momentum, acceptance, refreshes)
    end
    return LQCDConfig(
        lattice,
        gauge,
        target_gauge_action,
        fermions,
        update,
    )
end

function test_exact_fermion_restart(config::LQCDConfig)
    schedule = SimulationSchedule(0, 2)
    spec = SimulationSpec(config, schedule)
    reference = build_simulation(spec, RESTART_TEST_ENVIRONMENT)
    step!(reference)

    mktempdir() do directory
        path = joinpath(directory, "restart.jld2")
        save_checkpoint(path, reference)
        expected = step!(reference)

        restored = build_simulation(spec, RESTART_TEST_ENVIRONMENT)
        load_checkpoint!(restored, path)
        actual = step!(restored)

        @test actual.phase == expected.phase
        @test typeof(actual.update) == typeof(expected.update)
        for field in fieldnames(typeof(expected.update))
            @test getfield(actual.update, field) ==
                  getfield(expected.update, field)
        end
        @test restored.simulation.state.trajectory == 2
        @test restored.simulation.state.accepted ==
              reference.simulation.state.accepted
        for direction in eachindex(reference.simulation.configuration.gauge)
            expected_values = Array(
                reference.simulation.configuration.gauge[direction].U.A,
            )
            actual_values = Array(
                restored.simulation.configuration.gauge[direction].U.A,
            )
            @test actual_values == expected_values
        end
    end
end

if get(ENV, "LQCD_DEFINE_RESTART_HELPERS_ONLY", "0") != "1"
@testset "Exact restart matrix for dynamical fermions" begin
    cases = (
        (
            "Wilson clover Nf=2",
            restart_test_input(WilsonCloverDiracConfig(0.05, 1.0, 1.0)),
        ),
        (
            "staggered RHMC Nf=2",
            restart_test_input(StaggeredDiracConfig(0.5); flavors=2),
        ),
        (
            "HISQ RHMC Nf=2",
            restart_test_input(HISQDiracConfig(0.5, -0.083); flavors=2),
        ),
        (
            "domain-wall Nf=2",
            restart_test_input(DomainwallDiracConfig(0.1, -1.0, 2)),
        ),
        (
            "Mobius domain-wall Nf=2",
            restart_test_input(
                MobiusDomainwallDiracConfig(0.1, -1.0, 2, 2.0, 1.0),
            ),
        ),
        (
            "one-layer stout Wilson Nf=2",
            restart_test_input(
                WilsonDiracConfig(0.05, 1.0);
                smearing=StoutFermionSmearingConfig(
                    1,
                    (0.1,),
                    ("plaquette",),
                ),
            ),
        ),
        (
            "fermionic SLHMC",
            restart_test_input(WilsonDiracConfig(0.05, 1.0); slhmc=true),
        ),
    )

    requested_case = get(ENV, "LQCD_RESTART_CASE", "")
    selected_cases = isempty(requested_case) ? cases : filter(
        case -> occursin(lowercase(requested_case), lowercase(first(case))),
        cases,
    )
    isempty(selected_cases) && error(
        "LQCD_RESTART_CASE=$(repr(requested_case)) matched no restart case",
    )
    for (name, config) in selected_cases
        @testset "$name" begin
            test_exact_fermion_restart(config)
        end
    end
end
end
