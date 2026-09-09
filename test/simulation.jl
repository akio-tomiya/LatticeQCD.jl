using LatticeQCD
using Test
import Gaugefields

function test_hmc_input(;
    initial="hot",
    ordering=QPQConfig(),
    lattice_size=(2, 2, 2, 2),
    colors=2,
    halo=1,
    momentum_sigma=sqrt(2.0),
)
    lattice = LatticeConfig(lattice_size)
    gauge = if initial isa AbstractGaugeInitializationConfig
        GaugeConfig(colors, halo, initial)
    else
        GaugeConfig(colors, halo, initial, nothing, 0x1234)
    end
    action = GaugeActionConfig(
        GaugeActionTermConfig(
            :gauge_plaquette,
            "plaquette",
            1.9,
        ),
    )
    integrator = LeapfrogConfig(ordering, ForceGroupConfig(:gauge))
    md = MDConfig(0.02, 2, integrator)
    momentum = GaussianMomentumConfig(
        momentum_sigma,
        RandomStreamConfig(0x5678, :momentum),
    )
    acceptance = RankZeroMetropolisConfig(
        RandomStreamConfig(0x9abc, :metropolis),
    )
    update = HMCConfig(md, momentum, acceptance)
    return LQCDConfig(lattice, gauge, action, update)
end

function test_heatbath_input(;
    initial="cold",
    lattice_size=(2, 2, 2, 2),
    colors=2,
    beta=colors == 2 ? 1.9 : 5.7,
    even_odd=true,
    overrelaxation_steps=0,
    max_iterations=100_000,
)
    lattice = LatticeConfig(lattice_size)
    gauge = GaugeConfig(colors, 1, initial, nothing, 0x1234)
    action = GaugeActionConfig(
        GaugeActionTermConfig(
            :gauge_plaquette,
            "plaquette",
            beta,
        ),
    )
    update = HeatbathConfig(
        even_odd,
        max_iterations,
        overrelaxation_steps,
        RandomStreamConfig(0x5678, :heatbath),
    )
    return LQCDConfig(lattice, gauge, action, update)
end

function test_same_bridge_configuration(left, right)
    mktempdir() do directory
        left_path = joinpath(directory, "left.txt")
        right_path = joinpath(directory, "right.txt")
        Gaugefields.save_configuration(
            left_path,
            left.gauge;
            format=:bridge,
        )
        Gaugefields.save_configuration(
            right_path,
            right.gauge;
            format=:bridge,
        )
        @test read(left_path) == read(right_path)
    end
    return nothing
end

@testset "Gaugefields v1 embedded-instanton initialization" begin
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        element_type=ComplexF64,
        verbose=0,
    )
    initialization = EmbeddedInstantonConfig(
        center=(1.5, 1.5, 1.5, 1.5),
        radius=1.0,
        sign=-1,
        block=(2, 3),
    )
    input = test_hmc_input(
        initial=initialization,
        colors=3,
    )
    configuration = build_configuration(input, environment)
    gauge = configuration.gauge

    @test Gaugefields.gauge_backend(gauge) isa
          Gaugefields.LatticeMatricesBackend
    @test Gaugefields.gauge_lattice_size(gauge) == (2, 2, 2, 2)
    @test Gaugefields.gauge_num_colors(gauge) == 3
    @test Gaugefields.gauge_halo_width(gauge) == 1
    @test Gaugefields.gauge_process_grid(gauge) == (1, 1, 1, 1)
    @test isfinite(Gaugefields.measure_plaquette(gauge))

    default_configuration = build_configuration(
        test_hmc_input(initial="embedded instanton", colors=3),
        environment,
    )
    @test isfinite(Gaugefields.measure_plaquette(
        default_configuration.gauge,
    ))

    legacy_configuration = build_configuration(
        input,
        GaugefieldsEnvironment(
            backend=Gaugefields.LegacyBackend(),
            verbose=0,
        ),
    )
    @test Gaugefields.gauge_backend(legacy_configuration.gauge) isa
          Gaugefields.LegacyBackend
    @test isfinite(Gaugefields.measure_plaquette(legacy_configuration.gauge))

    @test_throws ArgumentError build_configuration(
        test_hmc_input(
            initial=initialization,
            colors=3,
            lattice_size=(2, 2, 2),
        ),
        GaugefieldsEnvironment(process_grid=(1, 1, 1)),
    )
    @test_throws ArgumentError build_configuration(
        test_hmc_input(
            initial=EmbeddedInstantonConfig(block=(1, 4)),
            colors=3,
        ),
        environment,
    )
    @test_throws ArgumentError build_configuration(
        input,
        GaugefieldsEnvironment(
            process_grid=(1, 1, 1, 1),
            element_type=Float64,
        ),
    )
end

@testset "Gaugefields v1 one-instanton initialization" begin
    input = test_hmc_input(initial="one instanton")
    lm_environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        element_type=ComplexF64,
        verbose=0,
    )
    lm_configuration = build_configuration(input, lm_environment)
    lm_gauge = lm_configuration.gauge

    @test Gaugefields.gauge_backend(lm_gauge) isa
          Gaugefields.LatticeMatricesBackend
    @test Gaugefields.gauge_lattice_size(lm_gauge) == (2, 2, 2, 2)
    @test Gaugefields.gauge_num_colors(lm_gauge) == 2
    @test Gaugefields.gauge_halo_width(lm_gauge) == 1
    @test Gaugefields.gauge_process_grid(lm_gauge) == (1, 1, 1, 1)
    @test isfinite(Gaugefields.measure_plaquette(lm_gauge))

    legacy_environment = GaugefieldsEnvironment(
        backend=Gaugefields.LegacyBackend(),
        verbose=0,
    )
    legacy_configuration = build_configuration(input, legacy_environment)
    @test Gaugefields.gauge_backend(legacy_configuration.gauge) isa
          Gaugefields.LegacyBackend
    @test isfinite(Gaugefields.measure_plaquette(legacy_configuration.gauge))

    @test_throws ArgumentError build_configuration(
        test_hmc_input(initial="one instanton", colors=3),
        lm_environment,
    )
    @test_throws ArgumentError build_configuration(
        test_hmc_input(
            initial="one instanton",
            lattice_size=(2, 2, 2),
        ),
        GaugefieldsEnvironment(process_grid=(1, 1, 1)),
    )
    @test_throws ArgumentError build_configuration(
        input,
        GaugefieldsEnvironment(
            process_grid=(1, 1, 1, 1),
            element_type=Float64,
        ),
    )
end

@testset "Gaugefields v1 pure-gauge simulation" begin
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        element_type=ComplexF64,
        verbose=0,
    )
    input = test_hmc_input()
    simulation = build_simulation(input, environment)
    gauge = simulation.configuration.gauge

    @test simulation.configuration isa GaugeConfiguration
    @test Gaugefields.gauge_backend(gauge) isa
          Gaugefields.LatticeMatricesBackend
    @test Gaugefields.gauge_lattice_size(gauge) == (2, 2, 2, 2)
    @test Gaugefields.gauge_num_colors(gauge) == 2
    @test Gaugefields.gauge_halo_width(gauge) == 1
    @test Gaugefields.gauge_process_grid(gauge) == (1, 1, 1, 1)
    @test simulation.action.dataset[1].β == 1.9 / 2
    @test simulation.updater.md_driver.integrator isa Gaugefields.QPQ
    @test simulation.updater.md_driver.steps == 2
    @test simulation.updater.md_driver.trajectory_length == 0.04
    @test simulation.updater.md_driver.momentum_denominator == 2.0
    @test fieldtype(typeof(simulation), :configuration) ===
          typeof(simulation.configuration)
    @test fieldtype(typeof(simulation), :updater) === typeof(simulation.updater)
    @test all(isconcretetype, fieldtypes(typeof(simulation)))
    @test occursin("updater=HMCUpdater", sprint(show, simulation))

    initial_plaquette = Gaugefields.measure_plaquette(gauge)
    result = update!(simulation)
    @test result isa HMCUpdateResult
    @test result.trajectory == 0
    @test isfinite(result.initial_hamiltonian)
    @test isfinite(result.final_hamiltonian)
    @test isfinite(result.delta_hamiltonian)
    @test simulation.state.trajectory == 1
    @test simulation.state.accepted == Int(result.accepted)
    @test Gaugefields.measure_plaquette(gauge) != initial_plaquette ||
          !result.accepted

    repeated = build_simulation(input, environment)
    repeated_result = update!(repeated)
    @test repeated_result.accepted == result.accepted
    @test repeated_result.delta_hamiltonian ≈ result.delta_hamiltonian
    @test Gaugefields.measure_plaquette(repeated.configuration.gauge) ≈
          Gaugefields.measure_plaquette(gauge)

    pqp = build_simulation(test_hmc_input(ordering=PQPConfig()), environment)
    @test pqp.updater.md_driver.integrator isa Gaugefields.PQP

    ltk = build_simulation(
        test_hmc_input(momentum_sigma=1.0),
        environment,
    )
    @test ltk.updater.md_driver.momentum_denominator == 1.0
end

@testset "HMC acceptance and device-local rollback" begin
    environment = GaugefieldsEnvironment(process_grid=(1, 1, 1, 1))
    cold = build_configuration(test_hmc_input(initial="cold"), environment)
    hot = build_configuration(test_hmc_input(initial="hot"), environment)
    backup = GaugeConfiguration([similar(link) for link in cold.gauge])
    copy_configuration!(backup, cold)

    @test Gaugefields.measure_plaquette(cold.gauge) ≈ 1
    copy_configuration!(cold, hot)
    @test !(Gaugefields.measure_plaquette(cold.gauge) ≈ 1)
    @test !apply_metropolis_decision!(cold, backup, false)
    @test Gaugefields.measure_plaquette(cold.gauge) ≈ 1

    @test metropolis_rule(-1.0, 0.999)
    @test metropolis_rule(1.0, 0.1)
    @test !metropolis_rule(1.0, 0.9)
    @test !metropolis_rule(Inf, 0.1)
    @test_throws ArgumentError metropolis_rule(1.0, 0.0)
end

@testset "Gaugefields v1 pure-gauge heatbath simulation" begin
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        element_type=ComplexF64,
        verbose=0,
    )

    for colors in (2, 3), overrelaxation_steps in (0, 3)
        input = test_heatbath_input(
            colors=colors,
            overrelaxation_steps=overrelaxation_steps,
        )
        simulation = build_simulation(input, environment)
        reference = build_configuration(input, environment)
        beta = only(input.gauge_action.terms).coupling
        seed = LatticeQCD.Simulation_module.random_stream_seed(
            input.update.random,
        )
        reference_kernel = Gaugefields.heatbath_updater(
            reference.gauge;
            beta,
            ITERATION_MAX=input.update.max_iterations,
            seed,
            sweep=0,
            overrelaxation_sweep=0,
        )

        @test simulation.updater isa HeatbathUpdater
        @test simulation.state isa HeatbathState
        @test simulation.updater.overrelaxation_steps ==
              overrelaxation_steps
        @test fieldtype(typeof(simulation), :updater) ===
              typeof(simulation.updater)
        @test fieldtype(typeof(simulation), :state) ===
              typeof(simulation.state)
        @test all(isconcretetype, fieldtypes(typeof(simulation)))
        @test all(isconcretetype, fieldtypes(typeof(simulation.updater)))
        @test all(isconcretetype, fieldtypes(typeof(simulation.state)))

        for trajectory in 0:1
            Gaugefields.heatbath!(reference.gauge, reference_kernel)
            for _ in 1:overrelaxation_steps
                Gaugefields.overrelaxation!(
                    reference.gauge,
                    reference_kernel,
                )
            end

            result = update!(simulation)
            @test result isa HeatbathUpdateResult
            @test result.trajectory == trajectory
            @test result.heatbath_sweep == trajectory + 1
            @test result.overrelaxation_sweep ==
                  (trajectory + 1) * overrelaxation_steps
            @test simulation.state.trajectory == trajectory + 1
            @test simulation.state.heatbath_sweep == reference_kernel.sweep
            @test simulation.state.overrelaxation_sweep ==
                  reference_kernel.overrelaxation_sweep
            @test Gaugefields.measure_plaquette(
                simulation.configuration.gauge,
            ) ≈ Gaugefields.measure_plaquette(reference.gauge)
            test_same_bridge_configuration(
                simulation.configuration,
                reference,
            )
        end
    end
end

@testset "Gaugefields v1 general-action heatbath dispatch" begin
    environment = GaugefieldsEnvironment(process_grid=(1, 1, 1, 1))
    input = test_heatbath_input(
        colors=2,
        even_odd=false,
        overrelaxation_steps=1,
    )
    simulation = build_simulation(input, environment)
    reference = build_configuration(input, environment)
    reference_action = build_gauge_action(input.gauge_action, reference)
    seed = LatticeQCD.Simulation_module.random_stream_seed(
        input.update.random,
    )
    reference_kernel = Gaugefields.heatbath_updater(
        reference.gauge,
        reference_action;
        ITERATION_MAX=input.update.max_iterations,
        seed,
        sweep=0,
        overrelaxation_sweep=0,
    )

    Gaugefields.heatbath!(reference.gauge, reference_kernel)
    Gaugefields.overrelaxation!(reference.gauge, reference_kernel)
    result = update!(simulation)

    @test result.heatbath_sweep == 1
    @test result.overrelaxation_sweep == 1
    @test simulation.updater.kernel.colorings !== nothing
    test_same_bridge_configuration(simulation.configuration, reference)
end

@testset "Heatbath validation and failure counters" begin
    environment = GaugefieldsEnvironment(process_grid=(1, 1, 1, 1))
    invalid_action = GaugeActionConfig(
        GaugeActionTermConfig(:plaquette, "plaquette", 1.9),
        GaugeActionTermConfig(:rectangle, "rectangular", -0.1),
    )
    invalid = LQCDConfig(
        LatticeConfig((2, 2, 2, 2)),
        GaugeConfig(2, 1, "cold", nothing),
        invalid_action,
        HeatbathConfig(
            true,
            100_000,
            0,
            RandomStreamConfig(1, :heatbath),
        ),
    )
    @test_throws ArgumentError build_simulation(invalid, environment)

    failing = test_heatbath_input(
        colors=2,
        beta=1e-12,
        max_iterations=1,
    )
    simulation = build_simulation(failing, environment)
    @test_throws ErrorException update!(simulation)
    @test simulation.state.trajectory == 0
    @test simulation.state.heatbath_sweep == 0
    @test simulation.state.overrelaxation_sweep == 0
    @test simulation.updater.kernel.sweep == 0
    @test simulation.updater.kernel.overrelaxation_sweep == 0
end
