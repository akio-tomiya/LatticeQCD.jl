using LatticeQCD
using Test
using TOML
import Gaugefields

function fermion_params(filename, directory; remove_halo=false)
    parameters = TOML.parsefile(joinpath(@__DIR__, filename))
    remove_halo && delete!(parameters["Physical setting"], "Nwing")
    control = parameters["System Control"]
    relative_directory = relpath(directory, pwd())
    control["log_dir"] = joinpath(relative_directory, "logs")
    control["logfile"] = "fermion-config.log"
    control["measurement_basedir"] = joinpath(
        relative_directory,
        "measurements",
    )
    control["measurement_dir"] = "fermion"
    return LatticeQCD.Parameters_TOML.construct_Params_from_TOML(parameters)
end

@testset "Typed dynamical-fermion config from Params" begin
    cases = (
        ("test_wilson.toml", WilsonDiracConfig, true),
        ("test_wilsonclover.toml", WilsonCloverDiracConfig, true),
        ("test_staggered.toml", StaggeredDiracConfig, false),
        (
            "wizard_fermion/su3_hisq_nf4_thin_leapfrog.toml",
            HISQDiracConfig,
            false,
        ),
        ("test_domainwallhmc.toml", DomainwallDiracConfig, false),
        (
            "wizard_fermion/su3_mobius_domainwall_thin_leapfrog.toml",
            MobiusDomainwallDiracConfig,
            false,
        ),
    )

    mktempdir() do directory
        for (filename, operator_type, uses_sw) in cases
            parameters = fermion_params(filename, directory)
            try
                config = LQCDConfig(parameters)
                @test length(config.fermions) == 1
                fermion = only(config.fermions)
                @test fermion.name === :fermion_1
                @test fermion.operator isa operator_type
                @test fermion.solver.tolerance == parameters.eps
                @test fermion.solver.max_steps == parameters.MaxCGstep
                @test fermion.solver.boundary_conditions ==
                      Tuple(parameters.BoundaryCondition)
                @test fermion.smearing isa NoFermionSmearingConfig
                @test all(isconcretetype, fieldtypes(typeof(fermion)))
                @test all(isconcretetype, fieldtypes(typeof(fermion.solver)))

                @test length(config.update.pseudofermions) == 1
                refresh = only(config.update.pseudofermions)
                @test refresh.action_name === :fermion_1
                @test refresh.random.name === :pseudofermion
                @test refresh.subgroup == 1

                if uses_sw
                    @test config.update.md.integrator isa
                          SextonWeingartenConfig
                    @test config.update.md.integrator.slow_forces.names ==
                          (:fermion_1,)
                    @test config.update.md.integrator.fast_integrator.forces.names ==
                          (:gauge,)
                    @test config.update.md.integrator.fast_steps ==
                          parameters.N_SextonWeingargten ÷ 2
                else
                    @test config.update.md.integrator isa LeapfrogConfig
                    @test config.update.md.integrator.forces.names ==
                          (:gauge, :fermion_1)
                end

                output = sprint(show_config, config)
                @test occursin("fermion_1", output)
                @test occursin(String(nameof(operator_type)), output)
                @test occursin("pseudofermion refreshes", output)
            finally
                isopen(parameters.load_fp) && close(parameters.load_fp)
            end
        end

        legacy_parameters = fermion_params(
            "test_wilson.toml",
            directory;
            remove_halo=true,
        )
        try
            @test legacy_parameters.Nwing == 0
            legacy_config = @test_logs (
                :warn,
                r"Nwing=0.*using one gauge-field halo layer",
            ) LQCDConfig(legacy_parameters)
            @test legacy_config.gauge.halo == 1
        finally
            isopen(legacy_parameters.load_fp) && close(legacy_parameters.load_fp)
        end
    end
end

function dynamical_test_input(
    operator;
    colors=2,
    flavors=2,
    sexton_weingarten=false,
    smearing=NoFermionSmearingConfig(),
)
    uses_hisq = operator isa HISQDiracConfig
    colors = uses_hisq ? 3 : colors
    lattice = LatticeConfig(uses_hisq ? (4, 4, 4, 4) : (2, 2, 2, 2))
    gauge = GaugeConfig(colors, uses_hisq ? 3 : 1, "cold", nothing, 0x1234)
    gauge_action = GaugeActionConfig(
        GaugeActionTermConfig(:gauge_plaquette, "plaquette", 1.9),
    )
    solver = FermionSolverConfig(
        1e-10,
        2_000,
        0,
        (1, 1, 1, -1),
    )
    fermions = (
        FermionActionConfig(
            :fermion_1,
            operator,
            flavors,
            solver,
            smearing,
        ),
    )
    integrator = if sexton_weingarten
        SextonWeingartenConfig(
            QPQConfig(),
            ForceGroupConfig(:fermion_1),
            LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge)),
            2,
        )
    else
        LeapfrogConfig(
            QPQConfig(),
            ForceGroupConfig(:gauge, :fermion_1),
        )
    end
    md = MDConfig(0.001, 1, integrator)
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
    update = HMCConfig(md, momentum, acceptance, refreshes)
    return LQCDConfig(
        lattice,
        gauge,
        gauge_action,
        fermions,
        update,
    )
end

@testset "LDO pseudofermion actions in the Gaugefields MD driver" begin
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        communicator=Gaugefields.SerialCommunicator(),
        element_type=ComplexF64,
        verbose=0,
    )
    cases = (
        (
            WilsonDiracConfig(0.05, 1.0),
            2,
            2,
            NoFermionSmearingConfig(),
        ),
        (
            WilsonCloverDiracConfig(0.05, 1.0, 1.0),
            2,
            2,
            NoFermionSmearingConfig(),
        ),
        (
            StaggeredDiracConfig(0.5),
            4,
            2,
            NoFermionSmearingConfig(),
        ),
        (
            HISQDiracConfig(0.5, -0.083),
            4,
            3,
            NoFermionSmearingConfig(),
        ),
        (
            DomainwallDiracConfig(1.0, -1.0, 2),
            2,
            2,
            NoFermionSmearingConfig(),
        ),
        (
            MobiusDomainwallDiracConfig(0.1, -1.0, 2, 2.0, 1.0),
            2,
            2,
            NoFermionSmearingConfig(),
        ),
        (
            WilsonDiracConfig(0.05, 1.0),
            2,
            2,
            StoutFermionSmearingConfig(1, (0.1,), ("plaquette",)),
        ),
        (
            WilsonDiracConfig(0.05, 1.0),
            2,
            2,
            StoutFermionSmearingConfig(4, (0.025,), ("plaquette",)),
        ),
    )

    for (operator, flavors, colors, smearing) in cases
        input = dynamical_test_input(
            operator;
            colors,
            flavors,
            smearing,
        )
        simulation = build_simulation(input, environment)
        @test simulation.configuration isa GaugeFermionConfiguration
        @test keys(simulation.configuration.fermions) == (:fermion_1,)
        if operator isa Union{
            DomainwallDiracConfig,
            MobiusDomainwallDiracConfig,
        }
            field = simulation.configuration.fermions.fermion_1
            @test Tuple(field.f.phases) == (1, 1, 1, -1, 1)
        end
        @test simulation.action isa Gaugefields.MDActionSet
        @test keys(simulation.action.terms) == (:gauge, :fermion_1)
        provider = simulation.action.terms.fermion_1
        if smearing isa NoFermionSmearingConfig
            @test provider.smearing === nothing
        else
            @test provider.smearing isa Gaugefields.CovNeuralnet
            @test length(provider.smearing.layers) == smearing.layers
        end
        @test simulation.updater.pseudofermion_refreshes isa Tuple
        @test length(simulation.updater.pseudofermion_refreshes) == 1
        @test simulation.updater.md_driver.integrator isa Gaugefields.QPQ
        @test all(isconcretetype, fieldtypes(typeof(simulation)))

        result = update!(simulation)
        @test result isa HMCUpdateResult
        @test isfinite(result.initial_hamiltonian)
        @test isfinite(result.final_hamiltonian)
        @test isfinite(result.delta_hamiltonian)
        @test simulation.state.trajectory == 1
    end

    valid_hisq = dynamical_test_input(HISQDiracConfig(0.5, -0.083);
        flavors=4)
    bad_color_hisq = LQCDConfig(
        valid_hisq.lattice,
        GaugeConfig(2, 3, "cold", nothing, 0x1234),
        valid_hisq.gauge_action,
        valid_hisq.fermions,
        valid_hisq.update,
    )
    @test_throws ArgumentError build_simulation(bad_color_hisq, environment)
    bad_halo_hisq = LQCDConfig(
        valid_hisq.lattice,
        GaugeConfig(3, 2, "cold", nothing, 0x1234),
        valid_hisq.gauge_action,
        valid_hisq.fermions,
        valid_hisq.update,
    )
    @test_throws ArgumentError build_simulation(bad_halo_hisq, environment)

    hisq_issues = validate(SimulationSpec(
        bad_color_hisq,
        SimulationSchedule(0, 1),
    ))
    @test getproperty.(hisq_issues, :code) == [:must_equal_three]
    hisq_issues = validate(SimulationSpec(
        bad_halo_hisq,
        SimulationSchedule(0, 1),
    ))
    @test getproperty.(hisq_issues, :code) == [:minimum]

    sw_input = dynamical_test_input(
        WilsonDiracConfig(0.05, 1.0);
        sexton_weingarten=true,
    )
    sw_simulation = build_simulation(sw_input, environment)
    @test sw_simulation.updater.md_driver.integrator isa
          Gaugefields.SextonWeingarten
    @test sw_simulation.updater.md_driver.integrator.n_fast == 2
    sw_result = update!(sw_simulation)
    @test isfinite(sw_result.delta_hamiltonian)

    repeated = build_simulation(sw_input, environment)
    repeated_result = update!(repeated)
    @test repeated_result.accepted == sw_result.accepted
    @test repeated_result.delta_hamiltonian ≈ sw_result.delta_hamiltonian
    @test Gaugefields.measure_plaquette(repeated.configuration.gauge) ≈
          Gaugefields.measure_plaquette(sw_simulation.configuration.gauge)

    solver = only(sw_input.fermions).solver
    multiple_fermions = (
        FermionActionConfig(
            :light,
            WilsonDiracConfig(0.05, 1.0),
            2,
            solver,
        ),
        FermionActionConfig(
            :heavy,
            WilsonDiracConfig(0.03, 1.0),
            2,
            solver,
        ),
    )
    multiple_integrator = SextonWeingartenConfig(
        QPQConfig(),
        ForceGroupConfig(:light, :heavy),
        LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge)),
        2,
    )
    multiple_refreshes = (
        PseudofermionRefreshConfig(
            :light,
            RandomStreamConfig(0xdef0, :pseudofermion),
            1,
        ),
        PseudofermionRefreshConfig(
            :heavy,
            RandomStreamConfig(0xdef0, :pseudofermion),
            2,
        ),
    )
    multiple_update = HMCConfig(
        MDConfig(0.001, 1, multiple_integrator),
        sw_input.update.momentum,
        sw_input.update.acceptance,
        multiple_refreshes,
    )
    multiple_input = LQCDConfig(
        sw_input.lattice,
        sw_input.gauge,
        sw_input.gauge_action,
        multiple_fermions,
        multiple_update,
    )
    multiple = build_simulation(multiple_input, environment)
    @test keys(multiple.configuration.fermions) == (:light, :heavy)
    @test keys(multiple.action.terms) == (:gauge, :light, :heavy)
    @test length(multiple.updater.pseudofermion_refreshes) == 2
    multiple_result = update!(multiple)
    @test isfinite(multiple_result.delta_hamiltonian)

    dynamical_measurements = MeasurementProgram(MeasurementPlan(
        ScheduledMeasurementConfig(
            PlaquetteObservableConfig(),
            PeriodicSchedule(1),
        ),
    ))
    dynamical_spec = SimulationSpec(
        dynamical_test_input(WilsonDiracConfig(0.05, 1.0)),
        SimulationSchedule(
            0,
            1;
            measurements=dynamical_measurements,
        ),
    )
    dynamical_events = RecordingSimulationEventSink()
    dynamical_session = build_simulation(
        dynamical_spec,
        environment;
        sink=dynamical_events,
    )
    @test dynamical_session.simulation.configuration isa
          GaugeFermionConfiguration
    @test dynamical_session.measurements isa MeasurementProgramRuntime
    dynamical_step = step!(dynamical_session)
    @test dynamical_step.phase === :production
    @test length(dynamical_step.measurements) == 1
    @test only(dynamical_step.measurements).name === :plaquette
    @test count(
        event -> event isa MeasurementFinished,
        dynamical_events.events,
    ) == 1
end
