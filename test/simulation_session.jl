using LatticeQCD
using Test
using TOML
import Gaugefields

function session_hmc_input()
    lattice = LatticeConfig((2, 2, 2, 2))
    gauge = GaugeConfig(2, 1, "cold", nothing, 0x1020)
    gauge_action = GaugeActionConfig(
        GaugeActionTermConfig(:gauge_plaquette, "plaquette", 1.9),
    )
    integrator = LeapfrogConfig(
        QPQConfig(),
        ForceGroupConfig(:gauge),
    )
    md = MDConfig(0.001, 1, integrator)
    momentum = GaussianMomentumConfig(
        1.0,
        RandomStreamConfig(0x3040, :momentum),
    )
    acceptance = RankZeroMetropolisConfig(
        RandomStreamConfig(0x5060, :metropolis),
    )
    return LQCDConfig(
        lattice,
        gauge,
        gauge_action,
        HMCConfig(md, momentum, acceptance),
    )
end

function session_wilson_hmc_input()
    lattice = LatticeConfig((2, 2, 2, 2))
    gauge = GaugeConfig(2, 1, "cold", nothing, 0x1234)
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
            WilsonDiracConfig(0.05, 1.0),
            2,
            solver,
        ),
    )
    integrator = LeapfrogConfig(
        QPQConfig(),
        ForceGroupConfig(:gauge, :fermion_1),
    )
    update = HMCConfig(
        MDConfig(0.001, 1, integrator),
        GaussianMomentumConfig(
            1.0,
            RandomStreamConfig(0x5678, :momentum),
        ),
        RankZeroMetropolisConfig(
            RandomStreamConfig(0x9abc, :metropolis),
        ),
        (
            PseudofermionRefreshConfig(
                :fermion_1,
                RandomStreamConfig(0xdef0, :pseudofermion),
                1,
            ),
        ),
    )
    return LQCDConfig(
        lattice,
        gauge,
        gauge_action,
        fermions,
        update,
    )
end

function session_measurement_program()
    measurement = ScheduledMeasurementConfig(
        PlaquetteObservableConfig(),
        PeriodicSchedule(1),
    )
    return MeasurementProgram(MeasurementPlan(measurement))
end

function quietly_construct_session_params(dictionary, directory)
    control = dictionary["System Control"]
    relative_directory = relpath(directory, pwd())
    control["log_dir"] = joinpath(relative_directory, "logs")
    control["logfile"] = "simulation-session.log"
    control["measurement_basedir"] = joinpath(
        relative_directory,
        "measurements",
    )
    control["measurement_dir"] = "simulation-session"
    quiet_display = Base.Multimedia.TextDisplay(devnull)
    Base.Multimedia.pushdisplay(quiet_display)
    try
        return redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                LatticeQCD.Parameters_TOML.construct_Params_from_TOML(
                    dictionary,
                )
            end
        end
    finally
        Base.Multimedia.popdisplay(quiet_display)
    end
end

@testset "GUI-neutral simulation specification" begin
    program = session_measurement_program()
    schedule = SimulationSchedule(
        2,
        3;
        initial_trajectory=7,
        measurements=program,
    )
    output = OutputConfig(JLD2ConfigurationOutput(
        "confs";
        prefix="configuration_",
        every=5,
        width=6,
    ))
    spec = SimulationSpec(session_hmc_input(), schedule, output)

    @test schedule.thermalization_steps == 2
    @test schedule.production_steps == 3
    @test schedule.initial_trajectory == 7
    @test output.configurations.every == 5
    @test configuration_output_path(output.configurations, 12) ==
          joinpath("confs", "configuration_000012.jld2")
    @test isempty(validate(spec))
    @test all(isconcretetype, fieldtypes(typeof(spec)))
    @test all(isconcretetype, fieldtypes(typeof(schedule)))
    @test all(isconcretetype, fieldtypes(typeof(output)))
    @test all(isconcretetype, fieldtypes(typeof(output.configurations)))

    text = sprint(show_config, spec)
    @test occursin("SimulationSpec", text)
    @test occursin("thermalization steps: 2", text)
    @test occursin("portable JLD2", text)

    @test_throws ArgumentError SimulationSchedule(-1, 1)
    @test_throws ArgumentError SimulationSchedule(0, -1)
    @test_throws ArgumentError SimulationSchedule(
        0,
        1;
        initial_trajectory=-1,
    )
    @test_throws ArgumentError JLD2ConfigurationOutput("")
    @test_throws ArgumentError JLD2ConfigurationOutput("confs"; every=0)

    sequence_config = LQCDConfig(
        LatticeConfig((2, 2, 2, 2)),
        GaugeConfig(2, 1, "cold", nothing),
        GaugeActionConfig(
            GaugeActionTermConfig(:gauge, "plaquette", 1.9),
        ),
        ConfigurationSequenceConfig(
            DirectorySourceConfig("confs"),
            :jld2,
        ),
    )
    invalid_sequence = SimulationSpec(
        sequence_config,
        SimulationSchedule(
            1,
            2;
            initial_trajectory=3,
            measurements=program,
        ),
        output,
    )
    issues = validate(invalid_sequence)
    @test length(issues) == 3
    @test getproperty.(issues, :code) == [
        :not_supported,
        :must_be_zero,
        :not_supported,
    ]
    @test getproperty.(issues, :path) == [
        ["schedule", "thermalization_steps"],
        ["schedule", "initial_trajectory"],
        ["output", "configurations"],
    ]
    @test all(issue -> all(
        isconcretetype,
        fieldtypes(typeof(issue)),
    ), issues)
    @test_throws ArgumentError build_simulation(invalid_sequence)
end

@testset "Legacy Params to SimulationSpec" begin
    mktempdir() do directory
        dictionary = TOML.parsefile(joinpath(
            @__DIR__,
            "wizard_fermion",
            "su2_wilson_thin_leapfrog.toml",
        ))
        physical = dictionary["Physical setting"]
        physical["initialtrj"] = 1
        physical["Nthermalization"] = 2
        physical["Nsteps"] = 4
        control = dictionary["System Control"]
        control["saveU_format"] = "JLD"
        control["saveU_every"] = 2
        control["saveU_dir"] = joinpath(
            relpath(directory, pwd()),
            "confs",
        )

        parameters = quietly_construct_session_params(dictionary, directory)
        try
            spec = SimulationSpec(parameters)
            @test spec.schedule.initial_trajectory == 0
            @test spec.schedule.thermalization_steps == 1
            @test spec.schedule.production_steps == 3
            @test length(spec.schedule.measurements.direct.measurements) == 1
            @test spec.output.configurations isa JLD2ConfigurationOutput
            @test spec.output.configurations.every == 2
            @test spec.config.fermions[1].operator isa WilsonDiracConfig
            @test isempty(validate(spec))
        finally
            isopen(parameters.load_fp) && close(parameters.load_fp)
        end
    end
end

@testset "Headless session events, measurements, and portable JLD2" begin
    mktempdir() do directory
        schedule = SimulationSchedule(
            1,
            3;
            measurements=session_measurement_program(),
        )
        output = OutputConfig(JLD2ConfigurationOutput(
            directory;
            every=2,
        ))
        spec = SimulationSpec(session_hmc_input(), schedule, output)
        sink = RecordingSimulationEventSink()
        session = build_simulation(
            spec,
            GaugefieldsEnvironment(
                process_grid=(1, 1, 1, 1),
                element_type=ComplexF64,
                verbose=0,
            );
            sink,
        )

        @test session isa SimulationSession
        @test all(isconcretetype, fieldtypes(typeof(session)))
        @test all(isconcretetype, fieldtypes(typeof(session.state)))
        @test !is_finished(session)
        @test !is_running(session)

        compact_display = sprint(show, session)
        @test startswith(compact_display, "SimulationSession(Simulation(")
        @test occursin("thermalization=0/1", compact_display)
        @test occursin("production=0/3", compact_display)
        @test occursin("status=ready", compact_display)
        @test ncodeunits(compact_display) < 500

        plain_display = sprint(show, MIME"text/plain"(), session)
        @test startswith(plain_display, "SimulationSession\n")
        @test occursin("  thermalization: 0 / 1", plain_display)
        @test occursin("  production: 0 / 3", plain_display)
        @test occursin("  status: ready", plain_display)
        @test ncodeunits(plain_display) < 1_000

        thermalization = step!(session)
        @test thermalization.phase === :thermalization
        @test isempty(thermalization.measurements)
        @test thermalization.saved_path === nothing
        @test session.state.thermalization_completed == 1
        @test session.state.production_completed == 0
        @test count(event -> event isa ThermalizationStepFinished, sink.events) == 1

        console = IOBuffer()
        summary = run!(session; io=console)
        console_output = String(take!(console))
        @test occursin("# run started:", console_output)
        @test occursin("# trajectory=2", console_output)
        @test occursin("accepted=", console_output)
        @test occursin("# measurement: trajectory=2 plaquette =", console_output)
        @test occursin("# configuration saved: trajectory=2", console_output)
        @test occursin("# run finished:", console_output)
        @test summary == run_summary(session)
        @test summary.initial_trajectory == 0
        @test summary.final_trajectory == 4
        @test summary.thermalization_completed == 1
        @test summary.production_completed == 3
        @test summary.saved_configurations == 2
        @test !summary.stopped
        @test occursin("production=3", sprint(show, summary))
        @test occursin(
            "production completed: 3",
            sprint(show, MIME"text/plain"(), summary),
        )
        @test is_finished(session)
        @test !is_running(session)
        @test_throws EOFError step!(session)

        @test count(event -> event isa RunStarted, sink.events) == 1
        @test count(event -> event isa TrajectoryFinished, sink.events) == 3
        @test count(event -> event isa MeasurementFinished, sink.events) == 3
        @test count(event -> event isa ConfigurationSaved, sink.events) == 2
        @test count(event -> event isa RunFinished, sink.events) == 1
        @test count(event -> event isa RunStopped, sink.events) == 0
        @test count(event -> event isa RunFailed, sink.events) == 0
        @test all(event -> all(
            isconcretetype,
            fieldtypes(typeof(event)),
        ), sink.events)

        paths = sort(filter(
            path -> endswith(path, ".jld2"),
            readdir(directory; join=true),
        ))
        @test basename.(paths) == [
            "conf_00000002.jld2",
            "conf_00000004.jld2",
        ]
        saved_events = filter(
            event -> event isa ConfigurationSaved,
            sink.events,
        )
        @test getproperty.(saved_events, :path) == paths
        @test getproperty.(saved_events, :trajectory) == [2, 4]

        loaded = Gaugefields.load_configuration(
            last(paths);
            process_grid=(1, 1, 1, 1),
            halo=1,
            verbose=0,
        )
        @test Gaugefields.measure_plaquette(loaded) ≈
              Gaugefields.measure_plaquette(
            session.simulation.configuration.gauge,
        )
    end
end

@testset "Cooperative GUI stop and resume" begin
    events = Any[]
    session_slot = Ref{Any}()
    stop_once = Ref(true)
    sink = FunctionSimulationEventSink() do event
        push!(events, event)
        if event isa TrajectoryFinished && stop_once[]
            stop_once[] = false
            request_stop!(session_slot[])
        end
    end
    spec = SimulationSpec(
        session_hmc_input(),
        SimulationSchedule(0, 4),
    )
    session = build_simulation(
        spec,
        GaugefieldsEnvironment(
            process_grid=(1, 1, 1, 1),
            element_type=ComplexF64,
            verbose=0,
        );
        sink,
    )
    session_slot[] = session

    silent_output = IOBuffer()
    stopped = run!(session; verbose=false, io=silent_output)
    @test isempty(take!(silent_output))
    @test stopped.stopped
    @test stopped.production_completed == 1
    @test !is_finished(session)
    @test count(event -> event isa RunStopped, events) == 1

    finished = run!(session; verbose=false)
    @test !finished.stopped
    @test finished.production_completed == 4
    @test is_finished(session)
    @test count(event -> event isa RunStarted, events) == 2
    @test count(event -> event isa RunFinished, events) == 1
end

@testset "Safe periodic HMC restart checkpoints" begin
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        element_type=ComplexF64,
        verbose=0,
    )
    schedule = SimulationSchedule(1, 4)

    reference_events = Any[]
    reference = build_simulation(
        SimulationSpec(session_hmc_input(), schedule),
        environment;
        sink=FunctionSimulationEventSink(event -> push!(reference_events, event)),
    )
    run!(reference; verbose=false)

    mktempdir() do directory
        checkpoints = JLD2CheckpointOutput(
            directory;
            every=2,
            prefix="restart_",
        )
        output = OutputConfig(; checkpoints)
        spec = SimulationSpec(session_hmc_input(), schedule, output)
        interrupted = build_simulation(spec, environment)

        first = step!(interrupted)
        second = step!(interrupted)
        checkpoint_path = checkpoint_output_path(checkpoints, 2)
        @test first.checkpoint_path === nothing
        @test second.checkpoint_path == checkpoint_path
        @test isfile(checkpoint_path)
        @test !ispath(checkpoint_path * ".pending")
        @test interrupted.state.saved_checkpoints == 1

        restored_events = Any[]
        restored = build_simulation(
            spec,
            environment;
            sink=FunctionSimulationEventSink(
                event -> push!(restored_events, event),
            ),
        )
        load_checkpoint!(restored, checkpoint_path)
        @test restored.simulation.state.trajectory == 2
        @test restored.simulation.state.accepted ==
              interrupted.simulation.state.accepted
        @test restored.state.thermalization_completed == 1
        @test restored.state.production_completed == 1
        @test restored.state.saved_checkpoints == 1
        @test count(event -> event isa CheckpointLoaded, restored_events) == 1

        summary = run!(restored; verbose=false)
        @test !summary.stopped
        @test is_finished(restored)
        @test restored.state.saved_checkpoints == 2
        @test summary.saved_checkpoints == 2
        @test isfile(checkpoint_output_path(checkpoints, 4))
        @test !ispath(checkpoint_output_path(checkpoints, 4) * ".pending")

        for direction in eachindex(reference.simulation.configuration.gauge)
            reference_values = Array(
                reference.simulation.configuration.gauge[direction].U.A,
            )
            restored_values = Array(
                restored.simulation.configuration.gauge[direction].U.A,
            )
            @test restored_values == reference_values
        end
        @test restored.simulation.state.accepted ==
              reference.simulation.state.accepted

        reference_updates = [
            event.update for event in reference_events
            if event isa TrajectoryFinished
        ]
        restored_updates = [
            event.update for event in restored_events
            if event isa TrajectoryFinished
        ]
        @test getproperty.(restored_updates, :accepted) ==
              getproperty.(reference_updates[2:4], :accepted)
        @test getproperty.(restored_updates, :delta_hamiltonian) ==
              getproperty.(reference_updates[2:4], :delta_hamiltonian)
    end

    @test_throws ArgumentError JLD2CheckpointOutput("")
    @test_throws ArgumentError JLD2CheckpointOutput("restart"; every=0)
end

@testset "Dynamical Wilson HMC restart is exact" begin
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        communicator=Gaugefields.SerialCommunicator(),
        element_type=ComplexF64,
        verbose=0,
    )
    config = session_wilson_hmc_input()
    schedule = SimulationSchedule(0, 3)

    reference_events = Any[]
    reference = build_simulation(
        SimulationSpec(config, schedule),
        environment;
        sink=FunctionSimulationEventSink(
            event -> push!(reference_events, event),
        ),
    )
    run!(reference; verbose=false)

    mktempdir() do directory
        checkpoints = JLD2CheckpointOutput(directory; every=1)
        spec = SimulationSpec(
            config,
            schedule,
            OutputConfig(; checkpoints),
        )
        interrupted = build_simulation(spec, environment)
        first = step!(interrupted)
        @test first.checkpoint_path == checkpoint_output_path(checkpoints, 1)

        checkpoint_metadata =
            LatticeQCD.SimulationSession_module.read_checkpoint_metadata(
                first.checkpoint_path,
            )
        @test checkpoint_metadata.input_signature_version == 1
        @test length(checkpoint_metadata.input_fingerprint) == 64
        @test checkpoint_metadata.package_versions.latticeqcd ==
              string(pkgversion(LatticeQCD))
        @test checkpoint_metadata.package_versions.gaugefields ==
              string(pkgversion(Gaugefields))
        independently_built = build_simulation(spec, environment)
        @test LatticeQCD.SimulationSession_module.checkpoint_input_fingerprint(
            LatticeQCD.SimulationSession_module.checkpoint_input_snapshot(
                independently_built,
            ),
        ) == checkpoint_metadata.input_fingerprint

        mismatched_gauge_action = GaugeActionConfig(
            GaugeActionTermConfig(:gauge_plaquette, "plaquette", 2.1),
        )
        mismatched_beta_config = LQCDConfig(
            config.lattice,
            config.gauge,
            mismatched_gauge_action,
            config.fermions,
            config.update,
        )
        mismatched_beta = build_simulation(
            SimulationSpec(mismatched_beta_config, schedule),
            environment,
        )
        @test_throws ArgumentError load_checkpoint!(
            mismatched_beta,
            first.checkpoint_path,
        )

        original_fermion = only(config.fermions)
        mismatched_fermion = FermionActionConfig(
            original_fermion.name,
            WilsonDiracConfig(0.06, 1.0),
            original_fermion.flavors,
            original_fermion.solver,
            original_fermion.smearing,
        )
        mismatched_kappa_config = LQCDConfig(
            config.lattice,
            config.gauge,
            config.gauge_action,
            (mismatched_fermion,),
            config.update,
        )
        mismatched_kappa = build_simulation(
            SimulationSpec(mismatched_kappa_config, schedule),
            environment,
        )
        @test_throws ArgumentError load_checkpoint!(
            mismatched_kappa,
            first.checkpoint_path,
        )

        fake_versions = merge(
            checkpoint_metadata.package_versions,
            (latticeqcd="0.0.0",),
        )
        version_mismatch = merge(
            checkpoint_metadata,
            (package_versions=fake_versions,),
        )
        @test_logs (:warn, r"package versions") (
            LatticeQCD.SimulationSession_module.validate_checkpoint_metadata(
                independently_built,
                version_mismatch,
            )
        )
        @test_throws ArgumentError (
            LatticeQCD.SimulationSession_module.validate_checkpoint_metadata(
                independently_built,
                version_mismatch;
                strict_versions=true,
            )
        )

        restored_events = Any[]
        restored = build_simulation(
            spec,
            environment;
            sink=FunctionSimulationEventSink(
                event -> push!(restored_events, event),
            ),
        )
        load_checkpoint!(restored, first.checkpoint_path)
        @test restored.simulation.configuration isa GaugeFermionConfiguration
        run!(restored; verbose=false)

        @test restored.simulation.state.accepted ==
              reference.simulation.state.accepted
        for direction in eachindex(reference.simulation.configuration.gauge)
            difference = restored.simulation.configuration.gauge[
                direction
            ].U.A .- reference.simulation.configuration.gauge[direction].U.A
            @test maximum(abs, difference) == 0
        end

        reference_updates = [
            event.update for event in reference_events
            if event isa TrajectoryFinished
        ]
        restored_updates = [
            event.update for event in restored_events
            if event isa TrajectoryFinished
        ]
        @test getproperty.(restored_updates, :accepted) ==
              getproperty.(reference_updates[2:3], :accepted)
        @test getproperty.(restored_updates, :delta_hamiltonian) ==
              getproperty.(reference_updates[2:3], :delta_hamiltonian)
    end
end
