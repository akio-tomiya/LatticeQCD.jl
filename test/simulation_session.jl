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
