using LatticeQCD
using Test
using TOML
import Gaugefields

const WIZARD_FERMION_IO_DIRECTORY = joinpath(
    @__DIR__,
    "wizard_fermion_io",
)

const WIZARD_FERMION_IO_FORMATS = (
    (id="jld2", value="JLD", format=:jld2, extension=".jld2"),
    (id="ildg", value="ILDG", format=:ildg, extension=".ildg"),
    (id="bridge", value="BridgeText", format=:bridge, extension=".txt"),
)

const WIZARD_FERMION_IO_SAVE_FORMATS = (
    (id="none", value="nothing", format=nothing, extension=nothing),
    WIZARD_FERMION_IO_FORMATS...,
)

function expected_wizard_fermion_io_cases()
    cases = NamedTuple[]
    for colors in (2, 3)
        for load_format in WIZARD_FERMION_IO_FORMATS
            for save_format in WIZARD_FERMION_IO_SAVE_FORMATS
                filename = join((
                    "su$colors",
                    "load_$(load_format.id)",
                    "save_$(save_format.id)",
                ), "_") * ".toml"
                push!(cases, (;
                    filename,
                    colors,
                    load_format,
                    save_format,
                ))
            end
        end
    end
    return cases
end

function wizard_fermion_io_quiet(function_to_run)
    quiet_display = Base.Multimedia.TextDisplay(devnull)
    Base.Multimedia.pushdisplay(quiet_display)
    try
        return redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                function_to_run()
            end
        end
    finally
        Base.Multimedia.popdisplay(quiet_display)
    end
end

function wizard_fermion_io_source_input(colors)
    lattice = LatticeConfig((2, 2, 2, 2))
    gauge = GaugeConfig(colors, 1, "hot", nothing, 0x77112233)
    gauge_action = GaugeActionConfig(GaugeActionTermConfig(
        :gauge_plaquette,
        "plaquette",
        colors == 2 ? 2.7 : 5.7,
    ))
    solver = FermionSolverConfig(
        1e-10,
        2_000,
        0,
        (1, 1, 1, -1),
    )
    fermions = (
        FermionActionConfig(
            :fermion_1,
            WilsonDiracConfig(0.141139, 1.0),
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
            RandomStreamConfig(0x44556677, :momentum),
        ),
        RankZeroMetropolisConfig(
            RandomStreamConfig(0x8899aabb, :metropolis),
        ),
        (
            PseudofermionRefreshConfig(
                :fermion_1,
                RandomStreamConfig(0xccddeeff, :pseudofermion),
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

function wizard_fermion_io_bridge_bytes(configuration)
    mktempdir() do directory
        path = joinpath(directory, "configuration.txt")
        save_configuration(path, configuration; format=:bridge)
        return read(path)
    end
end

function make_wizard_fermion_io_sources(root, colors, environment)
    simulation = wizard_fermion_io_quiet() do
        build_simulation(
            wizard_fermion_io_source_input(colors),
            environment,
        )
    end
    wizard_fermion_io_quiet() do
        update!(simulation)
    end

    directory = joinpath(root, "sources", "su$colors")
    mkpath(directory)
    paths = Dict{String,String}()
    for format in WIZARD_FERMION_IO_FORMATS
        Sys.iswindows() && format.format === :ildg && continue
        path = joinpath(directory, "wizard_fermion_source$(format.extension)")
        save_configuration(
            path,
            simulation.configuration;
            format=format.format,
        )
        paths[format.id] = path
    end
    return (;
        paths,
        bytes=wizard_fermion_io_bridge_bytes(simulation.configuration),
    )
end

function prepare_wizard_fermion_io_dictionary!(
    dictionary,
    case,
    root,
    source,
)
    case_id = splitext(case.filename)[1]
    physical = dictionary["Physical setting"]
    control = dictionary["System Control"]
    hmc = dictionary["HMC related"]

    physical["initial"] = source.paths[case.load_format.id]
    physical["initialtrj"] = 1
    physical["Nthermalization"] = 0
    physical["Nsteps"] = 1
    hmc["MDsteps"] = 1
    hmc["Δτ"] = 0.001

    log_directory = joinpath(root, "logs")
    measurement_directory = joinpath(root, "measurements")
    save_directory = joinpath(root, "saved", case_id)
    mkpath(log_directory)
    mkpath(measurement_directory)
    control["log_dir"] = relpath(log_directory, pwd())
    control["logfile"] = "$case_id.log"
    control["measurement_basedir"] = relpath(measurement_directory, pwd())
    control["measurement_dir"] = case_id
    if case.save_format.format === nothing
        control["saveU_format"] = "nothing"
        control["saveU_dir"] = ""
    else
        mkpath(save_directory)
        control["saveU_format"] = case.save_format.value
        control["saveU_every"] = 1
        control["saveU_dir"] = relpath(save_directory, pwd())
    end
    return dictionary
end

function construct_wizard_fermion_io_params(dictionary)
    return wizard_fermion_io_quiet() do
        LatticeQCD.Parameters_TOML.construct_Params_from_TOML(dictionary)
    end
end

function test_wizard_fermion_io_update_equal(reference, candidate)
    @test typeof(candidate) === typeof(reference)
    @test candidate.trajectory == reference.trajectory
    @test candidate.accepted == reference.accepted
    @test candidate.initial_hamiltonian ≈ reference.initial_hamiltonian rtol=1e-12 atol=1e-12
    @test candidate.final_hamiltonian ≈ reference.final_hamiltonian rtol=1e-12 atol=1e-12
    @test candidate.delta_hamiltonian ≈ reference.delta_hamiltonian rtol=1e-12 atol=1e-12
end

function test_wizard_fermion_io_output(
    case,
    session,
    result,
    spec,
    environment,
)
    if case.save_format.format === nothing
        @test session.output.configurations isa NoConfigurationOutput
        @test result.saved_path === nothing
        return
    end

    @test isfile(result.saved_path)
    @test endswith(result.saved_path, case.save_format.extension)
    @test count(
        event -> event isa ConfigurationSaved,
        session.sink.events,
    ) == 1

    reloaded = wizard_fermion_io_quiet() do
        build_configuration(spec.config, environment)
    end
    @test reloaded isa GaugeFermionConfiguration
    returned = load_configuration!(
        reloaded,
        result.saved_path;
        format=case.save_format.value,
    )
    @test returned === reloaded
    @test wizard_fermion_io_bridge_bytes(reloaded) ==
          wizard_fermion_io_bridge_bytes(session.simulation.configuration)
end

@testset "Wizard dynamical-fermion configuration I/O" begin
    cases = expected_wizard_fermion_io_cases()
    existing = sort(filter(
        filename -> endswith(filename, ".toml"),
        readdir(WIZARD_FERMION_IO_DIRECTORY),
    ))
    @test length(cases) == 24
    @test existing == sort(getproperty.(cases, :filename))

    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        communicator=Gaugefields.SerialCommunicator(),
        element_type=ComplexF64,
        verbose=0,
    )
    mktempdir() do root
        sources = Dict(
            colors => make_wizard_fermion_io_sources(
                root,
                colors,
                environment,
            ) for colors in (2, 3)
        )
        executed = 0
        for case in cases
            unsupported_on_windows = Sys.iswindows() && (
                case.load_format.format === :ildg ||
                case.save_format.format === :ildg
            )
            unsupported_on_windows && continue

            path = joinpath(WIZARD_FERMION_IO_DIRECTORY, case.filename)
            dictionary = prepare_wizard_fermion_io_dictionary!(
                TOML.parsefile(path),
                case,
                root,
                sources[case.colors],
            )
            parameters = construct_wizard_fermion_io_params(dictionary)
            try
                @test parameters.initial ==
                      sources[case.colors].paths[case.load_format.id]
                @test parameters.loadU_format == case.load_format.value
                @test parameters.Dirac_operator == "Wilson"
                @test !parameters.quench

                spec = wizard_fermion_io_quiet() do
                    SimulationSpec(parameters)
                end
                @test isempty(validate(spec))
                @test spec.config.gauge.initialization isa FileStartConfig

                reference = wizard_fermion_io_quiet() do
                    build_simulation(spec.config, environment)
                end
                @test reference.configuration isa GaugeFermionConfiguration
                @test wizard_fermion_io_bridge_bytes(reference.configuration) ==
                      sources[case.colors].bytes

                sink = RecordingSimulationEventSink()
                session = wizard_fermion_io_quiet() do
                    build_simulation(spec, environment; sink)
                end
                @test wizard_fermion_io_bridge_bytes(
                    session.simulation.configuration,
                ) == sources[case.colors].bytes

                reference_result = wizard_fermion_io_quiet() do
                    update!(reference)
                end
                result = wizard_fermion_io_quiet() do
                    step!(session)
                end
                test_wizard_fermion_io_update_equal(
                    reference_result,
                    result.update,
                )
                @test wizard_fermion_io_bridge_bytes(reference.configuration) ==
                      wizard_fermion_io_bridge_bytes(
                    session.simulation.configuration,
                )
                test_wizard_fermion_io_output(
                    case,
                    session,
                    result,
                    spec,
                    environment,
                )
                executed += 1
            finally
                isopen(parameters.load_fp) && close(parameters.load_fp)
            end
            GC.gc(false)
        end
        @test executed == (Sys.iswindows() ? 12 : 24)
    end
end
