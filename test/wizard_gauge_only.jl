using LatticeQCD
using Test
using TOML
import Gaugefields

const WIZARD_GAUGE_DIRECTORY = joinpath(@__DIR__, "wizard_gauge_only")

@testset "Gauge-only Wizard parameter matrix" begin
    directory = WIZARD_GAUGE_DIRECTORY
    files = sort(filter(
        path -> endswith(path, ".toml"),
        readdir(directory; join=true),
    ))

    @test length(files) == 256
    @test count(path -> occursin("/su3_fileloading_", path), files) == 24
    @test count(path -> occursin("/su2_fileloading_", path), files) == 24
    @test count(path -> occursin(
        r"/su3_(cold|hot|file_|instanton|embedded_instanton)",
        path,
    ), files) == 96
    @test count(path -> occursin(
        r"/su2_(cold|hot|file_|instanton|embedded_instanton)",
        path,
    ), files) == 112

    for path in files
        parameters = TOML.parsefile(path)
        @test parameters["Physical setting(fermions)"]["Dirac_operator"] ==
              "nothing"
        measurements = parameters["Measurement set"]["measurement_methods"]
        @test Set(keys(measurements)) == Set(["Plaquette"])
        @test measurements["Plaquette"]["methodname"] == "Plaquette"
        gradient = parameters["gradientflow_measurements"]
        @test !gradient["hasgradientflow"]
        @test isempty(gradient["measurements_for_flow"])
    end
end


function wizard_gauge_quiet(function_to_run)
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


function wizard_gauge_source_input(colors)
    lattice = LatticeConfig((2, 2, 2, 2))
    gauge = GaugeConfig(colors, 1, "hot", nothing, 0x7711)
    action = GaugeActionConfig(GaugeActionTermConfig(
        :gauge_plaquette,
        "plaquette",
        colors == 2 ? 2.7 : 5.7,
    ))
    update = HeatbathConfig(
        true,
        100_000,
        0,
        RandomStreamConfig(0x8822, :heatbath),
    )
    return LQCDConfig(lattice, gauge, action, update)
end


function make_wizard_gauge_sources(root, colors, environment)
    directory = joinpath(root, "su$colors", "confs")
    mkpath(directory)
    simulation = build_simulation(
        wizard_gauge_source_input(colors),
        environment,
    )
    plaquettes = Float64[]
    for index in 1:2
        update!(simulation)
        push!(
            plaquettes,
            Gaugefields.measure_plaquette(simulation.configuration.gauge),
        )
        stem = "conf_" * lpad(index, 8, '0')
        Gaugefields.save_configuration(
            joinpath(directory, stem * ".jld2"),
            simulation.configuration.gauge;
            format=:jld2,
        )
        Gaugefields.save_configuration(
            joinpath(directory, stem * ".txt"),
            simulation.configuration.gauge;
            format=:bridge,
        )
        if !Sys.iswindows()
            Gaugefields.save_configuration(
                joinpath(directory, stem * ".ildg"),
                simulation.configuration.gauge;
                format=:ildg,
            )
        end
    end
    return (; directory, plaquettes)
end


function wizard_gauge_format(control)
    value = control["loadU_format"]
    value == "JLD" && return (:jld2, ".jld2")
    value == "ILDG" && return (:ildg, ".ildg")
    value == "BridgeText" && return (:bridge, ".txt")
    error("unexpected Wizard load format $value")
end


function prepare_wizard_gauge_dictionary!(
    dictionary,
    path,
    root,
    source,
)
    case_id = splitext(basename(path))[1]
    physical = dictionary["Physical setting"]
    control = dictionary["System Control"]
    hmc = dictionary["HMC related"]

    physical["L"] = [2, 2, 2, 2]
    physical["Nsteps"] = 1
    physical["Nthermalization"] = 0
    physical["initialtrj"] = 1
    control["verboselevel"] = 0

    log_directory = joinpath(root, "logs")
    measurement_directory = joinpath(root, "measurements")
    save_directory = joinpath(root, "saved", case_id)
    mkpath(log_directory)
    mkpath(measurement_directory)
    mkpath(dirname(save_directory))
    control["log_dir"] = relpath(log_directory, pwd())
    control["logfile"] = "$case_id.log"
    control["measurement_basedir"] = relpath(measurement_directory, pwd())
    control["measurement_dir"] = case_id

    if physical["update_method"] == "HMC"
        hmc["MDsteps"] = 1
        hmc["Δτ"] = 0.001
    end

    save_format = get(control, "saveU_format", "nothing")
    if save_format != "nothing"
        control["saveU_every"] = 1
        control["saveU_dir"] = relpath(save_directory, pwd())
    end

    if occursin("_file_", case_id) &&
       physical["update_method"] != "Fileloading"
        _, extension = wizard_gauge_format(control)
        physical["initial"] = joinpath(
            source.directory,
            "conf_00000001$extension",
        )
    elseif physical["update_method"] == "Fileloading"
        control["loadU_dir"] = source.directory
        _, extension = wizard_gauge_format(control)
        if get(control, "loadU_fromfile", false)
            control["loadU_filename"] = "filelist.txt"
            open(joinpath(source.directory, "filelist.txt"), "w") do io
                println(io, "conf_00000002", extension)
                println(io, "conf_00000001", extension)
            end
        end
    end
    return dictionary
end


function construct_wizard_gauge_params(dictionary)
    return wizard_gauge_quiet() do
        LatticeQCD.Parameters_TOML.construct_Params_from_TOML(dictionary)
    end
end


function wizard_gauge_bridge_bytes(configuration)
    mktempdir() do directory
        path = joinpath(directory, "configuration.txt")
        Gaugefields.save_configuration(
            path,
            configuration.gauge;
            format=:bridge,
        )
        return read(path)
    end
end


function test_wizard_update_result_equal(reference, candidate)
    @test typeof(candidate) === typeof(reference)
    for field in fieldnames(typeof(reference))
        left = getfield(reference, field)
        right = getfield(candidate, field)
        if left isa AbstractFloat
            @test right ≈ left rtol=1e-13 atol=1e-13
        else
            @test right == left
        end
    end
end


function wizard_saved_format(path)
    endswith(path, ".jld2") && return :jld2
    endswith(path, ".ildg") && return :ildg
    return :bridge
end


function test_wizard_saved_configuration(session, result)
    output = session.output.configurations
    if output isa NoConfigurationOutput
        @test result.saved_path === nothing
        return
    end
    @test isfile(result.saved_path)
    loaded = GaugeConfiguration([
        similar(link) for link in session.simulation.configuration.gauge
    ])
    Gaugefields.load_configuration!(
        loaded.gauge,
        result.saved_path;
        format=wizard_saved_format(result.saved_path),
    )
    @test wizard_gauge_bridge_bytes(loaded) == wizard_gauge_bridge_bytes(
        session.simulation.configuration,
    )
end


@testset "All gauge-only Wizard files through typed sessions" begin
    files = sort(filter(
        path -> endswith(path, ".toml"),
        readdir(WIZARD_GAUGE_DIRECTORY; join=true),
    ))
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        element_type=ComplexF64,
        verbose=0,
    )

    mktempdir() do root
        sources = Dict(
            colors => make_wizard_gauge_sources(root, colors, environment)
            for colors in (2, 3)
        )
        executed = 0
        supported = 0
        for (index, path) in enumerate(files)
            colors = startswith(basename(path), "su2_") ? 2 : 3
            dictionary = prepare_wizard_gauge_dictionary!(
                TOML.parsefile(path),
                path,
                root,
                sources[colors],
            )
            control = dictionary["System Control"]
            unsupported_on_windows = Sys.iswindows() && (
                get(control, "loadU_format", nothing) == "ILDG" ||
                (
                    dictionary["Physical setting"]["update_method"] !=
                    "Fileloading" &&
                    get(control, "saveU_format", nothing) == "ILDG"
                )
            )
            if unsupported_on_windows
                continue
            end
            supported += 1

            parameters = construct_wizard_gauge_params(dictionary)
            try
                spec = wizard_gauge_quiet() do
                    SimulationSpec(parameters)
                end
                @test isempty(validate(spec))
                sink = RecordingSimulationEventSink()
                session = wizard_gauge_quiet() do
                    build_simulation(spec, environment; sink)
                end

                if parameters.update_method == "Fileloading"
                    summary = wizard_gauge_quiet() do
                        run!(session)
                    end
                    list_mode = parameters.loadU_fromfile
                    expected = list_mode ? reverse(sources[colors].plaquettes) :
                               sources[colors].plaquettes
                    records = [
                        event.record for event in sink.events
                        if event isa MeasurementFinished
                    ]
                    @test getproperty.(records, :value) ≈ expected
                    @test getproperty.(getproperty.(records, :point), :trajectory) ==
                          [0, 1]
                    @test summary.production_completed == 2
                    @test summary.final_trajectory == 1
                    @test summary.saved_configurations == 0
                    @test count(
                        event -> event isa ConfigurationLoaded,
                        sink.events,
                    ) == 2
                    @test count(
                        event -> event isa ConfigurationSaved,
                        sink.events,
                    ) == 0
                else
                    reference = wizard_gauge_quiet() do
                        build_simulation(spec.config, environment)
                    end
                    reference_result = wizard_gauge_quiet() do
                        update!(reference)
                    end
                    result = wizard_gauge_quiet() do
                        step!(session)
                    end
                    test_wizard_update_result_equal(
                        reference_result,
                        result.update,
                    )
                    @test wizard_gauge_bridge_bytes(reference.configuration) ==
                          wizard_gauge_bridge_bytes(
                        session.simulation.configuration,
                    )
                    @test length(result.measurements) == 1
                    @test only(result.measurements).name === :plaquette
                    @test only(result.measurements).value ≈
                          Gaugefields.measure_plaquette(
                        reference.configuration.gauge,
                    )
                    test_wizard_saved_configuration(session, result)
                end
                executed += 1
            finally
                isopen(parameters.load_fp) && close(parameters.load_fp)
            end
            index % 8 == 0 && GC.gc(false)
        end
        @test executed == supported
        @test supported == (Sys.iswindows() ? 164 : 256)
    end
end
