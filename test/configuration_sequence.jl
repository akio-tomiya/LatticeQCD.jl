using LatticeQCD
using Test
using TOML
import Gaugefields

const SEQUENCE_FIXTURE_DIRECTORY = joinpath(@__DIR__, "wizard_gauge_only")

function construct_sequence_parameters(path)
    parameters = TOML.parsefile(path)
    control = parameters["System Control"]
    control["log_dir"] = "./logs"
    control["logfile"] = "configuration-sequence.log"
    control["measurement_basedir"] = "./measurements"
    control["measurement_dir"] = "configuration-sequence"

    display_backend = TextDisplay(devnull)
    pushdisplay(display_backend)
    try
        return redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                LatticeQCD.Parameters_TOML.construct_Params_from_TOML(
                    parameters,
                )
            end
        end
    finally
        popdisplay(display_backend)
    end
end

function save_sequence_source(directory, gauge, index)
    stem = "conf_" * lpad(index, 8, '0')
    Gaugefields.save_configuration(
        joinpath(directory, stem * ".jld2"),
        gauge;
        format=:jld2,
    )
    Gaugefields.save_configuration(
        joinpath(directory, stem * ".txt"),
        gauge;
        format=:bridge,
    )
    if !Sys.iswindows()
        Gaugefields.save_configuration(
            joinpath(directory, stem * ".ildg"),
            gauge;
            format=:ildg,
            tempfile1=joinpath(dirname(directory), stem * "-payload.dat"),
            tempfile2=joinpath(dirname(directory), stem * "-filelist.dat"),
        )
    end
    return nothing
end

function make_sequence_sources(root, colors, environment)
    work_directory = joinpath(root, "su$colors")
    configuration_directory = joinpath(work_directory, "confs")
    mkpath(configuration_directory)
    input = test_hmc_input(
        initial="hot",
        lattice_size=(4, 4, 4, 4),
        colors=colors,
        halo=0,
    )
    simulation = build_simulation(input, environment)
    plaquettes = Float64[]
    bridge_data = Vector{UInt8}[]
    for index in 1:2
        result = update!(simulation)
        isfinite(result.delta_hamiltonian) || error(
            "non-finite HMC source trajectory for SU($colors)",
        )
        push!(
            plaquettes,
            Gaugefields.measure_plaquette(simulation.configuration.gauge),
        )
        save_sequence_source(
            configuration_directory,
            simulation.configuration.gauge,
            index,
        )
        push!(
            bridge_data,
            read(joinpath(
                configuration_directory,
                "conf_" * lpad(index, 8, '0') * ".txt",
            )),
        )
    end
    return (;
        work_directory,
        configuration_directory,
        plaquettes,
        bridge_data,
    )
end

function sequence_extension(format)
    format == "JLD" && return ".jld2"
    format == "ILDG" && return ".ildg"
    format == "BridgeText" && return ".txt"
    error("unexpected sequence format $format")
end

function write_sequence_manifest(directory, format)
    extension = sequence_extension(format)
    manifest = joinpath(directory, "filelist.txt")
    open(manifest, "w") do io
        println(io, "# reverse order checks manifest ordering")
        println(io, "conf_00000002", extension)
        println(io)
        println(io, "conf_00000001", extension, " # inline comment")
    end
    return manifest
end

function sequence_bridge_data(configuration)
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

@testset "Wizard configuration sequence runtime" begin
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        element_type=ComplexF64,
        verbose=0,
    )

    mktempdir() do root
        sources = Dict(
            colors => make_sequence_sources(root, colors, environment)
            for colors in (2, 3)
        )
        files = sort(filter(
            path -> occursin("_fileloading_", basename(path)),
            readdir(SEQUENCE_FIXTURE_DIRECTORY; join=true),
        ))
        @test length(files) == 48
        executed = 0

        for path in files
            colors = startswith(basename(path), "su2_") ? 2 : 3
            source = sources[colors]
            raw = TOML.parsefile(path)
            control = raw["System Control"]
            format = control["loadU_format"]
            list_mode = get(control, "loadU_fromfile", false)
            if Sys.iswindows() && format == "ILDG"
                continue
            end

            manifest = joinpath(source.configuration_directory, "filelist.txt")
            if list_mode
                write_sequence_manifest(source.configuration_directory, format)
                expected_indices = (2, 1)
            else
                isfile(manifest) && rm(manifest; force=true)
                expected_indices = (1, 2)
            end

            cd(source.work_directory) do
                parameters = construct_sequence_parameters(path)
                try
                    input = LQCDConfig(parameters)
                    simulation = build_simulation(input, environment)

                    @test input.update isa ConfigurationSequenceConfig
                    @test simulation.updater isa ConfigurationSequenceUpdater
                    @test simulation.state isa ConfigurationSequenceState
                    @test all(
                        isconcretetype,
                        fieldtypes(typeof(simulation.updater)),
                    )
                    @test simulation.state.trajectory == 0
                    @test simulation.state.current_index == 1
                    @test has_next_configuration(simulation)
                    @test basename(current_configuration_path(simulation)) ==
                          "conf_" * lpad(expected_indices[1], 8, '0') *
                          sequence_extension(format)
                    @test Gaugefields.measure_plaquette(
                        simulation.configuration.gauge,
                    ) ≈ source.plaquettes[expected_indices[1]]
                    @test sequence_bridge_data(simulation.configuration) ==
                          source.bridge_data[expected_indices[1]]

                    result = update!(simulation)
                    @test result isa ConfigurationSequenceUpdateResult
                    @test result.trajectory == 0
                    @test result.current_index == 2
                    @test simulation.state.trajectory == 1
                    @test simulation.state.current_index == 2
                    @test !has_next_configuration(simulation)
                    @test basename(current_configuration_path(simulation)) ==
                          "conf_" * lpad(expected_indices[2], 8, '0') *
                          sequence_extension(format)
                    @test Gaugefields.measure_plaquette(
                        simulation.configuration.gauge,
                    ) ≈ source.plaquettes[expected_indices[2]]
                    @test sequence_bridge_data(simulation.configuration) ==
                          source.bridge_data[expected_indices[2]]
                    @test_throws EOFError update!(simulation)
                    @test simulation.state.trajectory == 1
                    @test simulation.state.current_index == 2
                    executed += 1
                finally
                    isopen(parameters.load_fp) && close(parameters.load_fp)
                end
            end
        end

        @test executed == (Sys.iswindows() ? 32 : 48)
    end
end

@testset "Configuration sequence input failures" begin
    environment = GaugefieldsEnvironment(process_grid=(1, 1, 1, 1))
    mktempdir() do directory
        empty_source = DirectorySourceConfig(directory)
        input = LQCDConfig(
            LatticeConfig((2, 2, 2, 2)),
            GaugeConfig(2, 0, "cold", nothing),
            GaugeActionConfig(
                GaugeActionTermConfig(:gauge, "plaquette", 1.9),
            ),
            ConfigurationSequenceConfig(empty_source, :bridge),
        )
        @test_throws ArgumentError build_simulation(input, environment)

        open(joinpath(directory, "filelist.txt"), "w") do io
            println(io, "missing.txt")
        end
        manifest_source = ManifestSourceConfig(directory, "filelist.txt")
        manifest_input = LQCDConfig(
            input.lattice,
            input.gauge,
            input.gauge_action,
            ConfigurationSequenceConfig(manifest_source, :bridge),
        )
        @test_throws ArgumentError build_simulation(manifest_input, environment)

        valid_path = joinpath(directory, "valid.txt")
        invalid_path = joinpath(directory, "invalid.txt")
        source_configuration = build_configuration(
            test_hmc_input(
                initial="hot",
                lattice_size=(2, 2, 2, 2),
                colors=2,
                halo=0,
            ),
            environment,
        )
        Gaugefields.save_configuration(
            valid_path,
            source_configuration.gauge;
            format=:bridge,
        )
        open(invalid_path, "w") do io
            println(io, "not a Bridge configuration")
        end
        open(joinpath(directory, "filelist.txt"), "w") do io
            println(io, basename(valid_path))
            println(io, basename(invalid_path))
        end

        failing_simulation = build_simulation(manifest_input, environment)
        @test_throws Exception update!(failing_simulation)
        @test failing_simulation.state.trajectory == 0
        @test failing_simulation.state.current_index == 1
        @test basename(current_configuration_path(failing_simulation)) ==
              basename(valid_path)
    end
end
