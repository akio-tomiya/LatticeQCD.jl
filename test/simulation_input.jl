using LatticeQCD
using Test
import TOML
import Logging
import Gaugefields

function simulation_input_structurally_equal(left, right)
    typeof(left) == typeof(right) || return false
    if left isa Union{Number,AbstractString,Symbol,Nothing,Bool}
        return isequal(left, right)
    elseif left isa Union{Tuple,AbstractArray}
        return length(left) == length(right) && all(
            simulation_input_structurally_equal(l, r) for
            (l, r) in zip(left, right)
        )
    end
    return all(
        simulation_input_structurally_equal(
            getfield(left, name),
            getfield(right, name),
        ) for name in fieldnames(typeof(left))
    )
end

function legacy_params_for_input_test(document, directory, index)
    values = deepcopy(document)
    control = values["System Control"]
    relative_directory = relpath(directory, pwd())
    control["log_dir"] = joinpath(relative_directory, "logs_$index")
    control["logfile"] = "legacy.log"
    control["measurement_basedir"] = joinpath(
        relative_directory,
        "measurements_$index",
    )
    control["measurement_dir"] = "run"
    pushdisplay(TextDisplay(devnull))
    try
        return redirect_stdout(devnull) do
            Logging.with_logger(Logging.NullLogger()) do
                LatticeQCD.Parameters_TOML.construct_Params_from_TOML(values)
            end
        end
    finally
        popdisplay()
    end
end

function wizard_input_files()
    roots = (
        joinpath(@__DIR__, "wizard_gauge_only"),
        joinpath(@__DIR__, "wizard_fermion"),
        joinpath(@__DIR__, "wizard_fermion_io"),
    )
    return sort(vcat([
        filter(path -> endswith(path, ".toml"), readdir(root; join=true))
        for root in roots
    ]...))
end

function historical_input_files()
    return sort(filter(
        path -> endswith(path, ".toml"),
        readdir(@__DIR__; join=true),
    ))
end

@testset "Param-free SimulationSpec TOML input" begin
    @testset "legacy defaults remain identical" begin
        mktempdir(pwd()) do directory
            relative_directory = relpath(directory, pwd())
            document = Dict{String,Any}(
                # The legacy Params constructor treats L as a mandatory key:
                # its struct default is a Vector although Params stores a
                # Tuple. Keep that historical requirement in this oracle and
                # test every other omitted value against the old path.
                "Physical setting" => Dict{String,Any}(
                    "L" => [4, 4, 4, 4],
                ),
                "Physical setting(fermions)" => Dict{String,Any}(),
                "System Control" => Dict{String,Any}(
                    "log_dir" => joinpath(relative_directory, "logs"),
                    "logfile" => "legacy.log",
                    "measurement_basedir" => joinpath(
                        relative_directory,
                        "measurements",
                    ),
                    "measurement_dir" => "run",
                ),
                "HMC related" => Dict{String,Any}(),
                "Measurement set" => Dict{String,Any}(
                    "measurement_methods" => Dict{String,Any}(),
                ),
                "gradientflow_measurements" => Dict{String,Any}(
                    "measurements_for_flow" => Dict{String,Any}(),
                ),
            )
            parameters = legacy_params_for_input_test(document, directory, 0)
            old_spec = Logging.with_logger(Logging.NullLogger()) do
                SimulationSpec(parameters)
            end
            close(parameters.load_fp)
            new_spec = Logging.with_logger(Logging.NullLogger()) do
                simulation_spec_from_legacy_toml(document)
            end
            @test simulation_input_structurally_equal(old_spec, new_spec)
            @test new_spec.config.lattice.L == (4, 4, 4, 4)
            @test new_spec.config.gauge.NC == 3
            @test new_spec.config.gauge.halo == 0
            @test new_spec.config.gauge.initialization isa ColdStartConfig
            @test new_spec.config.gauge_action.terms[1].coupling == 5.7
            @test isempty(new_spec.config.fermions)
            @test new_spec.config.update.md.step_size == 0.05
            @test new_spec.config.update.md.steps == 20
            @test new_spec.schedule.thermalization_steps == 0
            @test new_spec.schedule.production_steps == 100
            @test new_spec.output.configurations isa NoConfigurationOutput
        end
    end

    @testset "all current Wizard files match Params" begin
        files = wizard_input_files()
        @test length(files) == 464
        mktempdir(pwd()) do directory
            for (index, filename) in enumerate(files)
                document = TOML.parsefile(filename)
                parameters = legacy_params_for_input_test(
                    document,
                    directory,
                    index,
                )
                old_spec = Logging.with_logger(Logging.NullLogger()) do
                    SimulationSpec(parameters)
                end
                close(parameters.load_fp)
                new_spec = Logging.with_logger(Logging.NullLogger()) do
                    load_simulation_spec(filename)
                end
                @test simulation_input_structurally_equal(old_spec, new_spec)
            end
        end
    end

    @testset "historical parameter files match Params" begin
        files = historical_input_files()
        @test length(files) == 14
        mktempdir(pwd()) do directory
            for (index, filename) in enumerate(files)
                document = TOML.parsefile(filename)
                parameters = legacy_params_for_input_test(
                    document,
                    directory,
                    1_000 + index,
                )
                old_spec = Logging.with_logger(Logging.NullLogger()) do
                    SimulationSpec(parameters)
                end
                close(parameters.load_fp)
                new_spec = Logging.with_logger(Logging.NullLogger()) do
                    load_simulation_spec(filename)
                end
                @test simulation_input_structurally_equal(old_spec, new_spec)
            end
        end
    end

    @testset "canonical TOML round trip" begin
        representatives = (
            joinpath(
                @__DIR__,
                "wizard_fermion",
                "su3_wilson_clover_thin_leapfrog.toml",
            ),
            joinpath(
                @__DIR__,
                "wizard_fermion",
                "su3_hisq_nf4_thin_leapfrog.toml",
            ),
            joinpath(
                @__DIR__,
                "wizard_fermion",
                "su3_domainwall_stout_sw.toml",
            ),
            joinpath(
                @__DIR__,
                "wizard_fermion",
                "su3_mobius_domainwall_thin_leapfrog.toml",
            ),
            joinpath(
                @__DIR__,
                "wizard_fermion",
                "su2_slhmc_staggered_nf3_stout_leapfrog.toml",
            ),
            joinpath(
                @__DIR__,
                "wizard_gauge_only",
                "su3_embedded_instanton_heatbath_or_save_jld2.toml",
            ),
            joinpath(
                @__DIR__,
                "wizard_gauge_only",
                "su2_fileloading_bridge_list_save_none.toml",
            ),
        )
        mktempdir() do directory
            for (index, filename) in enumerate(representatives)
                spec = load_simulation_spec(filename)
                output = joinpath(directory, "spec_$index.toml")
                @test write_simulation_spec(output, spec) == output
                document = TOML.parsefile(output)
                @test document["format"] == SIMULATION_SPEC_FORMAT
                @test document["schema_version"] ==
                      SIMULATION_SPEC_SCHEMA_VERSION
                restored = load_simulation_spec(output)
                @test simulation_input_structurally_equal(spec, restored)
            end

            base = load_simulation_spec(first(representatives))
            checkpoint_spec = SimulationSpec(
                base.config,
                base.schedule,
                OutputConfig(
                    base.output.configurations,
                    JLD2CheckpointOutput(
                        "restart";
                        prefix="checkpoint_",
                        every=17,
                        width=6,
                    ),
                ),
            )
            checkpoint_toml = joinpath(directory, "checkpoint_spec.toml")
            write_simulation_spec(checkpoint_toml, checkpoint_spec)
            checkpoint_document = TOML.parsefile(checkpoint_toml)
            checkpoint_values = checkpoint_document["output"]["checkpoints"]
            @test checkpoint_values["format"] == "jld2"
            @test checkpoint_values["every"] == 17
            @test checkpoint_values["prefix"] == "checkpoint_"
            restored_checkpoint_spec = load_simulation_spec(checkpoint_toml)
            @test simulation_input_structurally_equal(
                checkpoint_spec,
                restored_checkpoint_spec,
            )
        end
    end

    @testset "one HMC trajectory matches Params path" begin
        filename = joinpath(
            @__DIR__,
            "wizard_gauge_only",
            "su2_hot_hmc_no_sw_save_none.toml",
        )
        document = TOML.parsefile(filename)
        document["Physical setting"]["L"] = [2, 2, 2, 2]
        document["Physical setting"]["Nsteps"] = 1
        document["HMC related"]["Δτ"] = 0.001
        document["HMC related"]["MDsteps"] = 1
        mktempdir(pwd()) do directory
            parameters = legacy_params_for_input_test(document, directory, 500)
            @test parameters isa Params
            old_spec = SimulationSpec(parameters)
            close(parameters.load_fp)
            new_spec = simulation_spec_from_legacy_toml(document)
            environment = GaugefieldsEnvironment(
                communicator=Gaugefields.SerialCommunicator(),
                process_grid=(1, 1, 1, 1),
                element_type=ComplexF64,
                verbose=0,
            )
            old_simulation = build_simulation(old_spec.config, environment)
            new_simulation = build_simulation(new_spec.config, environment)
            old_result = update!(old_simulation)
            new_result = update!(new_simulation)
            @test old_result.accepted == new_result.accepted
            @test old_result.delta_hamiltonian == new_result.delta_hamiltonian
            @test old_result.initial_hamiltonian ==
                  new_result.initial_hamiltonian
            @test old_result.final_hamiltonian == new_result.final_hamiltonian
            maximum_difference = maximum(
                maximum(
                    abs,
                    old_simulation.configuration.gauge[direction].U.A .-
                    new_simulation.configuration.gauge[direction].U.A,
                ) for direction in eachindex(
                    old_simulation.configuration.gauge,
                )
            )
            @test maximum_difference == 0
        end
    end

    @testset "direct legacy read has no filesystem side effects" begin
        filename = joinpath(
            @__DIR__,
            "wizard_gauge_only",
            "su2_hot_hmc_no_sw_save_jld2.toml",
        )
        mktempdir() do directory
            before = readdir(directory)
            cd(directory) do
                spec = load_simulation_spec(filename)
                @test spec isa SimulationSpec
            end
            @test readdir(directory) == before
        end
    end
end
