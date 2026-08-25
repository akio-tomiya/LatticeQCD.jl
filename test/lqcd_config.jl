using LatticeQCD
using Test
using TOML

function gauge_only_params(filename, directory)
    parameters = TOML.parsefile(joinpath(@__DIR__, filename))
    control = parameters["System Control"]
    relative_directory = relpath(directory, pwd())
    control["log_dir"] = joinpath(relative_directory, "logs")
    control["logfile"] = "lqcd-config.log"
    control["measurement_basedir"] = joinpath(
        relative_directory,
        "measurements",
    )
    control["measurement_dir"] = "gauge-only"
    return LatticeQCD.Parameters_TOML.construct_Params_from_TOML(parameters)
end

@testset "Typed gauge-only LQCDConfig" begin
    mktempdir() do directory
        parameters = gauge_only_params("test02.toml", directory)
        try
            config = LQCDConfig(parameters)

            @test config.lattice.L == (4, 4, 4, 4)
            @test config.gauge.NC == 2
            @test config.gauge.halo == 1
            @test config.gauge.initialization isa FileStartConfig
            @test config.gauge.initialization.path == parameters.initial
            @test config.gauge.initialization.format == "BridgeText"
            @test length(config.gauge_action.terms) == 1
            @test config.fermions == ()
            gauge_term = only(config.gauge_action.terms)
            @test gauge_term.name == :gauge_plaquette
            @test gauge_term.loop == "plaquette"
            @test gauge_term.coupling == 1.9

            @test typeof(config.lattice) === LatticeConfig{4,Int}
            @test typeof(config.gauge) ===
                  GaugeConfig{FileStartConfig{String,String}}
            gauge_action_type = GaugeActionConfig{
                Tuple{GaugeActionTermConfig{Symbol,String,Float64}}
            }
            @test typeof(config.gauge_action) === gauge_action_type

            @test config.update isa HMCConfig
            @test config.update.md.step_size == parameters.Δτ
            @test config.update.md.steps == parameters.MDsteps
            @test md_trajectory_length(config.update.md) ==
                  parameters.Δτ * parameters.MDsteps
            @test config.update.md.integrator isa LeapfrogConfig
            @test config.update.md.integrator.ordering isa QPQConfig
            @test config.update.md.integrator.forces.names == (:gauge,)
            @test config.update.momentum.sigma == 1.0
            @test config.update.momentum.random.seed == parameters.randomseed
            @test config.update.momentum.random.name == :momentum
            @test config.update.acceptance.random.seed == parameters.randomseed
            @test config.update.acceptance.random.name == :metropolis
            @test config.update.pseudofermions == ()

            @test fieldtype(typeof(config), :lattice) === typeof(config.lattice)
            @test fieldtype(typeof(config), :gauge) === typeof(config.gauge)
            @test fieldtype(typeof(config), :gauge_action) ===
                  typeof(config.gauge_action)
            @test fieldtype(typeof(config), :fermions) === Tuple{}
            @test fieldtype(typeof(config), :update) === typeof(config.update)
            for component in (
                config,
                config.lattice,
                config.gauge,
                config.gauge.initialization,
                config.gauge_action,
                gauge_term,
                config.update,
                config.update.md,
                config.update.md.integrator,
                config.update.md.integrator.forces,
                config.update.momentum,
                config.update.momentum.random,
                config.update.acceptance,
                config.update.acceptance.random,
            )
                @test all(isconcretetype, fieldtypes(typeof(component)))
            end

            output = sprint(show_config, config)
            @test occursin("LQCDConfig", output)
            @test occursin("gauge_plaquette", output)
            @test occursin("file start", output)
            @test occursin("format: \"BridgeText\"", output)
            @test occursin("include adjoint=true", output)
            @test occursin("fermion actions:\n    none", output)
            @test occursin("update:", output)
            @test occursin("HMC", output)
            @test occursin("scheme: QPQ", output)
            @test occursin("MD steps: 15", output)

            compact_output = sprint(show, config)
            @test occursin("update=HMCConfig", compact_output)
        finally
            isopen(parameters.load_fp) && close(parameters.load_fp)
        end
    end

    cold = GaugeConfig(3, 1, "cold", nothing)
    @test typeof(cold) === GaugeConfig{ColdStartConfig}
    @test cold.initialization isa ColdStartConfig

    hot = GaugeConfig(3, 1, "hot", nothing, 1234)
    @test hot.initialization isa HotStartConfig{Int}
    @test hot.initialization.seed == 1234

    instanton = GaugeConfig(2, 1, "one instanton", nothing)
    @test instanton.initialization isa InstantonConfig

    embedded = GaugeConfig(
        3,
        1,
        EmbeddedInstantonConfig(
            center=(2.5, 2.5, 2.5, 2.5),
            radius=2.0,
            sign=-1,
            block=(2, 3),
        ),
    )
    @test embedded.initialization isa EmbeddedInstantonConfig
    @test embedded.initialization.sign == -1
    @test embedded.initialization.block == (2, 3)
    @test all(isconcretetype, fieldtypes(typeof(embedded)))
    @test all(isconcretetype, fieldtypes(typeof(embedded.initialization)))
    @test_throws ArgumentError EmbeddedInstantonConfig(sign=0)
    @test_throws ArgumentError EmbeddedInstantonConfig(radius=0)
    @test_throws ArgumentError EmbeddedInstantonConfig(block=(1, 2, 3))

    improved = GaugeActionConfig(
        GaugeActionTermConfig(:plaquette, "plaquette", 5.5),
        GaugeActionTermConfig(:rectangle, "rectangle", -0.1),
    )
    @test length(improved.terms) == 2
    @test improved.terms[2].loop == "rectangle"

    nested = SextonWeingartenConfig(
        QPQConfig(),
        ForceGroupConfig(:fermion_1),
        LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge)),
        2,
    )
    @test nested.slow_forces.names == (:fermion_1,)
    @test nested.fast_integrator.forces.names == (:gauge,)
    @test nested.fast_steps == 2
    odd_nested = SextonWeingartenConfig(
        QPQConfig(),
        ForceGroupConfig(:fermion_1),
        LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge)),
        3,
    )
    @test odd_nested.fast_steps == 3
    @test_throws ArgumentError SextonWeingartenConfig(
        QPQConfig(),
        ForceGroupConfig(:fermion_1),
        LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge)),
        0,
    )

    hisq = HISQDiracConfig(0.17, -0.083)
    @test hisq.mass == 0.17
    @test hisq.naik_epsilon == -0.083
    @test all(isconcretetype, fieldtypes(typeof(hisq)))
    @test_throws ArgumentError HISQDiracConfig(0.0, -0.083)
    @test_throws ArgumentError HISQDiracConfig(0.17, Inf)

    mobius = MobiusDomainwallDiracConfig(0.1, -1.0, 4, 2.0, 1.0)
    @test mobius.mass == 0.1
    @test mobius.domainwall_height == -1.0
    @test mobius.fifth_dimension == 4
    @test mobius.b == 2.0
    @test mobius.c == 1.0
    @test all(isconcretetype, fieldtypes(typeof(mobius)))
    @test_throws ArgumentError MobiusDomainwallDiracConfig(
        0.1, -1.0, 0, 2.0, 1.0)
    @test_throws ArgumentError MobiusDomainwallDiracConfig(
        0.1, -1.0, 4, Inf, 1.0)
end

@testset "Typed gauge-only Heatbath config from Params" begin
    mktempdir() do directory
        parameters = gauge_only_params("test02-hb.toml", directory)
        try
            config = LQCDConfig(parameters)
            update = config.update

            @test update isa HeatbathConfig
            @test update.even_odd == parameters.isevenodd
            @test update.max_iterations == parameters.ITERATION_MAX
            @test update.overrelaxation_steps == 0
            @test update.random.seed == parameters.randomseed
            @test update.random.name == :heatbath
            @test fieldtype(typeof(config), :update) === typeof(update)
            @test all(isconcretetype, fieldtypes(typeof(update)))

            output = sprint(show_config, config)
            @test occursin("Heatbath", output)
            @test occursin("overrelaxation steps: 0", output)
        finally
            isopen(parameters.load_fp) && close(parameters.load_fp)
        end
    end
end


@testset "Typed configuration sequence config from Params" begin
    cases = (
        (
            "wizard_gauge_only/su3_fileloading_jld2_all_save_none.toml",
            DirectorySourceConfig,
            :jld2,
        ),
        (
            "wizard_gauge_only/su2_fileloading_bridge_list_save_none.toml",
            ManifestSourceConfig,
            :bridge,
        ),
    )

    mktempdir() do directory
        for (filename, source_type, format) in cases
            parameters = gauge_only_params(filename, directory)
            try
                config = LQCDConfig(parameters)
                sequence = config.update

                @test sequence isa ConfigurationSequenceConfig
                @test sequence.source isa source_type
                @test sequence.source.directory == parameters.loadU_dir
                @test sequence.format === format
                @test fieldtype(typeof(config), :update) === typeof(sequence)
                @test all(isconcretetype, fieldtypes(typeof(sequence)))
                @test all(isconcretetype, fieldtypes(typeof(sequence.source)))

                output = sprint(show_config, config)
                @test occursin("Configuration sequence", output)
                @test occursin("format: $format", output)
            finally
                isopen(parameters.load_fp) && close(parameters.load_fp)
            end
        end
    end

    directory = DirectorySourceConfig("confs")
    manifest = ManifestSourceConfig("confs", "filelist.txt")
    @test ConfigurationSequenceConfig(directory, "JLD").format === :jld2
    @test ConfigurationSequenceConfig(directory, "ILDG").format === :ildg
    @test ConfigurationSequenceConfig(manifest, "BridgeText").format === :bridge
    @test_throws ArgumentError DirectorySourceConfig("")
    @test_throws ArgumentError ManifestSourceConfig("confs", "")
    @test_throws ArgumentError ConfigurationSequenceConfig(directory, "raw")
end


@testset "Gauge-only Sexton-Weingarten compatibility" begin
    mktempdir() do directory
        sw_parameters = gauge_only_params(
            "wizard_gauge_only/su3_hot_hmc_sw_save_none.toml",
            directory,
        )
        leapfrog_parameters = gauge_only_params(
            "wizard_gauge_only/su3_hot_hmc_no_sw_save_none.toml",
            directory,
        )
        try
            sw_config = @test_logs (
                :warn,
                r"Sexton-Weingarten is ignored for a gauge-only simulation",
            ) LQCDConfig(sw_parameters)
            leapfrog_config = LQCDConfig(leapfrog_parameters)

            @test sw_config.update isa HMCConfig
            @test sw_config.update.md.integrator isa LeapfrogConfig
            @test sw_config.update.md.step_size ==
                  leapfrog_config.update.md.step_size
            @test sw_config.update.md.steps == leapfrog_config.update.md.steps
            @test typeof(sw_config.update.md.integrator.ordering) ===
                  typeof(leapfrog_config.update.md.integrator.ordering)
            @test sw_config.update.md.integrator.forces.names == (:gauge,)
        finally
            isopen(sw_parameters.load_fp) && close(sw_parameters.load_fp)
            isopen(leapfrog_parameters.load_fp) &&
                close(leapfrog_parameters.load_fp)
        end
    end
end
