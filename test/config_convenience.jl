using LatticeQCD, Test, TOML
import Logging

@testset "show_config stdout convenience methods" begin
    config = LQCDConfig(
        LatticeConfig((2, 2, 2, 2)), GaugeConfig(2, 0, "cold", nothing, 11),
        GaugeActionConfig(GaugeActionTermConfig(:gauge, "plaquette", 1.9)),
        HeatbathConfig(true, 100_000, 0, RandomStreamConfig(11, :heatbath)),
    )
    spec = SimulationSpec(config, SimulationSchedule(0, 1), OutputConfig())
    session = build_simulation(spec)
    for value in (spec, session, spec.schedule, spec.output, config, config.update)
        mktemp() do path, io
            result = redirect_stdout(io) do
                show_config(value)
            end
            @test result === nothing
            flush(io)
            @test read(path, String) == sprint(show_config, value)
        end
    end
end

@testset "Automatic legacy TOML halos" begin
    for (fixture, required) in (
        ("su2_wilson_thin_leapfrog.toml", 1),
        ("su2_wilson_clover_thin_leapfrog.toml", 1),
        ("su2_staggered_nf4_thin_leapfrog.toml", 1),
        ("su2_domainwall_thin_leapfrog.toml", 1),
        ("su2_mobius_domainwall_thin_leapfrog.toml", 1),
        ("su2_hisq_nf4_thin_leapfrog.toml", 3),
    )
        document = TOML.parsefile(joinpath(@__DIR__, "wizard_fermion", fixture))
        physical = document["Physical setting"]
        delete!(physical, "Nwing")
        original = deepcopy(document)
        spec = @test_logs min_level=Logging.Warn simulation_spec_from_legacy_toml(document)
        @test spec.config.gauge.halo == required
        @test document == original
        @test isempty(validate(spec))
        for explicit in 0:required-1
            physical["Nwing"] = explicit
            adjusted = @test_logs (:warn, r"Nwing=.*legacy dynamical-fermion") simulation_spec_from_legacy_toml(document)
            @test adjusted.config.gauge.halo == required
        end
        physical["Nwing"] = required + 1
        preserved = @test_logs min_level=Logging.Warn simulation_spec_from_legacy_toml(document)
        @test preserved.config.gauge.halo == required + 1
    end

    document = TOML.parsefile(joinpath(@__DIR__, "wizard_gauge_only",
                                     "su2_cold_heatbath_no_or_save_none.toml"))
    physical = document["Physical setting"]
    delete!(physical, "Nwing")
    gauge_only = @test_logs min_level=Logging.Warn simulation_spec_from_legacy_toml(document)
    @test gauge_only.config.gauge.halo == 0
    pion = Dict("methodname" => "Pion_correlator", "fermiontype" => "Wilson", "hop" => 0.164)
    for flowed in (false, true)
        valence = deepcopy(document)
        if flowed
            valence["gradientflow_measurements"] = Dict(
                "hasgradientflow" => true,
                "measurements_for_flow" => Dict("pion" => pion),
            )
        else
            valence["Measurement set"]["measurement_methods"] = Dict("pion" => pion)
        end
        automatic = @test_logs min_level=Logging.Warn simulation_spec_from_legacy_toml(valence)
        @test automatic.config.gauge.halo == 1
        @test isempty(validate(automatic))
        valence["Physical setting"]["Nwing"] = 0
        adjusted = @test_logs (:warn, r"Nwing=0 is too small for fermion measurements") simulation_spec_from_legacy_toml(valence)
        @test adjusted.config.gauge.halo == 1
    end
end
