using LatticeQCD
using Test
using TOML
import Gaugefields

function slhmc_gauge_hmc_input()
    lattice = LatticeConfig((2, 2, 2, 2))
    gauge = GaugeConfig(2, 1, "hot", nothing, 0x1234)
    action = GaugeActionConfig(
        GaugeActionTermConfig(:gauge_plaquette, "plaquette", 1.9),
    )
    md = MDConfig(
        0.02,
        2,
        LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge)),
    )
    update = HMCConfig(
        md,
        GaussianMomentumConfig(
            sqrt(2.0),
            RandomStreamConfig(0x5678, :momentum),
        ),
        RankZeroMetropolisConfig(
            RandomStreamConfig(0x9abc, :metropolis),
        ),
    )
    return LQCDConfig(lattice, gauge, action, update)
end

function slhmc_dynamical_hmc_input()
    gauge_input = slhmc_gauge_hmc_input()
    gauge = GaugeConfig(2, 1, "cold", nothing, 0x1234)
    solver = FermionSolverConfig(1e-10, 2_000, 0, (1, 1, 1, -1))
    fermions = (
        FermionActionConfig(
            :fermion_1,
            WilsonDiracConfig(0.05, 1.0),
            2,
            solver,
        ),
    )
    md = MDConfig(
        0.001,
        1,
        LeapfrogConfig(
            QPQConfig(),
            ForceGroupConfig(:gauge, :fermion_1),
        ),
    )
    update = HMCConfig(
        md,
        gauge_input.update.momentum,
        gauge_input.update.acceptance,
        (
            PseudofermionRefreshConfig(
                :fermion_1,
                RandomStreamConfig(0xdef0, :pseudofermion),
                1,
            ),
        ),
    )
    return LQCDConfig(
        gauge_input.lattice,
        gauge,
        gauge_input.gauge_action,
        fermions,
        update,
    )
end

function slhmc_from_hmc(input::LQCDConfig, md_gauge_action; md_fermions=input.fermions)
    hmc = input.update
    update = SLHMCConfig(
        hmc.md,
        hmc.momentum,
        hmc.acceptance,
        hmc.pseudofermions,
        md_gauge_action,
        md_fermions,
    )
    return LQCDConfig(
        input.lattice,
        input.gauge,
        input.gauge_action,
        input.fermions,
        update,
    )
end

@testset "SLHMC target and MD actions" begin
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        communicator=Gaugefields.SerialCommunicator(),
        element_type=ComplexF64,
        verbose=0,
    )

    hmc_input = slhmc_gauge_hmc_input()
    identical_input = slhmc_from_hmc(
        hmc_input,
        hmc_input.gauge_action;
        md_fermions=(),
    )
    hmc_simulation = build_simulation(hmc_input, environment)
    identical_simulation = build_simulation(identical_input, environment)

    @test identical_simulation.updater isa SLHMCUpdater
    @test identical_simulation.updater.md_driver.momentum_denominator ==
          2.0
    @test identical_simulation.updater.target_driver.momentum_denominator ==
          2.0
    @test identical_simulation.state isa HMCState
    @test all(isconcretetype, fieldtypes(typeof(identical_input.update)))
    @test all(isconcretetype, fieldtypes(typeof(identical_simulation.updater)))
    @test occursin("target action", sprint(show_config, identical_input))

    hmc_result = update!(hmc_simulation)
    identical_result = update!(identical_simulation)
    @test identical_result isa SLHMCUpdateResult
    @test identical_result.accepted == hmc_result.accepted
    @test identical_result.initial_hamiltonian ≈ hmc_result.initial_hamiltonian
    @test identical_result.final_hamiltonian ≈ hmc_result.final_hamiltonian
    @test identical_result.delta_hamiltonian ≈ hmc_result.delta_hamiltonian
    @test identical_result.md_initial_hamiltonian ≈
          identical_result.initial_hamiltonian
    @test identical_result.md_final_hamiltonian ≈
          identical_result.final_hamiltonian
    @test Gaugefields.measure_plaquette(
        identical_simulation.configuration.gauge,
    ) ≈ Gaugefields.measure_plaquette(hmc_simulation.configuration.gauge)

    approximate_gauge_action = GaugeActionConfig(
        GaugeActionTermConfig(:gauge_plaquette, "plaquette", 0.7),
    )
    approximate_simulation = build_simulation(
        slhmc_from_hmc(
            hmc_input,
            approximate_gauge_action;
            md_fermions=(),
        ),
        environment,
    )
    approximate_result = update!(approximate_simulation)
    @test isfinite(approximate_result.delta_hamiltonian)
    @test isfinite(approximate_result.md_delta_hamiltonian)
    @test approximate_result.delta_hamiltonian ≈
          approximate_result.final_hamiltonian -
          approximate_result.initial_hamiltonian
    @test approximate_result.md_delta_hamiltonian ≈
          approximate_result.md_final_hamiltonian -
          approximate_result.md_initial_hamiltonian
    @test !(approximate_result.delta_hamiltonian ≈
            approximate_result.md_delta_hamiltonian)

    dynamical_hmc = slhmc_dynamical_hmc_input()
    target_fermion = only(dynamical_hmc.fermions)
    md_fermions = (
        FermionActionConfig(
            :fermion_1,
            WilsonDiracConfig(0.04, 1.0),
            target_fermion.flavors,
            target_fermion.solver,
            target_fermion.smearing,
        ),
    )
    dynamical_slhmc = build_simulation(
        slhmc_from_hmc(
            dynamical_hmc,
            dynamical_hmc.gauge_action;
            md_fermions,
        ),
        environment,
    )
    dynamical_result = update!(dynamical_slhmc)
    @test dynamical_result isa SLHMCUpdateResult
    @test isfinite(dynamical_result.delta_hamiltonian)
    @test length(dynamical_slhmc.updater.pseudofermion_refreshes) == 1

    wrong_name = (
        FermionActionConfig(
            :wrong_name,
            target_fermion.operator,
            target_fermion.flavors,
            target_fermion.solver,
            target_fermion.smearing,
        ),
    )
    @test_throws ArgumentError build_simulation(
        slhmc_from_hmc(
            dynamical_hmc,
            dynamical_hmc.gauge_action;
            md_fermions=wrong_name,
        ),
        environment,
    )
end

@testset "Legacy SLMC Params convert to typed SLHMC" begin
    mktempdir() do directory
        parameters = TOML.parsefile(joinpath(@__DIR__, "test06_slmc_ks.toml"))
        parameters["SLHMC related"] = Dict("βeff" => 4.2)
        parameters["Physical setting(fermions)"]["quench"] = false
        control = parameters["System Control"]
        relative_directory = relpath(directory, pwd())
        control["log_dir"] = joinpath(relative_directory, "logs")
        control["measurement_basedir"] = joinpath(
            relative_directory,
            "measurements",
        )
        params = LatticeQCD.Parameters_TOML.construct_Params_from_TOML(parameters)
        try
            @test params.update_method == "SLMC"
            @test !params.quench
            input = LQCDConfig(params)
            @test input.update isa SLHMCConfig
            @test only(input.update.md_gauge_action.terms).coupling == 4.2
            @test input.update.md_fermions == input.fermions
        finally
            isopen(params.load_fp) && close(params.load_fp)
        end
    end
end
