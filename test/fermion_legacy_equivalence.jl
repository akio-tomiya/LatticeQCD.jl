using LatticeQCD
using Test
import Gaugefields
import LatticeDiracOperators

const LEGACY_MD = LatticeQCD.AbstractMD_module

function legacy_fermion_action(config::FermionActionConfig, gauge, field)
    solver = config.solver
    parameters = Dict{String,Any}(
        "eps_CG" => solver.tolerance,
        "MaxCGstep" => solver.max_steps,
        "verbose_level" => solver.verbose,
        "boundarycondition" => collect(solver.boundary_conditions),
    )
    operator = config.operator
    if operator isa WilsonDiracConfig
        parameters["Dirac_operator"] = "Wilson"
        parameters["κ"] = operator.hopping_parameter
        parameters["r"] = operator.wilson_parameter
        parameters["faster version"] = true
    elseif operator isa WilsonCloverDiracConfig
        parameters["Dirac_operator"] = "WilsonClover"
        parameters["κ"] = operator.hopping_parameter
        parameters["r"] = operator.wilson_parameter
        parameters["cSW"] = operator.clover_coefficient
    elseif operator isa StaggeredDiracConfig
        parameters["Dirac_operator"] = "staggered"
        parameters["mass"] = operator.mass
    elseif operator isa HISQDiracConfig
        parameters["Dirac_operator"] = "HISQ"
        parameters["mass"] = operator.mass
        parameters["naik_epsilon"] = operator.naik_epsilon
    elseif operator isa DomainwallDiracConfig
        parameters["Dirac_operator"] = "Domainwall"
        parameters["mass"] = operator.mass
        parameters["M"] = operator.domainwall_height
        parameters["L5"] = operator.fifth_dimension
    elseif operator isa MobiusDomainwallDiracConfig
        parameters["Dirac_operator"] = "MobiusDomainwall"
        parameters["mass"] = operator.mass
        parameters["M"] = operator.domainwall_height
        parameters["L5"] = operator.fifth_dimension
        parameters["b"] = operator.b
        parameters["c"] = operator.c
    else
        error("unsupported legacy comparison operator $(typeof(operator))")
    end

    dirac = LatticeDiracOperators.Dirac_operator(gauge, field, parameters)
    action_parameters = operator isa Union{StaggeredDiracConfig,HISQDiracConfig} ?
        Dict{String,Any}("Nf" => config.flavors) : Dict{String,Any}()
    return LatticeDiracOperators.FermiAction(dirac, action_parameters)
end

function maximum_link_difference(left, right)
    return maximum(
        maximum(abs, left[direction].U.A .- right[direction].U.A)
        for direction in eachindex(left)
    )
end

function maximum_momentum_difference(left, right)
    return maximum(
        maximum(abs, left[direction].a.A .- right[direction].a.A)
        for direction in eachindex(left)
    )
end

function copy_momenta!(destination, source)
    for direction in eachindex(destination, source)
        copyto!(destination[direction].a.A, source[direction].a.A)
    end
    return destination
end

legacy_fermion_smearing(::NoFermionSmearingConfig, gauge) = nothing

function legacy_fermion_smearing(
    config::StoutFermionSmearingConfig,
    gauge,
)
    smearing = Gaugefields.CovNeuralnet(gauge)
    loops = collect(config.loops)
    coefficients = collect(config.coefficients)
    for _ in 1:config.layers
        push!(
            smearing,
            Gaugefields.STOUT_Layer(loops, coefficients, gauge),
        )
    end
    return smearing
end

@testset "Legacy and typed fermion action/force equivalence" begin
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        communicator=Gaugefields.SerialCommunicator(),
        element_type=ComplexF64,
        verbose=0,
    )
    cases = (
        (WilsonDiracConfig(0.05, 1.0), 2),
        (WilsonCloverDiracConfig(0.05, 1.0, 1.0), 2),
        (StaggeredDiracConfig(0.5), 1),
        (StaggeredDiracConfig(0.5), 2),
        (StaggeredDiracConfig(0.5), 3),
        (StaggeredDiracConfig(0.5), 4),
        (StaggeredDiracConfig(0.5), 8),
        (HISQDiracConfig(0.5, -0.083), 4),
        (DomainwallDiracConfig(0.1, -1.0, 2), 2),
        (MobiusDomainwallDiracConfig(0.1, -1.0, 2, 2.0, 1.0), 2),
    )

    for (operator, flavors) in cases
        input = dynamical_test_input(operator; flavors)
        configuration = build_configuration(input, environment)
        gauge = configuration.gauge
        field = only(values(configuration.fermions))
        config = only(input.fermions)

        new_action = redirect_stdout(devnull) do
            build_fermion_action(config, gauge, field)
        end
        old_action = redirect_stdout(devnull) do
            legacy_fermion_action(config, gauge, field)
        end
        provider = LatticeDiracOperators.PseudofermionMDAction(
            new_action,
            field,
        )
        noise = similar(field)
        LatticeDiracOperators.refresh_pseudofermion!(
            provider,
            gauge,
            noise;
            seed=0x314159,
            sweep=7,
            subgroup=1,
        )

        new_value = LatticeDiracOperators.evaluate_FermiAction(
            new_action,
            gauge,
            field,
        )
        old_value = LatticeDiracOperators.evaluate_FermiAction(
            old_action,
            gauge,
            field,
        )
        @test isfinite(new_value)
        @test new_value ≈ old_value rtol=2e-11 atol=2e-11

        new_force = LatticeDiracOperators.calc_UdSfdU(
            new_action,
            gauge,
            field,
        )
        old_force = LatticeDiracOperators.calc_UdSfdU(
            old_action,
            gauge,
            field,
        )
        @test maximum_link_difference(new_force, old_force) < 2e-10
    end
end

@testset "Legacy StandardMD and Gaugefields driver equivalence" begin
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        communicator=Gaugefields.SerialCommunicator(),
        element_type=ComplexF64,
        verbose=0,
    )
    wizard_stout_loops = Tuple(
        LatticeQCD.Parameter_structs.kindsof_loops,
    )
    smearings = (
        NoFermionSmearingConfig(),
        (
            StoutFermionSmearingConfig(
                1,
                (loop == "plaquette" ? 0.1 : 0.02,),
                (loop,),
            ) for loop in wizard_stout_loops
        )...,
    )
    @test length(smearings) == 1 + length(wizard_stout_loops) == 8
    for smearing in smearings, sexton_weingarten in (false, true)
        input = dynamical_test_input(
            WilsonDiracConfig(0.05, 1.0);
            sexton_weingarten,
            smearing,
        )
        new_configuration = build_configuration(input, environment)
        old_configuration = build_configuration(input, environment)
        new_gauge = new_configuration.gauge
        old_gauge = old_configuration.gauge
        new_field = only(values(new_configuration.fermions))
        old_field = only(values(old_configuration.fermions))
        fermion_config = only(input.fermions)

        new_fermion_action = build_fermion_action(
            fermion_config,
            new_gauge,
            new_field,
        )
        old_fermion_action = legacy_fermion_action(
            fermion_config,
            old_gauge,
            old_field,
        )
        old_smearing = legacy_fermion_smearing(
            fermion_config.smearing,
            old_gauge,
        )
        new_provider = LatticeDiracOperators.PseudofermionMDAction(
            new_fermion_action,
            new_field,
        )
        if smearing isa NoFermionSmearingConfig
            @test new_provider.smearing === nothing
        else
            @test new_provider.smearing isa Gaugefields.CovNeuralnet
            @test length(new_provider.smearing.layers) == smearing.layers
        end
        noise = similar(new_field)
        LatticeDiracOperators.refresh_pseudofermion!(
            new_provider,
            new_gauge,
            noise;
            seed=0xabcdef,
            sweep=3,
            subgroup=1,
        )

        new_gauge_action = build_gauge_action(
            input.gauge_action,
            new_configuration,
        )
        old_gauge_action = build_gauge_action(
            input.gauge_action,
            old_configuration,
        )
        old_md = LEGACY_MD.MD(
            old_gauge,
            old_gauge_action,
            false,
            input.update.md.step_size,
            input.update.md.steps,
            old_fermion_action,
            old_smearing;
            QPQ=true,
            SextonWeingargten=sexton_weingarten,
            Nsw=2,
        )
        LatticeDiracOperators.substitute_fermion!(old_md.η, new_field)

        new_momenta = Gaugefields.gaussian_momenta(
            new_gauge;
            seed=UInt64(0x55667788),
            sweep=2,
        )
        copy_momenta!(old_md.p, new_momenta)
        new_actions = Gaugefields.MDActionSet(;
            gauge=new_gauge_action,
            fermion=new_provider,
        )
        old_provider = LatticeDiracOperators.PseudofermionMDAction(
            old_fermion_action,
            old_md.η,
            old_smearing,
        )
        old_actions = Gaugefields.MDActionSet(;
            gauge=old_gauge_action,
            fermion=old_provider,
        )
        integrator = sexton_weingarten ? Gaugefields.SextonWeingarten(;
            slow=:fermion,
            fast=:gauge,
            n_fast=1,
        ) : Gaugefields.QPQ()
        trajectory_length = input.update.md.step_size * input.update.md.steps
        new_driver = Gaugefields.md_driver(
            new_gauge,
            new_actions;
            steps=input.update.md.steps,
            trajectory_length,
            integrator,
        )
        old_diagnostics_driver = Gaugefields.md_driver(
            old_gauge,
            old_actions;
            steps=input.update.md.steps,
            trajectory_length,
            integrator,
        )
        old_initial_hamiltonian = Gaugefields.md_hamiltonian(
            old_gauge,
            old_md.p,
            old_diagnostics_driver,
        )

        new_result = Gaugefields.md_trajectory!(
            new_gauge,
            new_momenta,
            new_driver,
        )
        LEGACY_MD.runMD!(old_gauge, old_md)
        old_final_hamiltonian = Gaugefields.md_hamiltonian(
            old_gauge,
            old_md.p,
            old_diagnostics_driver,
        )
        old_delta = old_final_hamiltonian - old_initial_hamiltonian

        @test maximum_link_difference(new_gauge, old_gauge) < 2e-11
        @test maximum_momentum_difference(new_momenta, old_md.p) < 2e-11
        @test new_result.initial_hamiltonian ≈ old_initial_hamiltonian rtol=2e-11
        @test new_result.final_hamiltonian ≈ old_final_hamiltonian rtol=2e-11
        @test new_result.delta_hamiltonian ≈ old_delta rtol=2e-10 atol=2e-11
        @test metropolis_rule(new_result.delta_hamiltonian, 0.5) ==
              metropolis_rule(old_delta, 0.5)
    end
end

@testset "Legacy and typed SLHMC proposal equivalence" begin
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        communicator=Gaugefields.SerialCommunicator(),
        element_type=ComplexF64,
        verbose=0,
    )
    smearings = (
        NoFermionSmearingConfig(),
        StoutFermionSmearingConfig(1, (0.1,), ("plaquette",)),
    )

    for smearing in smearings, sexton_weingarten in (false, true)
        target_input = dynamical_test_input(
            WilsonDiracConfig(0.05, 1.0);
            sexton_weingarten,
            smearing,
        )
        hmc = target_input.update
        md_integrator = if sexton_weingarten
            SextonWeingartenConfig(
                QPQConfig(),
                ForceGroupConfig(:fermion_1),
                LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge)),
                1,
            )
        else
            hmc.md.integrator
        end
        md = MDConfig(hmc.md.step_size, hmc.md.steps, md_integrator)
        target_fermion = only(target_input.fermions)
        md_fermions = (
            FermionActionConfig(
                :fermion_1,
                WilsonDiracConfig(0.04, 1.0),
                target_fermion.flavors,
                target_fermion.solver,
                target_fermion.smearing,
            ),
        )
        md_gauge_config = GaugeActionConfig(
            GaugeActionTermConfig(:gauge_plaquette, "plaquette", 1.2),
        )
        slhmc_update = SLHMCConfig(
            md,
            hmc.momentum,
            hmc.acceptance,
            hmc.pseudofermions,
            md_gauge_config,
            md_fermions,
        )
        slhmc_input = LQCDConfig(
            target_input.lattice,
            target_input.gauge,
            target_input.gauge_action,
            target_input.fermions,
            slhmc_update,
        )

        new_simulation = build_simulation(slhmc_input, environment)
        old_configuration = build_configuration(target_input, environment)
        new_gauge = new_simulation.configuration.gauge
        old_gauge = old_configuration.gauge
        new_field = only(values(new_simulation.configuration.fermions))
        old_field = only(values(old_configuration.fermions))

        target_provider = new_simulation.action.terms.fermion_1
        noise = similar(new_field)
        LatticeDiracOperators.refresh_pseudofermion!(
            target_provider,
            new_gauge,
            noise;
            seed=0x13579bdf,
            sweep=4,
            subgroup=1,
        )

        old_md_fermion_action = legacy_fermion_action(
            only(md_fermions),
            old_gauge,
            old_field,
        )
        old_smearing = legacy_fermion_smearing(smearing, old_gauge)
        old_md_gauge_action = build_gauge_action(
            md_gauge_config,
            old_configuration,
        )
        old_md = LEGACY_MD.MD(
            old_gauge,
            old_md_gauge_action,
            false,
            md.step_size,
            md.steps,
            old_md_fermion_action,
            old_smearing;
            QPQ=true,
            SextonWeingargten=sexton_weingarten,
            Nsw=2,
        )
        LatticeDiracOperators.substitute_fermion!(old_md.η, new_field)

        momenta = Gaugefields.gaussian_momenta(
            new_gauge;
            seed=UInt64(0x2468ace0),
            sweep=3,
        )
        copy_momenta!(old_md.p, momenta)

        old_target_gauge_action = build_gauge_action(
            target_input.gauge_action,
            old_configuration,
        )
        old_target_fermion_action = legacy_fermion_action(
            target_fermion,
            old_gauge,
            old_md.η,
        )
        old_target_provider = LatticeDiracOperators.PseudofermionMDAction(
            old_target_fermion_action,
            old_md.η,
            old_smearing,
        )
        old_target_actions = Gaugefields.MDActionSet(;
            gauge=old_target_gauge_action,
            fermion_1=old_target_provider,
        )
        integrator = sexton_weingarten ? Gaugefields.SextonWeingarten(;
            slow=:fermion_1,
            fast=:gauge,
            n_fast=1,
        ) : Gaugefields.QPQ()
        trajectory_length = md.step_size * md.steps
        old_target_driver = Gaugefields.md_driver(
            old_gauge,
            old_target_actions;
            steps=md.steps,
            trajectory_length,
            integrator,
        )

        old_md_provider = LatticeDiracOperators.PseudofermionMDAction(
            old_md_fermion_action,
            old_md.η,
            old_smearing,
        )
        old_md_actions = Gaugefields.MDActionSet(;
            gauge=old_md_gauge_action,
            fermion_1=old_md_provider,
        )
        old_md_driver = Gaugefields.md_driver(
            old_gauge,
            old_md_actions;
            steps=md.steps,
            trajectory_length,
            integrator,
        )

        new_target_initial = Gaugefields.md_hamiltonian(
            new_gauge,
            momenta,
            new_simulation.updater.target_driver,
        )
        old_target_initial = Gaugefields.md_hamiltonian(
            old_gauge,
            old_md.p,
            old_target_driver,
        )
        old_md_initial = Gaugefields.md_hamiltonian(
            old_gauge,
            old_md.p,
            old_md_driver,
        )

        new_md_result = Gaugefields.md_trajectory!(
            new_gauge,
            momenta,
            new_simulation.updater.md_driver,
        )
        LEGACY_MD.runMD!(old_gauge, old_md)

        new_target_final = Gaugefields.md_hamiltonian(
            new_gauge,
            momenta,
            new_simulation.updater.target_driver,
        )
        old_target_final = Gaugefields.md_hamiltonian(
            old_gauge,
            old_md.p,
            old_target_driver,
        )
        old_md_final = Gaugefields.md_hamiltonian(
            old_gauge,
            old_md.p,
            old_md_driver,
        )
        new_target_delta = new_target_final - new_target_initial
        old_target_delta = old_target_final - old_target_initial

        @test maximum_link_difference(new_gauge, old_gauge) < 2e-11
        @test maximum_momentum_difference(momenta, old_md.p) < 2e-11
        @test new_target_initial ≈ old_target_initial rtol=2e-11
        @test new_target_final ≈ old_target_final rtol=2e-11
        @test new_target_delta ≈ old_target_delta rtol=2e-10 atol=2e-11
        @test new_md_result.initial_hamiltonian ≈ old_md_initial rtol=2e-11
        @test new_md_result.final_hamiltonian ≈ old_md_final rtol=2e-11
        @test metropolis_rule(new_target_delta, 0.5) ==
              metropolis_rule(old_target_delta, 0.5)
    end
end
