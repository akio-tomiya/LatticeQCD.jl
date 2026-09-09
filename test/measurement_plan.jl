using LatticeQCD
using Test
import Gaugefields
import QCDMeasurements
import TOML
import Random

function measurement_test_params(filename, directory)
    parameters = TOML.parsefile(joinpath(@__DIR__, filename))
    control = parameters["System Control"]
    relative_directory = relpath(directory, pwd())
    control["log_dir"] = joinpath(relative_directory, "logs")
    control["logfile"] = "measurement-plan.log"
    control["measurement_basedir"] = joinpath(
        relative_directory,
        "measurements",
    )
    control["measurement_dir"] = "measurement-plan"
    return LatticeQCD.Parameters_TOML.construct_Params_from_TOML(parameters)
end

function measurement_test_hmc_input(; initial="cold")
    lattice = LatticeConfig((2, 2, 2, 2))
    gauge = GaugeConfig(2, 1, initial, nothing, 0x1234)
    action = GaugeActionConfig(
        GaugeActionTermConfig(:gauge_plaquette, "plaquette", 1.9),
    )
    integrator = LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge))
    md = MDConfig(0.02, 2, integrator)
    momentum = GaussianMomentumConfig(
        1.0,
        RandomStreamConfig(0x5678, :momentum),
    )
    acceptance = RankZeroMetropolisConfig(
        RandomStreamConfig(0x9abc, :metropolis),
    )
    return LQCDConfig(
        lattice,
        gauge,
        action,
        HMCConfig(md, momentum, acceptance),
    )
end

function wizard_measurement_dictionaries()
    return Dict[
        Dict(
            "methodname" => "Plaquette",
            "measure_every" => 1,
            "printvalues" => true,
        ),
        Dict(
            "methodname" => "Polyakov_loop",
            "measure_every" => 2,
            "printvalues" => true,
        ),
        Dict(
            "methodname" => "Topological_charge",
            "measure_every" => 3,
            "kinds_of_topological_charge" => ["plaquette", "clover"],
            "improved_topological_charge_definition" => "alexandrou",
            "printvalues" => true,
        ),
        Dict(
            "methodname" => "Chiral_condensate",
            "measure_every" => 4,
            "fermiontype" => "Staggered",
            "mass" => 1.0,
            "Nf" => 4,
            "Nr" => 1,
            "eps" => 1.0e-10,
            "MaxCGstep" => 1_000,
            "smearing_for_fermion" => "nothing",
            "printvalues" => true,
        ),
        Dict(
            "methodname" => "Pion_correlator",
            "measure_every" => 5,
            "fermiontype" => "Staggered",
            "mass" => 1.0,
            "Nf" => 4,
            "eps" => 1.0e-10,
            "MaxCGstep" => 1_000,
            "method_CG" => "bicg",
            "smearing_for_fermion" => "nothing",
            "printvalues" => true,
        ),
        Dict(
            "methodname" => "Wilson_loop",
            "measure_every" => 6,
            "Tmax" => 1,
            "Rmax" => 1,
            "printvalues" => true,
        ),
        Dict(
            "methodname" => "Energy_density",
            "measure_every" => 7,
            "printvalues" => true,
        ),
    ]
end

@testset "Typed Wizard measurement plan" begin
    plan = measurement_plan(wizard_measurement_dictionaries(); start=3)
    @test length(plan.measurements) == 7

    observables = map(measurement -> measurement.observable, plan.measurements)
    @test observables[1] isa PlaquetteObservableConfig
    @test observables[2] isa PolyakovLoopObservableConfig
    @test observables[3] isa TopologicalChargeObservableConfig
    @test observables[4] isa ChiralCondensateObservableConfig
    @test observables[5] isa PionCorrelatorObservableConfig
    @test observables[6] isa WilsonLoopObservableConfig
    @test observables[7] isa EnergyDensityObservableConfig

    @test observables[3].kinds == (:plaquette, :clover)
    @test observables[3].improved_definition == :alexandrou
    @test observables[4].fermion isa StaggeredMeasurementFermionConfig
    @test observables[4].noise_vectors == 1
    @test observables[5].fermion isa StaggeredMeasurementFermionConfig
    @test observables[5].solver.method == :bicg
    @test observables[6].Tmax == 1
    @test observables[6].Rmax == 1

    @test map(observable_name, observables) == (
        :plaquette,
        :polyakov_loop,
        :topological_charge,
        :chiral_condensate,
        :pion_correlator,
        :wilson_loop,
        :energy_density,
    )
    @test map(measurement -> measurement.schedule.every, plan.measurements) ==
          (1, 2, 3, 4, 5, 6, 7)
    @test all(
        measurement -> measurement.schedule.start == 3,
        plan.measurements,
    )

    @test !is_due(PeriodicSchedule(2, 3), 2)
    @test !is_due(PeriodicSchedule(2, 3), 3)
    @test is_due(PeriodicSchedule(2, 3), 4)
    @test_throws ArgumentError PeriodicSchedule(0)
    @test_throws ArgumentError TopologicalChargeObservableConfig((:unknown,))
    @test WilsonLoopObservableConfig(0, 1).Tmax == 0
    @test_throws ArgumentError WilsonLoopObservableConfig(-1, 1)

    stout_chiral = copy(wizard_measurement_dictionaries()[4])
    stout_chiral["smearing_for_fermion"] = "stout"
    stout_chiral["stout_numlayers"] = 1
    stout_chiral["stout_ρ"] = [0.1]
    stout_chiral["stout_loops"] = ["plaquette"]
    stout_config = measurement_config(stout_chiral)
    @test stout_config.observable.smearing isa
          StoutMeasurementSmearingConfig
    @test stout_config.observable.smearing.coefficients == (0.1,)
    @test stout_config.observable.smearing.loops == ("plaquette",)

    wilson_pion = copy(wizard_measurement_dictionaries()[5])
    wilson_pion["fermiontype"] = "Wilson"
    delete!(wilson_pion, "mass")
    delete!(wilson_pion, "Nf")
    wilson_pion["hop"] = 0.141139
    wilson_pion["r"] = 1.0
    wilson_config = measurement_config(wilson_pion)
    @test wilson_config.observable.fermion isa
          WilsonMeasurementFermionConfig
    @test !wilson_config.observable.fermion.clover

    for component in (
        plan,
        plan.measurements...,
        observables...,
        (measurement.schedule for measurement in plan.measurements)...,
        observables[4].fermion,
        observables[4].solver,
        observables[4].smearing,
        observables[5].fermion,
        observables[5].solver,
        observables[5].smearing,
    )
        @test all(isconcretetype, fieldtypes(typeof(component)))
    end

    output = sprint(show_config, plan)
    @test occursin("MeasurementPlan", output)
    @test occursin("plaquette every=1 start=3", output)
    @test occursin("pion_correlator fermion=Staggered", output)
    @test occursin("wilson_loop Tmax=1 Rmax=1", output)
end

@testset "Wizard measurements through Params" begin
    mktempdir() do directory
        parameters = measurement_test_params("test02.toml", directory)
        try
            empty!(parameters.measurement_methods)
            append!(
                parameters.measurement_methods,
                wizard_measurement_dictionaries(),
            )
            plan = measurement_plan(parameters)
            @test length(plan.measurements) == 7
            @test all(
                measurement ->
                    measurement.schedule.start == parameters.Nthermalization,
                plan.measurements,
            )

            program = measurement_program(parameters)
            @test program.direct isa MeasurementPlan
            @test program.gradient_flow isa GradientFlowMeasurementConfig
            @test program.gradient_flow.step_size == parameters.eps_flow
            @test program.gradient_flow.samples == parameters.numflow
            @test program.gradient_flow.integration_steps == parameters.Nflow
            @test length(program.gradient_flow.measurements.measurements) == 1
            flow_measurement = only(
                program.gradient_flow.measurements.measurements,
            )
            @test flow_measurement.observable isa
                  TopologicalChargeObservableConfig
            @test flow_measurement.schedule.every == 10
            @test flow_measurement.schedule.start == 0
            @test all(isconcretetype, fieldtypes(typeof(program)))
            @test all(
                isconcretetype,
                fieldtypes(typeof(program.gradient_flow)),
            )

            output = sprint(show_config, program)
            @test occursin("MeasurementProgram", output)
            @test occursin("GradientFlowMeasurements", output)
            @test occursin("every_flow_sample=10", output)

            empty!(parameters.measurement_methods)
            push!(
                parameters.measurement_methods,
                wizard_measurement_dictionaries()[1],
            )
            empty!(parameters.measurements_for_flow)
            push!(
                parameters.measurements_for_flow,
                wizard_measurement_dictionaries()[1],
            )
            simulation = build_simulation(
                measurement_test_hmc_input(initial="cold"),
                GaugefieldsEnvironment(
                    backend=Gaugefields.LegacyBackend(),
                    verbose=0,
                ),
            )
            params_runner = build_simulation_runner(simulation, parameters)
            @test params_runner.measurements isa MeasurementProgramRuntime
            @test params_runner.measurements.gradient_flow isa
                  GradientFlowMeasurementRuntime
        finally
            isopen(parameters.load_fp) && close(parameters.load_fp)
        end
    end

    empty_plan = measurement_plan(Dict[])
    @test isempty(empty_plan.measurements)
    @test occursin("measurements: none", sprint(show_config, empty_plan))

    no_flow = MeasurementProgram(empty_plan)
    @test no_flow.gradient_flow isa NoGradientFlowMeasurementConfig
    @test occursin("disabled", sprint(show_config, no_flow))

    one_measurement = measurement_plan([
        Dict("methodname" => "Plaquette", "measure_every" => 1),
    ])
    @test_throws ArgumentError GradientFlowMeasurementConfig(
        0.0,
        1,
        1,
        one_measurement,
    )
    @test_throws ArgumentError GradientFlowMeasurementConfig(
        0.01,
        0,
        1,
        one_measurement,
    )
    @test_throws ArgumentError GradientFlowMeasurementConfig(
        0.01,
        1,
        0,
        one_measurement,
    )
    @test_throws ArgumentError GradientFlowMeasurementConfig(
        0.01,
        1,
        1,
        empty_plan,
    )
end

function gauge_observable_plan()
    values = wizard_measurement_dictionaries()
    gauge_values = (values[1], values[2], values[3], values[6], values[7])
    return measurement_plan(collect(gauge_values))
end

@testset "Gauge-observable measurement runtime" begin
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        element_type=ComplexF64,
        verbose=0,
    )
    simulation = build_simulation(
        measurement_test_hmc_input(initial="cold"),
        environment,
    )
    plan = gauge_observable_plan()
    runtime = build_measurement_runtime(plan, simulation)

    @test length(runtime.measurements) == 5
    @test runtime.measurements[1].measurement isa
          QCDMeasurements.PlaquetteMeasurement
    @test runtime.measurements[2].measurement isa
          QCDMeasurements.PolyakovMeasurement
    @test runtime.measurements[3].measurement isa
          QCDMeasurements.TopologicalChargeMeasurement
    @test runtime.measurements[4].measurement isa
          QCDMeasurements.WilsonLoopMeasurement
    @test runtime.measurements[5].measurement isa
          QCDMeasurements.EnergyDensityMeasurement
    @test all(isconcretetype, fieldtypes(typeof(runtime)))
    @test all(
        measurement -> all(
            isconcretetype,
            fieldtypes(typeof(measurement)),
        ),
        runtime.measurements,
    )

    records = measure_now!(runtime, simulation)
    @test map(record -> record.name, records) == (
        :plaquette,
        :polyakov_loop,
        :topological_charge,
        :wilson_loop,
        :energy_density,
    )
    @test records[1].value ≈ Gaugefields.measure_plaquette(
        simulation.configuration.gauge,
    )
    @test isfinite(records[1].value)
    @test isfinite(real(records[2].value))
    @test all(isfinite, values(records[3].value))
    @test size(records[4].value) == (1, 1)
    @test all(isfinite, records[4].value)
    @test isfinite(records[5].value)
    @test all(record -> record.point.trajectory == 0, records)
    @test all(record -> record.point.flow_time === nothing, records)

    due = measure_due!(runtime, simulation; trajectory=1)
    @test length(due) == 1
    @test only(due).name == :plaquette
    @test only(due).point.trajectory == 1

    seen = Symbol[]
    runner = build_simulation_runner(
        simulation,
        measurement_plan([wizard_measurement_dictionaries()[1]]);
        sink=FunctionMeasurementSink(record -> push!(seen, record.name)),
    )
    result = step!(runner)
    @test result.update isa HMCUpdateResult
    @test length(result.measurements) == 1
    @test only(result.measurements).name == :plaquette
    @test only(result.measurements).point.trajectory == 1
    @test seen == [:plaquette]
    @test all(isconcretetype, fieldtypes(typeof(runner)))

    run!(runner; steps=1)
    @test runner.simulation.state.trajectory == 2
    @test seen == [:plaquette, :plaquette]

    flow_values = wizard_measurement_dictionaries()
    flowed_topology = copy(flow_values[3])
    flowed_topology["measure_every"] = 1
    gpu_flow_program = MeasurementProgram(
        MeasurementPlan(),
        GradientFlowMeasurementConfig(
            0.01,
            1,
            1,
            measurement_plan([flow_values[1], flowed_topology]),
        ),
    )
    gpu_flow_runtime = build_measurement_runtime(
        gpu_flow_program,
        simulation,
    )
    source_plaquette = Gaugefields.measure_plaquette(
        simulation.configuration.gauge,
    )
    gpu_flow_records = measure_due!(
        gpu_flow_runtime,
        simulation;
        trajectory=2,
    )
    @test map(record -> record.name, gpu_flow_records) ==
          [:plaquette, :topological_charge]
    @test all(record -> record.point.flow_time == 0.01, gpu_flow_records)
    @test isfinite(gpu_flow_records[1].value)
    @test all(isfinite, values(gpu_flow_records[2].value))
    @test Gaugefields.measure_plaquette(simulation.configuration.gauge) ==
          source_plaquette
end

@testset "Gradient-flow measurement program" begin
    base_input = measurement_test_hmc_input(initial="cold")
    input = LQCDConfig(
        base_input.lattice,
        GaugeConfig(2, 0, HotStartConfig()),
        base_input.gauge_action,
        base_input.update,
    )
    simulation = build_simulation(
        input,
        GaugefieldsEnvironment(
            backend=Gaugefields.LegacyBackend(),
            verbose=0,
        ),
    )
    direct = measurement_plan([
        Dict("methodname" => "Plaquette", "measure_every" => 1),
    ])
    flowed = measurement_plan([
        Dict("methodname" => "Plaquette", "measure_every" => 1),
        Dict("methodname" => "Energy_density", "measure_every" => 2),
    ])
    flow_config = GradientFlowMeasurementConfig(0.01, 2, 2, flowed)
    program = MeasurementProgram(direct, flow_config)
    seen = Tuple{Symbol,Any}[]
    runner = build_simulation_runner(
        simulation,
        program;
        sink=FunctionMeasurementSink(record -> push!(
            seen,
            (record.name, record.point.flow_time),
        )),
    )

    runtime = runner.measurements
    @test runtime isa MeasurementProgramRuntime
    @test runtime.gradient_flow isa GradientFlowMeasurementRuntime
    @test runtime.gradient_flow.flow.Nflow == 2
    @test runtime.gradient_flow.flow.eps == 0.01
    @test runtime.gradient_flow.configuration.gauge !==
          simulation.configuration.gauge
    @test all(isconcretetype, fieldtypes(typeof(runtime)))
    @test all(
        isconcretetype,
        fieldtypes(typeof(runtime.gradient_flow)),
    )

    original_plaquette = Gaugefields.measure_plaquette(
        simulation.configuration.gauge,
    )
    records = measure_due!(runner; trajectory=3)
    @test map(record -> record.name, records) == [
        :plaquette,
        :plaquette,
        :plaquette,
        :energy_density,
    ]
    @test map(record -> record.point.flow_time, records) ==
          Any[nothing, 0.02, 0.04, 0.04]
    @test all(record -> record.point.trajectory == 3, records)
    @test all(record -> begin
        value = record.value
        value isa Number ? isfinite(value) : true
    end, records)
    @test Gaugefields.measure_plaquette(simulation.configuration.gauge) ==
          original_plaquette
    @test seen == map(
        record -> (record.name, record.point.flow_time),
        records,
    )

    all_records = measure_now!(runner; trajectory=4)
    @test map(record -> record.name, all_records) == [
        :plaquette,
        :plaquette,
        :energy_density,
        :plaquette,
        :energy_density,
    ]
    @test map(record -> record.point.flow_time, all_records) ==
          Any[nothing, 0.02, 0.02, 0.04, 0.04]
    @test_throws ArgumentError measure_due!(
        runtime,
        simulation;
        trajectory=4,
        flow_time=0.1,
    )

    direct_runner = build_simulation_runner(
        simulation,
        MeasurementProgram(direct),
    )
    @test direct_runner.measurements.gradient_flow isa
          NoGradientFlowMeasurementRuntime
    direct_records = measure_due!(direct_runner; trajectory=5)
    @test direct_records isa Tuple
    @test length(direct_records) == 1
    @test only(direct_records).point.flow_time === nothing
end

@testset "QCDMeasurements owns fermionic measurement runtime" begin
    gpu_environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        element_type=ComplexF64,
        verbose=0,
    )
    gpu_simulation = build_simulation(
        measurement_test_hmc_input(initial="cold"),
        gpu_environment,
    )
    values = wizard_measurement_dictionaries()
    plan = measurement_plan([values[4], values[5]])
    runtime = build_measurement_runtime(plan, gpu_simulation)

    @test runtime.measurements[1].measurement isa
          QCDMeasurements.ChiralCondensateMeasurement
    @test runtime.measurements[2].measurement isa
          QCDMeasurements.PionCorrelatorMeasurement
    @test gpu_simulation.configuration isa GaugeConfiguration
    @test fieldnames(typeof(gpu_simulation.configuration)) == (:gauge,)

    stout_values = map(index -> copy(values[index]), (4, 5))
    for stout_value in stout_values
        stout_value["smearing_for_fermion"] = "stout"
        stout_value["stout_numlayers"] = 1
        stout_value["stout_ρ"] = [0.1]
        stout_value["stout_loops"] = ["plaquette"]
    end
    stout_runtime = build_measurement_runtime(
        measurement_plan(collect(stout_values)),
        gpu_simulation,
    )
    @test stout_runtime.measurements[1].measurement isa
          QCDMeasurements.ChiralCondensateMeasurement
    @test stout_runtime.measurements[2].measurement isa
          QCDMeasurements.PionCorrelatorMeasurement

    cpu_simulation = build_simulation(
        measurement_test_hmc_input(initial="cold"),
        GaugefieldsEnvironment(
            backend=Gaugefields.LegacyBackend(),
            verbose=0,
        ),
    )
    cpu_runtime = build_measurement_runtime(plan, cpu_simulation)
    records = measure_now!(cpu_runtime, cpu_simulation)
    @test map(record -> record.name, records) ==
          (:chiral_condensate, :pion_correlator)
    @test isfinite(records[1].value)
    @test all(isfinite, records[2].value)

    flowed_fermion_values = map(index -> copy(values[index]), (4, 5))
    foreach(value -> value["measure_every"] = 1, flowed_fermion_values)
    flowed_fermion_program = MeasurementProgram(
        MeasurementPlan(),
        GradientFlowMeasurementConfig(
            0.01,
            1,
            1,
            measurement_plan(collect(flowed_fermion_values)),
        ),
    )
    flowed_fermion_runtime = build_measurement_runtime(
        flowed_fermion_program,
        cpu_simulation,
    )
    flowed_fermion_records = measure_due!(
        flowed_fermion_runtime,
        cpu_simulation;
        trajectory=1,
    )
    @test map(record -> record.name, flowed_fermion_records) ==
          [:chiral_condensate, :pion_correlator]
    @test all(
        record -> record.point.flow_time == 0.01,
        flowed_fermion_records,
    )
    @test isfinite(flowed_fermion_records[1].value)
    @test all(isfinite, flowed_fermion_records[2].value)
end


function wizard_measurement_value_equal(left, right)
    if left isa Number && right isa Number
        return isapprox(left, right; rtol=1e-11, atol=1e-11)
    elseif left isa AbstractArray && right isa AbstractArray
        size(left) == size(right) || return false
        return all(
            wizard_measurement_value_equal(l, r)
            for (l, r) in zip(left, right)
        )
    elseif left isa AbstractDict && right isa AbstractDict
        Set(keys(left)) == Set(keys(right)) || return false
        return all(
            wizard_measurement_value_equal(left[key], right[key])
            for key in keys(left)
        )
    elseif left isa Tuple && right isa Tuple
        length(left) == length(right) || return false
        return all(
            wizard_measurement_value_equal(l, r)
            for (l, r) in zip(left, right)
        )
    end
    return isequal(left, right)
end


@testset "All Wizard measurements through SimulationSession" begin
    direct_values = deepcopy(wizard_measurement_dictionaries())
    flow_values = deepcopy(wizard_measurement_dictionaries())
    foreach(value -> value["measure_every"] = 1, direct_values)
    foreach(value -> value["measure_every"] = 1, flow_values)
    program = MeasurementProgram(
        measurement_plan(direct_values),
        GradientFlowMeasurementConfig(
            0.01,
            1,
            1,
            measurement_plan(flow_values),
        ),
    )
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        element_type=ComplexF64,
        verbose=0,
    )
    input = measurement_test_hmc_input(initial="cold")

    Random.seed!(0x123456)
    runner = build_simulation_runner(
        build_simulation(input, environment),
        program,
    )
    Random.seed!(0x123456)
    events = RecordingSimulationEventSink()
    session = build_simulation(
        SimulationSpec(input, SimulationSchedule(
            0,
            1;
            measurements=program,
        )),
        environment;
        sink=events,
    )

    Random.seed!(0xabcdef)
    reference = step!(runner)
    Random.seed!(0xabcdef)
    candidate = step!(session)

    @test candidate.update.accepted == reference.update.accepted
    @test candidate.update.delta_hamiltonian ≈
          reference.update.delta_hamiltonian rtol=1e-12 atol=1e-12
    @test length(reference.measurements) == 14
    @test length(candidate.measurements) == 14
    @test getproperty.(reference.measurements, :name) ==
          getproperty.(candidate.measurements, :name)
    @test getproperty.(reference.measurements, :name) == [
        :plaquette,
        :polyakov_loop,
        :topological_charge,
        :chiral_condensate,
        :pion_correlator,
        :wilson_loop,
        :energy_density,
        :plaquette,
        :polyakov_loop,
        :topological_charge,
        :chiral_condensate,
        :pion_correlator,
        :wilson_loop,
        :energy_density,
    ]
    @test all(
        record -> record.point.flow_time === nothing,
        reference.measurements[1:7],
    )
    @test all(
        record -> record.point.flow_time == 0.01,
        reference.measurements[8:14],
    )
    for (left, right) in zip(
        reference.measurements,
        candidate.measurements,
    )
        @test left.point.trajectory == right.point.trajectory == 1
        @test left.point.flow_time == right.point.flow_time
        @test wizard_measurement_value_equal(left.value, right.value)
    end
    @test count(
        event -> event isa MeasurementFinished,
        events.events,
    ) == 14
end


@testset "Wizard fermionic measurement variants through session" begin
    values = wizard_measurement_dictionaries()
    stout_chiral = deepcopy(values[4])
    stout_staggered_pion = deepcopy(values[5])
    for measurement in (stout_chiral, stout_staggered_pion)
        measurement["measure_every"] = 1
        measurement["smearing_for_fermion"] = "stout"
        measurement["stout_numlayers"] = 1
        measurement["stout_ρ"] = [0.1]
        measurement["stout_loops"] = ["plaquette"]
    end

    wilson_pion = deepcopy(values[5])
    wilson_pion["measure_every"] = 1
    wilson_pion["fermiontype"] = "Wilson"
    delete!(wilson_pion, "mass")
    delete!(wilson_pion, "Nf")
    wilson_pion["hop"] = 0.141139
    wilson_pion["r"] = 1.0
    stout_wilson_pion = deepcopy(wilson_pion)
    stout_wilson_pion["smearing_for_fermion"] = "stout"
    stout_wilson_pion["stout_numlayers"] = 1
    stout_wilson_pion["stout_ρ"] = [0.1]
    stout_wilson_pion["stout_loops"] = ["plaquette"]

    program = MeasurementProgram(measurement_plan([
        stout_chiral,
        wilson_pion,
        stout_wilson_pion,
        stout_staggered_pion,
    ]))
    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        element_type=ComplexF64,
        verbose=0,
    )
    input = measurement_test_hmc_input(initial="cold")

    Random.seed!(0x2468)
    runner = build_simulation_runner(
        build_simulation(input, environment),
        program,
    )
    Random.seed!(0x2468)
    session = build_simulation(
        SimulationSpec(input, SimulationSchedule(
            0,
            1;
            measurements=program,
        )),
        environment,
    )
    Random.seed!(0x1357)
    reference = step!(runner)
    Random.seed!(0x1357)
    candidate = step!(session)

    @test getproperty.(reference.measurements, :name) == (
        :chiral_condensate,
        :pion_correlator,
        :pion_correlator,
        :pion_correlator,
    )
    @test getproperty.(candidate.measurements, :name) ==
          getproperty.(reference.measurements, :name)
    for (left, right) in zip(
        reference.measurements,
        candidate.measurements,
    )
        @test wizard_measurement_value_equal(left.value, right.value)
    end
end
