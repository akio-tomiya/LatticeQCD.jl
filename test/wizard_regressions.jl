using LatticeQCD, Test, TOML
import Logging

const RegressionWizard = LatticeQCD.Wizard

function regression_draft()
    w = RegressionWizard
    draft = w.WizardV2Draft()
    draft.mode = w.expert
    draft.physicalparams.L = [4, 4, 4, 4]
    draft.physicalparams.update_method = "Heatbath"
    draft.fermionparams.quench = true
    draft.fermionparams.Dirac_operator = nothing
    draft.fermion_parameters = w.Quench_parameters()
    return draft
end

function regression_pion(kind; stout=false)
    w = RegressionWizard
    method = w.Pion_parameters()
    answers = kind == :Wilson ? Any[5, 1, 0.164, 1e-12, 1700] :
                               Any[5, 2, 0.037, 4, 1e-12, 1700]
    append!(answers, stout ? Any[2, Set([1]), 0.12] : Any[1])
    ui = w.ScriptedWizardV2UI(answers)
    @test w.configure_wizard_v2_measurement!(ui, method, Int(w.Pion_correlator), [4,4,4,4]; isnew=true) === nothing
    @test isempty(ui.answers)
    return method
end

@testset "Wizard measurement parameters survive save and reload" begin
    w = RegressionWizard
    for kind in (:Wilson, :Staggered), flowed in (false, true), stout in (false, true)
        draft = regression_draft()
        method = regression_pion(kind; stout)
        set = flowed ? draft.measurement_gradientflow : draft.measurement
        set.measurement_methods = w.Measurement_parameters[method]
        draft.gradient_params.hasgradientflow = flowed
        dict = w.wizard_v2_parameter_dictionary(draft)
        document = TOML.parse(sprint(io -> TOML.print(io, dict)))
        spec = simulation_spec_from_legacy_toml(document)
        @test spec.config.gauge.halo == 1
        plan = flowed ? spec.schedule.measurements.gradient_flow.measurements : spec.schedule.measurements.direct
        observable = only(plan.measurements).observable
        if kind == :Wilson
            @test observable.fermion.hopping == 0.164
        else
            @test observable.fermion.mass == 0.037
            @test observable.fermion.flavors == 8
        end
        @test observable.solver.tolerance == 1e-12
        @test observable.solver.maximum_iterations == 1700
        @test only(plan.measurements).schedule.every == 5
        if stout
            @test observable.smearing.coefficients == (0.12,)
        end
        canonical = simulation_spec_dictionary(spec)
        @test simulation_spec_dictionary(simulation_spec_from_toml(canonical)) == canonical
        # Old saved input without a halo must also be made runnable.
        delete!(document["Physical setting"], "Nwing")
        upgraded = Logging.with_logger(Logging.NullLogger()) do
            simulation_spec_from_legacy_toml(document)
        end
        @test upgraded.config.gauge.halo == 1
    end

    # Revisiting the form and accepting defaults must preserve prior values.
    for kind in (:Wilson, :Staggered)
        method = regression_pion(kind)
        count = kind == :Wilson ? 6 : 7
        ui = w.ScriptedWizardV2UI(fill(:default, count))
        w.configure_wizard_v2_measurement!(ui, method, Int(w.Pion_correlator), [4,4,4,4])
        @test isempty(ui.answers)
        @test method.measure_every == 5
        @test method.eps == 1e-12
        @test method.MaxCGstep == 1700
        if kind == :Wilson
            @test method.fermion_parameters.hop == 0.164
        else
            @test method.fermion_parameters.mass == 0.037
            @test method.fermion_parameters.Nf == 8
        end
    end

    draft = regression_draft()
    chiral = w.ChiralCondensate_parameters()
    ui = w.ScriptedWizardV2UI(Any[7, 0.023, 1e-11, 1300, 1])
    w.configure_wizard_v2_measurement!(ui, chiral, Int(w.Chiral_condensate), [4,4,4,4])
    draft.measurement.measurement_methods = w.Measurement_parameters[chiral]
    spec = simulation_spec_from_legacy_toml(w.wizard_v2_parameter_dictionary(draft))
    obs = only(spec.schedule.measurements.direct.measurements).observable
    @test spec.config.gauge.halo == 1
    @test obs.fermion.mass == 0.023
    @test obs.solver.tolerance == 1e-11
    @test obs.solver.maximum_iterations == 1300
end

@testset "Wizard trajectory counts and legacy thresholds" begin
    w = RegressionWizard
    for first in (1, 11), warmup in (0, 1, 50)
        draft = regression_draft()
        draft.physicalparams.initialtrj = first
        ui = w.ScriptedWizardV2UI(Any[1, 1, 3, warmup, 300])
        w.edit_wizard_v2_update!(ui, draft)
        dict = w.wizard_v2_parameter_dictionary(draft)
        spec = simulation_spec_from_legacy_toml(dict)
        @test spec.schedule.thermalization_steps == warmup
        @test spec.schedule.production_steps == 300
        @test spec.schedule.initial_trajectory == first - 1
        @test dict["Physical setting"]["Nsteps"] == first - 1 + warmup + 300
        # An existing legacy input must keep its historical thresholds.
        dict["Physical setting"]["initialtrj"] = 1
        dict["Physical setting"]["Nthermalization"] = 50
        dict["Physical setting"]["Nsteps"] = 300
        legacy = simulation_spec_from_legacy_toml(dict)
        @test legacy.schedule.thermalization_steps == 49
        @test legacy.schedule.production_steps == 251
    end
end

@testset "Nested measurement input never silently drops explicit values" begin
    normalize = LatticeQCD.SimulationInput_module.normalize_legacy_measurement
    input = Dict{String,Any}("methodname" => "Pion_correlator", "hop" => 0.164,
        "fermion_parameters" => Dict{String,Any}("Dirac_operator" => "Wilson"))
    @test normalize(input)["hop"] == 0.164
    input["fermion_parameters"]["hop"] = 0.12
    @test_throws ArgumentError normalize(input)
    delete!(input["fermion_parameters"], "hop")
    input["fermion_parameters"]["hopp"] = 0.164
    @test_throws ArgumentError normalize(input)
    delete!(input["fermion_parameters"], "hopp")
    input["fermiontype"] = "Staggered"
    @test_throws ArgumentError normalize(input)
end

@testset "Typed text output, quiet execution, and resume" begin
    w = RegressionWizard
    mktempdir() do directory
        draft = regression_draft()
        draft.filename = joinpath(directory, "input.toml")
        draft.physicalparams.Nthermalization = 0
        draft.physicalparams.Nsteps = 2
        pion = regression_pion(:Wilson)
        pion.measure_every = 1
        draft.measurement.measurement_methods = w.Measurement_parameters[w.Plaq_parameters(), pion]
        draft.gradient_params.hasgradientflow = true
        draft.gradient_params.numflow = 1
        draft.measurement_gradientflow.measurement_methods =
            w.Measurement_parameters[deepcopy(pion)]
        draft.controlparams.log_dir = joinpath(directory, "nested", "logs")
        draft.controlparams.logfile = "run.log"
        draft.controlparams.measurement_basedir = joinpath(directory, "measurements")
        draft.controlparams.measurement_dir = "run"
        spec = w.wizard_v2_write_and_build_spec(draft)
        @test !isdir(draft.controlparams.log_dir)
        @test !isdir(draft.controlparams.measurement_basedir)
        canonical = simulation_spec_dictionary(spec)
        @test simulation_spec_dictionary(simulation_spec_from_toml(canonical)) == canonical
        sink = RecordingSimulationEventSink()
        session = build_simulation(spec; sink)
        step!(session)
        summary = run!(session; verbose=false)
        @test summary.production_completed == 2
        @test isfile(spec.output.text.logfile)
        @test occursin("run finished", read(spec.output.text.logfile, String))
        records = [e.record for e in sink.events if e isa MeasurementFinished]
        for (name, flowed) in ((:plaquette, false), (:pion_correlator, false), (:pion_correlator, true))
            path = flowed ? joinpath(spec.output.text.measurement_directory, "flow", "$name.txt") :
                            joinpath(spec.output.text.measurement_directory, "$name.txt")
            data = filter(l -> !startswith(l, "#"), readlines(path))
            expected = filter(r -> r.name == name && (r.point.flow_time !== nothing) == flowed, records)
            @test length(data) == length(expected)
            for (line, record) in zip(data, expected)
                columns = split(line, '\t')
                @test parse(Int, columns[1]) == record.point.trajectory
                @test columns[3] == sprint(io -> show(IOContext(io, :compact => true, :limit => false), record.value))
            end
        end
        # Saving an array must never inherit the console's truncation limit.
        long_values = collect(1.0:128.0)
        long_record = MeasurementRecord(:long_vector, long_values, MeasurementPoint(2))
        LatticeQCD.SimulationSession_module.emit_event!(session, MeasurementFinished(long_record))
        line = last(readlines(joinpath(spec.output.text.measurement_directory, "long_vector.txt")))
        @test split(line, '\t')[3] == sprint(io -> show(IOContext(io, :compact => true, :limit => false), long_values))
        @test !occursin("…", line)
        # Explicit typed input should fail validation before allocating fields.
        bad = simulation_spec_dictionary(spec)
        bad["gauge"]["halo"] = 0
        @test !isempty(validate(simulation_spec_from_toml(bad)))
    end
end

@testset "Nondefault dynamical fermion settings survive both TOML formats" begin
    w = RegressionWizard
    cases = (
        (Any[2, 0.137, 1e-11, 1400, 1], Dict("hopping_parameter" => 0.137)),
        (Any[3, 0.139, 1.27, 1e-11, 1400, 1], Dict("hopping_parameter" => 0.139, "clover_coefficient" => 1.27)),
        (Any[4, 0.039, 4, 1e-11, 1400, 1], Dict("mass" => 0.039)),
        (Any[5, 0.041, -0.07, 4, 1e-11, 1400], Dict("mass" => 0.041, "naik_epsilon" => -0.07)),
        (Any[6, 6, -1.3, 0.17, 1e-11, 1400, 1], Dict("fifth_dimension" => 6, "domainwall_height" => -1.3, "mass" => 0.17)),
        (Any[7, 6, -1.3, 0.17, 2.3, 0.7, 1e-11, 1400, 1], Dict("fifth_dimension" => 6, "domainwall_height" => -1.3, "mass" => 0.17, "b" => 2.3, "c" => 0.7)),
    )
    for (answers, expected) in cases
        draft = regression_draft()
        draft.physicalparams.update_method = "HMC"
        ui = w.ScriptedWizardV2UI(answers)
        w.edit_wizard_v2_fermion!(ui, draft)
        @test isempty(ui.answers)
        legacy = TOML.parse(sprint(io -> TOML.print(io, w.wizard_v2_parameter_dictionary(draft))))
        spec = simulation_spec_from_legacy_toml(legacy)
        canonical = simulation_spec_dictionary(spec)
        fermion = only(canonical["fermions"])
        for (key, value) in expected
            @test fermion["operator"][key] == value
        end
        @test fermion["solver"]["tolerance"] == 1e-11
        @test fermion["solver"]["max_steps"] == 1400
        if answers[1] in (4, 5)
            @test fermion["flavors"] == 8
        end
        @test simulation_spec_dictionary(simulation_spec_from_toml(canonical)) == canonical
    end
end
