using LatticeQCD, Test, TOML

const AuditWizard = LatticeQCD.Wizard

function audit_draft()
    d = AuditWizard.WizardV2Draft()
    d.mode = AuditWizard.expert
    d.physicalparams.update_method = "Heatbath"
    return d
end

# Exercise the same page-completion bookkeeping as run_wizard, then revisit
# every prompt with Enter. Compare persisted input, not mutable object identity.
function audit_page!(d, page, answers)
    w = AuditWizard
    ui = w.ScriptedWizardV2UI(answers)
    result = w.edit_wizard_v2_page!(ui, page, d)
    @test result.action == w.WizardV2Next
    @test isempty(ui.answers)
    push!(d.completed_pages, page)
    before = w.wizard_v2_parameter_dictionary(d)
    again = w.ScriptedWizardV2UI(fill(:default, length(ui.prompts)))
    @test w.edit_wizard_v2_page!(again, page, d).action == w.WizardV2Next
    @test isempty(again.answers)
    @test w.wizard_v2_parameter_dictionary(d) == before
    return d
end

# Use real files for both legacy Wizard output and canonical typed input.
function audit_reload(d)
    mktempdir() do dir
        d.filename = joinpath(dir, "wizard-custom.toml")
        spec = AuditWizard.wizard_v2_write_and_build_spec(d)
        @test isfile(d.filename)
        expected = simulation_spec_dictionary(spec)
        @test simulation_spec_dictionary(load_simulation_spec(d.filename)) == expected
        canonical_path = joinpath(dir, "canonical.toml")
        write_simulation_spec(canonical_path, spec)
        @test simulation_spec_dictionary(load_simulation_spec(canonical_path)) == expected
        return expected
    end
end

@testset "Wizard parameter audit" begin
    @testset "Wizard parameter audit: lattice and initialization" begin
        w = AuditWizard
        for (choice, nc) in ((1, 3), (2, 2), (3, 4))
            d = audit_draft()
            answers = Any[6, 8, 10, 12, choice]
            choice == 3 && push!(answers, nc)
            append!(answers, [793, 3, 6.12])
            audit_page!(d, w.WizardV2LatticePage, answers)
            audit_page!(d, w.WizardV2SourcePage, Any[1, 2])
            c = audit_reload(d)
            @test c["lattice"]["size"] == [6, 8, 10, 12]
            @test c["gauge"]["colors"] == nc
            @test c["gauge_action"]["terms"][1]["coupling"] == 6.12
            @test c["gauge"]["initialization"]["kind"] == "hot"
            @test c["gauge"]["initialization"]["seed"] == 793
            @test c["update"]["random"]["seed"] == 793
        end
        for (index, format) in enumerate(("JLD", "ILDG", "BridgeText"))
            d = audit_draft()
            audit_page!(d, w.WizardV2SourcePage, Any[1, 3, index, "custom.conf", 13])
            audit_page!(d, w.WizardV2UpdatePage, Any[1, 2, 7, 23])
            c = audit_reload(d)
            @test c["gauge"]["initialization"]["path"] == "custom.conf"
            @test c["gauge"]["initialization"]["format"] == format
            @test c["schedule"]["initial_trajectory"] == 12
            @test c["schedule"]["thermalization_steps"] == 7
            @test c["schedule"]["production_steps"] == 23
        end
        for (nc, choice, kind) in ((2, 1, "cold"), (2, 4, "instanton"),
                                   (2, 5, "embedded_instanton"), (3, 4, "embedded_instanton"))
            d = audit_draft()
            d.physicalparams.NC = nc
            audit_page!(d, w.WizardV2SourcePage, Any[1, choice])
            @test audit_reload(d)["gauge"]["initialization"]["kind"] == kind
        end
        d = audit_draft()
        @test_throws ErrorException w.edit_wizard_v2_source!(
            w.ScriptedWizardV2UI(Any[1, 3, 1, "custom.conf", 0]), d)
    end

    @testset "Wizard parameter audit: file sequences" begin
        w = AuditWizard
        for (index, format) in enumerate(("jld2", "ildg", "bridge")), manifest in (false, true)
            d = audit_draft()
            answers = Any[2, index, "./custom-input", manifest ? 2 : 1]
            manifest && push!(answers, "custom-list.txt")
            audit_page!(d, w.WizardV2SourcePage, answers)
            c = audit_reload(d)
            @test c["update"]["kind"] == "configuration_sequence"
            @test c["update"]["format"] == format
            @test c["update"]["source"]["directory"] == "./custom-input"
            @test c["update"]["source"]["kind"] == (manifest ? "manifest" : "directory")
            manifest && @test c["update"]["source"]["manifest"] == "custom-list.txt"
        end
    end

    @testset "Wizard parameter audit: fermions, solvers, smearing and MD" begin
        w = AuditWizard
        cases = (
            (Any[2, 0.137, 1e-11, 1400], Dict("kind"=>"wilson", "hopping_parameter"=>0.137)),
            (Any[3, 0.139, 1.27, 1e-11, 1400], Dict("kind"=>"wilson_clover", "hopping_parameter"=>0.139, "clover_coefficient"=>1.27)),
            (Any[4, 0.039, 4, 1e-11, 1400], Dict("kind"=>"staggered", "mass"=>0.039)),
            (Any[5, 0.041, -0.07, 4, 1e-11, 1400], Dict("kind"=>"hisq", "mass"=>0.041, "naik_epsilon"=>-0.07)),
            (Any[6, 6, -1.3, 0.17, 1e-11, 1400], Dict("kind"=>"domainwall", "fifth_dimension"=>6, "domainwall_height"=>-1.3, "mass"=>0.17)),
            (Any[7, 6, -1.3, 0.17, 2.3, 0.7, 1e-11, 1400], Dict("kind"=>"mobius_domainwall", "fifth_dimension"=>6, "domainwall_height"=>-1.3, "mass"=>0.17, "b"=>2.3, "c"=>0.7)),
        )
        for (input, expected) in cases, slhmc in (false, true)
            d = audit_draft()
            d.controlparams.verboselevel = 3
            d.controlparams.randomseed = 793
            answers = copy(input)
            input[1] != 5 && append!(answers, Any[2, Set([1, 2]), 0.12, 0.03])
            audit_page!(d, w.WizardV2FermionPage, answers)
            update_answers = slhmc ? Any[2, 4.91] : Any[1]
            append!(update_answers, Any[7, 23, 11, 0.031, 2, 6])
            audit_page!(d, w.WizardV2UpdatePage, update_answers)
            c = audit_reload(d)
            f = only(c["fermions"])
            for (key, value) in expected
                @test f["operator"][key] == value
            end
            @test f["flavors"] == (input[1] in (4, 5) ? 8 : 2)
            @test f["solver"]["tolerance"] == 1e-11
            @test f["solver"]["max_steps"] == 1400
            @test f["solver"]["verbose"] == 3
            if input[1] != 5
                @test f["smearing"]["coefficients"] == [0.12, 0.03]
                @test f["smearing"]["loops"] == collect(w.kindsof_loops)[[1, 2]]
            end
            # Check the values handed to LatticeDiracOperators, beyond TOML.
            action = only(simulation_spec_from_toml(c).config.fermions)
            dirac = LatticeQCD.Simulation_module.dirac_parameters(action)
            @test dirac["eps_CG"] == 1e-11
            @test dirac["MaxCGstep"] == 1400
            for (external, key) in (("κ", "hopping_parameter"), ("cSW", "clover_coefficient"),
                                    ("mass", "mass"), ("naik_epsilon", "naik_epsilon"),
                                    ("L5", "fifth_dimension"), ("M", "domainwall_height"),
                                    ("b", "b"), ("c", "c"))
                haskey(expected, key) && @test dirac[external] == expected[key]
            end
            u = c["update"]
            @test u["kind"] == (slhmc ? "SLHMC" : "HMC")
            @test u["md"]["steps"] == 11
            @test u["md"]["step_size"] == 0.031
            @test u["md"]["integrator"]["fast_steps"] == 3
            @test u["momentum"]["random"]["seed"] == 793
            @test u["acceptance"]["seed"] == 793
            @test only(u["pseudofermions"])["random"]["seed"] == 793
            if slhmc
                @test u["md_gauge_action"]["terms"][1]["coupling"] == 4.91
                @test u["md_fermions"] == c["fermions"]
            end
            @test c["schedule"]["thermalization_steps"] == 7
            @test c["schedule"]["production_steps"] == 23
        end
        # Every staggered/HISQ flavor choice, not just the default and Nf=8.
        for kind in (4, 5), (index, nf) in enumerate((2, 3, 4, 8, 1))
            d = audit_draft()
            answers = kind == 4 ? Any[4, 0.037, index, 1e-10, 901, 1] :
                                  Any[5, 0.037, -0.03, index, 1e-10, 901]
            audit_page!(d, w.WizardV2FermionPage, answers)
            d.physicalparams.update_method = "HMC"
            @test only(audit_reload(d)["fermions"])["flavors"] == nf
        end
        for overrelax in (false, true)
            d = audit_draft()
            answers = Any[1, overrelax ? 1 : 2]
            overrelax && push!(answers, 5)
            append!(answers, [7, 23])
            audit_page!(d, w.WizardV2UpdatePage, answers)
            @test audit_reload(d)["update"]["overrelaxation_steps"] == (overrelax ? 5 : 0)
        end
    end

    @testset "Wizard parameter audit: all direct and flowed measurements" begin
        w = AuditWizard
        for flowed in (false, true), pion_kind in (1, 2)
            d = audit_draft()
            # Distinct intervals detect accidental reuse across measurements.
            answers = Any[Set(1:7), 2, 3, 4,
                5, 0.023, 1e-11, 1301, 2, Set([1, 2]), 0.13, 0.04,
                6, pion_kind]
            append!(answers, pion_kind == 1 ? Any[0.164, 1e-12, 1701] :
                                              Any[0.037, 5, 1e-12, 1701])
            append!(answers, Any[2, Set([1, 2]), 0.14, 0.05, 7, 3, 5, 8])
            if flowed
                prepend!(answers, Any[2, 0.027, 9])
                audit_page!(d, w.WizardV2GradientFlowPage, answers)
            else
                audit_page!(d, w.WizardV2MeasurementPage, answers)
            end
            c = audit_reload(d)["schedule"]["measurements"]
            if flowed
                @test c["gradient_flow"]["step_size"] == 0.027
                @test c["gradient_flow"]["samples"] == 9
            end
            measurements = flowed ? c["gradient_flow"]["measurements"] : c["direct"]
            byname = Dict(m["methodname"] => m for m in measurements)
            @test length(byname) == 7
            for (index, name) in enumerate(w.WIZARD_V2_MEASUREMENT_NAMES)
                @test byname[name]["every"] == index + 1
            end
            chiral = byname["Chiral_condensate"]
            @test chiral["mass"] == 0.023
            @test chiral["eps"] == 1e-11
            @test chiral["MaxCGstep"] == 1301
            @test chiral["stout_ρ"] == [0.13, 0.04]
            pion = byname["Pion_correlator"]
            @test pion["eps"] == 1e-12
            @test pion["MaxCGstep"] == 1701
            @test pion["stout_ρ"] == [0.14, 0.05]
            if pion_kind == 1
                @test pion["hop"] == 0.164
            else
                @test pion["mass"] == 0.037
                @test pion["Nf"] == 1
            end
            @test byname["Wilson_loop"]["Rmax"] == 3
            @test byname["Wilson_loop"]["Tmax"] == 5
        end
    end

    @testset "Wizard parameter audit: outputs and flow-only paths" begin
        w = AuditWizard
        for direct in (false, true), (index, format) in enumerate(("none", "jld2", "ildg", "bridge"))
            d = audit_draft()
            d.physicalparams.update_method = "HMC"
            if direct
                d.measurement.measurement_methods = w.Measurement_parameters[w.Plaq_parameters()]
            else
                d.gradient_params.hasgradientflow = true
                d.measurement_gradientflow.measurement_methods = w.Measurement_parameters[w.Plaq_parameters()]
            end
            answers = Any["./custom-measurements", "custom-run", "./custom-logs", "custom.log", index]
            index > 1 && append!(answers, Any[13, "./custom-confs"])
            append!(answers, Any[2, 17, "./custom-restart"])
            audit_page!(d, w.WizardV2OutputPage, answers)
            c = audit_reload(d)["output"]
            @test c["text"]["measurement_directory"] == abspath("custom-measurements", "custom-run")
            @test c["text"]["logfile"] == abspath("custom-logs", "custom.log")
            @test c["configurations"]["format"] == format
            if index > 1
                @test c["configurations"]["every"] == 13
                @test c["configurations"]["directory"] == "./custom-confs"
            end
            @test c["checkpoints"]["every"] == 17
            @test c["checkpoints"]["directory"] == "./custom-restart"
            # Deselect all measurements, revisit output, then re-enable flow.
            d.measurement.measurement_methods = w.Measurement_parameters[]
            d.gradient_params.hasgradientflow = false
            audit_page!(d, w.WizardV2OutputPage, answers[3:end])
            @test d.controlparams.measurement_basedir == "./custom-measurements"
            @test d.controlparams.measurement_dir == "custom-run"
        end
    end

    @testset "Wizard parameter audit: manual TOML fields and legacy reader" begin
        w = AuditWizard
        d = audit_draft()
        audit_page!(d, w.WizardV2FermionPage, Any[2, 0.137, 1e-11, 1400, 1])
        d.physicalparams.update_method = "HMC"
        d.physicalparams.Nwing = 2
        d.fermion_parameters.r = 0.91
        d.fermionparams.BoundaryCondition = [-1, 1, 1, 1]
        d.hmcparams.QPQ = false
        d.fermionparams.smearing_for_fermion = "stout"
        d.fermionparams.stout_numlayers = 3
        d.fermionparams.stout_loops = ["plaquette"]
        d.fermionparams.stout_ρ = [0.025]
        c = audit_reload(d)
        f = only(c["fermions"])
        @test c["gauge"]["halo"] == 2
        @test f["operator"]["wilson_parameter"] == 0.91
        @test f["solver"]["boundary_conditions"] == [-1, 1, 1, 1]
        @test f["smearing"]["layers"] == 3
        @test c["update"]["md"]["integrator"]["ordering"] == "PQP"
        d = audit_draft()
        d.controlparams.isevenodd = false
        @test !audit_reload(d)["update"]["even_odd"]

        d = audit_draft()
        chiral = w.ChiralCondensate_parameters()
        chiral.Nr = 7
        pion = w.Pion_parameters()
        pion.method_CG = "bicgstab"
        topology = w.TopologicalCharge_parameters()
        topology.kinds_of_topological_charge = ["clover"]
        topology.improved_topological_charge_definition = "both"
        d.measurement.measurement_methods = w.Measurement_parameters[chiral, pion, topology]
        byname = Dict(m["methodname"]=>m for m in audit_reload(d)["schedule"]["measurements"]["direct"])
        @test byname["Chiral_condensate"]["Nr"] == 7
        @test byname["Pion_correlator"]["method_CG"] == "bicgstab"
        @test byname["Topological_charge"]["kinds_of_topological_charge"] == ["clover"]
        @test byname["Topological_charge"]["improved_topological_charge_definition"] == "both"

        readers = (LatticeQCD.SimulationInput_module.normalize_legacy_measurement,
            values -> only(LatticeQCD.Parameters_TOML.construct_measurement_dir(Dict("pion"=>values))))
        for reader in readers
            input = Dict{String,Any}("methodname"=>"Pion_correlator", "hop"=>0.164,
                "fermion_parameters"=>Dict{String,Any}("Dirac_operator"=>"Wilson"))
            @test reader(input)["hop"] == 0.164
            input["fermion_parameters"]["hop"] = 0.12
            @test_throws ArgumentError reader(input)
            delete!(input["fermion_parameters"], "hop")
            input["fermion_parameters"]["hopp"] = 0.164
            @test_throws ArgumentError reader(input)
        end
    end

    @testset "Wizard parameter audit: flow-only runtime and output cadence" begin
        w = AuditWizard
        mktempdir() do dir
            d = audit_draft()
            d.physicalparams.L = [2, 2, 2, 2]
            d.physicalparams.NC = 2
            audit_page!(d, w.WizardV2UpdatePage, Any[2, 1, 3, 2, 0.003, 1])
            audit_page!(d, w.WizardV2GradientFlowPage, Any[2, 0.017, 3, Set([1]), 2])
            d.gradient_params.Nflow = 3
            audit_page!(d, w.WizardV2OutputPage,
                Any[joinpath(dir, "meas"), "flow-only", joinpath(dir, "logs"), "audit.log",
                    2, 2, joinpath(dir, "confs"), 2, 3, joinpath(dir, "restart")])
            d.filename = joinpath(dir, "audit.toml")
            spec = w.wizard_v2_write_and_build_spec(d)
            sink = RecordingSimulationEventSink()
            session = build_simulation(spec; sink)
            @test isempty(session.measurements.direct.measurements)
            @test session.measurements.gradient_flow.config.integration_steps == 3
            @test session.simulation.input.update.md.step_size == 0.003
            @test session.simulation.input.update.md.steps == 2
            summary = run!(session; verbose=false)
            @test summary.thermalization_completed == 1
            @test summary.production_completed == 3
            @test summary.saved_configurations == 2
            @test summary.saved_checkpoints == 1
            records = [event.record for event in sink.events if event isa MeasurementFinished]
            @test [r.point.trajectory for r in records] == [2, 3, 4]
            @test all(r -> r.point.flow_time ≈ 0.102, records)
            output = joinpath(dir, "meas", "flow-only", "flow", "plaquette.txt")
            @test isfile(output)
            lines = filter(line -> !startswith(line, "#"), readlines(output))
            @test length(lines) == 3
            @test parse.(Int, first.(split.(lines, '\t'))) == [2, 3, 4]
            @test isfile(joinpath(dir, "logs", "audit.log"))
            @test length(readdir(joinpath(dir, "confs"))) == 2
            @test length(readdir(joinpath(dir, "restart"))) == 1
        end
    end
    @testset "Wizard parameter audit: simple mode and invalid inputs" begin
        w = AuditWizard
        d = w.WizardV2Draft()
        audit_page!(d, w.WizardV2ModePage, Any[1, "simple-custom.toml"])
        @test d.filename == "simple-custom.toml"
        audit_page!(d, w.WizardV2LatticePage, Any[6, 12, 6.12])
        audit_page!(d, w.WizardV2SourcePage, Any[1, 2])
        audit_page!(d, w.WizardV2FermionPage, Any[0.137, 1])
        audit_page!(d, w.WizardV2UpdatePage, Any[23])
        audit_page!(d, w.WizardV2MeasurementPage, Any[2, 3, 5, 1, 0.164, 1e-12, 1601, 1])
        audit_page!(d, w.WizardV2OutputPage, Any["./simple-meas", "run", "./simple-logs", "run.log", 1])
        c = audit_reload(d)
        @test c["lattice"]["size"] == [6, 6, 6, 12]
        @test c["gauge_action"]["terms"][1]["coupling"] == 6.12
        @test c["gauge"]["initialization"]["kind"] == "hot"
        @test only(c["fermions"])["operator"]["hopping_parameter"] == 0.137
        @test c["schedule"]["production_steps"] == 23
        @test c["schedule"]["thermalization_steps"] == 0
        pion = only(filter(m -> m["methodname"] == "Pion_correlator", c["schedule"]["measurements"]["direct"]))
        @test pion["hop"] == 0.164
        @test pion["eps"] == 1e-12
        @test pion["MaxCGstep"] == 1601
        @test c["output"]["text"]["measurement_directory"] == abspath("simple-meas", "run")
        @test c["output"]["text"]["logfile"] == abspath("simple-logs", "run.log")
        @test c["output"]["checkpoints"]["format"] == "none"

        d = audit_draft()
        @test_throws ErrorException w.edit_wizard_v2_lattice!(
            w.ScriptedWizardV2UI(Any[4, 4, 4, 4, 1, -1]), d)
        audit_page!(d, w.WizardV2FermionPage, Any[2, 0.137, 1e-11, 1400, 1])
        @test_throws ErrorException w.edit_wizard_v2_update!(
            w.ScriptedWizardV2UI(Any[1, 7, 23, 11, 0.031, 2, 3]), d)
    end
end
