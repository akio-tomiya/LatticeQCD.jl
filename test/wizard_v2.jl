using TOML

const WizardV2 = LatticeQCD.Wizard

@testset "Wizard v2 navigation" begin
    @test LatticeQCD.run_wizard === WizardV2.run_wizard
    @test LatticeQCD.run_wizardv2 === LatticeQCD.run_wizard
    @test WizardV2.run_wizardv2 === WizardV2.run_wizard
    @test LatticeQCD.run_wizard_legacy === WizardV2.run_wizard_legacy
    @test occursin(
        "two-flavor Wilson-fermion HMC",
        WizardV2.wizard_v2_simple_mode_description(),
    )
    @test occursin(
        "plaquette, Polyakov loop, and pion correlator",
        WizardV2.wizard_v2_simple_mode_description(),
    )
    @test occursin("previous section", WizardV2.WIZARD_V2_BACK_OPTION)
    @test occursin("discard edits", WizardV2.WIZARD_V2_BACK_OPTION)
    @test WizardV2.wizard_v2_navigation_command("back") isa
          WizardV2.WizardV2Back
    @test WizardV2.wizard_v2_navigation_command(" :BACK ") isa
          WizardV2.WizardV2Back
    @test WizardV2.wizard_v2_navigation_command("quit") isa
          WizardV2.WizardV2Quit
    @test WizardV2.wizard_v2_navigation_command("configuration.jld2") ===
          nothing
    @test occursin(
        "Lattice and gauge settings' to 'Mode and filename",
        WizardV2.wizard_v2_back_message(
            WizardV2.WizardV2LatticePage,
            WizardV2.WizardV2ModePage,
        ),
    )
    @test occursin(
        "Restarting 'Mode and filename'",
        WizardV2.wizard_v2_back_message(
            WizardV2.WizardV2ModePage,
            WizardV2.WizardV2ModePage,
        ),
    )
    v2_welcome_io = IOBuffer()
    WizardV2.print_wizard_logo(v2_welcome_io; navigation=:v2)
    v2_welcome = String(take!(v2_welcome_io))
    @test occursin("Back to previous section", v2_welcome)
    @test occursin("type `back`", v2_welcome)
    @test occursin("Quit wizard without saving", v2_welcome)
    @test occursin("type `quit`", v2_welcome)

    legacy_welcome_io = IOBuffer()
    WizardV2.print_wizard_logo(legacy_welcome_io)
    legacy_welcome = String(take!(legacy_welcome_io))
    @test occursin("To exit, press Ctrl + c.", legacy_welcome)
    @test !occursin("Back to previous section", legacy_welcome)
    @test all(
        isconcretetype,
        fieldtypes(typeof(WizardV2.WizardV2Draft())),
    )

    @testset "Back discards unfinished page edits" begin
        draft = WizardV2.WizardV2Draft()
        ui = WizardV2.ScriptedWizardV2UI(Any[8, :back])

        result = WizardV2.edit_wizard_v2_lattice!(ui, draft)

        @test result.action == WizardV2.WizardV2BackAction
        @test draft.physicalparams.L == [4, 4, 4, 4]
        @test isempty(ui.answers)
    end

    @testset "File-loading route skips update pages" begin
        draft = WizardV2.WizardV2Draft()
        draft.isfileloading = true

        route = WizardV2.wizard_v2_route(draft)

        @test WizardV2.WizardV2FermionPage ∉ route
        @test WizardV2.WizardV2UpdatePage ∉ route
        @test route[end] == WizardV2.WizardV2ReviewPage
    end

    @testset "File-loading output never offers checkpoint saving" begin
        draft = WizardV2.WizardV2Draft()
        draft.mode = WizardV2.expert
        draft.isfileloading = true
        draft.physicalparams.update_method = "Fileloading"
        draft.controlparams.saveU_format = "ILDG"
        draft.controlparams.saveU_dir = "./stale-output"
        draft.controlparams.saveU_every = 17
        draft.controlparams.checkpoint_dir = "./stale-restart"
        draft.controlparams.checkpoint_every = 19
        ui = WizardV2.ScriptedWizardV2UI(Any[
            "./logs",
            "fileloading.log",
        ])

        result = WizardV2.edit_wizard_v2_output!(ui, draft)

        @test result.action == WizardV2.WizardV2Next
        @test isempty(ui.answers)
        @test !any(
            occursin("Configuration format for saving"),
            ui.prompts,
        )
        @test draft.controlparams.saveU_format === nothing
        @test isempty(draft.controlparams.saveU_dir)
        @test draft.controlparams.saveU_every == 1
        @test isempty(draft.controlparams.checkpoint_dir)
        @test draft.controlparams.checkpoint_every == 0
        @test !any(
            prompt -> occursin("restart checkpoint", lowercase(prompt)),
            ui.prompts,
        )
    end

    @testset "Expert HMC configures an independent checkpoint interval" begin
        draft = WizardV2.WizardV2Draft()
        draft.mode = WizardV2.expert
        draft.physicalparams.update_method = "HMC"
        ui = WizardV2.ScriptedWizardV2UI(Any[
            "./logs",
            "hmc.log",
            1,              # no ordinary configuration output
            2,              # portable JLD2 restart checkpoints
            17,
            "./restart",
        ])

        result = WizardV2.edit_wizard_v2_output!(ui, draft)

        @test result.action == WizardV2.WizardV2Next
        @test isempty(ui.answers)
        @test draft.controlparams.saveU_format == "nothing"
        @test draft.controlparams.checkpoint_every == 17
        @test draft.controlparams.checkpoint_dir == "./restart"
        document = WizardV2.wizard_v2_parameter_dictionary(draft)
        spec = simulation_spec_from_legacy_toml(document)
        @test spec.output.configurations isa NoConfigurationOutput
        @test spec.output.checkpoints isa JLD2CheckpointOutput
        @test spec.output.checkpoints.every == 17
    end

    @testset "Embedded instanton is an initialization choice" begin
        draft = WizardV2.WizardV2Draft()
        draft.physicalparams.NC = 3
        ui = WizardV2.ScriptedWizardV2UI(Any[
            1,  # not Fileloading mode
            4,  # embedded instanton
        ])

        result = WizardV2.edit_wizard_v2_source!(ui, draft)

        @test result.action == WizardV2.WizardV2Next
        @test isempty(ui.answers)
        @test draft.physicalparams.initial == "embedded instanton"
        @test draft.controlparams.loadU_format === nothing
    end

    @testset "Expert mode accepts general SU(N)" begin
        draft = WizardV2.WizardV2Draft()
        draft.mode = WizardV2.expert
        ui = WizardV2.ScriptedWizardV2UI(Any[
            4, 4, 4, 4,
            3,      # Other SU(N)
            4,      # SU(4)
            111,
            2,
            7.2,
        ])

        result = WizardV2.edit_wizard_v2_lattice!(ui, draft)

        @test result.action == WizardV2.WizardV2Next
        @test isempty(ui.answers)
        @test draft.physicalparams.NC == 4
        @test draft.physicalparams.β == 7.2
    end

    @testset "Expert fermion and stout branches" begin
        clover = WizardV2.WizardV2Draft()
        clover.mode = WizardV2.expert
        clover_ui = WizardV2.ScriptedWizardV2UI(Any[
            3,          # Wilson--clover
            0.13,       # kappa
            1.2,        # cSW
            1e-10,
            2_000,
            1,          # no smearing
        ])

        result = WizardV2.edit_wizard_v2_fermion!(clover_ui, clover)

        @test result.action == WizardV2.WizardV2Next
        @test isempty(clover_ui.answers)
        @test clover.fermionparams.Dirac_operator == "WilsonClover"
        @test clover.fermion_parameters.Dirac_operator == "WilsonClover"
        @test clover.fermion_parameters.hasclover
        @test clover.fermion_parameters.hop == 0.13
        @test clover.fermion_parameters.Clover_coefficient == 1.2
        clover_dictionary = WizardV2.wizard_v2_parameter_dictionary(clover)
        @test clover_dictionary["Physical setting(fermions)"][
            "Clover_coefficient"
        ] == 1.2

        hisq = WizardV2.WizardV2Draft()
        hisq.mode = WizardV2.expert
        hisq.physicalparams.NC = 3
        hisq_ui = WizardV2.ScriptedWizardV2UI(Any[
            5,          # HISQ
            0.2,        # mass
            -0.083,     # Naik correction
            3,          # Nf=4
            1e-10,
            2_000,
        ])

        result = WizardV2.edit_wizard_v2_fermion!(hisq_ui, hisq)

        @test result.action == WizardV2.WizardV2Next
        @test isempty(hisq_ui.answers)
        @test hisq.fermionparams.Dirac_operator == "HISQ"
        @test hisq.fermion_parameters isa WizardV2.HISQ_parameters
        @test hisq.fermion_parameters.mass == 0.2
        @test hisq.fermion_parameters.naik_epsilon == -0.083
        @test hisq.fermion_parameters.Nf == 4
        @test hisq.fermionparams.smearing_for_fermion == "nothing"
        @test !any(occursin("stout"), hisq_ui.prompts)
        hisq_dictionary = WizardV2.wizard_v2_parameter_dictionary(hisq)
        @test hisq_dictionary["Physical setting"]["Nwing"] == 3
        @test hisq_dictionary["Physical setting(fermions)"][
            "naik_epsilon"
        ] == -0.083

        for colors in (2, 4)
            generic_hisq = WizardV2.WizardV2Draft()
            generic_hisq.mode = WizardV2.expert
            generic_hisq.physicalparams.NC = colors
            generic_hisq_ui = WizardV2.ScriptedWizardV2UI(Any[
                5,          # HISQ for any SU(N)
                0.2,
                -0.083,
                3,          # Nf=4
                1e-10,
                2_000,
            ])

            result = WizardV2.edit_wizard_v2_fermion!(
                generic_hisq_ui,
                generic_hisq,
            )

            @test result.action == WizardV2.WizardV2Next
            @test isempty(generic_hisq_ui.answers)
            @test generic_hisq.fermionparams.Dirac_operator == "HISQ"
            @test WizardV2.wizard_v2_parameter_dictionary(generic_hisq)[
                "Physical setting"
            ]["NC"] == colors
        end

        staggered = WizardV2.WizardV2Draft()
        staggered.mode = WizardV2.expert
        staggered_ui = WizardV2.ScriptedWizardV2UI(Any[
            4,          # Staggered
            0.25,       # mass
            2,          # Nf=3
            1e-10,
            2_000,
            2,          # stout
            Set(1:7),   # every Wizard stout loop
            0.01,
            0.02,
            0.03,
            0.04,
            0.05,
            0.06,
            0.07,
        ])

        result = WizardV2.edit_wizard_v2_fermion!(staggered_ui, staggered)

        @test result.action == WizardV2.WizardV2Next
        @test isempty(staggered_ui.answers)
        @test staggered.fermionparams.Dirac_operator == "Staggered"
        @test staggered.fermion_parameters.Nf == 3
        @test staggered.fermionparams.smearing_for_fermion == "stout"
        @test staggered.fermionparams.stout_loops ==
              LatticeQCD.Parameter_structs.kindsof_loops
        @test staggered.fermionparams.stout_ρ == collect(0.01:0.01:0.07)

        domainwall = WizardV2.WizardV2Draft()
        domainwall.mode = WizardV2.expert
        domainwall_ui = WizardV2.ScriptedWizardV2UI(Any[
            6,      # Domain wall (HISQ is choice 5 for SU(3))
            6,
            -1.2,
            0.15,
            1e-10,
            2_000,
            1,      # no smearing
        ])

        result = WizardV2.edit_wizard_v2_fermion!(domainwall_ui, domainwall)

        @test result.action == WizardV2.WizardV2Next
        @test isempty(domainwall_ui.answers)
        @test domainwall.fermionparams.Dirac_operator == "Domainwall"
        @test domainwall.fermion_parameters.N5 == 6
        @test domainwall.fermion_parameters.M == -1.2
        @test domainwall.fermion_parameters.m == 0.15

        mobius = WizardV2.WizardV2Draft()
        mobius.mode = WizardV2.expert
        mobius_ui = WizardV2.ScriptedWizardV2UI(Any[
            7,      # Möbius domain wall
            6,
            -1.2,
            0.15,
            2.2,
            0.8,
            1e-10,
            2_000,
            1,      # no smearing
        ])

        result = WizardV2.edit_wizard_v2_fermion!(mobius_ui, mobius)

        @test result.action == WizardV2.WizardV2Next
        @test isempty(mobius_ui.answers)
        @test mobius.fermionparams.Dirac_operator == "MobiusDomainwall"
        @test mobius.fermion_parameters isa
              WizardV2.MobiusDomainwall_parameters
        @test mobius.fermion_parameters.N5 == 6
        @test mobius.fermion_parameters.M == -1.2
        @test mobius.fermion_parameters.m == 0.15
        @test mobius.fermion_parameters.b == 2.2
        @test mobius.fermion_parameters.c == 0.8
        mobius_dictionary = WizardV2.wizard_v2_parameter_dictionary(mobius)
        @test mobius_dictionary["Physical setting(fermions)"]["b"] == 2.2
        @test mobius_dictionary["Physical setting(fermions)"]["c"] == 0.8
    end

    @testset "Expert Sexton-Weingarten and gradient-flow branches" begin
        draft = WizardV2.WizardV2Draft()
        draft.mode = WizardV2.expert
        draft.fermionparams.quench = false
        draft.fermionparams.Dirac_operator = "Wilson"
        update_ui = WizardV2.ScriptedWizardV2UI(Any[
            1,      # HMC
            0,      # thermalization
            2,      # trajectories
            1,      # MD steps
            0.001,
            2,      # Sexton-Weingarten
            4,
        ])

        result = WizardV2.edit_wizard_v2_update!(update_ui, draft)

        @test result.action == WizardV2.WizardV2Next
        @test isempty(update_ui.answers)
        @test draft.physicalparams.update_method == "HMC"
        @test draft.hmcparams.SextonWeingargten
        @test draft.hmcparams.N_SextonWeingargten == 4

        flow_ui = WizardV2.ScriptedWizardV2UI(Any[
            2,          # enable flow
            0.01,
            1,
            Set(1:7),   # all flow measurements
            1,          # plaquette interval
            1,          # Polyakov interval
            1,          # topology interval
            1,          # chiral interval
            0.5,
            1e-10,
            1_000,
            1,          # no chiral smearing
            1,          # pion interval
            1,          # Wilson pion
            0.141139,
            1e-10,
            1_000,
            1,          # no pion smearing
            1,          # Wilson-loop interval
            1,
            1,
            1,          # energy-density interval
        ])

        result = WizardV2.edit_wizard_v2_gradient_flow!(flow_ui, draft)

        @test result.action == WizardV2.WizardV2Next
        @test isempty(flow_ui.answers)
        @test draft.gradient_params.hasgradientflow
        @test getproperty.(
            draft.measurement_gradientflow.measurement_methods,
            :methodname,
        ) == WizardV2.WIZARD_V2_MEASUREMENT_NAMES
    end

    @testset "Expert SLHMC writes its MD-action beta" begin
        draft = WizardV2.WizardV2Draft()
        draft.mode = WizardV2.expert
        draft.fermionparams.quench = false
        draft.fermionparams.Dirac_operator = "Wilson"
        ui = WizardV2.ScriptedWizardV2UI(Any[
            2,      # SLHMC
            4.8,    # effective beta
            0,      # thermalization
            2,      # trajectories
            1,      # MD steps
            0.001,
            1,      # no Sexton-Weingarten
        ])

        result = WizardV2.edit_wizard_v2_update!(ui, draft)

        @test result.action == WizardV2.WizardV2Next
        @test isempty(ui.answers)
        @test draft.physicalparams.update_method == "SLHMC"
        @test draft.slhmc_beta == 4.8
        dictionary = WizardV2.wizard_v2_parameter_dictionary(draft)
        @test dictionary["SLHMC related"]["βeff"] == 4.8
    end
end

@testset "Wizard v2 legacy-compatible output" begin
    mktempdir() do directory
        cd(directory) do
            filename = joinpath(directory, "parameters-without-forced-extension")
            header = "HMC_L04040404_beta5.7_Wilson_kappa0.141139"
            ui = WizardV2.ScriptedWizardV2UI(Any[
                1,                              # simple mode
                filename,
                4, 4, 5.7,                     # lattice and beta
                1, 1,                          # no directory loading, cold start
                0.141139, 1,                   # Wilson fermion, no smearing
                :back,                         # return from update to fermions
                0.141139, 1,                   # keep the fermion settings
                101,                           # trajectories
                1, 1,                          # plaquette and Polyakov intervals
                1, 1, 0.141139, 1e-19, 3000, 1, # pion measurement
                "./measurements", header,
                "./logs", "$header.txt",
                2, 17, "./restart",            # restart checkpoints
                1,                              # save at review
            ])

            spec = WizardV2.run_wizard(ui)

            @test spec isa LatticeQCD.SimulationSpec
            @test isempty(ui.answers)
            @test isfile(filename)
            @test !isfile("$filename.toml")
            @test !isdir("./logs")
            @test !isdir("./measurements")

            physical = WizardV2.Print_Physical_parameters(Nsteps=101)
            fermions = WizardV2.Print_Fermions_parameters(
                quench=false,
                Dirac_operator="Wilson",
            )
            fermion_parameters = WizardV2.Wilson_parameters(hop=0.141139)
            control = WizardV2.Print_System_control_parameters(
                log_dir="./logs",
                logfile="$header.txt",
                measurement_basedir="./measurements",
                measurement_dir=header,
                checkpoint_dir="./restart",
                checkpoint_every=17,
            )
            hmc = WizardV2.Print_HMCrelated_parameters()
            plaquette = WizardV2.Plaq_parameters(measure_every=1)
            polyakov = WizardV2.Poly_parameters(measure_every=1)
            pion = WizardV2.Pion_parameters(measure_every=1)
            measurement = WizardV2.Measurement_parameterset(
                measurement_methods=WizardV2.Measurement_parameters[
                    plaquette,
                    polyakov,
                    pion,
                ],
            )
            gradient = WizardV2.Print_Gradientflow_parameters()
            flow_measurement = WizardV2.Measurement_parameterset()
            expected = WizardV2.wizard_parameter_dictionary(
                physical,
                fermions,
                fermion_parameters,
                control,
                hmc,
                measurement,
                gradient,
                flow_measurement,
            )

            @test TOML.parsefile(filename) == expected
            @test spec.schedule.production_steps == 101
            @test spec.config.fermions[1].operator isa WilsonDiracConfig
            @test spec.config.gauge.halo == 1
            @test spec.output.checkpoints isa JLD2CheckpointOutput
            @test spec.output.checkpoints.every == 17
            @test length(spec.schedule.measurements.direct.measurements) == 3
            # Two fermion-page visits plus the pion-measurement prompt.
            @test count(
                ==("Hopping parameter kappa for the two-flavor Wilson fermion"),
                ui.prompts,
            ) == 2
            @test count(==("Hopping parameter kappa"), ui.prompts) == 1
            summary = WizardV2.wizard_v2_review_summary(
                let draft = WizardV2.WizardV2Draft()
                    draft.fermionparams.quench = false
                    draft.fermionparams.Dirac_operator = "Wilson"
                    draft.fermion_parameters =
                        WizardV2.Wilson_parameters(hop=0.141139)
                    draft.controlparams.checkpoint_dir = "./restart"
                    draft.controlparams.checkpoint_every = 17
                    draft
                end,
            )
            @test occursin("two-flavor Wilson", summary)
            @test occursin("portable JLD2 every 17", summary)
        end
    end
end

@testset "Wizard v2 expert heatbath output" begin
    mktempdir() do directory
        cd(directory) do
            filename = joinpath(directory, "expert.toml")
            header = "Heatbath_L04040404_beta5.7_quenched"
            ui = WizardV2.ScriptedWizardV2UI(Any[
                2, filename,                     # expert mode
                4, 4, 4, 4, 1, 111, 2, 5.7,    # lattice and gauge settings
                1, 1,                            # no directory loading, cold start
                1,                               # quenched approximation
                1, 1, 3, 0, 101,                # Heatbath and trajectory settings
                Int[],                           # no measurements
                1,                               # no gradient flow
                "./logs", "$header.txt",
                1,                               # do not save configurations
                1,                               # save at review
            ])

            # The compatibility name is the same typed, Param-free Wizard.
            spec = WizardV2.run_wizardv2(ui)
            generated = TOML.parsefile(filename)

            @test spec isa LatticeQCD.SimulationSpec
            @test spec.config.update isa HeatbathConfig
            @test spec.config.update.overrelaxation_steps == 3
            @test isempty(spec.schedule.measurements.direct.measurements)
            @test isempty(ui.answers)
            @test generated["Physical setting"]["update_method"] == "Heatbath"
            @test generated["Physical setting"]["useOR"]
            @test generated["Physical setting"]["numOR"] == 3
            @test generated["System Control"]["measurement_basedir"] ==
                  "./measurements"
            @test isempty(generated["Measurement set"]["measurement_methods"])
            @test !haskey(generated["HMC related"], "MDsteps")
            @test !haskey(generated["HMC related"], "Δτ")
            @test !isdir("./logs")
            @test !isdir("./measurements")
        end
    end
end

@testset "Wizard v2 file-loading output" begin
    mktempdir() do directory
        cd(directory) do
            filename = joinpath(directory, "file-loading.toml")
            header = "Fileloading_L04040404_beta5.7_nothing"
            ui = WizardV2.ScriptedWizardV2UI(Any[
                1, filename,                    # simple mode
                4, 4, 5.7,                     # lattice and beta
                2, 1, "./confs", 1,            # load every JLD file in a directory
                1, 1,                          # plaquette and Polyakov intervals
                1, 1, 0.141139, 1e-19, 3000, 1, # pion measurement
                "./measurements", header,
                "./logs", "$header.txt",
                1,                              # save at review
            ])

            spec = WizardV2.run_wizard(ui)
            generated = TOML.parsefile(filename)

            @test spec isa LatticeQCD.SimulationSpec
            @test isempty(ui.answers)
            @test spec.config.update isa ConfigurationSequenceConfig
            @test isempty(spec.config.fermions)
            @test generated["System Control"]["loadU_format"] == "JLD"
            @test generated["System Control"]["loadU_dir"] == "./confs"
            @test count(==("Hopping parameter kappa"), ui.prompts) == 1
            @test !any(
                prompt -> occursin("restart checkpoint", lowercase(prompt)),
                ui.prompts,
            )
            @test !isdir("./logs")
            @test !isdir("./measurements")
        end
    end
end
