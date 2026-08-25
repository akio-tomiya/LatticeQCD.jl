using LatticeQCD
using Test
using TOML
import Gaugefields

const WIZARD_FERMION_DIRECTORY = joinpath(@__DIR__, "wizard_fermion")
const WIZARD_FERMION_VARIANTS = (
    (id="wilson", operator="Wilson", flavors=2, type=WilsonDiracConfig),
    (id="wilson_clover", operator="WilsonClover", flavors=2, type=WilsonCloverDiracConfig),
    (id="staggered_nf1", operator="Staggered", flavors=1, type=StaggeredDiracConfig),
    (id="staggered_nf2", operator="Staggered", flavors=2, type=StaggeredDiracConfig),
    (id="staggered_nf3", operator="Staggered", flavors=3, type=StaggeredDiracConfig),
    (id="staggered_nf4", operator="Staggered", flavors=4, type=StaggeredDiracConfig),
    (id="staggered_nf8", operator="Staggered", flavors=8, type=StaggeredDiracConfig),
    (id="hisq_nf1", operator="HISQ", flavors=1, type=HISQDiracConfig),
    (id="hisq_nf2", operator="HISQ", flavors=2, type=HISQDiracConfig),
    (id="hisq_nf3", operator="HISQ", flavors=3, type=HISQDiracConfig),
    (id="hisq_nf4", operator="HISQ", flavors=4, type=HISQDiracConfig),
    (id="hisq_nf8", operator="HISQ", flavors=8, type=HISQDiracConfig),
    (id="domainwall", operator="Domainwall", flavors=2, type=DomainwallDiracConfig),
    (id="mobius_domainwall", operator="MobiusDomainwall", flavors=2, type=MobiusDomainwallDiracConfig),
)

function expected_wizard_fermion_cases()
    cases = NamedTuple[]
    for variant in WIZARD_FERMION_VARIANTS
        colors_values = variant.operator == "HISQ" ? (3,) : (2, 3)
        stout_values = variant.operator == "HISQ" ? (false,) : (false, true)
        for colors in colors_values
            for stout in stout_values
                for sexton in (false, true)
                    for update_method in ("HMC", "SLHMC")
                        parts = String["su$colors"]
                        update_method == "SLHMC" && push!(parts, "slhmc")
                        append!(parts, [
                            variant.id,
                            stout ? "stout" : "thin",
                            sexton ? "sw" : "leapfrog",
                        ])
                        filename = join(parts, "_") * ".toml"
                        push!(
                            cases,
                            (;
                                filename,
                                colors,
                                variant,
                                stout,
                                sexton,
                                update_method,
                            ),
                        )
                    end
                end
            end
        end
    end
    return cases
end

function wizard_fermion_params(path, directory)
    parameters = TOML.parsefile(path)
    control = parameters["System Control"]
    relative_directory = relpath(directory, pwd())
    control["log_dir"] = joinpath(relative_directory, "logs")
    control["logfile"] = basename(path) * ".log"
    control["measurement_basedir"] = joinpath(
        relative_directory,
        "measurements",
    )
    control["measurement_dir"] = splitext(basename(path))[1]
    return LatticeQCD.Parameters_TOML.construct_Params_from_TOML(parameters)
end

function suppress_wizard_fermion_output(function_to_run)
    quiet_display = Base.Multimedia.TextDisplay(devnull)
    Base.Multimedia.pushdisplay(quiet_display)
    try
        return redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                function_to_run()
            end
        end
    finally
        Base.Multimedia.popdisplay(quiet_display)
    end
end


function wizard_fermion_bridge_bytes(configuration)
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


function test_wizard_fermion_update_equal(reference, candidate)
    @test typeof(candidate) === typeof(reference)
    @test candidate.trajectory == reference.trajectory
    @test candidate.accepted == reference.accepted
    @test isapprox(
        candidate.initial_hamiltonian,
        reference.initial_hamiltonian;
        rtol=1e-12,
        atol=1e-12,
    )
    @test isapprox(
        candidate.final_hamiltonian,
        reference.final_hamiltonian;
        rtol=1e-12,
        atol=1e-12,
    )
    @test isapprox(
        candidate.delta_hamiltonian,
        reference.delta_hamiltonian;
        rtol=1e-12,
        atol=1e-12,
    )
end

@testset "Wizard dynamical-fermion parameter matrix" begin
    cases = expected_wizard_fermion_cases()
    existing = sort(filter(
        filename -> endswith(filename, ".toml"),
        readdir(WIZARD_FERMION_DIRECTORY),
    ))
    @test length(cases) == 164
    @test existing == sort(getproperty.(cases, :filename))

    # The complete HISQ trajectories are substantially more expensive than
    # the one-link cases.  Keep the default test exhaustive, while allowing
    # developers and CI jobs to split the runtime matrix without changing
    # which files and filenames are validated above.
    operator_filter = get(ENV, "LQCD_WIZARD_FERMION_FILTER", "")
    if !isempty(operator_filter)
        cases = filter(
            case -> lowercase(case.variant.operator) ==
                    lowercase(operator_filter),
            cases,
        )
    end
    shard_count = parse(Int, get(
        ENV,
        "LQCD_WIZARD_FERMION_SHARD_COUNT",
        "1",
    ))
    shard_index = parse(Int, get(
        ENV,
        "LQCD_WIZARD_FERMION_SHARD_INDEX",
        "1",
    ))
    shard_count > 0 || error("the Wizard fermion shard count must be positive")
    1 <= shard_index <= shard_count || error(
        "the Wizard fermion shard index must be between 1 and $shard_count",
    )
    cases = [
        case for (index, case) in enumerate(cases)
        if mod1(index, shard_count) == shard_index
    ]
    @test !isempty(cases)

    environment = GaugefieldsEnvironment(
        process_grid=(1, 1, 1, 1),
        communicator=Gaugefields.SerialCommunicator(),
        element_type=ComplexF64,
        verbose=0,
    )
    mktempdir() do directory
        for (index, case) in enumerate(cases)
            path = joinpath(WIZARD_FERMION_DIRECTORY, case.filename)
            dictionary = TOML.parsefile(path)
            measurements = dictionary["Measurement set"]["measurement_methods"]
            @test Set(keys(measurements)) == Set(["Plaquette"])
            @test dictionary["gradientflow_measurements"]["hasgradientflow"] == false
            @test isempty(
                dictionary["gradientflow_measurements"]["measurements_for_flow"],
            )

            parameters = suppress_wizard_fermion_output() do
                wizard_fermion_params(path, directory)
            end
            try
                @test parameters.L ==
                      (case.variant.operator == "HISQ" ?
                       (4, 4, 4, 4) : (2, 2, 2, 2))
                @test parameters.NC == case.colors
                @test parameters.initial == "hot"
                @test parameters.Nwing ==
                      (case.variant.operator == "HISQ" ? 3 : 1)
                @test parameters.update_method == case.update_method
                @test parameters.quench == false
                @test parameters.Dirac_operator == case.variant.operator
                if case.variant.operator in ("Staggered", "HISQ")
                    @test parameters.Nf == case.variant.flavors
                elseif case.variant.operator == "WilsonClover"
                    @test parameters.Clover_coefficient == 1.5612
                end
                if case.variant.operator == "HISQ"
                    @test parameters.naik_epsilon == -0.083
                    @test case.colors == 3
                    @test !case.stout
                end
                if case.variant.operator == "MobiusDomainwall"
                    @test parameters.b == 2.0
                    @test parameters.c == 1.0
                end
                @test parameters.SextonWeingargten == case.sexton
                @test (parameters.smearing_for_fermion == "stout") == case.stout

                config = LQCDConfig(parameters)
                fermion = only(config.fermions)
                @test fermion.operator isa case.variant.type
                @test fermion.flavors == case.variant.flavors
                @test (fermion.smearing isa StoutFermionSmearingConfig) ==
                      case.stout
                @test (config.update.md.integrator isa SextonWeingartenConfig) ==
                      case.sexton
                @test (config.update isa SLHMCConfig) ==
                      (case.update_method == "SLHMC")
                if config.update isa SLHMCConfig
                    expected_beta = case.colors == 3 ? 5.2 : 2.2
                    @test only(config.update.md_gauge_action.terms).coupling ==
                          expected_beta
                end

                simulation = suppress_wizard_fermion_output() do
                    build_simulation(config, environment)
                end
                @test simulation.configuration isa GaugeFermionConfiguration
                @test keys(simulation.configuration.fermions) == (:fermion_1,)
                @test simulation.action isa Gaugefields.MDActionSet
                @test keys(simulation.action.terms) == (:gauge, :fermion_1)
                @test (simulation.updater.md_driver.integrator isa
                       Gaugefields.SextonWeingarten) == case.sexton
                @test (simulation.updater isa SLHMCUpdater) ==
                      (case.update_method == "SLHMC")
                @test Gaugefields.measure_plaquette(
                    simulation.configuration.gauge,
                ) isa Real

                spec = suppress_wizard_fermion_output() do
                    SimulationSpec(parameters)
                end
                @test isempty(validate(spec))
                events = RecordingSimulationEventSink()
                session = suppress_wizard_fermion_output() do
                    build_simulation(spec, environment; sink=events)
                end
                reference_result = suppress_wizard_fermion_output() do
                    update!(simulation)
                end
                session_result = suppress_wizard_fermion_output() do
                    step!(session)
                end
                test_wizard_fermion_update_equal(
                    reference_result,
                    session_result.update,
                )
                if case.update_method == "SLHMC"
                    @test reference_result isa SLHMCUpdateResult
                    @test isfinite(reference_result.md_delta_hamiltonian)
                end
                @test wizard_fermion_bridge_bytes(
                    simulation.configuration,
                ) == wizard_fermion_bridge_bytes(
                    session.simulation.configuration,
                )
                @test length(session_result.measurements) == 1
                @test only(session_result.measurements).name === :plaquette
                @test only(session_result.measurements).value ≈
                      Gaugefields.measure_plaquette(
                    simulation.configuration.gauge,
                )
                @test count(
                    event -> event isa TrajectoryFinished,
                    events.events,
                ) == 1
                @test count(
                    event -> event isa MeasurementFinished,
                    events.events,
                ) == 1
            finally
                isopen(parameters.load_fp) && close(parameters.load_fp)
            end
            index % 8 == 0 && GC.gc(false)
        end
    end
end
