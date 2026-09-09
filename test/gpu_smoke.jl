# Optional direct-CUDA regression. Run this file from an environment that
# contains LatticeQCD, JACC, CUDA, and has the JACC backend set to "cuda".
using CUDA
CUDA.device!(parse(Int, get(ENV, "CUDA_TEST_DEVICE", "0")))

import JACC

JACC.@init_backend

using LatticeQCD
using Test
import Gaugefields
import LatticeDiracOperators

CUDA.functional() || error("CUDA is not functional")
JACC.backend == "cuda" || error("JACC CUDA backend is not active")

function gpu_hmc_input(operator=nothing)
    uses_hisq = operator isa HISQDiracConfig
    lattice = LatticeConfig(uses_hisq ? (4, 4, 4, 4) : (2, 2, 2, 2))
    colors = uses_hisq ? 3 : 2
    halo = uses_hisq ? 3 : 1
    gauge = GaugeConfig(colors, halo, "cold", nothing, 0x1234)
    gauge_action = GaugeActionConfig(
        GaugeActionTermConfig(:gauge_plaquette, "plaquette", 1.9),
    )

    if operator === nothing
        fermions = ()
        forces = ForceGroupConfig(:gauge)
        refreshes = ()
    else
        solver = FermionSolverConfig(1e-10, 2_000, 0, (1, 1, 1, -1))
        fermions = (
            FermionActionConfig(
                :fermion_1,
                operator,
                4,
                solver,
            ),
        )
        forces = ForceGroupConfig(:gauge, :fermion_1)
        refreshes = (
            PseudofermionRefreshConfig(
                :fermion_1,
                RandomStreamConfig(0xdef0, :pseudofermion),
                1,
            ),
        )
    end

    update = HMCConfig(
        MDConfig(0.001, 1, LeapfrogConfig(QPQConfig(), forces)),
        GaussianMomentumConfig(
            sqrt(2.0),
            RandomStreamConfig(0x5678, :momentum),
        ),
        RankZeroMetropolisConfig(
            RandomStreamConfig(0x9abc, :metropolis),
        ),
        refreshes,
    )
    return LQCDConfig(
        lattice,
        gauge,
        gauge_action,
        fermions,
        update,
    )
end

environment = GaugefieldsEnvironment(
    process_grid=(1, 1, 1, 1),
    communicator=Gaugefields.SerialCommunicator(),
    element_type=ComplexF64,
    verbose=0,
)

@testset "LatticeQCD typed CUDA smoke" begin
    requested_cases = Set(Symbol.(split(
        get(ENV, "LQCD_CUDA_CASES", "gauge,staggered,hisq"),
        ',',
    )))
    cases = filter(
        case -> first(case) in requested_cases,
        (
            (:gauge, nothing),
            (:staggered, StaggeredDiracConfig(0.5)),
            (:hisq, HISQDiracConfig(0.5, -0.083)),
        ),
    )
    isempty(cases) && error("LQCD_CUDA_CASES selected no GPU smoke cases")
    for (name, operator) in cases
        config = gpu_hmc_input(operator)
        simulation = build_simulation(config, environment)
        @test simulation.configuration.gauge[1].U.A isa CUDA.CuArray
        result = update!(simulation)
        CUDA.synchronize()
        @test isfinite(result.initial_hamiltonian)
        @test isfinite(result.final_hamiltonian)
        @test isfinite(result.delta_hamiltonian)
        plaquette = Gaugefields.measure_plaquette(
            simulation.configuration.gauge,
        )
        @test isfinite(plaquette)

        mktempdir() do directory
            path = joinpath(directory, "$(name).jld2")
            save_configuration(
                path,
                simulation.configuration;
                format=:jld2,
            )
            @test isfile(path)
            reloaded = build_configuration(config, environment)
            load_configuration!(reloaded, path; format=:jld2)
            @test Gaugefields.measure_plaquette(reloaded.gauge) ≈ plaquette
        end
        println(
            "CUDA_CASE name=", name,
            " plaquette=", plaquette,
            " delta_hamiltonian=", result.delta_hamiltonian,
        )
    end

    requested_restart = Symbol(get(
        ENV,
        "LQCD_CUDA_RESTART_CASE",
        "staggered",
    ))
    restart_matches = filter(case -> first(case) === requested_restart, cases)
    isempty(restart_matches) && error(
        "LQCD_CUDA_RESTART_CASE=$requested_restart is not selected by " *
        "LQCD_CUDA_CASES",
    )
    restart_name, restart_operator = only(restart_matches)
    @testset "$(restart_name) restart checkpoint remains on CUDA" begin
        config = gpu_hmc_input(restart_operator)
        schedule = SimulationSchedule(0, 2)
        reference = build_simulation(
            SimulationSpec(config, schedule),
            environment,
        )
        run!(reference; verbose=false)

        mktempdir() do directory
            checkpoints = JLD2CheckpointOutput(directory; every=1)
            spec = SimulationSpec(
                config,
                schedule,
                OutputConfig(; checkpoints),
            )
            interrupted = build_simulation(spec, environment)
            checkpoint = step!(interrupted).checkpoint_path
            restored = build_simulation(spec, environment)
            load_checkpoint!(restored, checkpoint)
            @test restored.simulation.configuration.gauge[1].U.A isa
                  CUDA.CuArray
            run!(restored; verbose=false)
            CUDA.synchronize()
            @test restored.simulation.state.accepted ==
                  reference.simulation.state.accepted
            for direction in eachindex(
                reference.simulation.configuration.gauge,
            )
                @test Array(
                    restored.simulation.configuration.gauge[direction].U.A,
                ) == Array(
                    reference.simulation.configuration.gauge[direction].U.A,
                )
            end
        end
    end
end

println("GPU=", CUDA.name(CUDA.device()))
println("JACC_BACKEND=", JACC.backend)
println("GF=", pkgversion(Gaugefields))
println("LM=", pkgversion(LatticeDiracOperators.LatticeMatrices))
