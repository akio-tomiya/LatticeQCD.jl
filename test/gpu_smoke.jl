# Optional direct-CUDA regression. Run this file from an environment that
# contains LatticeQCD, JACC, CUDA, and has the JACC backend set to "cuda".
using CUDA
import JACC

JACC.@init_backend

using LatticeQCD
using Test
import Gaugefields
import LatticeDiracOperators

CUDA.device!(parse(Int, get(ENV, "CUDA_TEST_DEVICE", "0")))
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
            1.0,
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
    for (name, operator) in (
        (:gauge, nothing),
        (:staggered, StaggeredDiracConfig(0.5)),
        (:hisq, HISQDiracConfig(0.5, -0.083)),
    )
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
end

println("GPU=", CUDA.name(CUDA.device()))
println("JACC_BACKEND=", JACC.backend)
println("GF=", pkgversion(Gaugefields))
println("LM=", pkgversion(LatticeDiracOperators.LatticeMatrices))
