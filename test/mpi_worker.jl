using MPI
MPI.Initialized() || MPI.Init()

using LatticeQCD
using Test
import Gaugefields
import LatticeDiracOperators

const LatticeMatrices = LatticeDiracOperators.LatticeMatrices

function mpi_hmc_input()
    lattice = LatticeConfig((4, 2, 2, 2))
    gauge = GaugeConfig(2, 1, "hot", nothing, 0x1234)
    action = GaugeActionConfig(
        GaugeActionTermConfig(:gauge_plaquette, "plaquette", 1.9),
    )
    integrator = LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge))
    update = HMCConfig(
        MDConfig(0.02, 2, integrator),
        GaussianMomentumConfig(
            1.0,
            RandomStreamConfig(0x5678, :momentum),
        ),
        RankZeroMetropolisConfig(
            RandomStreamConfig(0x9abc, :metropolis),
        ),
    )
    return LQCDConfig(lattice, gauge, action, update)
end

@testset "LatticeQCD MPI extension worker" begin
    communicator = MPI.COMM_WORLD
    rank = MPI.Comm_rank(communicator)
    size = MPI.Comm_size(communicator)
    size in (1, 2) || error("the MPI test supports one or two ranks; got $size")

    @test Base.get_extension(LatticeQCD, :LatticeQCDMPIExt) !== nothing
    @test Base.get_extension(Gaugefields, :GaugefieldsMPIExt) !== nothing
    @test Base.get_extension(
        LatticeDiracOperators,
        :LatticeDiracOperatorsMPIExt,
    ) !== nothing
    @test Base.get_extension(LatticeMatrices, :LatticeMatricesMPIExt) !== nothing

    communication = LatticeQCD.LQCDCommunication
    @test communication.default_communicator() == communicator
    @test communication.communicator_ready(communicator)
    @test communication.comm_size(communicator) == size
    @test communication.comm_rank(communicator) == rank
    @test communication.is_distributed(communicator) == (size > 1)
    @test communication.is_root(communicator) == (rank == 0)
    @test get_myrank() == rank
    @test get_nprocs() == size

    root_value = Ref(rank == 0)
    communication.broadcast!(root_value, 0, communicator)
    @test root_value[]

    environment = GaugefieldsEnvironment(
        process_grid=(size, 1, 1, 1),
        communicator=communicator,
        element_type=ComplexF64,
        verbose=0,
    )
    simulation = build_simulation(mpi_hmc_input(), environment)
    @test Gaugefields.gauge_process_grid(
        simulation.configuration.gauge,
    ) == (size, 1, 1, 1)

    result = update!(simulation)
    @test isfinite(result.initial_hamiltonian)
    @test isfinite(result.final_hamiltonian)
    @test isfinite(result.delta_hamiltonian)

    minimum_delta = MPI.Allreduce(
        result.delta_hamiltonian,
        MPI.MIN,
        communicator,
    )
    maximum_delta = MPI.Allreduce(
        result.delta_hamiltonian,
        MPI.MAX,
        communicator,
    )
    @test minimum_delta ≈ maximum_delta rtol=0 atol=1e-13

    root_accepted = MPI.bcast(result.accepted, 0, communicator)
    @test result.accepted == root_accepted
end
