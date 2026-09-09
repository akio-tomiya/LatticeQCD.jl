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
            sqrt(2.0),
            RandomStreamConfig(0x5678, :momentum),
        ),
        RankZeroMetropolisConfig(
            RandomStreamConfig(0x9abc, :metropolis),
        ),
    )
    return LQCDConfig(lattice, gauge, action, update)
end

function mpi_wilson_hmc_input()
    lattice = LatticeConfig((4, 2, 2, 2))
    gauge = GaugeConfig(2, 1, "cold", nothing, 0x1234)
    action = GaugeActionConfig(
        GaugeActionTermConfig(:gauge_plaquette, "plaquette", 1.9),
    )
    solver = FermionSolverConfig(
        1e-10,
        2_000,
        0,
        (1, 1, 1, -1),
    )
    fermions = (
        FermionActionConfig(
            :fermion_1,
            WilsonDiracConfig(0.05, 1.0),
            2,
            solver,
        ),
    )
    integrator = LeapfrogConfig(
        QPQConfig(),
        ForceGroupConfig(:gauge, :fermion_1),
    )
    update = HMCConfig(
        MDConfig(0.001, 1, integrator),
        GaussianMomentumConfig(
            sqrt(2.0),
            RandomStreamConfig(0x5678, :momentum),
        ),
        RankZeroMetropolisConfig(
            RandomStreamConfig(0x9abc, :metropolis),
        ),
        (
            PseudofermionRefreshConfig(
                :fermion_1,
                RandomStreamConfig(0xdef0, :pseudofermion),
                1,
            ),
        ),
    )
    return LQCDConfig(lattice, gauge, action, fermions, update)
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

    schedule = SimulationSchedule(0, 3)
    reference = build_simulation(
        SimulationSpec(mpi_hmc_input(), schedule),
        environment,
    )
    run!(reference; verbose=false)

    checkpoint_directory = MPI.bcast(
        rank == 0 ? mktempdir() : "",
        0,
        communicator,
    )
    checkpoints = JLD2CheckpointOutput(
        checkpoint_directory;
        every=1,
    )
    checkpoint_spec = SimulationSpec(
        mpi_hmc_input(),
        schedule,
        OutputConfig(; checkpoints),
    )
    interrupted = build_simulation(checkpoint_spec, environment)
    step_result = step!(interrupted)
    @test step_result.checkpoint_path ==
          checkpoint_output_path(checkpoints, 1)

    restored = build_simulation(checkpoint_spec, environment)
    load_checkpoint!(restored, step_result.checkpoint_path)
    run!(restored; verbose=false)
    @test restored.simulation.state.accepted ==
          reference.simulation.state.accepted
    for direction in eachindex(reference.simulation.configuration.gauge)
        @test restored.simulation.configuration.gauge[direction].U.A ==
              reference.simulation.configuration.gauge[direction].U.A
    end

    fermion_schedule = SimulationSchedule(0, 2)
    fermion_reference = build_simulation(
        SimulationSpec(mpi_wilson_hmc_input(), fermion_schedule),
        environment,
    )
    run!(fermion_reference; verbose=false)
    fermion_checkpoints = JLD2CheckpointOutput(
        checkpoint_directory;
        prefix="fermion_",
        every=1,
    )
    fermion_spec = SimulationSpec(
        mpi_wilson_hmc_input(),
        fermion_schedule,
        OutputConfig(; checkpoints=fermion_checkpoints),
    )
    fermion_interrupted = build_simulation(fermion_spec, environment)
    fermion_checkpoint = step!(fermion_interrupted).checkpoint_path
    fermion_restored = build_simulation(fermion_spec, environment)
    load_checkpoint!(fermion_restored, fermion_checkpoint)
    run!(fermion_restored; verbose=false)
    @test fermion_restored.simulation.state.accepted ==
          fermion_reference.simulation.state.accepted
    for direction in eachindex(
        fermion_reference.simulation.configuration.gauge,
    )
        @test fermion_restored.simulation.configuration.gauge[direction].U.A ==
              fermion_reference.simulation.configuration.gauge[direction].U.A
    end

    MPI.Barrier(communicator)
    rank == 0 && rm(checkpoint_directory; recursive=true)
end
