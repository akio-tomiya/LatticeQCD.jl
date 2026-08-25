using Test
import Gaugefields

@testset "MPI remains optional" begin
    @test Base.find_package("MPI") === nothing
    @test all(package_id.name != "MPI" for package_id in keys(Base.loaded_modules))
    @test Base.get_extension(LatticeQCD, :LatticeQCDMPIExt) === nothing
    @test Base.get_extension(Gaugefields, :GaugefieldsMPIExt) === nothing

    communication = LatticeQCD.LQCDCommunication
    serial = Gaugefields.SerialCommunicator()
    @test communication.default_communicator() === nothing
    @test communication.communicator_ready(serial)
    @test !communication.is_distributed(serial)
    @test communication.is_root(serial)
    @test communication.comm_size(serial) == 1
    @test communication.comm_rank(serial) == 0

    decision = Ref(true)
    @test communication.broadcast!(decision, 0, serial) === decision
    @test decision[]

    @test get_myrank() == 0
    @test get_nprocs() == 1
    original_grid = get_PEs()
    try
        @test set_PEs((1, 1, 1, 1)) == [1, 1, 1, 1]
        @test_throws DimensionMismatch set_PEs((1, 1))
        @test_throws ArgumentError set_PEs((1, 1, 0, 1))
    finally
        set_PEs(original_grid)
    end
end
