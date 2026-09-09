using LatticeQCD
using Test

@testset "MPI weak dependency lifecycle" begin
    @test Base.get_extension(LatticeQCD, :LatticeQCDMPIExt) === nothing
    @test all(package_id.name != "MPI" for package_id in keys(Base.loaded_modules))

    @eval using MPI

    @test Base.get_extension(LatticeQCD, :LatticeQCDMPIExt) !== nothing
    @test !MPI.Initialized()
    @test LatticeQCD.LQCDCommunication.default_communicator() === nothing

    MPI.Init()
    @test LatticeQCD.LQCDCommunication.default_communicator() == MPI.COMM_WORLD
    @test get_myrank() == 0
    @test get_nprocs() == 1
end
