using MPI
using Test

@testset "LatticeQCD MPI integration" begin
    project_directory = dirname(Base.active_project())
    lifecycle = joinpath(@__DIR__, "mpi_lifecycle.jl")
    worker = joinpath(@__DIR__, "mpi_worker.jl")

    lifecycle_command = `$(Base.julia_cmd()) --startup-file=no --project=$(project_directory) $lifecycle`
    lifecycle_process = run(ignorestatus(lifecycle_command))
    @test lifecycle_process.exitcode == 0

    for ranks in (1, 2)
        command = `$(MPI.mpiexec()) -n $ranks $(Base.julia_cmd()) --startup-file=no --project=$(project_directory) $worker`
        process = run(ignorestatus(command))
        @test process.exitcode == 0
    end
end
