module MPImodules
    using MPI
    import ..Simpleprint:println_rank0
    MPI.Init()
    const comm = MPI.COMM_WORLD
    const myrank = MPI.Comm_rank(comm)
    const nprocs = MPI.Comm_size(comm)

    struct MPI_lattice
        PEs::Vector{Int64}
    end

    const mpi_lattice = MPI_lattice([1,1,1,1])

    function get_myrank()
        return myrank
    end

    function get_nprocs()
        return nprocs
    end

    function println_rank0(jj...)
        if myrank == 0
            println(jj...)
        end
    end

    function set_PEs(PEs)
        for i=1:length(PEs)
            mpi_lattice.PEs[i] = PEs[i]
        end
        
    end

    function get_PEs()
        return mpi_lattice.PEs
    end


end