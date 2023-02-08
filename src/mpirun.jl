using MPI
include("./LatticeQCD.jl")
using .LatticeQCD
#using LatticeQCD

function main()
    
    myrank = get_myrank()
    nprocs = get_nprocs()
    println(myrank)

    println_rank0(nprocs)
    println_rank0(get_PEs()) 

    println_rank0(ARGS)

    PEs = parse.(Int64,ARGS[2:5])

    set_PEs(PEs)

    println_rank0(get_PEs())

    run_LQCD(ARGS[1];MPIparallel=true)
end

main()