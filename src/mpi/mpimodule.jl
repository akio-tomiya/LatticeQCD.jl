module MPImodules

import ..Simpleprint: println_rank0
import ..LQCDCommunication:
    get_myrank,
    get_nprocs,
    get_PEs,
    set_PEs

export get_myrank, get_nprocs, println_rank0, set_PEs, get_PEs

end
