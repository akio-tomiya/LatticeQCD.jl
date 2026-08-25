module LatticeQCDMPIExt

using LatticeQCD
using MPI

import LatticeQCD.LQCDCommunication:
    broadcast!,
    comm_rank,
    comm_size,
    communicator_ready

function default_communicator()
    return MPI.Initialized() && !MPI.Finalized() ? MPI.COMM_WORLD : nothing
end

@inline communicator_ready(::MPI.Comm) =
    MPI.Initialized() && !MPI.Finalized()
@inline comm_size(comm::MPI.Comm) = MPI.Comm_size(comm)
@inline comm_rank(comm::MPI.Comm) = MPI.Comm_rank(comm)
@inline broadcast!(value, root::Integer, comm::MPI.Comm) =
    MPI.Bcast!(value, root, comm)

end
