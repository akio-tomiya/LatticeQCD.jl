module LQCDCommunication

import Gaugefields: SerialCommunicator

"""Return the active optional communicator, or `nothing` in serial mode."""
function default_communicator()
    extension = Base.get_extension(
        parentmodule(@__MODULE__),
        :LatticeQCDMPIExt,
    )
    return extension === nothing ? nothing : extension.default_communicator()
end

@inline communicator_ready(::Nothing) = true
@inline comm_size(::Nothing) = 1
@inline comm_rank(::Nothing) = 0

@inline communicator_ready(::SerialCommunicator) = true
@inline comm_size(::SerialCommunicator) = 1
@inline comm_rank(::SerialCommunicator) = 0

@inline function broadcast!(value, root::Integer, ::Nothing)
    iszero(root) || throw(ArgumentError(
        "the serial communicator only contains root rank 0, got root $root",
    ))
    return value
end

@inline function broadcast!(value, root::Integer, ::SerialCommunicator)
    iszero(root) || throw(ArgumentError(
        "the serial communicator only contains root rank 0, got root $root",
    ))
    return value
end

function check_communicator(comm)
    communicator_ready(comm) || throw(ArgumentError(
        "the communicator is not ready; initialize it before starting a run",
    ))
    return comm
end

@inline is_distributed(comm) = comm_size(check_communicator(comm)) > 1
@inline is_root(comm) = comm_rank(check_communicator(comm)) == 0

get_myrank() = comm_rank(default_communicator())
get_nprocs() = comm_size(default_communicator())

const LEGACY_PROCESS_GRID = [1, 1, 1, 1]

function set_PEs(process_grid)
    length(process_grid) == length(LEGACY_PROCESS_GRID) || throw(
        DimensionMismatch(
            "the legacy process grid must have " *
            "$(length(LEGACY_PROCESS_GRID)) entries",
        ),
    )
    all(>(0), process_grid) || throw(ArgumentError(
        "all process-grid entries must be positive; got $process_grid",
    ))
    LEGACY_PROCESS_GRID .= Int.(process_grid)
    return get_PEs()
end

get_PEs() = copy(LEGACY_PROCESS_GRID)

end
