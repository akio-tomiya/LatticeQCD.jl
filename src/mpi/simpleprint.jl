module Simpleprint
import ..LQCDCommunication: default_communicator, is_root

function println_rank0(jj...)
    is_root(default_communicator()) && println(jj...)
end

end
