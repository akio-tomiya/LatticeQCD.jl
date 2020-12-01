module Parallel
    using Base.Threads
    function get_looprange(N)
        id = threadid()
        n = nthreads()
        dn = N % n
        nbun = div(N - dn,n)
        ista = (id-1)*nbun + 1 
        ista += ifelse(id <= dn,(id-1),dn)
        nbun += ifelse(id <= dn,1,0)
        iend = ista + nbun - 1
        return ista:iend
    end
end