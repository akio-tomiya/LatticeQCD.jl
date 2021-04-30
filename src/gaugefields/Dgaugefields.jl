
module DFields
    using Distributed
    using DistributedArrays
    using LinearAlgebra

    import ..Gaugefields:SUn,SU,Gauge


    struct DGaugeFields{T <: SUn}  <: Gauge{T}
        g::DArray{ComplexF64,6}
        NX::Int64
        NY::Int64
        NZ::Int64
        NT::Int64
        NC::Int64
        NV::Int64
        parallel::Array{Int64,1}

        function DGaugeFields(NC,NX,NY,NZ,NT,parallel)
            sutype = SU{NC}
            @assert nworkers() == prod(parallel) "There are $(workers()) but $parallel"
                
            NV = NX*NY*NZ*NT
            g = dzeros(ComplexF64,(NC,NC,NX,NY,NZ,NT),workers(),parallel)

            return new{sutype}(g,NX,NY,NZ,NT,NC,NV,parallel)
        end
    end
end