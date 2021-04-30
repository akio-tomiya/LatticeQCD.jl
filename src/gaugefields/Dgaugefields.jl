module DFields
    using Distributed
    using DistributedArrays
    using LinearAlgebra

    import ..Gaugefields:SUn,SU,Gauge,GaugeFields,Gauge_1d,Adjoint_Gauge_1d
    import ..Gaugefields


    struct DGaugeFields{T <: SUn}  <: Gauge{T}
        g::DArray{ComplexF64,6}
        NX::Int64
        NY::Int64
        NZ::Int64
        NT::Int64
        NC::Int64
        NV::Int64
        parallel::Array{Int64,1}
        g_workers::Array{Int64,1}

        function DGaugeFields(NC,NX,NY,NZ,NT,parallel)
            sutype = SU{NC}
            @assert nworkers() == prod(parallel) "There are $(workers()) but $parallel"
                
            NV = NX*NY*NZ*NT
            g_workers = workers()
            g = dzeros(ComplexF64,(NC,NC,NX,NY,NZ,NT),g_workers,parallel)

            return new{sutype}(g,NX,NY,NZ,NT,NC,NV,parallel,g_workers)
        end

        #g = zeros(ComplexF64,NC,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW)

        function DGaugeFields(x::GaugeFields{T},parallel) where T
            g_workers = workers()
            g = distribute(x.g[1:x.NC,1:x.NC,1+x.NDW:x.NX+x.NDW,1+x.NDW:x.NY+x.NDW,1+x.NDW:x.NZ+x.NDW,1+x.NDW:x.NT+x.NDW],procs = g_workers,dist = parallel)
            return new{T}(g,x.NX,x.NY,x.NZ,x.NT,x.NC,x.NV,parallel,g_workers)
        end
    end

    struct DGaugeFields_1d{T <: SUn} <: Gauge_1d{T}
        g::DArray{ComplexF64,3}
        NC::Int64
        NV::Int64
        numworkers::Int64
    
        function DGaugeFields_1d(NC,NX,NY,NZ,NT)
            sutype = SU{NC}

            NV = NX*NY*NZ*NT
            numworkers = nworkers()
            g = dzeros(ComplexF64,(NC,NC,NV),workers(),[1,1,numworkers])
            return new{sutype}(g,NC,NV,numworkers)
        end
    
        function DGaugeFields_1d(NC,NV) 
            sutype = SU{NC}
            g = dzeros(ComplexF64,(NC,NC,NV),workers,[1,1,numworkers])
            return new{sutype}(g,NC,NV,numworkers)
        end
    end

    struct Adjoint_DGaugeFields_1d{T <: SUn} <: Adjoint_Gauge_1d{T}
        parent::DGaugeFields_1d{T}
    end

    function Base.adjoint(x::DGaugeFields_1d{T}) where T <: SUn
        Adjoint_DGaugeFields_1d{T}(x)
    end

    function DistributedArrays.localindices(a::DGaugeFields)
        return localindices(a.g)
    end

    function DistributedArrays.localpart(a::DGaugeFields)
        return localpart(a.g)
    end

    function DistributedArrays.localpart(a::DGaugeFields_1d)
        return localpart(a.g)
    end

    function DistributedArrays.localindices(a::DGaugeFields_1d)
        return localindices(a.g)
    end

    function Base.getindex(x::DGaugeFields,i1,i2,i3,i4,i5,i6)
        #indices_g = localindices(x.g)
        return x.g[i1,i2,i3,i4,i5,i6]
    end

    function Base.getindex(x::DGaugeFields_1d,i1,i2,i3)
        #indices_g = localindices(x.g)
        return x.g[i1,i2,i3]
    end

    function Base.setindex!(x::DGaugeFields,v,i1,i2,i3,i4,i5,i6) 
        indices_g = localindices(x.g)
        #println(indices_g)
        localpart(x.g)[i1,i2,i3-indices_g[3][begin]+1,i4-indices_g[4][begin]+1,i5-indices_g[5][begin]+1,i6-indices_g[6][begin]+1] = v
        return
    end

    function loop_substitute!(a::DGaugeFields_1d{SU{NC}},b::DGaugeFields{SU{NC}}) where NC
        #println("U ",localindices(U))
        i_local = localindices(b)
        icum = 0

        for it in i_local[6]
            for iz in i_local[5]
                for iy in i_local[4]
                    for ix in i_local[3]
                        icum += 1
                        println("$ix $iy $iz $it")
                        for j=1:NC
                            for i=1:NC
                                localpart(a)[i,j,icum] = b[i,j,ix,iy,iz,it]
                            end
                        end
                    end
                end
            end
        end
    end


    function Gaugefields.substitute!(a::DGaugeFields_1d{SU{NC}},b::DGaugeFields{SU{NC}}) where NC

        f = []
        @async for id in workers()
            push!(f,remotecall(loop_substitute!,id,a,b))
        end

        for fi in f
            wait(fi)
        end

    end

    function loop_clear!(a::DGaugeFields_1d)
        localpart(a.g) .= 0
    end

    function Gaugefields.clear!(a::DGaugeFields_1d)
        f = []
        @async for id in workers()
            push!(f,remotecall(loop_clear!,id,a))
        end
        for fi in f
            wait(fi)
        end

        return 
    end

    function loop_mul!(c::DGaugeFields_1d{SU{NC}},a::DGaugeFields_1d{SU{NC}},b::DGaugeFields_1d{SU{NC}}) where NC
        i_local = localindices(c)
        icum = 0

        for i in i_local[3]
            icum += 1
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] = 0
                    @simd for k3=1:NC
                        localpart(c)[k1,k2,icum] += a[k1,k3,i]*b[k3,k2,i]
                    end
                end
            end
        end
    end

    function LinearAlgebra.mul!(c::DGaugeFields_1d{SU{NC}},a::DGaugeFields_1d{SU{NC}},b::DGaugeFields_1d{SU{NC}}) where NC
        f = []
        @async for id in workers()
            push!(f,remotecall(loop_mul!,id,c,a,b))
        end
        for fi in f
            wait(fi)
        end

    end

    function loop_mul!(c::DGaugeFields_1d{SU{NC}},a::DGaugeFields_1d{SU{NC}},b::Adjoint_DGaugeFields_1d{SU{NC}}) where NC
        i_local = localindices(c)
        icum = 0

        for i in i_local[3]
            icum += 1
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] = 0
                    @simd for k3=1:NC
                        localpart(c)[k1,k2,icum] += a[k1,k3,i]*conj(b.parent[k3,k2,i])
                    end
                end
            end
        end
        

    end

    function LinearAlgebra.mul!(c::DGaugeFields_1d{SU{NC}},a::DGaugeFields_1d{SU{NC}},b::Adjoint_DGaugeFields_1d{SU{NC}}) where NC

        f = []
        @async for id in workers()
            push!(f,remotecall(loop_mul!,id,c,a,b))
        end
        for fi in f
            wait(fi)
        end
        

    end

    function loop_gauge_shift!(a::DGaugeFields_1d{SU{NC}},ν::N,b::DGaugeFields{SU{NC}}) where {N <: Int,NC}
        i_local = localindices(a)
        icum = 0

        idel = zeros(Int8,4)
        if ν > 0
            idel[ν] = +1
        else
            idel[-ν] = -1
        end

        NT = b.NT
        NZ = b.NZ
        NY = b.NY
        NX = b.NX


        for it in i_local[6]
            it1 = it + idel[4]
            it1 += ifelse(it1 > NT,-NT,0) + ifelse(it1 < 1,+NT,0)
            for iz in i_local[5]
                iz1 = iz + idel[3]
                iz1 += ifelse(iz1 > NZ,-NZ,0) + ifelse(iz1 < 1,+NZ,0)

                for iy in i_local[4]
                    iy1 = iy+idel[2]
                    iy1 += ifelse(iy1 > NY,-NY,0) + ifelse(iy1 < 1,+NY,0)
                    for ix in i_local[3]
                        ix1 = ix+idel[1]
                        ix1 += ifelse(ix1 > NX,-NX,0) + ifelse(ix1 < 1,+NX,0)
                        icum += 1

                        for k2=1:NC
                            @simd for k1=1:NC
                                localpart(a)[k1,k2,icum] = b[k1,k2,ix1,iy1,iz1,it1]
                                #println(localpart(a)[k1,k2,icum])
                            end
                        end
                    end
                end
            end
        end
    end

    function Gaugefields.gauge_shift!(a::DGaugeFields_1d{SU{NC}},ν::N,b::DGaugeFields{SU{NC}}) where {N <: Int,NC}

        if ν == 0
            substitute!(a,b)
            return
        end

        f = []
        @async for id in workers()
            push!(f,remotecall(loop_gauge_shift!,id,c,a,b))
        end
        for fi in f
            wait(fi)
        end

        return
    end

    function loop_add!(c::DGaugeFields_1d{SU{NC}},a::DGaugeFields_1d{SU{NC}}) where NC
        i_local = localindices(c)
        icum = 0

        for i in i_local[3]
            icum += 1
            for k2=1:NC                            
                @simd  for k1=1:NC
                    localpart(c)[k1,k2,icum] += a[k1,k2,i] 
                end
            end
        end

    end

    function Gaugefields.add!(c::DGaugeFields_1d{SU{NC}},a::DGaugeFields_1d{SU{NC}}) where NC

        f = []
        @async for id in workers()
            push!(f,remotecall(loop_add!,id,c,a))
        end
        for fi in f
            wait(fi)
        end


    end

    function loop_tr(a::DGaugeFields_1d{SU{NC}}) where NC

        i_local = localindices(a)
        icum = 0

        s = 0
        for i in i_local[3]
            icum += 1
            @simd for k=1:NC
                
                s += a[k,k,i]
            end

        end
        return s

    end

    function LinearAlgebra.tr(a::DGaugeFields_1d{SU{NC}}) where NC

        f = []
        #s = 0
        @async for id in workers()
            push!(f,remotecall(loop_tr,id,a))

        end
        s = 0
        for fi in f
            wait(fi)
            s += fetch(fi)
        end
        #s = sum(fetch.(f))
        return s


    end

    function Gaugefields.calc_Plaq(U::Array{T,1}) where T <: DGaugeFields
        NC = U[1].NC
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT

        temp1 = DGaugeFields_1d(NC,NX,NY,NZ,NT)
        temp2 = DGaugeFields_1d(NC,NX,NY,NZ,NT)
        temp3 = DGaugeFields_1d(NC,NX,NY,NZ,NT)
        staple = DGaugeFields_1d(NC,NX,NY,NZ,NT)

        

        return Gaugefields.calc_Plaq!(U,temp1,temp2,temp3,staple)
    end


end