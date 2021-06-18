module DFields
    using Distributed
    using DistributedArrays
    using LinearAlgebra

    import ..Gaugefields:SUn,SU,Gauge,GaugeFields,Gauge_1d,Adjoint_Gauge_1d
    import ..Gaugefields


    struct DGaugeFields{T <: SUn}  <: Gauge{T}
        g::DArray{ComplexF64,6,Array{ComplexF64,6}}
        NX::Int64
        NY::Int64
        NZ::Int64
        NT::Int64
        NC::Int64
        NV::Int64
        parallel::Array{Int64,1}
        g_workers::Array{Int64,1}
        #localpart::Array{ComplexF64,6}
        

        function DGaugeFields(NC,NX,NY,NZ,NT,parallel)
            sutype = SU{NC}
            @assert nworkers() == prod(parallel) "There are $(workers()) but $parallel"
                
            NV = NX*NY*NZ*NT
            g_workers = workers()
            g = dzeros(ComplexF64,(NC,NC,NX,NY,NZ,NT),g_workers,parallel)
            #localpart = g.localpart

            return new{sutype}(g,NX,NY,NZ,NT,NC,NV,parallel,g_workers)
        end

        #g = zeros(ComplexF64,NC,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW)

        function DGaugeFields(x::GaugeFields{T},parallel) where T
            g_workers = workers()
            @assert nworkers() == prod(parallel) "There are $(workers()) but $parallel"

            g = distribute(x.g[1:x.NC,1:x.NC,1+x.NDW:x.NX+x.NDW,1+x.NDW:x.NY+x.NDW,1+x.NDW:x.NZ+x.NDW,1+x.NDW:x.NT+x.NDW],procs = g_workers,dist = parallel)
            #localpart = g.localpart
            return new{T}(g,x.NX,x.NY,x.NZ,x.NT,x.NC,x.NV,parallel,g_workers)
        end
    end

    struct DGaugeFields_1d{T <: SUn} <: Gauge_1d{T}
        g::DArray{ComplexF64,3,Array{ComplexF64,3}}
        NC::Int64
        NV::Int64
        numworkers::Int64
        glocalindices::NTuple{3,UnitRange{Int64}}
        #localpart::Array{ComplexF64,3}
    
        function DGaugeFields_1d(NC,NX,NY,NZ,NT)
            sutype = SU{NC}

            NV = NX*NY*NZ*NT
            numworkers = nworkers()
            g = dzeros(ComplexF64,(NC,NC,NV),workers(),[1,1,numworkers])
            glocalindices = localindices(g)

            #localpart = g.localpart
            return new{sutype}(g,NC,NV,numworkers,glocalindices)
        end
    
        function DGaugeFields_1d(NC,NV) 
            sutype = SU{NC}
            g = dzeros(ComplexF64,(NC,NC,NV),workers,[1,1,numworkers])
            glocalindices = localindices(g)

            #localpart = g.localpart
            return new{sutype}(g,NC,NV,numworkers,glocalindices)
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
        indices = localindices(x.g)
        lenind = 3
        isinside = true
        isinside *= i1 in indices[1]
        isinside *= i2 in indices[2]
        isinside *= i3 in indices[3]
        isinside *= i4 in indices[4]
        isinside *= i5 in indices[5]
        isinside *= i6 in indices[6]

        if isinside
            index = (i1-indices[1][1]+1,i2-indices[2][1]+1,i3-indices[3][1]+1,i4-indices[4][1]+1,i5-indices[5][1]+1,i6-indices[6][1]+1)
            v = localpart(x.g)[index...]
            return v
        else
            return x.g[i1,i2,i3,i4,i5,i6]
        end

        #=
        indices_g = localindices(x.g)
        #return x.g[i1,i2,i3,i4,i5,i6]
        if i3 in indices_g[3] && i4 in indices_g[4] && i5 in indices_g[5] && i6 in indices_g[6]
            return x.g.localpart[i1,i2,i3-indices_g[3][begin]+1,i4-indices_g[4][begin]+1,i5-indices_g[5][begin]+1,i6-indices_g[6][begin]+1]
        else
            return x.g[i1,i2,i3,i4,i5,i6]
        end
        =#
    end

    @inline @inbounds function getlocal(x::DArray,i...)
        return x.localpart[i...]
    end

    function Base.getindex(x::DGaugeFields_1d,i1,i2,i3)
        indices = localindices(x.g)
        lenind = 3
        isinside = true
        isinside *= i1 in indices[1]
        isinside *= i2 in indices[2]
        isinside *= i3 in indices[3]

        if isinside
            index = (i1-indices[1][1]+1,i2-indices[2][1]+1,i3-indices[3][1]+1)
            v = localpart(x.g)[index...]
            return v
        else
            return x.g[i1,i2,i3]
        end

    end

    function Base.setindex!(x::DGaugeFields,v,i1,i2,i3,i4,i5,i6) 
        indices_g = localindices(x.g)
        #println(indices_g)
        #localpart(x.g)[i1,i2,i3-indices_g[3][begin]+1,i4-indices_g[4][begin]+1,i5-indices_g[5][begin]+1,i6-indices_g[6][begin]+1] = v
        #x.g.localpart[i1,i2,i3-indices_g[3][begin]+1,i4-indices_g[4][begin]+1,i5-indices_g[5][begin]+1,i6-indices_g[6][begin]+1] = v
        x.localpart[i1,i2,i3-indices_g[3][begin]+1,i4-indices_g[4][begin]+1,i5-indices_g[5][begin]+1,i6-indices_g[6][begin]+1] = v
        return
    end

    function loop_substitute!(a::DGaugeFields_1d{SU{NC}},b::DGaugeFields{SU{NC}}) where NC
        #println("U ",localindices(U))
        i_local = localindices(b)
        icum = 0

        for it = 1:(i_local[6][end]-i_local[6][1]+1)
            for iz = 1:(i_local[5][end]-i_local[5][1]+1)
                for iy = 1:(i_local[4][end]-i_local[4][1]+1)
                    for ix = 1:(i_local[3][end]-i_local[3][1]+1)
        #for it in i_local[6]
        #    for iz in i_local[5]
        #        for iy in i_local[4]
        #            for ix in i_local[3]
                        icum += 1
                        #println("$ix $iy $iz $it")
                        for j=1:NC
                            for i=1:NC
                                #a.g.localpart[i,j,icum] = b[i,j,ix,iy,iz,it]
                                a.g.localpart[i,j,icum] = b.g.localpart[i,j,ix,iy,iz,it]
                                #a.localpart[i,j,icum] = b[i,j,ix,iy,iz,it]
                                #localpart(a)[i,j,icum] = b[i,j,ix,iy,iz,it]
                            end
                        end
                    end
                end
            end
        end
    end


    function Gaugefields.substitute!(a::DGaugeFields_1d{SU{NC}},b::DGaugeFields{SU{NC}}) where NC

        f = []
        for id in workers()
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
        for id in workers()
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

        for i = 1:(i_local[3][end]-i_local[3][1]+1)
        #for i in i_local[3]
            icum += 1
            for k2=1:NC                            
                for k1=1:NC
                    c.g.localpart[k1,k2,icum] = 0
                    #c.localpart[k1,k2,icum] = 0
                    #localpart(c)[k1,k2,icum] = 0
                    @simd for k3=1:NC
                        #=
                        println("a")
                        @time a1 = a[k1,k3,i]
                        
                        println("b")
                        @time b1 = b[k3,k2,i]
                        println("c")
                        #@time c.g.localpart[k1,k2,icum] += a1*b1
                        @time c.localpart[k1,k2,icum] += a1*b1
                        =#
                        #c.g.localpart[k1,k2,icum] += a[k1,k3,i]*b[k3,k2,i]
                        c.g.localpart[k1,k2,icum] += a.g.localpart[k1,k3,i]*b.g.localpart[k3,k2,i]
                        #c.localpart[k1,k2,icum] += a[k1,k3,i]*b[k3,k2,i]
                        #localpart(c)[k1,k2,icum] += a[k1,k3,i]*b[k3,k2,i]
                    end
                    #exit()
                end
            end
        end
    end

    function LinearAlgebra.mul!(c::DGaugeFields_1d{SU{NC}},a::DGaugeFields_1d{SU{NC}},b::DGaugeFields_1d{SU{NC}}) where NC
        f = []
        for id in workers()
            push!(f,remotecall(loop_mul!,id,c,a,b))
        end
        for fi in f
            wait(fi)
        end

    end

    function loop_mul!(c::DGaugeFields_1d{SU{NC}},a::DGaugeFields_1d{SU{NC}},b::Adjoint_DGaugeFields_1d{SU{NC}}) where NC
        i_local = localindices(c)
        icum = 0

        for i = 1:(i_local[3][end]-i_local[3][1]+1)
            icum += 1
            for k2=1:NC                            
                for k1=1:NC
                    c.g.localpart[k1,k2,icum] = 0
                    @simd for k3=1:NC
                        c.g.localpart[k1,k2,icum] += a.g.localpart[k1,k3,i]*conj(b.parent.g.localpart[k3,k2,i])
                    end
                end
            end
        end

        #=

        for i in i_local[3]
            icum += 1
            for k2=1:NC                            
                for k1=1:NC
                    #localpart(c)[k1,k2,icum] = 0
                    #c.localpart[k1,k2,icum] = 0
                    c.g.localpart[k1,k2,icum] = 0
                    @simd for k3=1:NC
                        #localpart(c)[k1,k2,icum] += a[k1,k3,i]*conj(b.parent[k3,k2,i])
                        c.g.localpart[k1,k2,icum] += a[k1,k3,i]*conj(b.parent[k3,k2,i])
                        #c.localpart[k1,k2,icum] += a[k1,k3,i]*conj(b.parent[k3,k2,i])
                        #println("a = ",a[k1,k3,i])
                        #println("b = ",conj(b.parent[k3,k2,i]))
                        #println(localpart(c)[k1,k2,icum])
                    end
                end
            end
        end
        =#
        

    end

    function LinearAlgebra.mul!(c::DGaugeFields_1d{SU{NC}},a::DGaugeFields_1d{SU{NC}},b::Adjoint_DGaugeFields_1d{SU{NC}}) where NC
        #println("mul!")
        f = []
        #println(workers())
        #w = workers()
        for id in workers()
            push!(f,remotecall(loop_mul!,id,c,a,b))
        end
        for fi in f
            wait(fi)
        end
        

    end

    function loop_gauge_shift!(a::DGaugeFields_1d{SU{NC}},ν::N,b::DGaugeFields{SU{NC}}) where {N <: Int,NC}
        i_local = localindices(b)
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
            #isinside = true
            it1 = it + idel[4]
            it1 += ifelse(it1 > NT,-NT,0) + ifelse(it1 < 1,+NT,0)
            isinside_t = it1 in i_local[6]
            for iz in i_local[5]
                iz1 = iz + idel[3]
                iz1 += ifelse(iz1 > NZ,-NZ,0) + ifelse(iz1 < 1,+NZ,0)
                isinside_z = iz1 in i_local[5]
                
                for iy in i_local[4]
                    iy1 = iy+idel[2]
                    iy1 += ifelse(iy1 > NY,-NY,0) + ifelse(iy1 < 1,+NY,0)
                    isinside_y = iy1 in i_local[4]
                    for ix in i_local[3]
                        ix1 = ix+idel[1]
                        ix1 += ifelse(ix1 > NX,-NX,0) + ifelse(ix1 < 1,+NX,0)
                        isinside_x = ix1 in i_local[3]

                        icum += 1
                        isinside = isinside_t*isinside_z*isinside_y*isinside_x

                        if isinside
                            index = (ix1-i_local[3][1]+1,iy1-i_local[4][1]+1,iz1-i_local[5][1]+1,it1-i_local[6][1]+1)
                            #v = localpart(x.g)[index...]
                            for k2=1:NC
                                @simd for k1=1:NC
                                    #localpart(a)[k1,k2,icum] = b[k1,k2,ix1,iy1,iz1,it1]
                                    a.g.localpart[k1,k2,icum] = b.g.localpart[k1,k2,index...]
                                    #a.localpart[k1,k2,icum] = b[k1,k2,ix1,iy1,iz1,it1]
                                    #println(localpart(a)[k1,k2,icum])
                                end
                            end

                        else
                            for k2=1:NC
                                @simd for k1=1:NC
                                    #localpart(a)[k1,k2,icum] = b[k1,k2,ix1,iy1,iz1,it1]
                                    a.g.localpart[k1,k2,icum] = b.g[k1,k2,ix1,iy1,iz1,it1]
                                    #a.localpart[k1,k2,icum] = b[k1,k2,ix1,iy1,iz1,it1]
                                    #println(localpart(a)[k1,k2,icum])
                                end
                            end

                            #a.g.localpart[k1,k2,icum] = b.g[k1,k2,ix1,iy1,iz1,it1]
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
        for id in workers()
            push!(f,remotecall(loop_gauge_shift!,id,a,ν,b))
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
                    #localpart(c)[k1,k2,icum] += a[k1,k2,i] 
                    c.g.localpart[k1,k2,icum] += a[k1,k2,i] 
                end
            end
        end

    end

    function Gaugefields.add!(c::DGaugeFields_1d{SU{NC}},a::DGaugeFields_1d{SU{NC}}) where NC
        #for id in workers()
        #    remotecall_fetch(loop_add!,id,c,a)
        #end
        #return

        f = []
        for id in workers()
            #remote_do(loop_add!,id,c,a)
            #push!(remote_do(loop_add!,id,c,a))
            push!(f,remotecall(loop_add!,id,c,a))
        end
        #return 


        for fi in f
            wait(fi)
        end


    end

    function loop_tr(a::DGaugeFields_1d{SU{NC}}) where NC

        i_local = localindices(a)
        icum = 0

        s = 0
        for i = 1:(i_local[3][end]-i_local[3][1]+1)
        #@time for i in i_local[3]
            #icum += 1
            @simd for k=1:NC
                s += a.g.localpart[k,k,i]
                #println(a[k,k,i])
            end
            
            #=
            isinside = true
            isinside *= i in i_local[3]
            if isinside
                @simd for k=1:NC
                    s += localpart(a)[k,k,i-i_local[3][1]+1]
                    #println(a[k,k,i])
                end
            else
                @simd for k=1:NC
                    s += a.g[k,k,i]
                    #println(a[k,k,i])
                end
            end
            =#

        end
        #println(typeof(s))
        #println("s = ",s)
        return s

    end

    function LinearAlgebra.tr(a::DGaugeFields_1d{SU{NC}}) where NC
        s = 0
        for id in workers()
            s += remotecall_fetch(loop_tr,id,a)
        end
        return s

        f = Future[]
        #s = 0
        for id in workers()
            fid = remotecall(loop_tr,id,a)
            #println(typeof(fid))
            push!(f,fid)
        end
        #println("workers, ",workers())

        s = 0
        for fi in f
            wait(fi)
            si = fetch(fi)
            #println("si ",si)
            s += si
            #println(s)
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