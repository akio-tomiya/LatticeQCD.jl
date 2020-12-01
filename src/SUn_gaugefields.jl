module Gaugefields
    using LinearAlgebra
    using Random
    using JLD
    import ..Actions:GaugeActionParam,GaugeActionParam_standard,GaugeActionParam_autogenerator
    import ..Wilsonloops:Wilson_loop_set,calc_coordinate,make_plaq_staple_prime,calc_shift,make_plaq
    

    abstract type SUn end

    abstract type SU2 <: SUn
    end

    abstract type SU3 <: SUn
    end

    struct GaugeFields{T <: SUn} 
        g::Array{ComplexF64,6}
        NX::Int64
        NY::Int64
        NZ::Int64
        NT::Int64
        NC::Int64
        NDW::Int64
        NV::Int64


        function GaugeFields(NC,NDW,NX,NY,NZ,NT)
            if NC == 3
                sutype = SU3
            elseif NC == 2
                sutype = SU2
            else
                error("NC = $NC")
            end 
                
            NV = NX*NY*NZ*NT
            g = zeros(ComplexF64,NC,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW)
            return new{sutype}(g,NX,NY,NZ,NT,NC,NDW,NV)
        end

    end



    const SU2GaugeFields  = GaugeFields{SU2}
    const SU3GaugeFields  = GaugeFields{SU3}







    struct Adjoint_GaugeFields{T <: SUn} 
        parent::GaugeFields{T}
    end

    const Adjoint_SU2GaugeFields  = Adjoint_GaugeFields{SU2}
    const Adjoint_SU3GaugeFields  = Adjoint_GaugeFields{SU3}

    
    function Base.adjoint(x::GaugeFields{T}) where T <: SUn
        Adjoint_GaugeFields{T}(x)
    end

    struct GaugeFields_1d{T <: SUn} 
        g::Array{ComplexF64,3}
        NC::Int64
        NV::Int64
    
        function GaugeFields_1d(NC,NX,NY,NZ,NT)
            if NC == 3
                sutype = SU3
            elseif NC == 2
                sutype = SU2
            end 

            NV = NX*NY*NZ*NT
            g = zeros(ComplexF64,NC,NC,NV)
            return new{sutype}(g,NC,NV)
        end
    
        function GaugeFields_1d(NC,NV) 
            if NC == 3
                sutype = SU3
            elseif NC == 2
                sutype = SU2
            end 

            g = zeros(ComplexF64,NC,NC,NV)
            return new{sutype}(g,NC,NV)
        end
    end

    const SU2GaugeFields_1d  = GaugeFields_1d{SU2}
    const SU3GaugeFields_1d  = GaugeFields_1d{SU3}

    struct Adjoint_GaugeFields_1d{T <: SUn} 
        parent::GaugeFields_1d{T}
    end

    const Adjoint_SU2GaugeFields_1d  = Adjoint_GaugeFields_1d{SU2}
    const Adjoint_SU3GaugeFields_1d  = Adjoint_GaugeFields_1d{SU3}

    function Base.adjoint(x::GaugeFields_1d{T}) where T <: SUn
        Adjoint_GaugeFields_1d{T}(x)
    end

    function Base.setindex!(x::GaugeFields,v,i1,i2,i3,i4,i5,i6) 
        x.g[i1,i2,i3 .+ x.NDW,i4 + x.NDW,i5 + x.NDW,i6 + x.NDW] = v
    end

    function Base.getindex(x::GaugeFields,i1,i2,i3,i4,i5,i6)
        return x.g[i1,i2,i3 .+ x.NDW,i4 .+ x.NDW,i5 .+ x.NDW,i6 .+ x.NDW]
    end



    function Base.setindex!(x::GaugeFields_1d,v,i1,i2,i3) 
        x.g[i1,i2,i3] = v
    end

    function Base.getindex(x::GaugeFields_1d,i1,i2,i3)
        return x.g[i1,i2,i3]
    end


    function evaluate_wilson_loops!(xout::T_1d,w::Wilson_loop_set,U::Array{GaugeFields{S},1},temps::Array{T_1d,1}) where {S <: SUn,T_1d <: GaugeFields_1d}

        num = length(w)
        clear!(xout)
        temp1 = temps[1]
        temp2 = temps[2]
        temp3 = temps[3]



        for i=1:num
            wi = w[i]
            numloops = length(wi)    
            #println("loops ",wi)
            #coordinates = calc_coordinate(wi)
            #println("positions ",coordinates)        
            shifts = calc_shift(wi)
            #println("shift ",shifts)
            
            loopk = wi[1]
            k = 1
            #println("k = $k shift: ",shifts[k])
            gauge_shift_all!(temp1,shifts[1],U[loopk[1]])


            
            loopk1_2 = loopk[2]
            for k=2:numloops
                loopk = wi[k]
                #println("k = $k shift: ",shifts[k])
                #println("gauge_shift!(temp2,$(shifts[k]),$(loopk[1]) )")
                gauge_shift_all!(temp2,shifts[k],U[loopk[1]])

                multiply_12!(temp3,temp1,temp2,k,loopk,loopk1_2)
                temp1,temp3 = temp3,temp1

                
            end
            add!(xout,temp1)
            
        end
    end

    function shift_xyzt(shift,ix,iy,iz,it)
        ix1 = ix + shift[1]
        iy1 = iy + shift[2]
        iz1 = iz + shift[3]
        it1 = it + shift[4]
        return ix1,iy1,iz1,it1
    end

    function multiply_12!(temp3,temp1,temp2,k,loopk,loopk1_2)
        if loopk[2] == 1
            if k==2
                if loopk1_2 == 1
                    mul!(temp3,temp1,temp2)
                else
                    mul!(temp3,temp1',temp2)
                end
            else
                mul!(temp3,temp1,temp2)
            end
        elseif loopk[2] == -1
            if k==2
                if loopk1_2 == 1
                    mul!(temp3,temp1,temp2')
                else
                    mul!(temp3,temp1',temp2')
                end
            else
                mul!(temp3,temp1,temp2')
            end
        else
            error("Second element should be 1 or -1 but now $(loopk)")
        end
        return
    end

    function evaluate_wilson_loops!(V,w::Wilson_loop_set,U::Array{T,1},ix,iy,iz,it) where T <: GaugeFields
        num = length(w)

        V .= 0 #zeros(ComplexF64,NC,NC)
        temp1 = zero(V)
        temp2 = zero(V)
        temp3 = zero(V)


        for i=1:num
            wi = w[i]
            numloops = length(wi)    
            coordinates = calc_coordinate(wi)
            #println("positions ",coordinates)        
            shifts = calc_shift(wi)
            #println("shift ",shifts)
            
            loopk = wi[1]
            ix1,iy1,iz1,it1 = shift_xyzt(shifts[1],ix,iy,iz,it)

            
            temp1[:,:] = U[loopk[1]][:,:,ix1,iy1,iz1,it1]
            loopk1_2 = loopk[2]

            #gauge_shift_all!(temp1,shifts[1],U[loopk[1]])
            for k=2:numloops
                loopk = wi[k]
                
                #gauge_shift_all!(temp2,shifts[k],U[loopk[1]])
                ix1,iy1,iz1,it1 = shift_xyzt(shifts[k],ix,iy,iz,it)
                temp2[:,:] = U[loopk[1]][:,:,ix1,iy1,iz1,it1]

                multiply_12!(temp3,temp1,temp2,k,loopk,loopk1_2)
                #=

                if loopk[2] == 1
                    if k==2
                        if loopk1_2 == 1
                            mul!(temp3,temp1,temp2)
                        else
                            mul!(temp3,temp1',temp2)
                        end
                    else
                        mul!(temp3,temp1,temp2)
                    end
                elseif loopk[2] == -1
                    if k==2
                        if loopk1_2 == 1
                            mul!(temp3,temp1,temp2')
                        else
                            mul!(temp3,temp1',temp2')
                        end
                    else
                        mul!(temp3,temp1,temp2')
                    end
                else
                    error("Second element should be 1 or -1 but now $(loopk)")
                end
                =#

                temp1,temp3 = temp3,temp1
            end
            @. V += temp1
            #add!(xout,temp1)
            
        end
    end





    function set_wing!(u::Array{T,1}) where T <: GaugeFields
        for μ=1:4
            set_wing!(u[μ])
        end
    end

    function set_wing!(u::GaugeFields{T}) where T <: SUn
        NT = u.NT
        NY = u.NY
        NZ = u.NZ
        NX = u.NX
        NDW = u.NDW

        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end
    
        #X direction 
        #Now we send data
    
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for id=1:NDW
                        for k2=1:NC
                            for k1=1:NC
                                u[k1,k2,-NDW+id,iy,iz,it] = u[k1,k2,NX+(id-NDW),iy,iz,it]
                            end
                        end
                    end
                end
            end
        end
    
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for id=1:NDW
                        for k2=1:NC
                            for k1=1:NC
                                u[k1,k2,NX+id,iy,iz,it] = u[k1,k2,id,iy,iz,it]
                            end
                        end
                    end
                end
            end
        end
    
    
        #Y direction 
        #Now we send data
        for it=1:NT
            for iz=1:NZ
                for ix=-NDW+1:NX+NDW
                    for id=1:NDW
                        for k1=1:NC
                            for k2=1:NC
                                u[k1,k2,ix,-NDW+id,iz,it] = u[k1,k2,ix,NY+(id-NDW),iz,it]
                            end
                        end
                    end
                end
            end
        end
    
        for it=1:NT
            for iz=1:NZ
                for ix=-NDW+1:NX+NDW
                    for id=1:NDW
                        for k1=1:NC
                            for k2=1:NC
                                u[k1,k2,ix,NY+id,iz,it] = u[k1,k2,ix,id,iz,it]
                            end
                        end
                    end
                end
            end
        end
    
        #Z direction 
        #Now we send data
        for id=1:NDW
            for it=1:NT
                for iy=-NDW+1:NY+NDW
                    for ix=-NDW+1:NX+NDW
                        for k1=1:NC
                            for k2=1:NC
                                u[k1,k2,ix,iy,id-NDW,it] = u[k1,k2,ix,iy,NZ+(id-NDW),it]
                                u[k1,k2,ix,iy,NZ+id,it] = u[k1,k2,ix,iy,id,it]
                            end
                        end
                    end
                end
            end
        end
    
    
        for id=1:NDW
            for iz=-NDW+1:NZ+NDW
                for iy=-NDW+1:NY+NDW
                    for ix=-NDW+1:NX+NDW
                        for k1=1:NC
                            for k2=1:NC
                                u[k1,k2,ix,iy,iz,id-NDW] = u[k1,k2,ix,iy,iz,NT+(id-NDW)]
                                u[k1,k2,ix,iy,iz,NT+id] = u[k1,k2,ix,iy,iz,id]
                            end
                        end
                    end
                end
            end
        end
    
        #display(u.g)
        #exit()
    
        return
    end
    
    function Base.display(x::GaugeFields{T}) where T <: SUn
        NX=x.NX
        NY=x.NY
        NZ=x.NZ
        NT=x.NT

        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end

        icum = 0
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        icum += 1
                        println("i = $icum")
                        println("ix,iy,iz,it = $ix $iy $iz $it ")
                        display(x[:,:,ix,iy,iz,it])
                        println("\t")
                    end
                end
            end
        end
    end

    function Base.display(a::GaugeFields_1d) 
        for i=1:a.NV
            println("i = $i")
            display(a[:,:,i])
            println("\t")
        end
    end

    function substitute!(a::Array{T,1},b::Array{T,1}) where T <: GaugeFields
        for μ=1:4
            substitute!(a[μ],b[μ])
        end
    end

    function substitute!(a::GaugeFields,b::GaugeFields)
        n1,n2,n3,n4,n5,n6 = size(a.g)
        #println(size(a.g))
        #println(size(b.g))
        for i6=1:n6
            for i5=1:n5
                for i4=1:n4
                    for i3=1:n3
                        for i2=1:n2
                            for i1=1:n1
                                a.g[i1,i2,i3,i4,i5,i6]= b.g[i1,i2,i3,i4,i5,i6]
                            end
                        end
                    end
                end
            end
        end
    end


    function substitute!(a::GaugeFields{T},b::GaugeFields_1d{T}) where T <: SUn
        NT = a.NT
        NY = a.NY
        NZ = a.NZ
        NX = a.NX
        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX               
                        icum = (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix 
                        for k2=1:NC
                            for k1=1:NC
                                a[k1,k2,ix,iy,iz,it] = b[k1,k2,icum]
                            end
                        end

                        #func!(a,b,icum,ix,iy,iz,it)

                    end
                end
            end
        end

    end

    function substitute!(a::GaugeFields_1d{T},b::GaugeFields{T}) where T <: SUn

        NT = b.NT
        NY = b.NY
        NZ = b.NZ
        NX = b.NX

        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX               
                        icum = (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix 


                        for j=1:NC
                            for i=1:NC
                                a[i,j,icum] = b[i,j,ix,iy,iz,it]
                            end
                        end

                    end
                end
            end
        end

    end

    function substitute!(a::GaugeFields_1d,α::Number,b::GaugeFields_1d)
        n1,n2,n3 = size(a.g)

        for i3=1:n3
            for i2=1:n2
                for i1=1:n1
                    a.g[i1,i2,i3]= α*b.g[i1,i2,i3]
                end
            end
        end

    end

    function clear!(a::GaugeFields)
        @. a.g = 0
        return 
    end

    function clear!(a::GaugeFields_1d)
        @. a.g = 0
        return 
    end


    function IdentityGauges(NC,NX,NY,NZ,NT,NDW)
        U = GaugeFields(NC,NDW,NX,NY,NZ,NT)
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for ic=1:NC
                            U[ic,ic,ix,iy,iz,it] = 1 
                        end
                    end
                end
            end
        end
        return U
    end

    
    function RandomGauges(NC,NX,NY,NZ,NT,NDW)
        U = GaugeFields(NC,NDW,NX,NY,NZ,NT)
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for j=1:NC
                            for i=1:NC
                                U[i,j,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
                            end
                        end
                    end
                end
            end
        end
        normalize!(U)
        return U
    end

    function Oneinstanton(NC,NX,NY,NZ,NT,NDW)
        U = Array{SU2GaugeFields,1}(undef,4)
        for μ = 1:4
            U[μ] = GaugeFields(NC,NDW,NX,NY,NZ,NT)
        end
        L = (NX,NY,NZ,NT)
        NV = prod(L)

        R = div(NX,2) # instanton radius

        println("# Starting from a instanton backgorund with radius R=$R ")
        inst_cent = [L[1]/2+0.5, L[2]/2+0.5, L[3]/2+0.5, L[4]/2+0.5]
        #eps = 1/10000000
        s1 = [
        0.0 1.0
        1.0 0.0]
        s2 = [
        0.0 -im*1.0
        im*1.0 0.0]
        s3 = [
        1.0 0.0
        0.0 -1.0]
        En = [
        1.0 0.0
        0.0 1.0]
        ss=[ im*s1,  im*s2,  im*s3, En]
        sd=[-im*s1, -im*s2, -im*s3, En]
        nn = 0
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX  
                        nn += 0
                        nv=[ix-1, iy-1, iz-1, it-1]-inst_cent
                        n2=nv⋅nv
                        for mu = 1:4
                            tau=[
                                0 0
                                0 0]
                            for nu = 1:4
                                smunu = sd[mu]*ss[nu]-sd[nu]*ss[mu]
                                tau+=smunu*nv[nu]
                            end
                            sq = sqrt(n2+R^2)
                            tau = exp(im*tau*(1/2)*(1/(n2))*(im*R^2/(n2+R^2)) ) #1b
                            for j=1:2
                                for i=1:2
                                    U[mu][i,j,ix,iy,iz,it] = tau[i,j]
                                end
                            end
                        
                        end
                    end
                end
            end
        end

        #normalize!(U)
        return U

    end

        



    function Base.similar(x::GaugeFields) 
        return GaugeFields(x.NC,x.NDW,x.NX,x.NY,x.NZ,x.NT)
    end

    function Base.similar(x::Array{T,1}) where T <: GaugeFields
        xout = Array{T,1}(undef,4)
        for μ=1:4
            xout[μ] = similar(x[μ])
        end
        return xout
    end

    function Base.copyto!(Unew::GaugeFields,U::GaugeFields)
        Unew = deepcopy(U)
        return
    end

    function elementwise_tr!(s,u::GaugeFields_1d{T},v::GaugeFields_1d{T}) where T <: SUn
        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end

        NV=u.NV

        for i=1:NV
            for k1=1:NC
                for k2=1:NC
                    s[i] += real(u[k1,k2,i]*conj(v[k1,k2,i]))
                end
            end
        end
        return

    end


    function LinearAlgebra.tr(a::GaugeFields{T}) where T <: SUn
        NX=a.NX
        NY=a.NY
        NZ=a.NZ
        NT=a.NT
        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end

        s = 0
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k=1:NC
                            s += a[k,k,ix,iy,iz,it]
                        end
                    end
                end
            end
        end
        return s

    end

    function LinearAlgebra.tr(a::GaugeFields_1d{T}) where T <: SUn
        NV=a.NV
        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end
        s = 0

        for i=1:NV

            for k=1:NC
                s += a[k,k,i]
            end

        end
        return s

    end

    function elementwise_apply!(a::GaugeFields{T},func!::Function) where T <: SUn
        NT = a.NT
        NY = a.NY
        NZ = a.NZ
        NX = a.NX
        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX   
                        for j=1:NC
                            for i=1:NC
                                func!(a,i,j,ix,iy,iz,it,b...)
                            end
                        end
                    end
                end
            end
        end
        return
    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::GaugeFields_1d{T},b::GaugeFields_1d{T}) where T <: SUn
        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end

        NV=a.NV
        


        for i=1:NV
            #mulabc!(a,b,c,i)
            
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] = 0
                    for k3=1:NC
                        c[k1,k2,i] += a[k1,k3,i]*b[k3,k2,i]
                    end
                end
            end
            
        end

    end



    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::GaugeFields_1d{T},b::GaugeFields_1d{T}) where T <: SU3
        NV=a.NV
        #NC=a.NC
        #mulabc! = NCmul(NC)

        for i=1:NV
            c[1,1,i] = a[1,1,i] * b[1,1,i] +
                a[1,2,i] * b[2,1,i] +
                a[1,3,i] * b[3,1,i]
            c[1,2,i] = a[1,1,i] * b[1,2,i] + 
                    a[1,2,i] * b[2,2,i] +
                    a[1,3,i] * b[3,2,i]
            c[1,3,i] = a[1,1,i] * b[1,3,i] +
                    a[1,2,i] * b[2,3,i] +
                    a[1,3,i] * b[3,3,i]

            c[2,1,i] = a[2,1,i] * b[1,1,i] +
                    a[2,2,i] * b[2,1,i] +
                    a[2,3,i] * b[3,1,i]
            c[2,2,i] = a[2,1,i] * b[1,2,i] + 
                    a[2,2,i] * b[2,2,i] +
                    a[2,3,i] * b[3,2,i]
            c[2,3,i] = a[2,1,i] * b[1,3,i] + 
                    a[2,2,i] * b[2,3,i] +
                    a[2,3,i] * b[3,3,i]

            c[3,1,i] = a[3,1,i] * b[1,1,i] + 
                    a[3,2,i] * b[2,1,i] +
                    a[3,3,i] * b[3,1,i]
            c[3,2,i] = a[3,1,i] * b[1,2,i] + 
                    a[3,2,i] * b[2,2,i] +
                    a[3,3,i] * b[3,2,i]
            c[3,3,i] = a[3,1,i] * b[1,3,i] + 
                    a[3,2,i] * b[2,3,i] + 
                    a[3,3,i] * b[3,3,i]
        end

    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::GaugeFields_1d{T},b::GaugeFields_1d{T}) where T <: SU2
        NV=a.NV
        #NC=a.NC
        #mulabc! = NCmul(NC)

        for i=1:NV
            c[1,1,i] = a[1,1,i] * b[1,1,i] +
                a[1,2,i] * b[2,1,i] 

            c[1,2,i] = a[1,1,i] * b[1,2,i] + 
                    a[1,2,i] * b[2,2,i] 


            c[2,1,i] = a[2,1,i] * b[1,1,i] +
                    a[2,2,i] * b[2,1,i] 
            c[2,2,i] = a[2,1,i] * b[1,2,i] + 
                    a[2,2,i] * b[2,2,i] 

        end

    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::GaugeFields_1d{T},b::Adjoint_GaugeFields_1d{T}) where T <: SUn
        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end
        NV=a.NV


        for i=1:NV
            #mulabc!(a,b,c,i)
            
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] = 0
                    for k3=1:NC
                        c[k1,k2,i] += a[k1,k3,i]*conj(b.parent[k2,k3,i])
                    end
                end
            end
            
                        
        end
        

    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::Number) where T <: SUn
        NV=c.NV
        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end

        #NC=a.NC
        #mulabc! = NCmul_aconjb(NC)

        for i=1:NV
            for k2= 1:NC
                for k1=1:NC
                    c[k1,k2,i] *= a
                end
            end                        
        end
        return
    
    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::GaugeFields_1d{T},b::Adjoint_GaugeFields_1d{T}) where T <: SU3
        
        NV=a.NV
        #NC=a.NC
        #mulabc! = NCmul_aconjb(NC)

        for i=1:NV
            c[1,1,i] = a[1,1,i] * conj(b.parent[1,1,i]) +
                a[1,2,i] * conj(b.parent[1,2,i]) +
                a[1,3,i] * conj(b.parent[1,3,i])
            c[1,2,i] = a[1,1,i] * conj(b.parent[2,1,i]) + 
                    a[1,2,i] * conj(b.parent[2,2,i]) +
                    a[1,3,i] * conj(b.parent[2,3,i])
            c[1,3,i] = a[1,1,i] * conj(b.parent[3,1,i])+
                    a[1,2,i] * conj(b.parent[3,2,i]) +
                    a[1,3,i] * conj(b.parent[3,3,i])

            c[2,1,i] = a[2,1,i] * conj(b.parent[1,1,i]) +
                    a[2,2,i] * conj(b.parent[1,2,i]) +
                    a[2,3,i] * conj(b.parent[1,3,i])
            c[2,2,i] = a[2,1,i] * conj(b.parent[2,1,i]) + 
                    a[2,2,i] * conj(b.parent[2,2,i]) +
                    a[2,3,i] * conj(b.parent[2,3,i])
            c[2,3,i] = a[2,1,i] * conj(b.parent[3,1,i]) + 
                    a[2,2,i] * conj(b.parent[3,2,i]) +
                    a[2,3,i] * conj(b.parent[3,3,i])

            c[3,1,i] = a[3,1,i] * conj(b.parent[1,1,i]) + 
                    a[3,2,i] * conj(b.parent[1,2,i]) +
                    a[3,3,i] * conj(b.parent[1,3,i])
            c[3,2,i] = a[3,1,i] * conj(b.parent[2,1,i]) + 
                    a[3,2,i] * conj(b.parent[2,2,i]) +
                    a[3,3,i] * conj(b.parent[2,3,i])
            c[3,3,i] = a[3,1,i] * conj(b.parent[3,1,i]) + 
                    a[3,2,i] * conj(b.parent[3,2,i]) + 
                    a[3,3,i] * conj(b.parent[3,3,i])

                        
        end
    
    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::GaugeFields_1d{T},b::Adjoint_GaugeFields_1d{T}) where T <: SU2
        
        NV=a.NV
        #NC=a.NC
        #mulabc! = NCmul_aconjb(NC)

        for i=1:NV
            c[1,1,i] = a[1,1,i] * conj(b.parent[1,1,i]) +
                a[1,2,i] * conj(b.parent[1,2,i]) 
            c[1,2,i] = a[1,1,i] * conj(b.parent[2,1,i]) + 
                    a[1,2,i] * conj(b.parent[2,2,i]) 

            c[2,1,i] = a[2,1,i] * conj(b.parent[1,1,i]) +
                    a[2,2,i] * conj(b.parent[1,2,i]) 
            c[2,2,i] = a[2,1,i] * conj(b.parent[2,1,i]) + 
                    a[2,2,i] * conj(b.parent[2,2,i]) 

                        
        end
    
    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::Adjoint_GaugeFields_1d{T},b::Adjoint_GaugeFields_1d{T}) where T <: SU2
        
        NV=a.parent.NV
        #NC=a.NC
        #mulabc! = NCmul_aconjb(NC)

        for i=1:NV
            c[1,1,i] = conj(a.parent[1,1,i]) * conj(b.parent[1,1,i]) +
                conj(a.parent[2,1,i] )* conj(b.parent[1,2,i]) 
            c[1,2,i] = conj(a.parent[1,1,i]) * conj(b.parent[2,1,i]) + 
                conj(a.parent[2,1,i]) * conj(b.parent[2,2,i]) 

            c[2,1,i] = conj(a.parent[1,2,i]) * conj(b.parent[1,1,i]) +
                conj(a.parent[2,2,i]) * conj(b.parent[1,2,i]) 
            c[2,2,i] = conj(a.parent[1,2,i]) * conj(b.parent[2,1,i]) + 
                conj(a.parent[2,2,i]) * conj(b.parent[2,2,i])        
        end
    
    end


    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::Adjoint_GaugeFields_1d{T},b::GaugeFields_1d{T}) where T <: SUn
        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end

        NV=c.NV
        #NC=c.NC


        for i=1:NV
            #mulabc!(a,b,c,i)
            
    
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] = 0
                    for k3=1:NC
                        c[k1,k2,i] += conj(a.parent[k3,k1,i])*b[k3,k2,i]
                    end
                end
            end
            
                        
        end

    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::Adjoint_GaugeFields_1d{T},b::GaugeFields_1d{T}) where T <: SU3
        NV=c.NV
        #NC=c.NC

        #mulabc! = NCmul_conjab(c.NC)

        for i=1:NV

            c[1,1,i] = conj(a.parent[1,1,i]) * b[1,1,i] +
                conj(a.parent[2,1,i]) * b[2,1,i] +
                conj(a.parent[3,1,i]) * b[3,1,i]
            c[1,2,i] = conj(a.parent[1,1,i]) * b[1,2,i] + 
                    conj(a.parent[2,1,i]) * b[2,2,i] +
                    conj(a.parent[3,1,i]) * b[3,2,i]
            c[1,3,i] = conj(a.parent[1,1,i]) * b[1,3,i] +
                    conj(a.parent[2,1,i]) * b[2,3,i] +
                    conj(a.parent[3,1,i]) * b[3,3,i]

            c[2,1,i] = conj(a.parent[1,2,i]) * b[1,1,i] +
                    conj(a.parent[2,2,i]) * b[2,1,i] +
                    conj(a.parent[3,2,i]) * b[3,1,i]
            c[2,2,i] = conj(a.parent[1,2,i]) * b[1,2,i] + 
                    conj(a.parent[2,2,i]) * b[2,2,i] +
                    conj(a.parent[3,2,i]) * b[3,2,i]
            c[2,3,i] = conj(a.parent[1,2,i]) * b[1,3,i] + 
                    conj(a.parent[2,2,i]) * b[2,3,i] +
                    conj(a.parent[3,2,i]) * b[3,3,i]

            c[3,1,i] = conj(a.parent[1,3,i]) * b[1,1,i] + 
                    conj(a.parent[2,3,i]) * b[2,1,i] +
                    conj(a.parent[3,3,i]) * b[3,1,i]
            c[3,2,i] = conj(a.parent[1,3,i]) * b[1,2,i] + 
                    conj(a.parent[2,3,i]) * b[2,2,i] +
                    conj(a.parent[3,3,i]) * b[3,2,i]
            c[3,3,i] = conj(a.parent[1,3,i]) * b[1,3,i] + 
                    conj(a.parent[2,3,i]) * b[2,3,i] + 
                    conj(a.parent[3,3,i]) * b[3,3,i]
            #=
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] = 0
                    for k3=1:NC
                        c[k1,k2,i] += conj(a.parent[k3,k1,i])*b[k3,k2,i]
                    end
                end
            end
            =#
                        
        end

    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::Adjoint_GaugeFields_1d{T},b::GaugeFields_1d{T}) where T <: SU2
        NV=c.NV
        #NC=c.NC

        #mulabc! = NCmul_conjab(c.NC)

        for i=1:NV

            c[1,1,i] = conj(a.parent[1,1,i]) * b[1,1,i] +
                conj(a.parent[2,1,i]) * b[2,1,i] 
            c[1,2,i] = conj(a.parent[1,1,i]) * b[1,2,i] + 
                    conj(a.parent[2,1,i]) * b[2,2,i] 


            c[2,1,i] = conj(a.parent[1,2,i]) * b[1,1,i] +
                    conj(a.parent[2,2,i]) * b[2,1,i] 
            c[2,2,i] = conj(a.parent[1,2,i]) * b[1,2,i] + 
                    conj(a.parent[2,2,i]) * b[2,2,i] 


                        
        end

    end

    function add!(c::GaugeFields_1d{T},a::GaugeFields_1d{T}) where T <: SUn
        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end
        NV=c.NV

        for i=1:NV
            #ncadd!(a,c,i)
            
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] += a[k1,k2,i] 
                end
            end
            
        end


    end

    function staggered_phase(μ,ix,iy,iz,it,NX,NY,NZ,NT)
        t = it-1
        t += ifelse(t<0,NT,0)
        t += ifelse(t ≥ NT,-NT,0)
        #boundary_factor_t = ifelse(t == NT -1,BoundaryCondition[4],1)
        z = iz-1
        z += ifelse(z<0,NZ,0)
        z += ifelse(z ≥ NZ,-NZ,0)
        #boundary_factor_z = ifelse(z == NZ -1,BoundaryCondition[3],1)
        y = iy-1
        y += ifelse(y<0,NY,0)
        y += ifelse(y ≥ NY,-NY,0)
        #boundary_factor_y = ifelse(y == NY -1,BoundaryCondition[2],1)
        x = ix-1
        x += ifelse(x<0,NX,0)
        x += ifelse(x ≥ NX,-NX,0)
        #boundary_factor_x = ifelse(x == NX -1,BoundaryCondition[1],1)
        if μ ==1
            η = 1
        elseif μ ==2
            #η = (-1.0)^(x)
            η = ifelse(x%2 == 0,1,-1)
        elseif μ ==3
            #η = (-1.0)^(x+y)
            η = ifelse((x+y)%2 == 0,1,-1)
        elseif μ ==4
            #η = (-1.0)^(x+y+z)
            η = ifelse((x+y+z)%2 == 0,1,-1)
        else
            error("η should be positive but η = $η")
        end
        return η#*boundary_factor_x*boundary_factor_y*boundary_factor_z*boundary_factor_t
    end

    function gauge_shift!(a::GaugeFields_1d{T},ν::N,b::GaugeFields{T}) where {N <: Int,T <: SUn}
        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end

        if ν == 0
            substitute!(a,b)
            return
        end

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
        #NC = b.NC


        
        for it=1:NT
            it1 = it + idel[4]
            for iz=1:NZ
                iz1 = iz + idel[3]
                for iy=1:NY
                    iy1 = iy+idel[2]
                    for ix=1:NX
                        ix1 = ix+idel[1]
                        icum = (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix
                        #func!(a,b,icum,ix1,iy1,iz1,it1)

                        
                        
                        for k2=1:NC
                            for k1=1:NC
                                a[k1,k2,icum] = b[k1,k2,ix1,iy1,iz1,it1]
                            end
                        end
                        
                        
                        
                    end
                end
            end
        end
        return
    end

    function gauge_shift!(a::GaugeFields_1d{T},dir::Tuple,b::GaugeFields{T}) where T <: SUn 
        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            error("NC >3 is not supported")
        end

        #....  Stop for the exceptional case   ...
        if dir[1] != 0 && dir[2] != 0
            if dir[1] == dir[2]
                println("Sorry this case is not yet considered!")
                println("dir : ",dir[1],"\t",dir[2])
                exit()
            end
        end

        if dir[1] == 0 && dir[2] == 0
            subst!(ibush,b,a)
            return
        end

        idel1 = zeros(Int64,4)
        idel2 = zeros(Int64,4)
        if dir[1] > 0
            idel1[+dir[1]] = +1
        elseif dir[1] < 0
            idel1[-dir[1]] = -1
        end

        if dir[2] > 0
            idel2[+dir[2]] = +1
        elseif dir[2] < 0
            idel2[-dir[2]] = -1
        end


        
        NT = b.NT
        NZ = b.NZ
        NY = b.NY
        NX = b.NX


        #func! = NCsubstitute_63(b.NC)
        #        NC = b.NC
        

        
        for it=1:NT
            it1 = it + idel1[4] + idel2[4]
            for iz=1:NZ
                iz1 = iz + idel1[3]+idel2[3]
                for iy=1:NY
                    iy1 = iy+idel1[2]+idel2[2]
                    for ix=1:NX
                        ix1 = ix+idel1[1] + idel2[1]
                        icum = (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix
                        #func!(a,b,icum,ix1,iy1,iz1,it1)

                        
                        for k2=1:NC
                            for k1=1:NC
                                a[k1,k2,icum] = b[k1,k2,ix1,iy1,iz1,it1]
                            end
                        end
                        
                    end
                end
            end
        end
        return
    end

    function gauge_shift_all!(a::GaugeFields_1d{T},dir,b::GaugeFields{T}) where T <: SUn 

        indx = Int64[]
        for mu = 1:4
            if dir[mu] != 0
                push!(indx,mu*sign(dir[mu]))
            end
        end
        #println(indx,"\t")

        if length(indx) == 1
            gauge_shift!(a,indx[1],b) 
        elseif length(indx) == 2
            gauge_shift!(a,Tuple(indx),b) 
        elseif length(indx) == 0
            substitute!(a,b)
        end
        return
    end

    function make_staple!(staple::GaugeFields_1d,U,μ,temp1,temp2,temp3)
        clear!(staple)

        idir2 = zeros(Int64,2)

        for ν=1:4
            if ν == μ
                continue
            end

                        #=
            c       x+nu temp2
            c        .---------.
            c        I         I
            c  temp1 I         I
            c        I         I
            c        .         .
            c        x        x+mu
            =#
            substitute!(temp1,U[ν])
            gauge_shift!(temp2,ν,U[μ])
            #mul!(temp3,U[ν],temp2)
            mul!(temp3,temp1,temp2)

            gauge_shift!(temp1,μ,U[ν])
            mul!(temp2,temp3,temp1')
            add!(staple,temp2)



        end
    end



    function make_staple_double!(staple::GaugeFields_1d,U,μ,temp1,temp2,temp3)
        clear!(staple)
        #println("mu = ",μ)

        #loops = make_plaq_staple(μ)
        loops = make_plaq_staple_prime(μ)
        evaluate_wilson_loops!(staple,loops,U,[temp1,temp2,temp3])
        return


        clear!(staple)
        
        for ν=1:4
            if ν == μ
                continue
            end

                        #=
            c       x+nu temp2
            c        .---------.
            c        I         I
            c  temp1 I         I
            c        I         I
            c        .         .
            c        x        x+mu
            =#
            substitute!(temp1,U[ν])
            #println("[]")
            gauge_shift!(temp2,ν,U[μ])

            #println(ν)
            
            mul!(temp3,temp1,temp2)

            gauge_shift!(temp1,μ,U[ν])

            #println(μ)
            mul!(temp2,temp3,temp1')



            add!(staple,temp2)

            #=
            c        x
            c        .         .
            c        I         I
            c  temp1 I         I
            c        I         I
            c        .---------.
            c      x-nu        x+mu-nu
            =#
            gauge_shift!(temp1,-ν,U[ν])
            #println(-ν)
            gauge_shift!(temp2,-ν,U[μ])
            #println(-ν)

            mul!(temp3,temp1',temp2)




            #idir2 = Ivec2(μ,-ν)
            
            #println((μ,-ν))
            gauge_shift!(temp1,(μ,-ν),U[ν])

            mul!(temp2,temp3,temp1)

            add!(staple,temp2)
            #println("mu")

            #if ν == 4
            #    display(staple)
            #    exit()

            #end



        end
        
    end

    

    
function m3complv!(a::SU3GaugeFields)
    aa = zeros(Float64,18)
    NX = u.NX
    NY = u.NY
    NZ = u.NZ
    NT = u.NT

    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX

                    aa[ 1] = real( a[1,1,ix,iy,iz,it])
                    aa[ 2] = imag(a[1,1,ix,iy,iz,it])
                    aa[ 3] = real( a[1,2,ix,iy,iz,it])
                    aa[ 4] = imag(a[1,2,ix,iy,iz,it])
                    aa[ 5] = real( a[1,3,ix,iy,iz,it])
                    aa[ 6] = imag(a[1,3,ix,iy,iz,it])
                    aa[ 7] = real( a[2,1,ix,iy,iz,it])
                    aa[ 8] = imag(a[2,1,ix,iy,iz,it])
                    aa[ 9] = real( a[2,2,ix,iy,iz,it])
                    aa[10] = imag(a[2,2,ix,iy,iz,it])
                    aa[11] = real( a[2,3,ix,iy,iz,it])
                    aa[12] = imag(a[2,3,ix,iy,iz,it])

                    aa[13] = aa[ 3]*aa[11] - aa[ 4]*aa[12] -
                                aa[ 5]*aa[ 9] + aa[ 6]*aa[10]
                    aa[14] = aa[ 5]*aa[10] + aa[ 6]*aa[ 9] -
                                aa[ 3]*aa[12] - aa[ 4]*aa[11]
                    aa[15] = aa[ 5]*aa[ 7] - aa[ 6]*aa[ 8] -
                                aa[ 1]*aa[11] + aa[ 2]*aa[12]
                    aa[16] = aa[ 1]*aa[12] + aa[ 2]*aa[11] -
                                aa[ 5]*aa[ 8] - aa[ 6]*aa[ 7]
                    aa[17] = aa[ 1]*aa[ 9] - aa[ 2]*aa[10] -
                                aa[ 3]*aa[ 7] + aa[ 4]*aa[ 8]
                    aa[18] = aa[ 3]*aa[ 8] + aa[ 4]*aa[ 7] -
                                aa[ 1]*aa[10] - aa[ 2]*aa[ 9]

                    a[3,1,ix,iy,iz,it] = aa[13]+im*aa[14]
                    a[3,2,ix,iy,iz,it] = aa[15]+im*aa[16]
                    a[3,3,ix,iy,iz,it] = aa[17]+im*aa[18]
                end
            end
        end
    end
end

    function TA(vin::T) where T <: GaugeFields_1d
        vout = deepcopy(vin)
        projlink!(vout,vin)
        return vout
    end

"""
    c-----------------------------------------------------c
c     !!!!!   vin and vout should be different vectors
c
c     Projectin of the etraceless antiermite part 
c     vout = x/2 - Tr(x)/6
c     wher   x = vin - Conjg(vin)      
c-----------------------------------------------------c
    """
    function projlink!(vout::SU3GaugeFields,vin::SU3GaugeFields)
        
        fac13 = 1/3
        NX = vin.NX
        NY = vin.NY
        NZ = vin.NZ
        NT = vin.NT

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        v11 = vin[1,1,ix,iy,iz,it]
                        v22 = vin[2,2,ix,iy,iz,it]
                        v33 = vin[3,3,ix,iy,iz,it]

                        tri = fac13*(imag(v11)+imag(v22)+imag(v33))

                        vout[1,1,ix,iy,iz,it] = (imag(v11)-tri)*im
                        vout[2,2,ix,iy,iz,it] = (imag(v22)-tri)*im
                        vout[3,3,ix,iy,iz,it] = (imag(v33)-tri)*im
                    end
                end
            end
        end

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX

                        v12 = vin[1,2,ix,iy,iz,it]
                        v13 = vin[1,3,ix,iy,iz,it]
                        v21 = vin[2,1,ix,iy,iz,it]
                        v23 = vin[2,3,ix,iy,iz,it]
                        v31 = vin[3,1,ix,iy,iz,it]
                        v32 = vin[3,2,ix,iy,iz,it]

                        x12 = v12 - conj(v21)
                        x13 = v13 - conj(v31)
                        x23 = v23 - conj(v32)
                    
                        x21 = - conj(x12)
                        x31 = - conj(x13)
                        x32 = - conj(x23)

                        vout[1,2,ix,iy,iz,it] = 0.5  * x12
                        vout[1,3,ix,iy,iz,it] = 0.5  * x13
                        vout[2,1,ix,iy,iz,it] = 0.5  * x21
                        vout[2,3,ix,iy,iz,it] = 0.5  * x23
                        vout[3,1,ix,iy,iz,it] = 0.5  * x31
                        vout[3,2,ix,iy,iz,it] = 0.5  * x32
                    end
                end
            end
        end

    end

    const sr3 = sqrt(3)
    const sr3i = 1/sr3

    function projlink!(vout::SU3GaugeFields_1d,vin::SU3GaugeFields_1d)
        fac13 = 1/3
        nv = vin.NV

        for i=1:nv
            v11 = vin[1,1,i]
            v22 = vin[2,2,i]
            v33 = vin[3,3,i]

            tri = fac13*(imag(v11)+imag(v22)+imag(v33))

            vout[1,1,i] = (imag(v11)-tri)*im
            vout[2,2,i] = (imag(v22)-tri)*im
            vout[3,3,i] = (imag(v33)-tri)*im

        end

        for i=1:nv
            v12 = vin[1,2,i]
            v13 = vin[1,3,i]
            v21 = vin[2,1,i]
            v23 = vin[2,3,i]
            v31 = vin[3,1,i]
            v32 = vin[3,2,i]

            x12 = v12 - conj(v21)
            x13 = v13 - conj(v31)
            x23 = v23 - conj(v32)
        
            x21 = - conj(x12)
            x31 = - conj(x13)
            x32 = - conj(x23)

            vout[1,2,i] = 0.5  * x12
            vout[1,3,i] = 0.5  * x13
            vout[2,1,i] = 0.5  * x21
            vout[2,3,i] = 0.5  * x23
            vout[3,1,i] = 0.5  * x31
            vout[3,2,i] = 0.5  * x32
            
        end
    end

    """
    c-----------------------------------------------------c
    c     !!!!!   vin and vout should be different vectors
    c
    c     Projectin of the etraceless antiermite part 
    c     vout = x/2 - Tr(x)/6
    c     wher   x = vin - Conjg(vin)      
    c-----------------------------------------------------c
    """
    function projlink!(vout::SU2GaugeFields,vin::SU2GaugeFields)
        
        fac12 = 1/2
        NX = vin.NX
        NY = vin.NY
        NZ = vin.NZ
        NT = vin.NT

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        v11 = vin[1,1,ix,iy,iz,it]
                        v22 = vin[2,2,ix,iy,iz,it]

                        tri = fac12*(imag(v11)+imag(v22))

                        vout[1,1,ix,iy,iz,it] = (imag(v11)-tri)*im
                        vout[2,2,ix,iy,iz,it] = (imag(v22)-tri)*im
                    end
                end
            end
        end

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX

                        v12 = vin[1,2,ix,iy,iz,it]
                        #v13 = vin[1,3,ix,iy,iz,it]
                        v21 = vin[2,1,ix,iy,iz,it]

                        x12 = v12 - conj(v21)

                        x21 = - conj(x12)

                        vout[1,2,ix,iy,iz,it] = 0.5  * x12

                        vout[2,1,ix,iy,iz,it] = 0.5  * x21


                    end
                end
            end
        end

    end



    function projlink!(vout::SU2GaugeFields_1d,vin::SU2GaugeFields_1d)
        fac12 = 1/2
        nv = vin.NV

        for i=1:nv
            v11 = vin[1,1,i]
            v22 = vin[2,2,i]

            tri = fac12*(imag(v11)+imag(v22))

            vout[1,1,i] = (imag(v11)-tri)*im
            vout[2,2,i] = (imag(v22)-tri)*im

        end

        for i=1:nv
            v12 = vin[1,2,i]
            v21 = vin[2,1,i]


            x12 = v12 - conj(v21)

        
            x21 = - conj(x12)

            vout[1,2,i] = 0.5  * x12
            vout[2,1,i] = 0.5  * x21

            
        end

    end

    const sr3ih = 0.5*sr3i
    const sqr3inv = sr3i


    function lambdamul(b::SU3GaugeFields_1d,a::SU3GaugeFields_1d,k)
        #=
        c----------------------------------------------------------------------c
        c     b = (lambda_k/2)*a
        C             lambda_k : GellMann matrices. k=1, 8 
        c----------------------------------------------------------------------c
        =#
        NV = a.NV


        if k==1
            for i=1:NV
                b[1,1,i] = 0.5 * a[2,1,i] 
                b[1,2,i] = 0.5 * a[2,2,i]
                b[1,3,i] = 0.5 * a[2,3,i]
                b[2,1,i] = 0.5 * a[1,1,i]
                b[2,2,i] = 0.5 * a[1,2,i]
                b[2,3,i] = 0.5 * a[1,3,i]
                b[3,1,i] = 0
                b[3,2,i] = 0
                b[3,3,i] = 0
            end
        elseif k==2
            for i=1:NV
                b[1,1,i] = -0.5*im * a[2,1,i] 
                b[1,2,i] = -0.5*im * a[2,2,i]
                b[1,3,i] = -0.5*im * a[2,3,i]
                b[2,1,i] =  0.5*im * a[1,1,i]
                b[2,2,i] =  0.5*im * a[1,2,i]
                b[2,3,i] =  0.5*im * a[1,3,i]
                b[3,1,i] = 0
                b[3,2,i] = 0
                b[3,3,i] = 0
            end
        elseif k==3
            for i=1:NV
                b[1,1,i] =  0.5 * a[1,1,i] 
                b[1,2,i] =  0.5 * a[1,2,i]
                b[1,3,i] =  0.5 * a[1,3,i]
                b[2,1,i] = -0.5 * a[2,1,i]
                b[2,2,i] = -0.5 * a[2,2,i]
                b[2,3,i] = -0.5 * a[2,3,i]
                b[3,1,i] = 0
                b[3,2,i] = 0
                b[3,3,i] = 0
            end
        elseif k==4
            for i=1:NV
                b[1,1,i] = 0.5 * a[3,1,i] 
                b[1,2,i] = 0.5 * a[3,2,i]
                b[1,3,i] = 0.5 * a[3,3,i]
                b[2,1,i] = 0
                b[2,2,i] = 0
                b[2,3,i] = 0
                b[3,1,i] = 0.5 * a[1,1,i]
                b[3,2,i] = 0.5 * a[1,2,i]
                b[3,3,i] = 0.5 * a[1,3,i]
            end
        elseif k==5
            for i=1:NV
                b[1,1,i] = -0.5*im * a[3,1,i] 
                b[1,2,i] = -0.5*im * a[3,2,i]
                b[1,3,i] = -0.5*im * a[3,3,i]
                b[2,1,i] = 0
                b[2,2,i] = 0
                b[2,3,i] = 0
                b[3,1,i] =  0.5*im * a[1,1,i]
                b[3,2,i] =  0.5*im * a[1,2,i]
                b[3,3,i] =  0.5*im * a[1,3,i]
            end
        elseif k==6
            for i=1:NV
                b[1,1,i] = 0
                b[1,2,i] = 0
                b[1,3,i] = 0
                b[2,1,i] = 0.5 * a[3,1,i] 
                b[2,2,i] = 0.5 * a[3,2,i]
                b[2,3,i] = 0.5 * a[3,3,i]
                b[3,1,i] = 0.5 * a[2,1,i]
                b[3,2,i] = 0.5 * a[2,2,i]
                b[3,3,i] = 0.5 * a[2,3,i]
            end
        elseif k==7
            for i=1:NV
                b[1,1,i] = 0
                b[1,2,i] = 0
                b[1,3,i] = 0
                b[2,1,i] = -0.5*im * a[3,1,i] 
                b[2,2,i] = -0.5*im * a[3,2,i]
                b[2,3,i] = -0.5*im * a[3,3,i]
                b[3,1,i] =  0.5*im * a[2,1,i]
                b[3,2,i] =  0.5*im * a[2,2,i]
                b[3,3,i] =  0.5*im * a[2,3,i]
            end
        elseif k==8
            for i=1:NV
                b[1,1,i] =  sr3ih * a[1,1,i] 
                b[1,2,i] =  sr3ih * a[1,2,i]
                b[1,3,i] =  sr3ih * a[1,3,i]
                b[2,1,i] =  sr3ih * a[2,1,i] 
                b[2,2,i] =  sr3ih * a[2,2,i]
                b[2,3,i] =  sr3ih * a[2,3,i]
                b[3,1,i] = -sqr3inv * a[3,1,i]
                b[3,2,i] = -sqr3inv * a[3,2,i]
                b[3,3,i] = -sqr3inv * a[3,3,i]
            end
        else
            error("k should be k <= 8 but k = $k")
        end

        return
    end

    function lambdamul(b::SU2GaugeFields_1d,a::SU2GaugeFields_1d,k)
        #=
        c----------------------------------------------------------------------c
        c     b = (lambda_k/2)*a
        C             lambda_k : GellMann matrices. k=1, 8 
        c----------------------------------------------------------------------c
        =#
        NV = a.NV
    
    
        if k==1
            for i=1:NV
                b[1,1,i] = -0.5*im* a[2,1,i]*im
                b[1,2,i] = -0.5*im * a[2,2,i]*im
    
                b[2,1,i] = -0.5*im * a[1,1,i]*im
                b[2,2,i] = -0.5*im * a[1,2,i]*im
    
    
            end
        elseif k==2
            for i=1:NV
                b[1,1,i] = -0.5 * a[2,1,i] *im
                b[1,2,i] = -0.5 * a[2,2,i]*im
    
                b[2,1,i] =  0.5 * a[1,1,i]*im
                b[2,2,i] =  0.5 * a[1,2,i]*im
    
            end
        elseif k==3
            for i=1:NV
                b[1,1,i] =  -0.5*im * a[1,1,i] *im
                b[1,2,i] =  -0.5*im * a[1,2,i]*im
    
                b[2,1,i] = 0.5*im * a[2,1,i]*im
                b[2,2,i] = 0.5*im * a[2,2,i]*im
    
            end
        else
            error("k should be k <= 3 but k = $k")
        end
    
        return
    end

        





    function calc_Polyakov(u::Array{T,1}) where T <: GaugeFields
        NX = u[1].NX
        NY = u[1].NY
        NZ = u[1].NZ
        NT = u[1].NT
        NC = u[1].NC

        Pol = zeros(ComplexF64,NC,NC,NX,NY,NZ)
        
        tmp1= zeros(ComplexF64,NC,NC)
        tmp2= zero(tmp1)

        set_wing!(u)

        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    for k2=1:NC
                        for k1=1:NC
                            Pol[k1,k2,ix,iy,iz] = u[4][k1,k2,ix,iy,iz,1]
                        end
                    end
                end
            end
        end

        for it=2:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k2=1:NC
                            for k1=1:NC
                                tmp1[k1,k2] = Pol[k1,k2,ix,iy,iz]
                                tmp2[k1,k2] = u[4][k1,k2,ix,iy,iz,it]
                            end
                        end
                        Pol[:,:,ix,iy,iz] = tmp1*tmp2
                    end
                end
            end
        end

        avePol = 0.0
        for iz=1:NX
            for iy=1:NY
                for ix=1:NX
                    for k1=1:NC
                        avePol += Pol[k1,k1,ix,iy,iz]
                    end
                end
            end
        end
        avePol /= NX*NY*NZ


        return avePol
    end

    function calc_Plaq_notrace_1d(U::Array{T,1},μ,ν,origin) where T <: GaugeFields

        NC = U[1].NC
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT

        temp1 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp2 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp3 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        plaq = GaugeFields_1d(NC,NX,NY,NZ,NT)

        loop = make_plaq(μ,ν,origin)

        evaluate_wilson_loops!(plaq,loop ,U,[temp1,temp2,temp3])

        return plaq
    end

    function calc_Plaq_notrace_1d(U::Array{T,1},μ,ν) where T <: GaugeFields

        NC = U[1].NC
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT

        temp1 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp2 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp3 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        plaq = GaugeFields_1d(NC,NX,NY,NZ,NT)

        loop = make_plaq(μ,ν)

        evaluate_wilson_loops!(plaq,loop ,U,[temp1,temp2,temp3])

        return plaq
    end

    function calc_Plaq_notrace_1d(U::Array{T,1}) where T <: GaugeFields

        NC = U[1].NC
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT

        temp1 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp2 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp3 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        plaq = GaugeFields_1d(NC,NX,NY,NZ,NT)

        loop = make_plaq()

        evaluate_wilson_loops!(plaq,loop ,U,[temp1,temp2,temp3])

        return plaq
    end


    function calc_Plaq(U::Array{T,1}) where T <: GaugeFields
        NC = U[1].NC
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT

        temp1 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp2 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp3 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        staple = GaugeFields_1d(NC,NX,NY,NZ,NT)

        return calc_Plaq!(U::Array{T,1},temp1,temp2,temp3,staple)
    end

    function calc_Plaq!(U::Array{T,1},temp1,temp2,temp3,staple) where T <: GaugeFields
        plaq = 0
        for μ=1:4
            clear!(staple)
            make_staple!(staple,U,μ,temp1,temp2,temp3)
            substitute!(temp1,U[μ])

            mul!(temp2,temp1,staple')
            #mul!(temp2,U[μ],staple')
            plaq += tr(temp2)
            #println(plaq)
            #exit()
        end
        #plaq /= 
        return plaq*0.5
    end

    function calc_Plaq!(U::Array{T,1},temps::Array{T_1d,1}) where {T <: GaugeFields,T_1d <: GaugeFields_1d}
        plaq = 0
        temp1 = temps[1]
        temp2 = temps[2]
        temp3 = temps[3]
        staple = temps[4]

        return calc_Plaq!(U,temp1,temp2,temp3,staple)
    end

    function calc_Plaq!(U::Array{T,1},temps::Array{T_1d,1}) where {T <: GaugeFields,T_1d <: GaugeFields}
        plaq = 0
        temp1 = temps[1]
        temp2 = temps[2]
        temp3 = temps[3]
        staple = temps[4]

        for μ=1:4
            clear!(staple)
            make_staple!(staple,U,μ,temp1,temp2,temp3)
            substitute!(temp1,U[μ])

            mul!(temp2,temp1,staple')
            #mul!(temp2,U[μ],staple')
            plaq += tr(temp2)
            #println(plaq)
            #exit()
        end
        #plaq /= 
        return plaq*0.5
    end

    function calc_GaugeAction(U::Array{T,1},gparam::GaugeActionParam_standard,temps::Array{T_1d,1}) where {T <: GaugeFields,T_1d <: GaugeFields_1d}
        plaq = calc_Plaq!(U,temps)

        Sg = -plaq*gparam.β/gparam.NTRACE
        return Sg,plaq
    end

    function calc_GaugeAction(U::Array{T,1},gparam::GaugeActionParam_standard,temps::Array{T_1d,1}) where {T <: GaugeFields,T_1d <: GaugeFields}
        plaq = calc_Plaq!(U,temps)

        Sg = -plaq*gparam.β/gparam.NTRACE
        return Sg,plaq
    end

    function calc_GaugeAction(U::Array{T,1},gparam::GaugeActionParam_autogenerator,temps::Array{T_1d,1}) where {T <: GaugeFields,T_1d <: GaugeFields_1d}
        Sg = 0
        loopaction = temps[4]        

        #=
        evaluate_wilson_loops!(loopaction,gparam.loops[1],U,temps[1:3])
        Sg += (-gparam.βs[1]/gparam.NTRACE)*tr(loopaction)/2
        println("Sg1: $Sg")

        evaluate_wilson_loops!(loopaction,gparam.loops[2],U,temps[1:3])
        Sg2 = (-gparam.βs[2]/gparam.NTRACE)*tr(loopaction)/2
        println("Sg2: $Sg2")
        println("dS12: $((Sg2-Sg)/Sg2)")
        Sg += Sg2
        return Sg,Sg/(-gparam.βs[1]/gparam.NTRACE)
        =#
        for i = 1:gparam.numactions

            evaluate_wilson_loops!(loopaction,gparam.loops[i],U,temps[1:3])
            Sg += (-gparam.βs[i]/gparam.NTRACE)*tr(loopaction)/2
        end

        
        return Sg,Sg/(-gparam.βs[1]/gparam.NTRACE)
        
    end



    """
    normalize!(u)
    c----------------------------------------------------------------------c
    c     normalizes links                                                 c
    c     input  x : arbitrary 3*3 complex matrix                          c
    c     output x : SU(3) matrix 
    """
    function normalize!(u::GaugeFields{T}) where T <: SU3
        NX = u.NX
        NY = u.NY
        NZ = u.NZ
        NT = u.NT

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        w1 = 0
                        w2 = 0
                        for ic=1:3
                            w1 += u[2,ic,ix,iy,iz,it]*conj(u[1,ic,ix,iy,iz,it])
                            w2 += u[1,ic,ix,iy,iz,it]*conj(u[1,ic,ix,iy,iz,it])
                        end
                        zerock2 = w2
                        if zerock2 == 0 
                            println("w2 is zero  !!  (in normlz)")
                            println("u[1,1),u[1,2),u[1,3) : ",u[1,1,ix,iy,iz,it], "\t",u[1,2,ix,iy,iz,it],"\t", u[1,3,ix,iy,iz,it])
                        end

                        w1 = -w1/w2

                        x4 = (u[2,1,ix,iy,iz,it]) + w1*u[1,1,ix,iy,iz,it]
                        x5 = (u[2,2,ix,iy,iz,it]) + w1*u[1,2,ix,iy,iz,it]
                        x6 = (u[2,3,ix,iy,iz,it]) + w1*u[1,3,ix,iy,iz,it]

                        w3 = x4*conj(x4) + x5*conj(x5) + x6*conj(x6)

                        zerock3 = w3
                        if zerock3 == 0
                            println("w3 is zero  !!  (in normlz)")
                            println("x4, x5, x6 : $x4, $x5, $x6")
                            exit()
                        end

                        u[2,1,ix,iy,iz,it] = x4
                        u[2,2,ix,iy,iz,it] = x5
                        u[2,3,ix,iy,iz,it] = x6

                        w3 = 1/sqrt(w3)
                        w2 = 1/sqrt(w2)

                        u[1,1,ix,iy,iz,it] = u[1,1,ix,iy,iz,it]*w2
                        u[1,2,ix,iy,iz,it] = u[1,2,ix,iy,iz,it]*w2
                        u[1,3,ix,iy,iz,it] = u[1,3,ix,iy,iz,it]*w2
                        u[2,1,ix,iy,iz,it] = u[2,1,ix,iy,iz,it]*w3
                        u[2,2,ix,iy,iz,it] = u[2,2,ix,iy,iz,it]*w3
                        u[2,3,ix,iy,iz,it] = u[2,3,ix,iy,iz,it]*w3

                        if zerock2*zerock3 == 0 
                            println("!! devided by zero !! (in normalize)")
                            println("w2 or w3 in normlz is zero !!")
                            println("w2, w3 : $w2, $w3   ")
                            exit()
                        end
                        #println(u[:,:,ix,iy,iz,it]'*u[:,:,ix,iy,iz,it])
                    end
                end
            end
        end
        m3complv!(u)
    end

    """
    normalize!(u)
    ----------------------------------------------------------------------c
         normalizes links                                                 c
         input  x : arbitrary 2*2 complex matrix                          c
         output x : SU(2) matrix 
    """
    function normalize!(u::GaugeFields{T}) where T <: SU2
        NX = u.NX
        NY = u.NY
        NZ = u.NZ
        NT = u.NT

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX

                        α = u[1,1,ix,iy,iz,it]
                        β = u[2,1,ix,iy,iz,it]
                        detU = abs(α)^2 + abs(β)^2
                        u[1,1,ix,iy,iz,it] = α/detU
                        u[2,1,ix,iy,iz,it] = β/detU
                        u[1,2,ix,iy,iz,it] = -conj(β)/detU
                        u[2,2,ix,iy,iz,it] = conj(α)/detU

                    end
                end
            end
        end
    end


    function m3complv!(a::SU3GaugeFields)
        aa = zeros(Float64,18)
        NX = a.NX
        NY = a.NY
        NZ = a.NZ
        NT = a.NT
    
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
    
                        aa[ 1] = real( a[1,1,ix,iy,iz,it])
                        aa[ 2] = imag(a[1,1,ix,iy,iz,it])
                        aa[ 3] = real( a[1,2,ix,iy,iz,it])
                        aa[ 4] = imag(a[1,2,ix,iy,iz,it])
                        aa[ 5] = real( a[1,3,ix,iy,iz,it])
                        aa[ 6] = imag(a[1,3,ix,iy,iz,it])
                        aa[ 7] = real( a[2,1,ix,iy,iz,it])
                        aa[ 8] = imag(a[2,1,ix,iy,iz,it])
                        aa[ 9] = real( a[2,2,ix,iy,iz,it])
                        aa[10] = imag(a[2,2,ix,iy,iz,it])
                        aa[11] = real( a[2,3,ix,iy,iz,it])
                        aa[12] = imag(a[2,3,ix,iy,iz,it])
    
                        aa[13] = aa[ 3]*aa[11] - aa[ 4]*aa[12] -
                                    aa[ 5]*aa[ 9] + aa[ 6]*aa[10]
                        aa[14] = aa[ 5]*aa[10] + aa[ 6]*aa[ 9] -
                                    aa[ 3]*aa[12] - aa[ 4]*aa[11]
                        aa[15] = aa[ 5]*aa[ 7] - aa[ 6]*aa[ 8] -
                                    aa[ 1]*aa[11] + aa[ 2]*aa[12]
                        aa[16] = aa[ 1]*aa[12] + aa[ 2]*aa[11] -
                                    aa[ 5]*aa[ 8] - aa[ 6]*aa[ 7]
                        aa[17] = aa[ 1]*aa[ 9] - aa[ 2]*aa[10] -
                                    aa[ 3]*aa[ 7] + aa[ 4]*aa[ 8]
                        aa[18] = aa[ 3]*aa[ 8] + aa[ 4]*aa[ 7] -
                                    aa[ 1]*aa[10] - aa[ 2]*aa[ 9]
    
                        a[3,1,ix,iy,iz,it] = aa[13]+im*aa[14]
                        a[3,2,ix,iy,iz,it] = aa[15]+im*aa[16]
                        a[3,3,ix,iy,iz,it] = aa[17]+im*aa[18]

                        #println(a[:,:,ix,iy,iz,it]'*a[:,:,ix,iy,iz,it] )
                    end
                end
            end
        end
    end



end