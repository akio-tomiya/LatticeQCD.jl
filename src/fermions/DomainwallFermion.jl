module DomainwallFermion_module
    using LinearAlgebra

    import ..Actions:FermiActionParam,FermiActionParam_Wilson,
                FermiActionParam_WilsonClover,FermiActionParam_Staggered
    import ..Gaugefields:GaugeFields,GaugeFields_1d,SU3GaugeFields,SU2GaugeFields,SU3GaugeFields_1d,SU2GaugeFields_1d,
                staggered_phase,SUn,SU2,SU3,SUNGaugeFields,SUNGaugeFields_1d,SU


    import ..AbstractFermion:FermionFields,
        Wx!,Wdagx!,clear!,substitute_fermion!,Dx!,fermion_shift!,fermion_shiftB!,add!,set_wing_fermi!,WdagWx!,apply_periodicity,
        gauss_distribution_fermi!,gauss_distribution_fermi_Z2!,Z4_distribution_fermi!,Ddagx!
    import ..WilsonFermion_module:WilsonFermion,mul_1minusγ5x_add!,mul_1plusγ5x_add!

    struct DomainwallFermion <: FermionFields
            NC::Int64
            NX::Int64
            NY::Int64
            NZ::Int64
            NT::Int64
            N5::Int64
            f::Array{WilsonFermion,1}
            r::Float64
            ωs::Array{Float64,1}
            b::Float64
            c::Float64
            bs::Array{Float64,1}
            cs::Array{Float64,1}
            M::Float64
            m::Float64

            eps::Float64
            Dirac_operator::String
            MaxCGstep::Int64
            BoundaryCondition::Array{Int8,1}
            sizeW::Int64
    end


    function DomainwallFermion(NC,NX,NY,NZ,NT,fparam::FermiActionParam,BoundaryCondition) 
        N5 = fparam.N5
        r = fparam.r
        M = fparam.M
        m = fparam.m
        eps = fparam.eps
        MaxCGstep = fparam.MaxCGstep
        b = fparam.b
        c = fparam.c
        ωs = fparam.ωs

        return DomainwallFermion(NC,NX,NY,NZ,NT,N5,r,M,m,eps,MaxCGstep,b,c,ωs,BoundaryCondition) 
    end

    function DomainwallFermion(NC,NX,NY,NZ,NT,N5,r,M,m,eps,MaxCGstep,b,c,ωs,BoundaryCondition) 
        hop = 1/(8r+2M)
        f = Array{WilsonFermion,1}(undef,N5)
        for i=1:N5
            f[i] = WilsonFermion(NC,NX,NY,NZ,NT,r,hop,eps,MaxCGstep,BoundaryCondition)
        end
        @assert b == 1 && c == 1 "Only standard domain wall fermion can be used. Others are not implemented, yet. Now we have b=$b c=$c. Put b=c=1."
        
        @assert length(ωs) == N5 "length(ωs) should be N5! now length(ωs) = $(length(ωs)) and N5 = $N5."

        for i=1:N5
            @assert ωs[i] == 1  "Only standard domain wall fermion can be used. Others are not implemented, yet. Now we have ωs[i] = $(ωs[i]). Put ωs[i]=1."
        end
        bs = b*ωs .+ c
        cs = b*ωs .- c
        Dirac_operator = "Domainwall"
        sizeW = NC*NX*NY*NZ*NT*4
        return DomainwallFermion(
            NC,#::Int64
            NX,#::Int64
            NY,#::Int64
            NZ,#::Int64
            NT,#::Int64
            N5,#::Int64
            f,#::Array{WilsonFermion,1}
            r,
            ωs,#::Array{Float64,1}
            b,#::Float64
            c,#::Float64
            bs,#::Array{Float64,1}
            cs,#::Array{Float64,1}
            M,#::Float64
            m,#::Float64

            eps,#::Float64
            Dirac_operator,#::String
            MaxCGstep,#::Int64
            BoundaryCondition,#::Array{Int8,1}  )  
            sizeW)
        #NC,NX,NY,NZ,NT,N5,f,ωs,b,c,bs,cs,M,m,Dirac_operator,eps,MaxCGstep,BoundaryCondition)
    end

    function Base.length(a::DomainwallFermion)
        return a.NC*a.NX*a.NY*a.NZ*a.NT*a.N5*4
    end

    function Base.setindex!(a::DomainwallFermion,v,i)
        # i = (ei5-1)*sizeW + ii
        ii = (i - 1) % a.sizeW + 1
        ei5 = (i-ii) ÷ a.sizeW + 1
        a.f[ei5][ii] = v
    end

    function Base.getindex(a::DomainwallFermion,i)
        # i = (ei5-1)*sizeW + ii
        ii = (i - 1) % a.sizeW + 1
        ei5 = (i-ii) ÷ a.sizeW + 1
        return a.f[ei5][ii]
    end


    function clear!(a::DomainwallFermion)
        n1,n2,n3,n4,n5,n6 = size(a.f[1].f)
        N5 = a.N5
        for ei5=1:N5
            for i6=1:n6
                for i5=1:n5
                    for i4=1:n4
                        for i3=1:n3
                            for i2=1:n2
                                @simd for i1=1:n1
                                    a.f[ei5].f[i1,i2,i3,i4,i5,i6]= 0
                                end
                            end
                        end
                    end
                end
            end
        end
    end


    function Base.similar(x::DomainwallFermion)
        return DomainwallFermion(x.NC,x.NX,x.NY,x.NZ,x.NT,x.N5,x.r,x.M,x.m,x.eps,x.MaxCGstep,x.b,x.c,x.ωs,x.BoundaryCondition) 
    end

    function Base.:*(a::DomainwallFermion,b::DomainwallFermion) #a^+ * b
        c = 0.0im
        for ei5=1:a.N5
            for α=1:4
                for it=1:a.NT
                    for iz=1:a.NZ
                        for iy=1:a.NY
                            for ix=1:a.NX
                                @simd for ic=1:a.NC
                                    c+= conj(a.f[ei5][ic,ix,iy,iz,it,α])*b.f[ei5][ic,ix,iy,iz,it,α]
                                end
                            end
                        end
                    end
                end            
            end
        end
        return c
    end

    function add!(c::DomainwallFermion,alpha::Number,a::DomainwallFermion,beta::Number,b::DomainwallFermion) #c = c + alpha*a + beta*b
        n1,n2,n3,n4,n5,n6 = size(a.f[1].f)

        for ei5=1:a.N5
            for i6=1:n6
                for i5=1:n5
                    for i4=1:n4
                        for i3=1:n3
                            for i2=1:n2
                                @simd for i1=1:n1
                                    #println(a.f[i1,i2,i3,i4,i5,i6],"\t",b.f[i1,i2,i3,i4,i5,i6] )
                                    c.f[ei5].f[i1,i2,i3,i4,i5,i6] += alpha*a.f[ei5].f[i1,i2,i3,i4,i5,i6] +beta*b.f[ei5].f[i1,i2,i3,i4,i5,i6] 
                                end
                            end
                        end
                    end
                end
            end
        end
        return
    end

    function add!(c::DomainwallFermion,alpha::Number,a::DomainwallFermion) #c = c + alpha*a 
        n1,n2,n3,n4,n5,n6 = size(a.f[1].f)

        for ei5=1:a.N5
            for i6=1:n6
                for i5=1:n5
                    for i4=1:n4
                        for i3=1:n3
                            for i2=1:n2
                                @simd for i1=1:n1
                                    #println(a.f[i1,i2,i3,i4,i5,i6],"\t",b.f[i1,i2,i3,i4,i5,i6] )
                                    c.f[ei5].f[i1,i2,i3,i4,i5,i6] += alpha*a.f[ei5].f[i1,i2,i3,i4,i5,i6] 
                                end
                            end
                        end
                    end
                end
            end
        end
        return
    end

    function add!(coeff::Number,c::DomainwallFermion,alpha::Number,a::DomainwallFermion) #c = coeff*c + alpha*a 
        n1,n2,n3,n4,n5,n6 = size(a.f[1].f)

        for ei5=1:a.N5
            for i6=1:n6
                for i5=1:n5
                    for i4=1:n4
                        for i3=1:n3
                            for i2=1:n2
                                @simd for i1=1:n1
                                    #println(a.f[i1,i2,i3,i4,i5,i6],"\t",b.f[i1,i2,i3,i4,i5,i6] )
                                    c.f[ei5].f[i1,i2,i3,i4,i5,i6] = coeff*c.f[ei5].f[i1,i2,i3,i4,i5,i6] + alpha*a.f[ei5].f[i1,i2,i3,i4,i5,i6] 
                                end
                            end
                        end
                    end
                end
            end
        end
        return
    end

    #add!(coeff::Number,c::FermionFields,alpha::Number,a::FermionFields) #c = coeff*c + alpha*a 


    function D5DWx!(xout::DomainwallFermion,U::Array{G,1},
        x::DomainwallFermion,m,temps::Array{TW,1},N5) where  {T <: DomainwallFermion,G <: GaugeFields,TW <:WilsonFermion}

        #temp = temps[4]
        #temp1 = temps[1]
        #temp2 = temps[2]
        clear!(xout)

        for i5=1:N5   
            j5=i5
            Dx!(xout.f[i5],U,x.f[j5],temps) #Dw*x
            #Wx!(xout.f[i5],U,x.f[j5],temps) #Dw*x
            set_wing_fermi!(xout.f[i5])
            add!(1,xout.f[i5],1,x.f[j5]) #D = x + Dw*x
            set_wing_fermi!(xout.f[i5])  

        
            j5=i5+1
            if 1 <= j5 <= N5
                #-P_-
                if N5 != 2
                    mul_1minusγ5x_add!(xout.f[i5],x.f[j5],-1) 
                    set_wing_fermi!(xout.f[i5])  
                end
            end

            j5=i5-1
            if 1 <= j5 <= N5
                #-P_+
                if N5 != 2
                    mul_1plusγ5x_add!(xout.f[i5],x.f[j5],-1) 
                    set_wing_fermi!(xout.f[i5])  
                end
            end

            if N5 != 1
                if i5==1
                    j5 = N5
                    mul_1plusγ5x_add!(xout.f[i5],x.f[j5],m) 
                    set_wing_fermi!(xout.f[i5])  
                end

                if i5== N5
                    j5 = 1
                    mul_1minusγ5x_add!(xout.f[i5],x.f[j5],m) 
                    set_wing_fermi!(xout.f[i5])  
                end
            end

        end  
        set_wing_fermi!(xout)   

        return
    end


    function D5DWdagx!(xout::DomainwallFermion,U::Array{G,1},
        x::DomainwallFermion,m,temps::Array{TW,1},N5) where  {T <: DomainwallFermion,G <: GaugeFields,TW <:WilsonFermion}

        #temp = temps[4]
        #temp1 = temps[1]
        #temp2 = temps[2]
        clear!(xout)

        for i5=1:N5   
            j5=i5
            Ddagx!(xout.f[i5],U,x.f[j5],temps) #Ddagw*x
            #Wdagx!(xout.f[i5],U,x.f[j5],temps) #Ddagw*x
            set_wing_fermi!(xout.f[i5])
            add!(1,xout.f[i5],1,x.f[j5]) #D = x + Ddagw*x
            set_wing_fermi!(xout.f[i5])  

        
            j5=i5+1
            if 1 <= j5 <= N5
                #-P_-
                if N5 != 2
                    mul_1plusγ5x_add!(xout.f[i5],x.f[j5],-1) 
                    set_wing_fermi!(xout.f[i5])  
                end
            end

            j5=i5-1
            if 1 <= j5 <= N5
                #-P_+
                if N5 != 2
                    mul_1minusγ5x_add!(xout.f[i5],x.f[j5],-1) 
                    set_wing_fermi!(xout.f[i5])  
                end
            end

            if N5 != 1
                if i5==1
                    j5 = N5
                    mul_1minusγ5x_add!(xout.f[i5],x.f[j5],m) 
                    set_wing_fermi!(xout.f[i5])  
                end

                if i5==N5
                    j5 = 1
                    mul_1plusγ5x_add!(xout.f[i5],x.f[j5],m) 
                    set_wing_fermi!(xout.f[i5])  
                end
            end

        end  
        set_wing_fermi!(xout)   

        return
    end


       


    """
c-------------------------------------------------c
c     Random number function for Gaussian  Noise
    with σ^2 = 1/2
c-------------------------------------------------c
    """
    function gauss_distribution_fermi!(x::DomainwallFermion)
        NC = x.NC
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        n6 = size(x.f[1])[6]
        σ = sqrt(1/2)

        for ei5=1:x.N5
            for ialpha = 1:n6
                for it=1:NT
                    for iz=1:NZ
                        for iy=1:NY
                            for ix=1:NX
                                for ic=1:NC
                                    
                                    x.f[ei5][ic,ix,iy,iz,it,ialpha] = σ*randn()+im*σ*randn()
                                end
                            end
                        end
                    end
                end
            end
        end

        set_wing_fermi!(x)

        return
    end

    """
c-------------------------------------------------c
c     Random number function for Gaussian  Noise
    with σ^2 = 1/2
c-------------------------------------------------c
    """
    function gauss_distribution_fermi!(x::DomainwallFermion,randomfunc,σ)
        NC = x.NC
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        n6 = size(x.f[1])[6]
        #σ = sqrt(1/2)

        for ei5=1:x.N5
            for mu = 1:n6
                for ic=1:NC
                    for it=1:NT
                        for iz=1:NZ
                            for iy=1:NY
                                for ix=1:NX
                                    v1 = sqrt(-log(randomfunc()+1e-10))
                                    v2 = 2pi*randomfunc()

                                    xr = v1*cos(v2)
                                    xi = v1 * sin(v2)

                                    x.f[ei5][ic,ix,iy,iz,it,mu] = σ*xr + σ*im*xi
                                end
                            end
                        end
                    end
                end
            end
        end

        set_wing_fermi!(x)

        return
    end

    function gauss_distribution_fermi!(x::DomainwallFermion,randomfunc)
        σ = 1
        gauss_distribution_fermi!(x,randomfunc,σ)
    end

    function gauss_distribution_fermi_Z2!(x::DomainwallFermion) 
        NC = x.NC
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        n6 = size(x.f[1])[6]
        #σ = sqrt(1/2)

        for ei5=1:x.N5
            for mu = 1:n6
                for it=1:NT
                    for iz=1:NZ
                        for iy=1:NY
                            for ix=1:NX
                                for ic=1:NC
                                    x.f[ei5][ic,ix,iy,iz,it,mu] = rand([-1,1])
                                end
                            end
                        end
                    end
                end
            end
        end

        set_wing_fermi!(x)

        return
    end

    """
c-------------------------------------------------c
c     Random number function Z4  Noise
c     https://arxiv.org/pdf/1611.01193.pdf
c-------------------------------------------------c
    """
    function Z4_distribution_fermi!(x::DomainwallFermion)
        NC = x.NC
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        n6 = size(x.f[1])[6]
        θ = 0.0
        N::Int32 = 4
        Ninv = Float64(1/N)
        for ei5=1:x.N5
            for ialpha = 1:n6
                for it=1:NT
                    for iz=1:NZ
                        for iy=1:NY
                            for ix=1:NX
                                for ic=1:NC
                                    θ = Float64(rand(0:N-1))*π*Ninv # r \in [0,π/4,2π/4,3π/4]
                                    x.f[ei5][ic,ix,iy,iz,it,ialpha] = cos(θ)+im*sin(θ) 
                                end
                            end
                        end
                    end
                end
            end
        end

        set_wing_fermi!(x)

        return
    end

    function set_wing_fermi!(a::DomainwallFermion)
        NT = a.NT
        NZ = a.NZ
        NY = a.NY
        NX = a.NX
        NC = a.NC
        for ei5=1:a.N5

            #!  X-direction
            for ialpha=1:4
                for it=1:NT
                    for iz = 1:NZ
                        for iy=1:NY
                            for k=1:NC
                                a.f[ei5][k,0,iy,iz,it,ialpha] = a.BoundaryCondition[1]*a.f[ei5][k,NX,iy,iz,it,ialpha]
                            end
                        end
                    end
                end
            end

            for ialpha=1:4
                for it=1:NT
                    for iz=1:NZ
                        for iy=1:NY
                            for k=1:NC
                                a.f[ei5][k,NX+1,iy,iz,it,ialpha] = a.BoundaryCondition[1]*a.f[ei5][k,1,iy,iz,it,ialpha]
                            end
                        end
                    end
                end
            end

            #Y-direction
            for ialpha = 1:4
                for it=1:NT
                    for iz=1:NZ
                        for ix=1:NX
                            for k=1:NC
                                a.f[ei5][k,ix,0,iz,it,ialpha] = a.BoundaryCondition[2]*a.f[ei5][k,ix,NY,iz,it,ialpha]
                            end
                        end
                    end
                end
            end

            for ialpha=1:4
                for it=1:NT
                    for iz=1:NZ
                        for ix=1:NX
                            for k=1:NC
                                a.f[ei5][k,ix,NY+1,iz,it,ialpha] = a.BoundaryCondition[2]*a.f[ei5][k,ix,1,iz,it,ialpha]
                            end
                        end
                    end
                end
            end

            
            for ialpha=1:4
                # Z-direction
                for it=1:NT
                    for iy=1:NY
                        for ix=1:NX
                            for k=1:NC
                                a.f[ei5][k,ix,iy,0,it,ialpha] = a.BoundaryCondition[3]*a.f[ei5][k,ix,iy,NZ,it,ialpha]
                                a.f[ei5][k,ix,iy,NZ+1,it,ialpha] = a.BoundaryCondition[3]*a.f[ei5][k,ix,iy,1,it,ialpha]

                            end
                        end
                    end
                end

                #T-direction
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            for k=1:NC
                                a.f[ei5][k,ix,iy,iz,0,ialpha] = a.BoundaryCondition[4]*a.f[ei5][k,ix,iy,iz,NT,ialpha]
                                a.f[ei5][k,ix,iy,iz,NT+1,ialpha] = a.BoundaryCondition[4]*a.f[ei5][k,ix,iy,iz,1,ialpha]
                            end
                        end
                    end
                end

            end

        end

        

        

    end


    function substitute_fermion!(H,j,x::DomainwallFermion)
        i = 0
        #println(x.N5)
        for ei5=1:x.N5
            for ialpha = 1:4
                for it=1:x.NT
                    for iz=1:x.NZ
                        for iy=1:x.NY
                            for ix=1:x.NX
                                @simd for ic=1:x.NC
                                    i += 1
                                    H[i,j] = x.f[ei5][ic,ix,iy,iz,it,ialpha]
                                end
                            end
                        end
                    end
                end
            end
        end
    end


end