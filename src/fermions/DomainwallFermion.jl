module DomainwallFermion_module
    using LinearAlgebra

    import ..Actions:FermiActionParam,FermiActionParam_Wilson,
                FermiActionParam_WilsonClover,FermiActionParam_Staggered
    import ..Gaugefields:GaugeFields,GaugeFields_1d,SU3GaugeFields,SU2GaugeFields,SU3GaugeFields_1d,SU2GaugeFields_1d,
                staggered_phase,SUn,SU2,SU3,SUNGaugeFields,SUNGaugeFields_1d,SU


    import ..AbstractFermion:FermionFields,
        Wx!,Wdagx!,clear!,substitute_fermion!,Dx!,fermion_shift!,fermion_shiftB!,add!,set_wing_fermi!,WdagWx!,apply_periodicity
    import ..WilsonFermion_module:WilsonFermion

    struct DomainwallFermion <: FermionFields
            NC::Int64
            NX::Int64
            NY::Int64
            NZ::Int64
            NT::Int64
            N5::Int64
            f::Array{WilsonFermion,1}
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
        
        for i=1:N5
            @assert ωs[i] == 1  "Only standard domain wall fermion can be used. Others are not implemented, yet. Now we have ωs[i] = $(ωs[i]). Put ωs[i]=1."
        end
        bs .= b*ωs + c
        cs .= b*ωs - c
        Dirac_operator = "Domainwall"
        return DomainwallFermion(NC,NX,NY,NZ,NT,N5,f,ωs,b,c,bs,cs,M,m,Dirac_operator,eps,MaxCGstep,BoundaryCondition)
    end


    function Base.similar(x::DomainwallFermion)
        return DomainwallFermion(x.NC,x.NX,x.NY,x.NZ,x.NT,x.N5,x.r,x.M,x.m,x.eps,x.MaxCGstep,x.b,x.c,x.ωs,x.BoundaryCondition) 
    end

    function Base.:*(a::DomainwallFermion,b::DomainwallFermion) #a^+ * b
        c = 0.0im
        for i5=1:a.N5
            for α=1:4
                for it=1:a.NT
                    for iz=1:a.NZ
                        for iy=1:a.NY
                            for ix=1:a.NX
                                @simd for ic=1:a.NC
                                    c+= conj(a.f[i5].f[ic,ix,iy,iz,it,α])*b.f[i5].f[ic,ix,iy,iz,it,α]
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

        for i5=1:a.N5
            for i6=1:n6
                for i5=1:n5
                    for i4=1:n4
                        for i3=1:n3
                            for i2=1:n2
                                @simd for i1=1:n1
                                    #println(a.f[i1,i2,i3,i4,i5,i6],"\t",b.f[i1,i2,i3,i4,i5,i6] )
                                    c.f[i5].f[i1,i2,i3,i4,i5,i6] += alpha*a.f[i5].f[i1,i2,i3,i4,i5,i6] +beta*b.f[i5].f[i1,i2,i3,i4,i5,i6] 
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

        for i5=1:a.N5
            for i6=1:n6
                for i5=1:n5
                    for i4=1:n4
                        for i3=1:n3
                            for i2=1:n2
                                @simd for i1=1:n1
                                    #println(a.f[i1,i2,i3,i4,i5,i6],"\t",b.f[i1,i2,i3,i4,i5,i6] )
                                    c.f[i5].f[i1,i2,i3,i4,i5,i6] += alpha*a.f[i5].f[i1,i2,i3,i4,i5,i6] 
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

        for i5=1:a.N5
            for i6=1:n6
                for i5=1:n5
                    for i4=1:n4
                        for i3=1:n3
                            for i2=1:n2
                                @simd for i1=1:n1
                                    #println(a.f[i1,i2,i3,i4,i5,i6],"\t",b.f[i1,i2,i3,i4,i5,i6] )
                                    c.f[i5].f[i1,i2,i3,i4,i5,i6] = coeff*c.f[i5].f[i1,i2,i3,i4,i5,i6] + alpha*a.f[i5].f[i1,i2,i3,i4,i5,i6] 
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
        x::DomainwallFermion,m,temps::Array{TW,1}) where  {T <: DomainwallFermion,G <: GaugeFields,TW <:WilsonFermion}

        #temp = temps[4]
        #temp1 = temps[1]
        #temp2 = temps[2]

        for i5=1:xout.N5   
            j5=i5
            Dx!(xout.f[i5],U,x.f[j5],temps) #Dw*x
            add!(1,xout.f[i5],1,x.f[j5]) #D = x + Dw*x

        
            j5=i5+1
            if 1 <= j5 <= xout.N5
                #-P_-
                mul_1minusγ5x_add!(xout.f[i5],x.f[j5],-1) 
            end

            j5=i5-1
            if 1 <= j5 <= xout.N5
                #-P_+
                mul_1plusγ5x_add!(xout.f[i5],x.f[j5],-1) 
            end

            if i5==1
                j5 = xout.N5
                mul_1plusγ5x_add!(xout.f[i5],x.f[j5],m) 
            end

            if i5==xout.N5
                j5 = 1
                mul_1minusγ5x_add!(xout.f[i5],x.f[j5],m) 
            end

        end     

        return
    end


    function D5DWdagx!(xout::DomainwallFermion,U::Array{G,1},
        x::DomainwallFermion,m,temps::Array{TW,1}) where  {T <: DomainwallFermion,G <: GaugeFields,TW <:WilsonFermion}


        for i5=1:xout.N5   
            j5=i5
            Ddagx!(xout.f[i5],U,x.f[j5],temps) #Dw*x
            add!(1,xout.f[i5],1,x.f[j5]) #D = x + Dw*x

        
            j5=i5+1
            if 1 <= j5 <= xout.N5
                #-P_-
                mul_1plusγ5x_add!(xout.f[i5],x.f[j5],-1) 
            end

            j5=i5-1
            if 1 <= j5 <= xout.N5
                #-P_+
                mul_1minusγ5x_add!(xout.f[i5],x.f[j5],-1) 
            end

            if i5==1
                j5 = xout.N5
                mul_1minusγ5x_add!(xout.f[i5],x.f[j5],m) 
            end

            if i5==xout.N5
                j5 = 1
                mul_1plusγ5x_add!(xout.f[i5],x.f[j5],m) 
            end

        end     

        return
    end
    

end