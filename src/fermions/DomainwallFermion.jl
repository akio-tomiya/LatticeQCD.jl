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
        hop = 1/(8r+2M)
        eps = fparam.eps
        MaxCGstep = fparam.MaxCGstep

        f = Array{WilsonFermion,1}(undef,N5)
        for i=1:N5
            f[i] = WilsonFermion(NC,NX,NY,NZ,NT,r,hop,eps,MaxCGstep,BoundaryCondition)
        end
        b = fparam.b
        c = fparam.c
        @assert b == 1 && c == 1 "Only standard domain wall fermion can be used. Others are not implemented, yet. Now we have b=$b c=$c. Put b=c=1."
        
        
        ωs = fparam.ωs
        for i=1:N5
            @assert ωs[i] == 1  "Only standard domain wall fermion can be used. Others are not implemented, yet. Now we have ωs[i] = $(ωs[i]). Put ωs[i]=1."
        end
        bs .= b*ωs + c
        cs .= b*ωs - c


        Dirac_operator = "Domainwall"
        return DomainwallFermion(NC,NX,NY,NZ,NT,N5,f,ωs,b,c,bs,cs,M,m,Dirac_operator,eps,MaxCGstep,BoundaryCondition)
    end


    #add!(coeff::Number,c::FermionFields,alpha::Number,a::FermionFields) #c = coeff*c + alpha*a 


    function D5DWx!(xout::DomainwallFermion,U::Array{G,1},
        x::DomainwallFermion,,m,temps::Array{T,1}) where  {T <: FermionFields,G <: GaugeFields}

        #temp = temps[4]
        #temp1 = temps[1]
        #temp2 = temps[2]

        for i5=1:xout.N5   
            j5=i5
            Dx!(xout.f[i5],U,x.f[j5],temps.f[i5]) #Dw*x
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
    

end