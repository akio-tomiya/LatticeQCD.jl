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
            b::Float64
            c::Float64
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
        hop = fparam.hop
        eps = fparam.eps
        MaxCGstep = fparam.MaxCGstep

        f = Array{WilsonFermion,1}(undef,N5)
        for i=1:N5
            f[i] = WilsonFermion(NC,NX,NY,NZ,NT,r,hop,eps,MaxCGstep,BoundaryCondition)
        end
        b = fparam.b
        c = fparam.c
        M = fparam.M
        m = fparam.m
        Dirac_operator = "Domainwall"
        return DomainwallFermion(NC,NX,NY,NZ,NT,N5,f,b,c,M,m,Dirac_operator,eps,MaxCGstep,BoundaryCondition)
    end

end