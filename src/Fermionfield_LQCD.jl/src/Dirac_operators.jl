module Dirac_operators_module
    import ..AbstractFermionfields_module:Fermionfields
    import ..Gaugefield:AbstractGaugefields
    import ..FermionAction_module:FermiActionParam,FermiActionParam_Wilson,FermiActionParam_WilsonClover,FermiActionParam_Staggered

    abstract type Operator end
        
    abstract type Dirac_operator  <: Operator
    end

    abstract type Adjoint_Dirac_operator <: Operator
    end

    abstract type DdagD_operator  <: Operator
    end

    function Dirac_operator(U::Array{T,1},x,fparam) where  T <: AbstractGaugefields
        if typeof(fparam) == FermiActionParam_Wilson
            W = Wilson_operator(U,x,fparam)
        elseif typeof(fparam) == FermiActionParam_WilsonClover
            W = WilsonClover_operator(U,x,fparam)
        elseif typeof(fparam) == FermiActionParam_Staggered
            W = Staggered_operator(U,x,fparam)
        #elseif typeof(fparam) == FermiActionParam_Domainwall
        #    W = Domainwall_operator(U,x,fparam)
        else
            error("unknown fparam: fparam is ",typeof(fparam))
        end
    end

    struct Staggered_operator{T,StaggeredFermion} <: Dirac_operator  where  T <: AbstractGaugefields
        U::Array{T,1}
        _temporal_fermi::Array{StaggeredFermion,1}
        mass::Float64

        function Staggered_operator(U::Array{T,1},x,fparam) where  T <: AbstractGaugefields
            num = 9
            StaggeredFermion = typeof(x)
            _temporal_fermi = Array{StaggeredFermion,1}(undef,num)
            for i=1:num
                _temporal_fermi[i] = similar(x)
            end
            return new{eltype(U),StaggeredFermion}(U,_temporal_fermi,fparam.mass)
        end
    end
end