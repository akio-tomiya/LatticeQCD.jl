module Diracoperators
    import ..Gaugefields:GaugeFields,GaugeFields_1d
    import ..Fermionfields:FermionFields,WilsonFermion,Wx!,Wdagx!,WdagWx!,
                            StaggeredFermion,set_wing_fermi!,Dx!,add!,clear!
    import ..Clover:Make_CloverFμν
    import ..Actions:FermiActionParam_Wilson,FermiActionParam_WilsonClover,
                FermiActionParam_Staggered
    using LinearAlgebra
    
    abstract type Dirac_operator 
    end

    abstract type Adjoint_Dirac_operator 
    end

    abstract type DdagD_operator 
    end





    function Dirac_operator(U::Array{T,1},x,fparam) where  T <: GaugeFields
        if typeof(fparam) == FermiActionParam_Wilson
            W = Wilson_operator(U,x,fparam)
        elseif typeof(fparam) == FermiActionParam_WilsonClover
            W = WilsonClover_operator(U,x,fparam)
        elseif typeof(fparam) == FermiActionParam_Staggered
            W = Staggered_operator(U,x,fparam)
        else
            error("unknown fparam: fparam is ",typeof(fparam))
        end
    end



    struct Wilson_operator{T} <: Dirac_operator  where  T <: GaugeFields
        U::Array{T,1}
        _temporal_fermi::Array{WilsonFermion,1}

        function Wilson_operator(U::Array{T,1},x,fparam) where  T <: GaugeFields
            num = 5
            _temporal_fermi = Array{WilsonFermion,1}(undef,num)
            for i=1:num
                _temporal_fermi[i] = similar(x)
            end
            return new{eltype(U)}(U,_temporal_fermi)
        end
    end

    struct DdagD_Wilson_operator <: DdagD_operator 
        dirac::Wilson_operator
        function DdagD_Wilson_operator(U::Array{T,1},x,fparam) where  T <: GaugeFields
            return new(Wilson_operator(U,x,fparam))
        end
    end

    struct WilsonClover_operator{T} <: Dirac_operator  where  T <: GaugeFields
        U::Array{T,1}
        _temporal_fermi::Array{WilsonFermion,1}
        CloverFμν::Array{ComplexF64,4}


        function WilsonClover_operator(U::Array{T,1},x,fparam) where  T <: GaugeFields
            num = 5
            _temporal_fermi = Array{WilsonFermion,1}(undef,num)
            for i=1:num
                _temporal_fermi[i] = similar(x)
            end

            NV = U[1].NV
            NX = U[1].NX
            NY = U[1].NY
            NZ = U[1].NZ
            NT = U[1].NT
            NC = U[1].NC

            temp1 = GaugeFields_1d(NC,NX,NY,NZ,NT)
            temps = Array{typeof(temp1),1}(undef,num)
            for i=1:num
                temps[i] = deepcopy(temp1)
            end


            CloverFμν = Make_CloverFμν(fparam,U,temps) 

            return new{eltype(U)}(U,_temporal_fermi,CloverFμν)
        end
    end

    struct DdagD_WilsonClover_operator <: DdagD_operator 
        dirac::WilsonClover_operator
        function DdagD_WilsonClover_operator(U::Array{T,1},x,fparam) where  T <: GaugeFields
            return new(WilsonClover_operator(U,x,fparam))
        end
    end

    struct Staggered_operator{T} <: Dirac_operator  where  T <: GaugeFields
        U::Array{T,1}
        _temporal_fermi::Array{StaggeredFermion,1}
        mass::Float64

        function Staggered_operator(U::Array{T,1},x,fparam) where  T <: GaugeFields
            num = 6
            _temporal_fermi = Array{StaggeredFermion,1}(undef,num)
            for i=1:num
                _temporal_fermi[i] = similar(x)
            end
            return new{eltype(U)}(U,_temporal_fermi,fparam.mass)
        end
    end

    struct DdagD_Staggered_operator <: DdagD_operator 
        dirac::Staggered_operator
        function DdagD_Staggered_operator(U::Array{T,1},x,fparam) where  T <: GaugeFields
            return new(Staggered_operator(U,x,fparam))
        end
    end



    function DdagD_operator(U::Array{T,1},x,fparam) where  T <: GaugeFields
        if typeof(fparam) == FermiActionParam_Wilson
            W = DdagD_Wilson_operator(U,x,fparam)
        elseif typeof(fparam) == FermiActionParam_WilsonClover
            W = DdagD_WilsonClover_operator(U,x,fparam)
        elseif typeof(fparam) == FermiActionParam_Staggered
            W = DdagD_Staggered_operator(U,x,fparam)
        else
            error("unknown fparam: fparam is ",typeof(fparam))
        end
    end

    struct Adjoint_Wilson_operator <: Adjoint_Dirac_operator 
        parent::Wilson_operator
    end

    struct Adjoint_WilsonClover_operator <: Adjoint_Dirac_operator 
        parent::WilsonClover_operator
    end

    struct Adjoint_Staggered_operator <: Adjoint_Dirac_operator 
        parent::Staggered_operator
    end
    
    function Base.adjoint(A::Wilson_operator)
        Adjoint_Wilson_operator(A)
    end

    function Base.adjoint(A::WilsonClover_operator)
        Adjoint_WilsonClover_operator(A)
    end

    function Base.adjoint(A::Staggered_operator)
        Adjoint_Staggered_operator(A)
    end

    Base.adjoint(A::Adjoint_Dirac_operator) = A.parent

    function LinearAlgebra.mul!(y::WilsonFermion,A::Wilson_operator,x::WilsonFermion) #y = A*x
        Wx!(y,A.U,x,A._temporal_fermi) 
        return
    end

    function LinearAlgebra.mul!(y::WilsonFermion,A::WilsonClover_operator,x::WilsonFermion) #y = A*x
        Wx!(y,A.U,x,A._temporal_fermi,A.CloverFμν) 
        return
    end

    function LinearAlgebra.mul!(y::FermionFields,A::T,x::FermionFields)  where T <: DdagD_operator #y = A*x
        temp = A.dirac._temporal_fermi[5]
        mul!(temp,A.dirac,x)
        mul!(y,A.dirac',temp)

        return
    end

    function LinearAlgebra.mul!(y::FermionFields,A::DdagD_operator,x::FermionFields,indices) #y = A*x        
        mul!(y,A,x)
        return
    end

    function LinearAlgebra.mul!(y::FermionFields,A::DdagD_Staggered_operator,x::FermionFields,indices) #y = A*x        
        WdagWx!(y,A.dirac.U,x,A.dirac._temporal_fermi,A.dirac.mass,indices)
        return
    end


    function LinearAlgebra.mul!(y::StaggeredFermion,A::Staggered_operator,x::StaggeredFermion) #y = A*x
        temps = A._temporal_fermi
        temp = temps[4]
        Dx!(temp,A.U,x,temps)
        clear!(y)
        add!(y,A.mass,x,1,temp)
        set_wing_fermi!(y)
        return
    end

    function Base.:*(A::Dirac_operator,x::FermionFields)
        y = similar(x)
        mul!(y,A,x)
        return y
    end

    function Base.:*(A::DdagD_operator,x::FermionFields)
        y = similar(x)
        mul!(y,A,x)
        return y
    end

    function LinearAlgebra.mul!(y::WilsonFermion,A::Adjoint_Wilson_operator,x::WilsonFermion) #y = A*x
        Wdagx!(y,A.parent.U,x,A.parent._temporal_fermi) 
    end

    function LinearAlgebra.mul!(y::WilsonFermion,A::Adjoint_WilsonClover_operator,x::WilsonFermion) #y = A*x
        Wdagx!(y,A.parent.U,x,A.parent._temporal_fermi,A.parent.CloverFμν) 
    end

    function LinearAlgebra.mul!(y::StaggeredFermion,A::Adjoint_Staggered_operator,x::StaggeredFermion) #y = A*x
        temps = A.parent._temporal_fermi
        temp = temps[4]
        Dx!(temp,A.parent.U,x,temps)
        clear!(y)
        add!(y,A.parent.mass,x,-1,temp)
        set_wing_fermi!(y)
        return
    end

    

    function Base.:*(A::Adjoint_Dirac_operator,x::WilsonFermion)
        y = similar(x)
        mul!(y,A,x)
        return y
    end
    
end