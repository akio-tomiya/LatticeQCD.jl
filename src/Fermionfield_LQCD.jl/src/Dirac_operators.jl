module Dirac_operators_module
    using LinearAlgebra
    import ..AbstractFermionfields_module:Fermionfields,AbstractFermionfields,clear_fermion!,set_wing_fermion!,
                                            Dx!,add_fermion!,Wx!,Wdagx!
    import ..Gaugefield:AbstractGaugefields
    import ..FermionAction_module:FermiActionParam,FermiActionParam_Wilson,FermiActionParam_WilsonClover,FermiActionParam_Staggered
    import ..Clover:Make_CloverFμν


    abstract type Operator end
        
    abstract type Dirac_operator  <: Operator
    end

    abstract type Adjoint_Dirac_operator <: Operator
    end

    abstract type DdagD_operator  <: Operator
    end

    function Dirac_operator(U::Array{T,1},x,fparam) where  T <: AbstractGaugefields
        if typeof(fparam) <: FermiActionParam_Wilson
            W = Wilson_operator(U,x,fparam)
        elseif typeof(fparam) <: FermiActionParam_WilsonClover
            W = WilsonClover_operator(U,x,fparam)
        elseif typeof(fparam) <: FermiActionParam_Staggered
            W = Staggered_operator(U,x,fparam)
        #elseif typeof(fparam) == FermiActionParam_Domainwall
        #    W = Domainwall_operator(U,x,fparam)
        else
            error("unknown fparam: fparam is ",typeof(fparam))
        end
    end

    function DdagD_operator(U::Array{T,1},x,fparam) where  T <: AbstractGaugefields
        if typeof(fparam) <: FermiActionParam_Wilson
            W = DdagD_Wilson_operator(U,x,fparam)
        elseif typeof(fparam) <: FermiActionParam_WilsonClover
            W = DdagD_WilsonClover_operator(U,x,fparam)
        elseif typeof(fparam) <: FermiActionParam_Staggered
            W = DdagD_Staggered_operator(U,x,fparam)
        #elseif typeof(fparam) == FermiActionParam_Domainwall
        #    W = DdagD_Domainwall_operator(U,x,fparam)
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

    struct DdagD_Staggered_operator <: DdagD_operator 
        dirac::Staggered_operator

        function DdagD_Staggered_operator(U::Array{T,1},x,fparam) where  T <: AbstractGaugefields
            return new(Staggered_operator(U,x,fparam))
        end
    end

    struct Wilson_operator{T,WilsonFermion} <: Dirac_operator  where T <: AbstractGaugefields
        U::Array{T,1}
        _temporal_fermi::Array{WilsonFermion,1}

        function Wilson_operator(U::Array{T,1},x,fparam) where  T <: AbstractGaugefields
            num = 5
            WilsonFermion = typeof(x)
            _temporal_fermi = Array{WilsonFermion,1}(undef,num)
            for i=1:num
                _temporal_fermi[i] = similar(x)
            end
            return new{eltype(U),WilsonFermion}(U,_temporal_fermi)
        end
    end

    struct DdagD_Wilson_operator <: DdagD_operator 
        dirac::Wilson_operator
        function DdagD_Wilson_operator(U::Array{T,1},x,fparam) where  T <: AbstractGaugefields
            return new(Wilson_operator(U,x,fparam))
        end
    end


    struct WilsonClover_operator{T,WilsonFermion} <: Dirac_operator  where  T <: AbstractGaugefields
        U::Array{T,1}
        _temporal_fermi::Array{WilsonFermion,1}
        CloverFμν::Array{ComplexF64,4}


        function WilsonClover_operator(U::Array{T,1},x,fparam) where  T <: AbstractGaugefields
            num = 5
            WilsonFermion = typeof(x)
            _temporal_fermi = Array{WilsonFermion,1}(undef,num)
            for i=1:num
                _temporal_fermi[i] = similar(x)
            end

            _,_,NN... = size(U[1])

            #NV = U[1].NV
            #NX = U[1].NX
            #NY = U[1].NY
            #NZ = U[1].NZ
            #NT = U[1].NT
            #NC = U[1].NC

            temp1 = similar(U[1])#GaugeFields_1d(NC,NX,NY,NZ,NT)
            temps = Array{typeof(temp1),1}(undef,num)
            for i=1:num
                temps[i] = deepcopy(temp1)
            end


            CloverFμν = Make_CloverFμν(fparam,U,temps) 

            

            return new{eltype(U),WilsonFermion}(U,_temporal_fermi,CloverFμν)
        end
    end

    struct DdagD_WilsonClover_operator <: DdagD_operator 
        dirac::WilsonClover_operator
        function DdagD_WilsonClover_operator(U::Array{T,1},x,fparam) where  T <: AbstractGaugefields
            return new(WilsonClover_operator(U,x,fparam))
        end
    end



    struct Adjoint_Staggered_operator{T} <: Adjoint_Dirac_operator 
        parent::T
    end

    struct Adjoint_Wilson_operator{T} <: Adjoint_Dirac_operator 
        parent::T
    end


    struct Adjoint_WilsonClover_operator{T} <: Adjoint_Dirac_operator 
        parent::T
    end


    function Base.adjoint(A::T) where T <: Staggered_operator
        Adjoint_Staggered_operator{typeof(A)}(A)
    end

    function Base.adjoint(A::T) where T <: Wilson_operator
        Adjoint_Wilson_operator{typeof(A)}(A)
    end

    function Base.adjoint(A::T) where T <: WilsonClover_operator
        Adjoint_WilsonClover_operator{typeof(A)}(A)
    end

    Base.adjoint(A::Adjoint_Dirac_operator) = A.parent



    function LinearAlgebra.mul!(y::T1,A::T2,x::T3) where {T1 <:AbstractFermionfields,T2 <:  Dirac_operator, T3 <:  AbstractFermionfields}
        error("LinearAlgebra.mul!(y,A,x) is not implemented in type y:$(typeof(y)),A:$(typeof(A)) and x:$(typeof(x))")
    end

    function LinearAlgebra.mul!(y::T1,A::T2,x::T3) where {T1 <:AbstractFermionfields,T2 <:  Operator , T3 <:  AbstractFermionfields}
        error("LinearAlgebra.mul!(y,A,x) is not implemented in type y:$(typeof(y)),A:$(typeof(A)) and x:$(typeof(x))")
    end


    function LinearAlgebra.mul!(y::T1,A::T,x::T2)  where {T <: DdagD_operator,T1 <: AbstractFermionfields,T2 <: AbstractFermionfields} #y = A*x
        temp = A.dirac._temporal_fermi[5]
        mul!(temp,A.dirac,x)
        mul!(y,A.dirac',temp)

        return
    end


    function LinearAlgebra.mul!(y::T1,A::T2,x::T3) where {T1 <:AbstractFermionfields,T2 <:  Staggered_operator, T3 <:  AbstractFermionfields}
        @assert typeof(A._temporal_fermi[1]) == typeof(x) "type of A._temporal_fermi[1] $(typeof(A._temporal_fermi[1])) should be type of x: $(typeof(x))"
        temps = A._temporal_fermi
        temp = temps[4]
        Dx!(temp,A.U,x,temps)
        clear_fermion!(y)
        add_fermion!(y,A.mass,x,1,temp)
        set_wing_fermion!(y)

        
        #error("xout")
        return
    end

    function LinearAlgebra.mul!(y::T1,A::T2,x::T3) where {T1 <:AbstractFermionfields,T2 <:  DdagD_Staggered_operator, T3 <:  AbstractFermionfields}
        @assert typeof(A.dirac._temporal_fermi[1]) == typeof(x) "type of A._temporal_fermi[1] $(typeof(A.dirac._temporal_fermi[1])) should be type of x: $(typeof(x))"
        temps = A.dirac._temporal_fermi
        temp = temps[5]
        temp2 = temps[6]
        Dx!(temp,A.dirac.U,x,temps)
        Dx!(temp2,A.dirac.U,temp,temps)

        clear_fermion!(y)
        add_fermion!(y,A.dirac.mass^2,x,-1,temp2)
        set_wing_fermion!(y)
        
        #error("xout")
        return
    end

    function LinearAlgebra.mul!(y::T1,A::T2,x::T3) where {T1 <:AbstractFermionfields,T2 <: Adjoint_Staggered_operator, T3 <:  AbstractFermionfields}
        temps = A.parent._temporal_fermi
        temp = temps[4]
        Dx!(temp,A.parent.U,x,temps)
        clear_fermion!(y)
        add_fermion!(y,A.parent.mass,x,-1,temp)
        set_wing_fermion!(y)
        #println(y[1,1,1,1,1,1])
        return
    end

    function LinearAlgebra.mul!(y::T1,A::T2,x::T3) where {T1 <:AbstractFermionfields,T2 <:  Wilson_operator, T3 <:  AbstractFermionfields}
        Wx!(y,A.U,x,A._temporal_fermi) 
        return
    end

    function LinearAlgebra.mul!(y::T1,A::T2,x::T3) where {T1 <:AbstractFermionfields,T2 <:  WilsonClover_operator, T3 <:  AbstractFermionfields}
        Wx!(y,A.U,x,A._temporal_fermi,A.CloverFμν) 
        return
    end

    function LinearAlgebra.mul!(y::T1,A::T2,x::T3) where {T1 <:AbstractFermionfields,T2 <:  Adjoint_Wilson_operator, T3 <:  AbstractFermionfields}
        Wdagx!(y,A.parent.U,x,A.parent._temporal_fermi) 
    end

    function LinearAlgebra.mul!(y::T1,A::T2,x::T3) where {T1 <:AbstractFermionfields,T2 <:  Adjoint_WilsonClover_operator, T3 <:  AbstractFermionfields}
        Wdagx!(y,A.parent.U,x,A.parent._temporal_fermi,A.parent.CloverFμν) 
        return
    end

end