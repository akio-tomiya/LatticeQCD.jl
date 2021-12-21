module AbstractFermionfields_module
    using LinearAlgebra
    #using InteractiveUtils
    import ..Gaugefield:AbstractGaugefields,Abstractfields,shift_U,Staggered_Gaugefields,staggered_U

    abstract type Abstractfermion end

    abstract type AbstractFermionfields{NC,Dim}  <: Abstractfermion
    end

    abstract type Adjoint_fermion{T} <: Abstractfermion
    end

    struct Adjoint_fermionfields{T} <: Adjoint_fermion{T}
        parent::T
    end

    function Base.adjoint(F::T) where T <: Abstractfermion
        Adjoint_fermionfields{T}(F)
    end

    abstract type Shifted_fermionfields{NC,Dim} <: Abstractfermion
    end






    include("./AbstractFermionfields_4D.jl")

    include("./StaggeredFermion/StaggeredFermion.jl")
    include("./WilsonFermion/WilsonFermion.jl")

    function Fermionfields(params,NC,Dirac_operator,NN...)
        if Dirac_operator == "Staggered"
            fermion = StaggeredFermion(params,NC,NN...)
        elseif Dirac_operator == "Wilson"
            fermion = WilsonFermion(params,NC,NN...)
        else
            error("$Dirac_operator is not supported")
        end
        return fermion
    end

    function clear_fermion!(F::T) where T <: AbstractFermionfields
        error("clear_fermion! is not implemented in type $(typeof(F)) ")
    end

    function Base.similar(F::T) where T <: AbstractFermionfields
        error("Base.similar is not implemented in type $(typeof(F)) ")
    end

    function Z4_distribution_fermion!(F::T) where T <: AbstractFermionfields
        error("Z4_distribution_fermi! is not implemented in type $(typeof(F)) ")
    end

    function set_wing_fermion!(F::T) where T <: AbstractFermionfields
        error("set_wing_fermion! is not implemented in type $(typeof(F)) ")
    end

    function Dx!(temp::T,U,x,temps)  where T <: AbstractFermionfields
        error("Dx! is not implemented in type $(typeof(temp)) ")
    end

    function LinearAlgebra.mul!(y::T1,A::T2,x::T3) where {T1 <:Abstractfermion,T2 <: Abstractfields, T3 <: Abstractfermion} 
        error("LinearAlgebra.mul!(y,A,x) is not implemented in type y:$(typeof(y)),A:$(typeof(A)) and x:$(typeof(x))")
    end

    function add_fermion!(coeff::Number,c::T1,alpha::Number,a::T2) where {T1 <: Abstractfermion, T2 <: Abstractfermion}
        error("add_fermion! (c <- c +  coeff*c +alpha*a) is not implemented in type c:$T1,a:$T2")
    end

    function LinearAlgebra.axpby!(a::Number, X::T1, b::Number, Y::T2) where {T1 <: Abstractfermion, T2 <: Abstractfermion}
        error("LinearAlgebra.axpby(! (Y <- a*X + b*Y ) is not implemented in type X:$T1,Y:$T2")
    end

    function LinearAlgebra.dot(a::T1,b::T2) where {T1 <: Abstractfermion, T2 <: Abstractfermion}
        error("LinearAlgebra.dot is not implemented in type X:$T1,Y:$T2")
    end

    
end
