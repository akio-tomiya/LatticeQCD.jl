module AbstractFermionfields_module
    using LinearAlgebra

    abstract type AbstractFermionfields{NC,Dim} 
    end

    include("./StaggeredFermion/StaggeredFermion.jl")

    function Fermionfields(params,NC,NDW,Dirac_operator,NN...)
        if Dirac_operator == "Staggered"
            fermion = StaggeredFermion(params,NC,NDW,NN...)
        else
            error("$Dirac_operator is not supported")
        end
    end

    function clear_F!(F::T) where T <: AbstractFermionfields
        error("clear_F! is not implemented in type $(typeof(F)) ")
    end

    function Base.similar(F::T) where T <: AbstractFermionfields
        error("Base.similar is not implemented in type $(typeof(F)) ")
    end

    
end
