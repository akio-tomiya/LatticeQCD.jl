#=
module Gaugefields_4D_module
    using LinearAlgebra
    import ..AbstractGaugefields_module:AbstractGaugefields,Shifted_Gaugefields,shift_U,
                        Adjoint_Gaugefields,set_wing_U!,Abstractfields,construct_staple!,clear_U!,
                        calculate_Plaquette
    import Base
    =#

    abstract type Gaugefields_4D{NC} <: AbstractGaugefields{NC,4}
    end




    include("./gaugefields_4D_wing.jl")
    include("./gaugefields_4D_mpi.jl")

    function Base.size(U::Gaugefields_4D{NC}) where NC
        return NC,NC,U.NX,U.NY,U.NZ,U.NT
    end




    function clear_U!(U::Array{T,1}) where T <: Gaugefields_4D
        for μ=1:4
            clear_U!(U[μ])
        end
        
    end



    #=
    function calculate_Plaquet(U::Array{T,1}) where T <: Gaugefields_4D
        error("calculate_Plaquet is not implemented in type $(typeof(U)) ")
    end
    =#



#end