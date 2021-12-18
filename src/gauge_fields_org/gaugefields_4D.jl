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

    function Base.size(U::Gaugefields_4D{NC}) where NC
        return NC,NC,U.NX,U.NY,U.NZ,U.NT
    end


    #lattice shift
    function shift_U(U::TU,ν::T) where {T <: Integer,TU <: Gaugefields_4D}
        if ν == 1
            shift = (1,0,0,0)
        elseif ν == 2
            shift = (0,1,0,0)
        elseif ν == 3
            shift = (0,0,1,0)
        elseif ν == 4
            shift = (0,0,0,1)
        elseif ν == -1
                shift = (-1,0,0,0)
        elseif ν == -2
                shift = (0,-1,0,0)
        elseif ν == -3
                shift = (0,0,-1,0)
        elseif ν == -4
                shift = (0,0,0,-1)
        end

        Shifted_Gaugefields(U,shift,4)
    end

    function shift_U(U::TU,shift::NTuple{Dim,T}) where {Dim,T <: Integer,TU <: Gaugefields_4D}
        Shifted_Gaugefields(U,shift,4)
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

    function clear_U!(Uμ::Gaugefields_4D{NC}) where NC
        NT = Uμ.NT
        NZ = Uμ.NZ
        NY = Uμ.NY
        NX = Uμ.NX
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k2=1:NC                            
                            for k1=1:NC
                                Uμ[k1,k2,ix,iy,iz,it] = 0
                            end
                        end
                    end
                end
            end
        end
        set_wing_U!(Uμ)

    end

#end