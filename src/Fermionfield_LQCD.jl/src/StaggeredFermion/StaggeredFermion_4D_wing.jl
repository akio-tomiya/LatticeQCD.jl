struct StaggeredFermion_4D_wing{NC} <: StaggeredFermion{NC,4}
    NC::Int64
    NX::Int64
    NY::Int64
    NZ::Int64
    NT::Int64
    NDW::Int64
    f::Array{ComplexF64,6}

    mass::Float64
    eps::Float64
    Dirac_operator::String
    MaxCGstep::Int64
    BoundaryCondition::Array{Int8,1}

    function StaggeredFermion_4D_wing(NC::T,NDW::T,NX::T,NY::T,NZ::T,NT::T,mass,eps,MaxCGstep,BoundaryCondition) where T<: Integer
        f = zeros(ComplexF64,NC,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW)
        Dirac_operator = "Staggered"
        return new{NC}(NC,NX,NY,NZ,NT,NDW,f,mass,eps,Dirac_operator,MaxCGstep,BoundaryCondition)
    end

    function StaggeredFermion_4D_wing(params,NC::T,NDW::T,NX::T,NY::T,NZ::T,NT::T) where T<: Integer
        mass = params["mass"]
        eps = params["eps"]
        MaxCGstep = params["MaxCGstep"]
        BoundaryCondition = params["BoundaryCondition"]
        f = zeros(ComplexF64,NC,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW)
        Dirac_operator = "Staggered"
        return new{NC}(NC,NX,NY,NZ,NT,NDW,f,mass,eps,Dirac_operator,MaxCGstep,BoundaryCondition)
    end

end

function Base.similar(x::T) where T <: StaggeredFermion_4D_wing
    return StaggeredFermion_4D_wing(x.NC,x.NDW,x.NX,x.NY,x.NZ,x.NT,
                x.mass,x.eps,x.MaxCGstep,x.BoundaryCondition)
end