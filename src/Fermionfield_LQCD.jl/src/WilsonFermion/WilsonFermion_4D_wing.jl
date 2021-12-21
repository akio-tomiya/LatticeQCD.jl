"""
Struct for WilsonFermion
"""
struct WilsonFermion_4D_wing{NC} <: AbstractFermionfields_4D{NC}
    NC::Int64
    NX::Int64
    NY::Int64
    NZ::Int64
    NT::Int64
    NG::Int64
    f::Array{ComplexF64,6}

    γ::Array{ComplexF64,3}
    rplusγ::Array{ComplexF64,3}
    rminusγ::Array{ComplexF64,3}
    hop::Float64 #Hopping parameter
    r::Float64 #Wilson term
    hopp::Array{ComplexF64,1}
    hopm::Array{ComplexF64,1}
    eps::Float64
    Dirac_operator::String
    MaxCGstep::Int64
    BoundaryCondition::Array{Int8,1}


    function WilsonFermion_4D_wing(NC::T,NX::T,NY::T,NZ::T,NT::T,r,hop,eps,MaxCGstep,BoundaryCondition) where T<: Integer
        γ,rplusγ,rminusγ = mk_gamma(r)
        NG = 4
        NDW = 1
        hopp = zeros(ComplexF64,4)
        hopm = zeros(ComplexF64,4)
        hopp .= hop
        hopm .= hop
        #@assert NDW == 1 "only NDW = 1 is supported. Now NDW = $NDW"
        f = zeros(ComplexF64,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW,NG)
        Dirac_operator = "Wilson"
        return new{NC}(NC,NX,NY,NZ,NT,NG,f,
                γ,rplusγ,rminusγ,hop,r,hopp,hopm,eps,Dirac_operator,MaxCGstep,BoundaryCondition)
    end


    function WilsonFermion_4D_wing(params,NC::T,NX::T,NY::T,NZ::T,NT::T) where T<: Integer
        eps = params["eps"]
        MaxCGstep = params["MaxCGstep"]
        BoundaryCondition = params["BoundaryCondition"]
        r = params["r"]
        hop = params["hop"]
        return WilsonFermion_4D_wing(NC,NX,NY,NZ,NT,r,hop,eps,MaxCGstep,BoundaryCondition)
    end



end
