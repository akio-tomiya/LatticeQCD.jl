
struct StaggeredFermion_4D_wing{NC} <: AbstractFermionfields_4D{NC}
    NC::Int64
    NX::Int64
    NY::Int64
    NZ::Int64
    NT::Int64
    NDW::Int64
    NG::Int64 #size of the Gamma matrix. In Staggered fermion, this is one. 
    f::Array{ComplexF64,6}

    mass::Float64
    eps::Float64
    Dirac_operator::String
    MaxCGstep::Int64
    BoundaryCondition::Array{Int8,1}

    function StaggeredFermion_4D_wing(NC::T,NX::T,NY::T,NZ::T,NT::T,mass,eps,MaxCGstep,BoundaryCondition) where T<: Integer
        NG = 1
        NDW = 1
        #@assert NDW == 1 "only NDW = 1 is supported. Now NDW = $NDW"
        f = zeros(ComplexF64,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW,NG)
        Dirac_operator = "Staggered"
        return new{NC}(NC,NX,NY,NZ,NT,NDW,NG,f,mass,eps,Dirac_operator,MaxCGstep,BoundaryCondition)
    end

    function StaggeredFermion_4D_wing(params,NC::T,NX::T,NY::T,NZ::T,NT::T) where T<: Integer
        NDW = 1
        #@assert NDW == 1 "only NDW = 1 is supported. Now NDW = $NDW"
        mass = params["mass"]
        eps = params["eps"]
        MaxCGstep = params["MaxCGstep"]
        BoundaryCondition = params["BoundaryCondition"]
        NG = 1
        f = zeros(ComplexF64,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW,NG)
        Dirac_operator = "Staggered"
        return new{NC}(NC,NX,NY,NZ,NT,NDW,NG,f,mass,eps,Dirac_operator,MaxCGstep,BoundaryCondition)
    end

end



function Base.similar(x::T) where T <: StaggeredFermion_4D_wing
    return StaggeredFermion_4D_wing(x.NC,x.NX,x.NY,x.NZ,x.NT,
                x.mass,x.eps,x.MaxCGstep,x.BoundaryCondition)
end

function Dx!(xout::T,U::Array{G,1},
    x::T,temps::Array{T,1}) where  {T <: StaggeredFermion_4D_wing,G <:AbstractGaugefields}
    #temp = temps[4]
    temp1 = temps[1]
    temp2 = temps[2]

    #clear!(temp)
    set_wing_fermion!(x)
    clear_fermion!(xout)
    for ν=1:4
        xplus = shift_fermion(x,ν)
        Us = staggered_U(U[ν],ν)
        mul!(temp1,Us,xplus)


        xminus = shift_fermion(x,-ν)
        Uminus = shift_U(U[ν],-ν)
        Uminus_s = staggered_U(Uminus,ν)
        mul!(temp2,Uminus_s',xminus)
        
        add_fermion!(xout,0.5,temp1,-0.5,temp2)

        #fermion_shift!(temp1,U,ν,x)
        #fermion_shift!(temp2,U,-ν,x)
        #add!(xout,0.5,temp1,-0.5,temp2)
        
    end

    
    set_wing_fermion!(xout)

    return
end


function clear_fermion!(x::StaggeredFermion_4D_wing{NC},evensite) where NC
    ibush = ifelse(evensite,0,1)
    for it=1:x.NT
        for iz=1:x.NZ
            for iy=1:x.NY
                xran =1+(1+ibush+iy+iz+it)%2:2:x.NX
                for ix in xran
                    @simd for ic=1:NC
                        x[ic,ix,iy,iz,it,1] = 0
                    end
                end
            end
        end
    end
    return
end





"""
-------------------------------------------------c
     Random number function Z4  Noise
     https://arxiv.org/pdf/1611.01193.pdf
-------------------------------------------------c
"""
function Z4_distribution_fermion!(x::StaggeredFermion_4D_wing{NC}) where NC 
    NX = x.NX
    NY = x.NY
    NZ = x.NZ
    NT = x.NT
    #n6 = size(x.f)[6]
    θ = 0.0
    N::Int32 = 4
    Ninv = Float64(1/N)
    for ialpha = 1:1
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for ic=1:NC
                            θ = Float64(rand(0:N-1))*π*Ninv # r \in [0,π/4,2π/4,3π/4]
                            x[ic,ix,iy,iz,it,ialpha] = cos(θ)+im*sin(θ) 
                        end
                    end
                end
            end
        end
    end

    set_wing_fermion!(x)

    return
end


function set_wing_fermion!(a::StaggeredFermion_4D_wing{NC}) where NC 
    NT = a.NT
    NZ = a.NZ
    NY = a.NY
    NX = a.NX

    #!  X-direction
    for ialpha=1:1
        for it=1:NT
            for iz = 1:NZ
                for iy=1:NY
                    @simd for k=1:NC
                        a[k,0,iy,iz,it,ialpha] = a.BoundaryCondition[1]*a[k,NX,iy,iz,it,ialpha]
                    end
                end
            end
        end
    end

    for ialpha=1:1
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for k=1:NC
                        a[k,NX+1,iy,iz,it,ialpha] = a.BoundaryCondition[1]*a[k,1,iy,iz,it,ialpha]
                    end
                end
            end
        end
    end

    #Y-direction
    for ialpha = 1:1
        for it=1:NT
            for iz=1:NZ
                for ix=1:NX
                    @simd for k=1:NC
                        a[k,ix,0,iz,it,ialpha] = a.BoundaryCondition[2]*a[k,ix,NY,iz,it,ialpha]
                    end
                end
            end
        end
    end

    for ialpha=1:1
        for it=1:NT
            for iz=1:NZ
                for ix=1:NX
                    @simd for k=1:NC
                        a[k,ix,NY+1,iz,it,ialpha] = a.BoundaryCondition[2]*a[k,ix,1,iz,it,ialpha]
                    end
                end
            end
        end
    end

    
    for ialpha=1:1
        # Z-direction
        for it=1:NT
            for iy=1:NY
                for ix=1:NX
                    @simd for k=1:NC
                        a[k,ix,iy,0,it,ialpha] = a.BoundaryCondition[3]*a[k,ix,iy,NZ,it,ialpha]
                        a[k,ix,iy,NZ+1,it,ialpha] = a.BoundaryCondition[3]*a[k,ix,iy,1,it,ialpha]

                    end
                end
            end
        end

        #T-direction
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    @simd for k=1:NC
                        a[k,ix,iy,iz,0,ialpha] = a.BoundaryCondition[4]*a[k,ix,iy,iz,NT,ialpha]
                        a[k,ix,iy,iz,NT+1,ialpha] = a.BoundaryCondition[4]*a[k,ix,iy,iz,1,ialpha]
                    end
                end
            end
        end

    end

end
