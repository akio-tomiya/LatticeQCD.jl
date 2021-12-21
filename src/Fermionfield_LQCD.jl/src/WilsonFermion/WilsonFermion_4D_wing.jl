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
    NDW::Int64
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
        return new{NC}(NC,NX,NY,NZ,NT,NG,NDW,f,
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

function Base.similar(x::T) where T <: WilsonFermion_4D_wing
    return WilsonFermion_4D_wing(x.NC,x.NX,x.NY,x.NZ,x.NT,
                x.r,x.hop,x.eps,x.MaxCGstep,x.BoundaryCondition)
end

function Wx!(xout::T,U::Array{G,1},
    x::T,temps::Array{T,1}) where  {T <: WilsonFermion_4D_wing,G <: AbstractGaugefields}
    temp = temps[4]
    temp1 = temps[1]
    temp2 = temps[2]

    clear_fermion!(temp)
    set_wing_fermion!(x)
    for ν=1:4
        xplus = shift_fermion(x,ν)
        mul!(temp1,U[ν],xplus)

        #fermion_shift!(temp1,U,ν,x)

        #... Dirac multiplication

        mul!(temp1,view(x.rminusγ,:,:,ν))

        xminus = shift_fermion(x,-ν)
        Uminus = shift_U(U[ν],-ν)

        mul!(temp2,Uminus',xminus)
     
        #
        #fermion_shift!(temp2,U,-ν,x)
        #mul!(temp2,view(x.rplusγ,:,:,ν),temp2)
        mul!(temp2,view(x.rplusγ,:,:,ν))

        add_fermion!(temp,x.hopp[ν],temp1,x.hopm[ν],temp2)
        
    end

    clear_fermion!(xout)
    add_fermion!(xout,1,x,-1,temp)

    #display(xout)
    #    exit()
    return
end

function Wdagx!(xout::T,U::Array{G,1},
    x::T,temps::Array{T,1}) where  {T <: WilsonFermion_4D_wing,G <: AbstractGaugefields}
    temp = temps[4]
    temp1 = temps[1]
    temp2 = temps[2]

    clear_fermion!(temp)
    set_wing_fermion!(x)
    for ν=1:4
        xplus = shift_fermion(x,ν)
        mul!(temp1,U[ν],xplus)

        #fermion_shift!(temp1,U,ν,x)

        #... Dirac multiplication
        #mul!(temp1,view(x.rminusγ,:,:,ν),temp1)
        mul!(temp1,view(x.rplusγ,:,:,ν))
        
        
        #
        xminus = shift_fermion(x,-ν)
        Uminus = shift_U(U[ν],-ν)

        mul!(temp2,Uminus',xminus)
        #fermion_shift!(temp2,U,-ν,x)
        #mul!(temp2,view(x.rminusγ,:,:,ν),temp2)
        mul!(temp2,view(x.rminusγ,:,:,ν))


        add_fermion!(temp,x.hopp[ν],temp1,x.hopm[ν],temp2)
        
        
        
    end

    clear_fermion!(xout)
    add_fermion!(xout,1,x,-1,temp)

    #display(xout)
    #    exit()
    return
end

function LinearAlgebra.mul!(x::WilsonFermion_4D_wing{NC},A::TA) where {TA <: AbstractMatrix, NC}
    NX = x.NX
    NY = x.NY
    NZ = x.NZ
    NT = x.NT

    #n6 = size(x.f)[6]
    #f = zeros(ComplexF64,4)
    #e = zeros(ComplexF64,4)
    
    for ic=1:NC
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX
                            e1 = x[ic,ix,iy,iz,it,1]
                            e2 = x[ic,ix,iy,iz,it,2]
                            e3 = x[ic,ix,iy,iz,it,3]
                            e4 = x[ic,ix,iy,iz,it,4]

                            x[ic,ix,iy,iz,it,1] = A[1,1]*e1+A[1,2]*e2+A[1,3]*e3+A[1,4]*e4
                            x[ic,ix,iy,iz,it,2] = A[2,1]*e1+A[2,2]*e2+A[2,3]*e3+A[2,4]*e4
                            x[ic,ix,iy,iz,it,3] = A[3,1]*e1+A[3,2]*e2+A[3,3]*e3+A[3,4]*e4
                            x[ic,ix,iy,iz,it,4] = A[4,1]*e1+A[4,2]*e2+A[4,3]*e3+A[4,4]*e4

                    end
                end
            end
        end
    end
end

function set_wing_fermion!(a::WilsonFermion_4D_wing{NC}) where NC 
    NT = a.NT
    NZ = a.NZ
    NY = a.NY
    NX = a.NX

    #!  X-direction
    for ialpha=1:4
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

    for ialpha=1:4
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
    for ialpha = 1:4
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

    for ialpha=1:4
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

    
    for ialpha=1:4
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

