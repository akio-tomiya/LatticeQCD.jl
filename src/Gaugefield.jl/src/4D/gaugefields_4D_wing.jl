#=
module Gaugefields_4D_wing_module
    using LinearAlgebra
    import ..AbstractGaugefields_module:AbstractGaugefields,Shifted_Gaugefields,shift_U,
                        Adjoint_Gaugefields,set_wing_U!,Abstractfields,construct_staple!,clear_U!,
                        calculate_Plaquette,substitute_U!,mul_skiplastindex!,partial_tr,add_U!,
                        Traceless_antihermitian!
    import Base
    import ..Gaugefields_4D_module:Gaugefields_4D
=#

    struct Gaugefields_4D_wing{NC} <: Gaugefields_4D{NC}
        U::Array{ComplexF64,6}
        NX::Int64
        NY::Int64
        NZ::Int64
        NT::Int64
        NDW::Int64
        NV::Int64
        NC::Int64

        function Gaugefields_4D_wing(NC::T,NDW::T,NX::T,NY::T,NZ::T,NT::T) where T<: Integer
            NV = NX*NY*NZ*NT
            U = zeros(ComplexF64,NC,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW)
            #U = Array{Array{ComplexF64,6}}(undef,4)
            #for μ=1:4
            #    U[μ] = zeros(ComplexF64,NC,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW)
            #end
            return new{NC}(U,NX,NY,NZ,NT,NDW,NV,NC)
        end
    end




    function Base.setindex!(x::Gaugefields_4D_wing,v,i1,i2,i3,i4,i5,i6) 
        @inbounds x.U[i1,i2,i3 + x.NDW,i4 + x.NDW,i5 + x.NDW,i6 + x.NDW] = v
    end

    @inline function Base.getindex(x::Gaugefields_4D_wing,i1,i2,i3,i4,i5,i6) 
        @inbounds return x.U[i1,i2,i3 .+ x.NDW,i4 .+ x.NDW,i5 .+ x.NDW,i6 .+ x.NDW]
    end

    @inline function get_latticeindex(i,NX,NY,NZ,NT)
        #i = (((it -1)*NZ+iz-1)*NY+iy-1)*NX + ix
        # i = (it-1)*NZ*NY*NX + (iz-1)*NY*NX + (iy-1)*NX + ix
        ix = (i-1) % NX + 1
        ii = div(i-ix,NX)
        #ii = ((it -1)*NZ+iz-1)*NY+iy-1)
        iy = ii % NY + 1
        ii = div(ii-(iy-1),NY)
        #ii = (it -1)*NZ+iz-1)
        iz = ii % NZ + 1
        it = div(ii -(iz-1),NZ)
        return ix,iy,iz,it        
    end

    function Base.setindex!(x::Gaugefields_4D_wing,v,i1,i2,ii) 
        ix,iy,iz,it =  get_latticeindex(ii,x.NX,x.NY,x.NZ,x.NT)
        @inbounds x.U[i1,i2,ix + x.NDW,iy + x.NDW,iz + x.NDW,it + x.NDW] = v
    end

    @inline function Base.getindex(x::Gaugefields_4D_wing,i1,i2,ii)
        ix,iy,iz,it =  get_latticeindex(ii,x.NX,x.NY,x.NZ,x.NT)
        @inbounds return x.U[i1,i2,ix + x.NDW,iy + x.NDW,iz + x.NDW,it + x.NDW]
    end


    @inline function Base.getindex(U::Adjoint_Gaugefields{T},i1,i2,i3,i4,i5,i6) where T <: Abstractfields #U'
        @inbounds return conj(U.parent[i2,i1,i3,i4,i5,i6])
    end

    function Base.setindex!(U::Adjoint_Gaugefields{T},v,i1,i2,i3,i4,i5,i6,μ)  where T <: Abstractfields
        error("type $(typeof(U)) has no setindex method. This type is read only.")
    end

    @inline function Base.getindex(U::Adjoint_Gaugefields{T},i1,i2,ii) where T <: Abstractfields #U'
        @inbounds return conj(U.parent[i2,i1,ii])
    end

    function substitute_U!(a::Array{T1,1},b::Array{T2,1}) where {T1 <: Gaugefields_4D_wing,T2 <: Gaugefields_4D_wing}
        for μ=1:4
            substitute_U!(a[μ],b[μ])
        end
    end

    function substitute_U!(a::Array{T1,1},b::Array{T2,1},iseven) where {T1 <: Gaugefields_4D_wing,T2 <: Gaugefields_4D_wing}
        for μ=1:4
            substitute_U!(a[μ],b[μ],iseven)
        end
    end

    function Base.similar(U::T) where T <: Gaugefields_4D_wing
        Uout = Gaugefields_4D_wing(U.NC,U.NDW,U.NX,U.NY,U.NZ,U.NT)
        #identityGaugefields_4D_wing(U.NC,U.NX,U.NY,U.NZ,U.NT,U.NDW)
        return Uout
    end

    function Base.similar(U::Array{T,1}) where T <: Gaugefields_4D_wing
        Uout = Array{T,1}(undef,4)
        for μ=1:4
            Uout[μ] = similar(U[μ]) 
        end
        return Uout
    end

    function substitute_U!(a::T,b::T) where T <: Gaugefields_4D_wing
        for i=1:length(a.U)
            a.U[i] = b.U[i]
        end
        return 
    end

    function substitute_U!(a::Gaugefields_4D_wing{NC},b::T2) where {NC, T2 <: Abstractfields}
        NT = a.NT
        NZ = a.NZ
        NY = a.NY
        NX = a.NX
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k2=1:NC                            
                            for k1=1:NC
                                @inbounds a[k1,k2,ix,iy,iz,it] = b[k1,k2,ix,iy,iz,it]
                            end
                        end
                    end
                end
            end
        end
        set_wing_U!(a)

    end


    function substitute_U!(a::Gaugefields_4D_wing{NC},b::T2,iseven::Bool) where {NC, T2 <: Abstractfields}
        NT = a.NT
        NZ = a.NZ
        NY = a.NY
        NX = a.NX
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        evenodd = ifelse( (ix+iy+iz+it) % 2 ==0, true,false)
                        if evenodd == iseven
                            for k2=1:NC                            
                                for k1=1:NC
                                    @inbounds a[k1,k2,ix,iy,iz,it] = b[k1,k2,ix,iy,iz,it]
                                end
                            end
                        end
                    end
                end
            end
        end
        set_wing_U!(a)

    end

    function randomGaugefields_4D_wing(NC,NX,NY,NZ,NT,NDW)
        U = Gaugefields_4D_wing(NC,NDW,NX,NY,NZ,NT)
    
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for j=1:NC
                            @simd for i=1:NC
                                U[i,j,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
                            end
                        end
                    end
                end
            end
        end
        normalize_U!(U)
        set_wing_U!(U)
        return U
    end


    function identityGaugefields_4D_wing(NC,NX,NY,NZ,NT,NDW)
        U = Gaugefields_4D_wing(NC,NDW,NX,NY,NZ,NT)
    
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        @simd for ic=1:NC
                            U[ic,ic,ix,iy,iz,it] = 1 
                        end
                    end
                end
            end
        end
        set_wing_U!(U)
        return U
    end

    function Oneinstanton_4D_wing(NC,NX,NY,NZ,NT,NDW)
        @assert NC == 2 "NC should be 2"
        u = Gaugefields_4D_wing(NC,NDW,NX,NY,NZ,NT)
        U = Array{typeof(u),1}(undef,4)
        U[1] = u
        for μ = 2:4
            U[μ] = similar(u)
        end
        L = (NX,NY,NZ,NT)
        NV = prod(L)

        R = div(NX,2) # instanton radius

        println("# Starting from a instanton backgorund with radius R=$R ")
        inst_cent = [L[1]/2+0.5, L[2]/2+0.5, L[3]/2+0.5, L[4]/2+0.5]
        #eps = 1/10000000
        s1 = [
        0.0 1.0
        1.0 0.0]
        s2 = [
        0.0 -im*1.0
        im*1.0 0.0]
        s3 = [
        1.0 0.0
        0.0 -1.0]
        En = [
        1.0 0.0
        0.0 1.0]
        ss=[ im*s1,  im*s2,  im*s3, En]
        sd=[-im*s1, -im*s2, -im*s3, En]
        nn = 0

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX  
                        nn += 0
                        nv=[ix-1, iy-1, iz-1, it-1]-inst_cent
                        n2=nv⋅nv
                        for mu = 1:4
                            tau=[
                                0 0
                                0 0]
                            for nu = 1:4
                                smunu = sd[mu]*ss[nu]-sd[nu]*ss[mu]
                                tau+=smunu*nv[nu]
                            end
                            sq = sqrt(n2+R^2)
                            tau = exp(im*tau*(1/2)*(1/(n2))*(im*R^2/(n2+R^2)) ) #1b
                            for j=1:2
                                for i=1:2
                                    U[mu][i,j,ix,iy,iz,it] = tau[i,j]
                                end
                            end
                        
                        end
                    end
                end
            end
        end


        return U



    end




    struct Shifted_Gaugefields_4D{NC,outside} <: Shifted_Gaugefields{NC,4} 
        parent::Gaugefields_4D_wing{NC}
        #parent::T
        shift::NTuple{4,Int8}
        NX::Int64
        NY::Int64
        NZ::Int64
        NT::Int64
        NDW::Int64

        #function Shifted_Gaugefields(U::T,shift,Dim) where {T <: AbstractGaugefields}
        function Shifted_Gaugefields_4D(U::Gaugefields_4D_wing{NC},shift) where NC
            outside = check_outside(U.NDW,shift)
            return new{NC,outside}(U,shift,U.NX,U.NY,U.NZ,U.NT,U.NDW)
        end
    end

    function check_outside(NDW,shift)
        outside = false
        for μ=1:4
            if outside == false
                outside = ifelse(abs(shift[μ]) > NDW,true,false)
            end
        end
        #if outside
        #    println(outside,"\t",shift)
        #end
        return outside        
    end

        #lattice shift
    function shift_U(U::Gaugefields_4D_wing{NC},ν::T) where {T <: Integer,NC}
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
    
            return Shifted_Gaugefields_4D(U,shift)
    end
    
    function shift_U(U::TU,shift::NTuple{Dim,T}) where {Dim,T <: Integer,TU <: Gaugefields_4D}
        return Shifted_Gaugefields_4D(U,shift)
    end



    #function Base.setindex!(U::Shifted_Gaugefields{T,4},v,i1,i2,i3,i4,i5,i6)  where T <: Gaugefields_4D_wing
    function Base.setindex!(U::Shifted_Gaugefields{NC,4},v,i1,i2,i3,i4,i5,i6)  where NC
        error("type $(typeof(U)) has no setindex method. This type is read only.")
    end

    #=
    function Base.getindex(U::Shifted_Gaugefields_4D{NC,false},i1,i2,i3,i4,i5,i6) where NC
    #function Base.getindex(U::Shifted_Gaugefields{T,4},i1,i2,i3,i4,i5,i6) where T <: Gaugefields_4D_wing
        return U.parent.U[i1,i2,i3 .+ U.parent.NDW .+ U.shift[1],i4 .+ U.parent.NDW .+ U.shift[2],i5 .+ U.parent.NDW .+ U.shift[3],i6 .+ U.parent.NDW .+ U.shift[4]]
    end
    =#

    @inline function Base.getindex(U::Shifted_Gaugefields_4D{NC,false},i1,i2,i3,i4,i5,i6) where NC
        #function Base.getindex(U::Shifted_Gaugefields{T,4},i1,i2,i3,i4,i5,i6) where T <: Gaugefields_4D_wing
        @inbounds return U.parent[i1,i2,i3 .+ U.shift[1],i4 .+ U.shift[2],i5 .+ U.shift[3],i6 .+ U.shift[4]]
    end

    @inline function Base.getindex(U::Shifted_Gaugefields_4D{NC,true},i1,i2,i3,i4,i5,i6) where NC
        i3_new = i3 + U.shift[1]
        i3_new += ifelse(i3_new > U.NX + U.NDW,-U.NX,0)
        i3_new += ifelse(i3_new < 1 - U.NDW,U.NX,0)
        i4_new = i4 + U.shift[2]
        i4_new += ifelse(i4_new > U.NY + U.NDW,-U.NY,0)
        i4_new += ifelse(i4_new < 1 - U.NDW,U.NY,0)
        i5_new = i5 + U.shift[3]
        i5_new += ifelse(i5_new > U.NZ + U.NDW,-U.NZ,0)
        i5_new += ifelse(i5_new < 1 - U.NDW,U.NZ,0)
        i6_new = i6 + U.shift[4]
        i6_new += ifelse(i6_new > U.NT + U.NDW,-U.NT,0)
        i6_new += ifelse(i6_new < 1 - U.NDW,U.NT,0)
        #function Base.getindex(U::Shifted_Gaugefields{T,4},i1,i2,i3,i4,i5,i6) where T <: Gaugefields_4D_wing
        @inbounds return U.parent[i1,i2,i3_new ,i4_new ,i5_new ,i6_new ]
    end
    
    #=

    function Base.getindex(U::Adjoint_Gaugefields{Shifted_Gaugefields_4D{NC,false}},i1,i2,i3,i4,i5,i6) where NC
    #function Base.getindex(U::Adjoint_Gaugefields{Shifted_Gaugefields{T,4}},i1,i2,i3,i4,i5,i6) where T <: Gaugefields_4D_wing
        return conj(U.parent[i2,i1,i3,i4,i5,i6])
    end

    function Base.setindex!(U::Adjoint_Gaugefields{Shifted_Gaugefields_4D{T}},v,i1,i2,i3,i4,i5,i6)  where T <: Gaugefields_4D_wing
        error("type $(typeof(U)) has no setindex method. This type is read only.")
    end

    =#


    function Base.getindex(u::Staggered_Gaugefields{T,μ},i1,i2,i3,i4,i5,i6) where {T <: Gaugefields_4D_wing,μ}
        NT = u.parent.NT
        NZ = u.parent.NZ
        NY = u.parent.NY
        NX = u.parent.NX

        t = i6-1
        t += ifelse(t<0,NT,0)
        t += ifelse(t ≥ NT,-NT,0)
        #boundary_factor_t = ifelse(t == NT -1,BoundaryCondition[4],1)
        z = i5-1
        z += ifelse(z<0,NZ,0)
        z += ifelse(z ≥ NZ,-NZ,0)
        #boundary_factor_z = ifelse(z == NZ -1,BoundaryCondition[3],1)
        y = i4-1
        y += ifelse(y<0,NY,0)
        y += ifelse(y ≥ NY,-NY,0)
        #boundary_factor_y = ifelse(y == NY -1,BoundaryCondition[2],1)
        x = i3-1
        x += ifelse(x<0,NX,0)
        x += ifelse(x ≥ NX,-NX,0)
        #boundary_factor_x = ifelse(x == NX -1,BoundaryCondition[1],1)
        if μ ==1
            η = 1
        elseif μ ==2
            #η = (-1.0)^(x)
            η = ifelse(x%2 == 0,1,-1)
        elseif μ ==3
            #η = (-1.0)^(x+y)
            η = ifelse((x+y)%2 == 0,1,-1)
        elseif μ ==4
            #η = (-1.0)^(x+y+z)
            η = ifelse((x+y+z)%2 == 0,1,-1)
        else
            error("η should be positive but η = $η")
        end

        @inbounds return η*u.parent[i1,i2,i3,i4,i5,i6]
    end

    function Base.getindex(u::Staggered_Gaugefields{T,μ},i1,i2,i3,i4,i5,i6) where {T <: Shifted_Gaugefields_4D,μ}
        NT = u.parent.NT
        NZ = u.parent.NZ
        NY = u.parent.NY
        NX = u.parent.NX

        t = i6-1 + u.parent.shift[4]
        t += ifelse(t<0,NT,0)
        t += ifelse(t ≥ NT,-NT,0)
        #boundary_factor_t = ifelse(t == NT -1,BoundaryCondition[4],1)
        z = i5-1+ u.parent.shift[3]
        z += ifelse(z<0,NZ,0)
        z += ifelse(z ≥ NZ,-NZ,0)
        #boundary_factor_z = ifelse(z == NZ -1,BoundaryCondition[3],1)
        y = i4-1+ u.parent.shift[2]
        y += ifelse(y<0,NY,0)
        y += ifelse(y ≥ NY,-NY,0)
        #boundary_factor_y = ifelse(y == NY -1,BoundaryCondition[2],1)
        x = i3-1+ u.parent.shift[1]
        x += ifelse(x<0,NX,0)
        x += ifelse(x ≥ NX,-NX,0)
        #boundary_factor_x = ifelse(x == NX -1,BoundaryCondition[1],1)
        if μ ==1
            η = 1
        elseif μ ==2
            #η = (-1.0)^(x)
            η = ifelse(x%2 == 0,1,-1)
        elseif μ ==3
            #η = (-1.0)^(x+y)
            η = ifelse((x+y)%2 == 0,1,-1)
        elseif μ ==4
            #η = (-1.0)^(x+y+z)
            η = ifelse((x+y+z)%2 == 0,1,-1)
        else
            error("η should be positive but η = $η")
        end

        @inbounds  return η*u.parent[i1,i2,i3,i4,i5,i6]
    end


    function map_U!(U::Gaugefields_4D_wing{NC},f!::Function,V::Gaugefields_4D_wing{NC},iseven::Bool) where {NC} 
        NT = U.NT
        NZ = U.NZ
        NY = U.NY
        NX = U.NX
        A = zeros(ComplexF64,NC,NC)
        B = zeros(ComplexF64,NC,NC)
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        evenodd = ifelse( (ix+iy+iz+it) % 2 ==0, true,false)
                        if evenodd == iseven     
                            for k2=1:NC                            
                                for k1=1:NC
                                    A[k1,k2] = V[k1,k2,ix,iy,iz,it]
                                    B[k1,k2] = U[k1,k2,ix,iy,iz,it]
                                end
                            end
                            f!(B,A)
                            for k2=1:NC                            
                                for k1=1:NC
                                    U[k1,k2,ix,iy,iz,it] = B[k1,k2]
                                end
                            end
                        end
                    end
                end
            end
        end
        set_wing_U!(U)
    end


    
    function clear_U!(Uμ::Gaugefields_4D_wing{NC}) where NC
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
                                @inbounds Uμ[k1,k2,ix,iy,iz,it] = 0
                            end
                        end
                    end
                end
            end
        end
        set_wing_U!(Uμ)

    end

    function clear_U!(Uμ::Gaugefields_4D_wing{NC},iseven::Bool) where NC
        NT = Uμ.NT
        NZ = Uμ.NZ
        NY = Uμ.NY
        NX = Uμ.NX
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        evenodd = ifelse( (ix+iy+iz+it) % 2 ==0, true,false)
                        if evenodd == iseven     
                            for k2=1:NC                            
                                for k1=1:NC
                                    @inbounds Uμ[k1,k2,ix,iy,iz,it] = 0
                                end
                            end
                        end
                    end
                end
            end
        end
        set_wing_U!(Uμ)

    end

        
    function unit_U!(Uμ::Gaugefields_4D_wing{NC}) where NC
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
                                @inbounds Uμ[k1,k2,ix,iy,iz,it] = ifelse(k1 == k2,1,0)
                            end
                        end
                    end
                end
            end
        end
        set_wing_U!(Uμ)

    end

    """
    M = (U*δ_prev) star (dexp(Q)/dQ)
    Λ = TA(M)
    """
    function construct_Λmatrix_forSTOUT!(Λ,δ_current::Gaugefields_4D_wing{NC},Q,u::Gaugefields_4D_wing{NC}) where NC
        NT = u.NT
        NY = u.NY
        NZ = u.NZ
        NX = u.NX
        Qn = zeros(ComplexF64,NC,NC)
        Un = zero(Qn)
        Mn = zero(Qn)
        Λn = zero(Qn)
        δn_current = zero(Qn)
        temp1 = similar(Qn)
        temp2 = similar(Qn)
        temp3 = similar(Qn)

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX

                        for jc=1:NC
                            for ic=1:NC
                                Un[ic,jc] = u[ic,jc,ix,iy,iz,it]
                                Qn[ic,jc] = Q[ic,jc,ix,iy,iz,it]#*im
                                #if (ix,iy,iz,it) == (1,1,1,1)
                                #    println("Qij $ic $jc ",Qn[ic,jc],"\t ",Q[ic,jc,ix,iy,iz,it])
                                #end
                                δn_current[ic,jc] = δ_current[ic,jc,ix,iy,iz,it]
                            end
                        end

                        #=
                        if (ix,iy,iz,it) == (1,1,1,1)
                            println("Qn = ",Qn[1:3,1:3])
                        end
                        =#

                        calc_Mmatrix!(Mn,δn_current,Qn,Un,u,[temp1,temp2,temp3])
                        #=
                        if (ix,iy,iz,it) == (1,1,1,1)
                            println(" Un ",  Un[1,1])
                            println("δn_current ", δn_current[1,1])
                            #println("Qn[1,1] = ",Qn[1,1])
                            println("M[1,1] = ",Mn[1,1])

                            #println("Qn = ",Qn[1:3,1:3])
                        end
                        =#
                        calc_Λmatrix!(Λn,Mn,NC)


                        for jc=1:NC
                            for ic=1:NC
                                Λ[ic,jc,ix,iy,iz,it] = Λn[ic,jc]
                            end
                        end

                    end
                end
            end
        end
        set_wing_U!(Λ)
    end


    function set_wing_U!(u::Array{Gaugefields_4D_wing{NC},1}) where NC
        for μ=1:4
            set_wing_U!(u[μ]) 
        end
    end
    

    function set_wing_U!(u::Gaugefields_4D_wing{NC}) where NC
        NT = u.NT
        NY = u.NY
        NZ = u.NZ
        NX = u.NX
        NDW = u.NDW
    
        #X direction 
        #Now we send data
    
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for id=1:NDW
                        for k2=1:NC
                            @simd for k1=1:NC
                                @inbounds u[k1,k2,-NDW+id,iy,iz,it] = u[k1,k2,NX+(id-NDW),iy,iz,it]
                            end
                        end
                    end
                end
            end
        end
    
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for id=1:NDW
                        for k2=1:NC
                            @simd for k1=1:NC
                                @inbounds u[k1,k2,NX+id,iy,iz,it] = u[k1,k2,id,iy,iz,it]
                            end
                        end
                    end
                end
            end
        end
    
    
        #Y direction 
        #Now we send data
        for it=1:NT
            for iz=1:NZ
                for ix=-NDW+1:NX+NDW
                    for id=1:NDW
                        for k1=1:NC
                            @simd for k2=1:NC
                                @inbounds u[k1,k2,ix,-NDW+id,iz,it] = u[k1,k2,ix,NY+(id-NDW),iz,it]
                            end
                        end
                    end
                end
            end
        end
    
        for it=1:NT
            for iz=1:NZ
                for ix=-NDW+1:NX+NDW
                    for id=1:NDW
                        for k1=1:NC
                            @simd for k2=1:NC
                                @inbounds u[k1,k2,ix,NY+id,iz,it] = u[k1,k2,ix,id,iz,it]
                            end
                        end
                    end
                end
            end
        end
    
        #Z direction 
        #Now we send data
        for id=1:NDW
            for it=1:NT
                for iy=-NDW+1:NY+NDW
                    for ix=-NDW+1:NX+NDW
                        for k1=1:NC
                            @simd for k2=1:NC
                                @inbounds u[k1,k2,ix,iy,id-NDW,it] = u[k1,k2,ix,iy,NZ+(id-NDW),it]
                                @inbounds u[k1,k2,ix,iy,NZ+id,it] = u[k1,k2,ix,iy,id,it]
                            end
                        end
                    end
                end
            end
        end
    
    
        for id=1:NDW
            for iz=-NDW+1:NZ+NDW
                for iy=-NDW+1:NY+NDW
                    for ix=-NDW+1:NX+NDW
                        for k1=1:NC
                            @simd for k2=1:NC
                                @inbounds u[k1,k2,ix,iy,iz,id-NDW] = u[k1,k2,ix,iy,iz,NT+(id-NDW)]
                                @inbounds u[k1,k2,ix,iy,iz,NT+id] = u[k1,k2,ix,iy,iz,id]
                            end
                        end
                    end
                end
            end
        end
    
        #display(u.g)
        #exit()
    
        return
    end

    function exptU!(uout::T,t::N,u::Gaugefields_4D_wing{NC},temps::Array{T,1}) where {N <: Number, T <: Gaugefields_4D_wing, NC} #uout = exp(t*u)
        @assert NC != 3 && NC != 2 "This function is for NC != 2,3"
        
        
        NT = u.NT
        NZ = u.NZ
        NY = u.NY
        NX = u.NX
        V0 = zeros(ComplexF64,NC,NC)
        V1 = zeros(ComplexF64,NC,NC)
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k2=1:NC
                            for k1=1:NC  
                                @inbounds V0[k1,k2] = im*t*u[k1,k2,ix,iy,iz,it]
                            end
                        end
                        V1[:,:] = exp(V0)
                        for k2=1:NC
                            for k1=1:NC  
                                @inbounds uout[k1,k2,ix,iy,iz,it] = V1[k1,k2]
                            end
                        end

                    end
                end
            end
        end
        #error("exptU! is not implemented in type $(typeof(u)) ")
    end

    const tinyvalue =1e-100
    const pi23 = 2pi/3
    const fac13 = 1/3

    # #=
    function exptU!(uout::T,t::N,v::Gaugefields_4D_wing{3},temps::Array{T,1}) where {N <: Number, T <: Gaugefields_4D_wing} #uout = exp(t*u)
        ww = temps[1]
        w = temps[2]
        

        NT = v.NT
        NZ = v.NZ
        NY = v.NY
        NX = v.NX
        #t = 1

        @inbounds for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        v11 = v[1,1,ix,iy,iz,it]
                        v22 = v[2,2,ix,iy,iz,it]
                        v33 = v[3,3,ix,iy,iz,it]

                        tri = fac13*(imag(v11)+imag(v22)+imag(v33))

                        #=
                        vout[1,1,ix,iy,iz,it] = (imag(v11)-tri)*im
                        vout[2,2,ix,iy,iz,it] = (imag(v22)-tri)*im
                        vout[3,3,ix,iy,iz,it] = (imag(v33)-tri)*im
                        =#
                        y11 = (imag(v11)-tri)*im
                        y22 = (imag(v22)-tri)*im
                        y33 = (imag(v33)-tri)*im

                        v12 = v[1,2,ix,iy,iz,it]
                        v13 = v[1,3,ix,iy,iz,it]
                        v21 = v[2,1,ix,iy,iz,it]
                        v23 = v[2,3,ix,iy,iz,it]
                        v31 = v[3,1,ix,iy,iz,it]
                        v32 = v[3,2,ix,iy,iz,it]

                        x12 = v12 - conj(v21)
                        x13 = v13 - conj(v31)
                        x23 = v23 - conj(v32)
                    
                        x21 = - conj(x12)
                        x31 = - conj(x13)
                        x32 = - conj(x23)

                        y12 = 0.5  * x12
                        y13 = 0.5  * x13
                        y21 = 0.5  * x21
                        y23 = 0.5  * x23
                        y31 = 0.5  * x31
                        y32 =  0.5  * x32

                        c1_0 = ( imag(y12) + imag(y21) )
                        c2_0 = ( real(y12) - real(y21) )
                        c3_0 = ( imag(y11) - imag(y22) )
                        c4_0= ( imag(y13) + imag(y31) )
                        c5_0 = ( real(y13) - real(y31) )
                        
                        c6_0 = ( imag(y23) + imag(y32) )
                        c7_0 = ( real(y23) - real(y32) )
                        c8_0 = sr3i *
                                ( imag(y11) + imag(y22) -
                                        2*imag(y33) )

                        c1 = t*c1_0 * 0.5 
                        c2 = t*c2_0 * 0.5
                        c3 = t*c3_0 * 0.5
                        c4 = t*c4_0 * 0.5
                        c5 = t*c5_0 * 0.5
                        c6 = t*c6_0 * 0.5
                        c7 = t*c7_0 * 0.5
                        c8 = t*c8_0 * 0.5
                        csum = c1+c2+c3+c4+c5+c6+c7+c8
                        if csum == 0
                            w[1,1,ix,iy,iz,it]=   1 
                            w[1,2,ix,iy,iz,it]=   0
                            w[1,3,ix,iy,iz,it]=   0 
                            w[2,1,ix,iy,iz,it]=  0
                            w[2,2,ix,iy,iz,it]=   1 
                            w[2,3,ix,iy,iz,it]=   0
                            w[3,1,ix,iy,iz,it]=   0 
                            w[3,2,ix,iy,iz,it]=   0  
                            w[3,3,ix,iy,iz,it]=   1  
                    
                            ww[1,1,ix,iy,iz,it]=   1
                            ww[1,2,ix,iy,iz,it]=   0
                            ww[1,3,ix,iy,iz,it]=   0
                            ww[2,1,ix,iy,iz,it]=   0
                            ww[2,2,ix,iy,iz,it]=   1 
                            ww[2,3,ix,iy,iz,it]=   0 
                            ww[3,1,ix,iy,iz,it]=   0
                            ww[3,2,ix,iy,iz,it]=   0  
                            ww[3,3,ix,iy,iz,it]=   1
                            continue
                        end


                        #x[1,1,icum] =  c3+sr3i*c8 +im*(  0.0 )
                        v1 = c3+sr3i*c8
                        v2 = 0.0
                        #x[1,2,icum] =  c1         +im*( -c2   )
                        v3 = c1
                        v4 = -c2
                        #x[1,3,icum] =  c4         +im*(-c5   )
                        v5 = c4
                        v6 = -c5
                    
                        #x[2,1,icum] =  c1         +im*(  c2   )
                        v7 = c1
                        v8 = c2

                        #x[2,2,icum] =  -c3+sr3i*c8+im*(  0.0 )
                        v9 =-c3+sr3i*c8
                        v10 = 0.0

                        #x[2,3,icum] =  c6         +im*( -c7   )
                        v11 = c6
                        v12 = -c7
                    
                        #x[3,1,icum] =  c4         +im*(  c5   )
                        v13 = c4
                        v14 = c5

                        #x[3,2,icum] =  c6         +im*(  c7   )
                        v15 = c6 
                        v16 = c7
                        #x[3,3,icum] =  -sr3i2*c8  +im*(  0.0 )
                        v17 = -sr3i2*c8
                        v18 = 0.0

                
                #c find eigenvalues of v
                        trv3 = (v1 + v9 + v17) / 3.0 
                        cofac = v1 *  v9  - v3^2  - v4^2 + 
                                v1 * v17 -  v5^2  - v6^2 + 
                                v9 * v17 - v11^2 - v12^2
                        det = v1 * v9 * v17 -
                                v1 * (v11^2 + v12^2) -
                                v9 * (v5^2 + v6^2) -
                                v17 * (v3^2 + v4^2) +
                                (v5 * (v3 * v11 - v4 * v12) +
                                v6 * (v3 * v12 + v4 * v11)) * 2.0
                        p3 = cofac / 3.0 - trv3^2
                        q = trv3 * cofac - det - 2.0 * trv3^3
                        x = sqrt(-4.0 * p3) + tinyvalue
                        arg = q / (x * p3)

                        arg = min(1, max(-1, arg))
                        theta = acos(arg) / 3.0
                        e1 = x * cos(theta) + trv3
                        theta = theta + pi23
                        e2 = x * cos(theta) + trv3
                #       theta = theta + pi23
                #       e3 = x * cos(theta) + trv3
                        e3 = 3.0 * trv3 - e1 - e2
                
                # solve for eigenvectors
                
                        w1 =   v5 * (v9 - e1) - v3 * v11 + v4 * v12
                        w2 = - v6 * (v9 - e1) + v4 * v11 + v3 * v12
                        w3 =   (v1 - e1) * v11 - v3 * v5 - v4 * v6
                        w4 = - (v1 - e1) * v12 - v4 * v5 + v3 * v6
                        w5 = - (v1 - e1) * (v9 - e1) +  v3^2 + v4^2
                        w6 = 0.0
                
                        coeff = 1.0 / sqrt(w1^2 + w2^2 + w3^2 + w4^2 + w5^2)
                        
                
                        w1 = w1 * coeff
                        w2 = w2 * coeff
                        w3 = w3 * coeff
                        w4 = w4 * coeff
                        w5 = w5 * coeff
                
                        w7 =   v5 * (v9 - e2) - v3 * v11 + v4 * v12
                        w8 = - v6 * (v9 - e2) + v4 * v11 + v3 * v12
                        w9 = (v1 - e2) * v11 - v3 * v5 - v4 * v6
                        w10 = - (v1 - e2) * v12 - v4 * v5 + v3 * v6
                        w11 = - (v1 - e2) * (v9 - e2) +  v3^2 + v4^2
                        w12 = 0.0
                
                        coeff = 1.0 / sqrt(w7^2  + w8^2 + w9^2+ w10^2 + w11^2)
                        
                        w7  = w7  * coeff
                        w8  = w8  * coeff
                        w9  = w9  * coeff
                        w10 = w10 * coeff
                        w11 = w11 * coeff
                
                        w13 =   v5 * (v9 - e3) - v3 * v11 + v4 * v12
                        w14 = - v6 * (v9 - e3) + v4 * v11 + v3 * v12
                        w15 =   (v1 - e3) * v11 - v3 * v5 - v4 * v6
                        w16 = - (v1 - e3) * v12 - v4 * v5 + v3 * v6
                        w17 = - (v1 - e3) * (v9 - e3) +  v3^2 + v4^2
                        w18 = 0.0
                
                        coeff = 1.0 / sqrt(w13^2 + w14^2 + w15^2+ w16^2 + w17^2)
                        w13 = w13 * coeff
                        w14 = w14 * coeff
                        w15 = w15 * coeff
                        w16 = w16 * coeff
                        w17 = w17 * coeff
                
                # construct the projection v
                        c1 = cos(e1)
                        s1 = sin(e1)
                        ww1  = w1  * c1 - w2 * s1
                        ww2  = w2  * c1 + w1 * s1
                        ww3  = w3  * c1 - w4 * s1
                        ww4  = w4  * c1 + w3 * s1
                        ww5  = w5  * c1 - w6 * s1
                        ww6  = w6  * c1 + w5 * s1
                
                        c2 = cos(e2)
                        s2 = sin(e2)
                        ww7  = w7  * c2 - w8 * s2
                        ww8  = w8  * c2 + w7 * s2
                        ww9  = w9  * c2 - w10 * s2
                        ww10 = w10 * c2 + w9 * s2
                        ww11 = w11 * c2 - w12 * s2
                        ww12 = w12 * c2 + w11 * s2
                
                        c3 = cos(e3)
                        s3 = sin(e3)
                        ww13 = w13 * c3 - w14 * s3
                        ww14 = w14 * c3 + w13 * s3
                        ww15 = w15 * c3 - w16 * s3
                        ww16 = w16 * c3 + w15 * s3
                        ww17 = w17 * c3 - w18 * s3
                        ww18 = w18 * c3 + w17 * s3
                
                        w[1,1,ix,iy,iz,it]=   w1+im*w2 
                        w[1,2,ix,iy,iz,it]=   w3+im*w4 
                        w[1,3,ix,iy,iz,it]=   w5+im*w6 
                        w[2,1,ix,iy,iz,it]=   w7+im*w8 
                        w[2,2,ix,iy,iz,it]=   w9+im*w10  
                        w[2,3,ix,iy,iz,it]=   w11+im*w12  
                        w[3,1,ix,iy,iz,it]=   w13+im*w14  
                        w[3,2,ix,iy,iz,it]=   w15+im*w16  
                        w[3,3,ix,iy,iz,it]=   w17+im*w18  
                
                        ww[1,1,ix,iy,iz,it]=   ww1+im*ww2 
                        ww[1,2,ix,iy,iz,it]=   ww3+im*ww4 
                        ww[1,3,ix,iy,iz,it]=   ww5+im*ww6 
                        ww[2,1,ix,iy,iz,it]=   ww7+im*ww8 
                        ww[2,2,ix,iy,iz,it]=   ww9+im*ww10  
                        ww[2,3,ix,iy,iz,it]=   ww11+im*ww12  
                        ww[3,1,ix,iy,iz,it]=   ww13+im*ww14  
                        ww[3,2,ix,iy,iz,it]=   ww15+im*ww16  
                        ww[3,3,ix,iy,iz,it]=   ww17+im*ww18  

                        #a = ww[:,:,ix,iy,iz,it]
                        #b = w[:,:,ix,iy,iz,it]
                        #println(b'*a)
                        #println(exp(im*t*v[:,:,ix,iy,iz,it]))
                        #error("d")
                    end
                end
            end
        end

        #mul!(v,w',ww)
        mul!(uout,w',ww)

        #error("exptU! is not implemented in type $(typeof(u)) ")
    end
    # =#


    """
-----------------------------------------------------c
     !!!!!   vin and vout should be different vectors

     Projectin of the etraceless antiermite part 
     vout = x/2 - Tr(x)/6
     wher   x = vin - Conjg(vin)      
-----------------------------------------------------c
    """

    #Q = -(1/2)*(Ω' - Ω) + (1/(2NC))*tr(Ω' - Ω)*I0_2
    #Omega' - Omega = -2i imag(Omega)
    function Traceless_antihermitian!(vout::Gaugefields_4D_wing{3},vin::Gaugefields_4D_wing{3}) 
        #error("Traceless_antihermitian! is not implemented in type $(typeof(vout)) ")
        fac13 = 1/3
        NX = vin.NX
        NY = vin.NY
        NZ = vin.NZ
        NT = vin.NT
    
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX
                        v11 = vin[1,1,ix,iy,iz,it]
                        v21 = vin[2,1,ix,iy,iz,it]
                        v31 = vin[3,1,ix,iy,iz,it]

                        v12 = vin[1,2,ix,iy,iz,it]
                        v22 = vin[2,2,ix,iy,iz,it]
                        v32 = vin[3,2,ix,iy,iz,it]

                        v13 = vin[1,3,ix,iy,iz,it]
                        v23 = vin[2,3,ix,iy,iz,it]
                        v33 = vin[3,3,ix,iy,iz,it]

    
                        tri = fac13*(imag(v11)+imag(v22)+imag(v33))
    
                        #=
                        vout[1,1,ix,iy,iz,it] = (imag(v11)-tri)*im
                        vout[2,2,ix,iy,iz,it] = (imag(v22)-tri)*im
                        vout[3,3,ix,iy,iz,it] = (imag(v33)-tri)*im
                        =#
                        y11 = (imag(v11)-tri)*im
                        y22 = (imag(v22)-tri)*im
                        y33 = (imag(v33)-tri)*im
    
                        
    
                        x12 = v12 - conj(v21)
                        x13 = v13 - conj(v31)
                        x23 = v23 - conj(v32)
                    
                        x21 = - conj(x12)
                        x31 = - conj(x13)
                        x32 = - conj(x23)
    
                        #=
                        vout[1,2,ix,iy,iz,it] = 0.5  * x12
                        vout[1,3,ix,iy,iz,it] = 0.5  * x13
                        vout[2,1,ix,iy,iz,it] = 0.5  * x21
                        vout[2,3,ix,iy,iz,it] = 0.5  * x23
                        vout[3,1,ix,iy,iz,it] = 0.5  * x31
                        vout[3,2,ix,iy,iz,it] = 0.5  * x32
                        =#
                        y12 = 0.5  * x12
                        y13 = 0.5  * x13
                        y21 = 0.5  * x21
                        y23 = 0.5  * x23
                        y31 = 0.5  * x31
                        y32 =  0.5  * x32

                        
                        
                        
                        
                        
                        

                        vout[1,1,ix,iy,iz,it] = y11
                        vout[2,1,ix,iy,iz,it] = y21
                        vout[3,1,ix,iy,iz,it] = y31

                        vout[1,2,ix,iy,iz,it] = y12
                        vout[2,2,ix,iy,iz,it] = y22
                        vout[3,2,ix,iy,iz,it] = y32

                        vout[1,3,ix,iy,iz,it] = y13
                        vout[2,3,ix,iy,iz,it] = y23
                        vout[3,3,ix,iy,iz,it] = y33


                    end
                end
            end
        end
    
    
    end


    #=
    function Traceless_antihermitian!(vout::Gaugefields_4D_wing{3},vin::Gaugefields_4D_wing{3})
        #error("Traceless_antihermitian! is not implemented in type $(typeof(vout)) ")
        fac13 = 1/3
        NX = vin.NX
        NY = vin.NY
        NZ = vin.NZ
        NT = vin.NT

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX
                        v11 = vin[1,1,ix,iy,iz,it]
                        v22 = vin[2,2,ix,iy,iz,it]
                        v33 = vin[3,3,ix,iy,iz,it]

                        tri = fac13*(imag(v11)+imag(v22)+imag(v33))

                        vout[1,1,ix,iy,iz,it] = (imag(v11)-tri)*im
                        vout[2,2,ix,iy,iz,it] = (imag(v22)-tri)*im
                        vout[3,3,ix,iy,iz,it] = (imag(v33)-tri)*im
                    end
                end
            end
        end

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX

                        v12 = vin[1,2,ix,iy,iz,it]
                        v13 = vin[1,3,ix,iy,iz,it]
                        v21 = vin[2,1,ix,iy,iz,it]
                        v23 = vin[2,3,ix,iy,iz,it]
                        v31 = vin[3,1,ix,iy,iz,it]
                        v32 = vin[3,2,ix,iy,iz,it]

                        x12 = v12 - conj(v21)
                        x13 = v13 - conj(v31)
                        x23 = v23 - conj(v32)
                    
                        x21 = - conj(x12)
                        x31 = - conj(x13)
                        x32 = - conj(x23)

                        vout[1,2,ix,iy,iz,it] = 0.5  * x12
                        vout[1,3,ix,iy,iz,it] = 0.5  * x13
                        vout[2,1,ix,iy,iz,it] = 0.5  * x21
                        vout[2,3,ix,iy,iz,it] = 0.5  * x23
                        vout[3,1,ix,iy,iz,it] = 0.5  * x31
                        vout[3,2,ix,iy,iz,it] = 0.5  * x32
                    end
                end
            end
        end

    end
    =#

    function Traceless_antihermitian!(vout::Gaugefields_4D_wing{2},vin::Gaugefields_4D_wing{2})
        #error("Traceless_antihermitian! is not implemented in type $(typeof(vout)) ")
        fac13 = 1/3
        NX = vin.NX
        NY = vin.NY
        NZ = vin.NZ
        NT = vin.NT


        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX

                        v11 = vin[1,1,ix,iy,iz,it]
                        v22 = vin[2,2,ix,iy,iz,it]

                        tri = (imag(v11)+imag(v22))*0.5

                        

                        v12 = vin[1,2,ix,iy,iz,it]
                        #v13 = vin[1,3,ix,iy,iz,it]
                        v21 = vin[2,1,ix,iy,iz,it]

                        x12 = v12 - conj(v21)

                        x21 = - conj(x12)

                        vout[1,1,ix,iy,iz,it] = (imag(v11)-tri)*im
                        vout[1,2,ix,iy,iz,it] = 0.5  * x12
                        vout[2,1,ix,iy,iz,it] = 0.5  * x21
                        vout[2,2,ix,iy,iz,it] = (imag(v22)-tri)*im
                    end
                end
            end
        end

    end


    function Traceless_antihermitian!(vout::Gaugefields_4D_wing{NC},vin::Gaugefields_4D_wing{NC}) where NC
        #NC = vout.NC
        fac1N = 1/NC
        nv = vin.NV

        NX = vin.NX
        NY = vin.NY
        NZ = vin.NZ
        NT = vin.NT

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX
                        tri = 0.0
                        @simd for k=1:NC
                            tri += imag(vin[k,k,ix,iy,iz,it])
                        end
                        tri *= fac1N
                        @simd for k=1:NC
                            vout[k,k,ix,iy,iz,it] = (imag(vin[k,k,ix,iy,iz,it])-tri)*im
                        end
                    end
                end
            end
        end


        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX
                        for k1=1:NC
                            @simd for k2=k1+1:NC
                                vv = 0.5*(vin[k1,k2,ix,iy,iz,it] - conj(vin[k2,k1,ix,iy,iz,it]))
                                vout[k1,k2,ix,iy,iz,it] = vv
                                vout[k2,k1,ix,iy,iz,it] = -conj(vv)
                            end
                        end
                    end
                end
            end
        end
            
        
    end

    function LinearAlgebra.tr(a::Gaugefields_4D_wing{NC},b::Gaugefields_4D_wing{NC}) where NC
        NX=a.NX
        NY=a.NY
        NZ=a.NZ
        NT=a.NT

        s = 0
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k=1:NC
                            for k2=1:NC
                                s += a[k,k2,ix,iy,iz,it]*b[k2,k,ix,iy,iz,it]
                            end
                        end
                    end
                end
            end
        end
        #println(3*NT*NZ*NY*NX*NC)
        return s
    end

    function LinearAlgebra.tr(a::Gaugefields_4D_wing{2},b::Gaugefields_4D_wing{2})
        NX=a.NX
        NY=a.NY
        NZ=a.NZ
        NT=a.NT

        s = 0
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        a11 = a[1,1,ix,iy,iz,it]
                        a21 = a[2,1,ix,iy,iz,it]
                        a12 = a[1,2,ix,iy,iz,it]
                        a22 = a[2,2,ix,iy,iz,it]

                        b11 = b[1,1,ix,iy,iz,it]
                        b21 = b[2,1,ix,iy,iz,it]
                        b12 = b[1,2,ix,iy,iz,it]
                        b22 = b[2,2,ix,iy,iz,it]

                        s += a11*b11+a12*b21+
                            a21*b12+a22*b22

                    end
                end
            end
        end
        #println(3*NT*NZ*NY*NX*NC)
        return s
    end

    function LinearAlgebra.tr(a::Gaugefields_4D_wing{3},b::Gaugefields_4D_wing{3}) 
        NX=a.NX
        NY=a.NY
        NZ=a.NZ
        NT=a.NT

        s = 0
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        a11 = a[1,1,ix,iy,iz,it]
                        a21 = a[2,1,ix,iy,iz,it]
                        a31 = a[3,1,ix,iy,iz,it]
                        a12 = a[1,2,ix,iy,iz,it]
                        a22 = a[2,2,ix,iy,iz,it]
                        a32 = a[3,2,ix,iy,iz,it]
                        a13 = a[1,3,ix,iy,iz,it]
                        a23 = a[2,3,ix,iy,iz,it]
                        a33 = a[3,3,ix,iy,iz,it]
                        b11 = b[1,1,ix,iy,iz,it]
                        b21 = b[2,1,ix,iy,iz,it]
                        b31 = b[3,1,ix,iy,iz,it]
                        b12 = b[1,2,ix,iy,iz,it]
                        b22 = b[2,2,ix,iy,iz,it]
                        b32 = b[3,2,ix,iy,iz,it]
                        b13 = b[1,3,ix,iy,iz,it]
                        b23 = b[2,3,ix,iy,iz,it]
                        b33 = b[3,3,ix,iy,iz,it]
                        s += a11*b11+a12*b21+a13*b31 + 
                            a21*b12+a22*b22+a23*b32 + 
                            a31*b13+a32*b23+a33*b33
                    end
                end
            end
        end
        #println(3*NT*NZ*NY*NX*NC)
        return s
    end



    function LinearAlgebra.tr(a::Gaugefields_4D_wing{NC}) where NC
        NX=a.NX
        NY=a.NY
        NZ=a.NZ
        NT=a.NT

        s = 0
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        @simd for k=1:NC
                            s += a[k,k,ix,iy,iz,it]
                            #println(a[k,k,ix,iy,iz,it])
                        end
                    end
                end
            end
        end
        #println(3*NT*NZ*NY*NX*NC)
        return s
    end

    function partial_tr(a::Gaugefields_4D_wing{NC},μ) where NC
        NX=a.NX
        NY=a.NY
        NZ=a.NZ
        NT=a.NT

        if μ == 1
            s = 0
            ix = 1
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        #for ix=1:NX
                            @simd for k=1:NC
                                s += a[k,k,ix,iy,iz,it]
                                #println(a[k,k,ix,iy,iz,it])
                            end
                        
                        #end
                    end
                end
            end
        elseif μ == 2
            s = 0
            iy =1
            for it=1:NT
                for iz=1:NZ
                    #for iy=1:NY
                        for ix=1:NX
                            @simd for k=1:NC
                                s += a[k,k,ix,iy,iz,it]
                                #println(a[k,k,ix,iy,iz,it])
                            end
                        end
                    #end
                end
            end
        elseif μ == 3
            s = 0
            iz = 1
            for it=1:NT
                #for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            @simd for k=1:NC
                                s += a[k,k,ix,iy,iz,it]
                                #println(a[k,k,ix,iy,iz,it])
                            end
                        end
                    end
                #end
            end
        else 
            s = 0
            it = 1
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            @simd for k=1:NC
                                s += a[k,k,ix,iy,iz,it]
                                #println(a[k,k,ix,iy,iz,it])
                            end
                        end
                    end
                end
            
        end
        


        #println(3*NT*NZ*NY*NX*NC)
        return s
    end

    function add_U!(c::Gaugefields_4D_wing{NC},a::T1) where {NC,T1 <: Abstractfields}
        @inbounds for i=1:length(c.U)
            c.U[i] += a.U[i]
        end
        return 


        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k2=1:NC                            
                            @simd for k1=1:NC
                                c[k1,k2,ix,iy,iz,it] += a[k1,k2,ix,iy,iz,it]
                            end
                        end
                    end
                end
            end
        end
    end

    function add_U!(c::Gaugefields_4D_wing{NC},a::T1,iseven::Bool) where {NC,T1 <: Abstractfields}
        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        evenodd = ifelse( (ix+iy+iz+it) % 2 ==0, true,false)
                        if evenodd == iseven
                            for k2=1:NC                            
                                @simd for k1=1:NC
                                    c[k1,k2,ix,iy,iz,it] += a[k1,k2,ix,iy,iz,it]
                                end
                            end
                        end
                    end
                end
            end
        end
    end


    function add_U!(c::Gaugefields_4D_wing{NC},α::N,a::T1) where {NC,T1 <: Abstractfields, N<:Number}
        #@inbounds for i=1:length(c.U)
        #    c.U[i] += α*a.U[i]
        #end
        #return 

        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        @inbounds for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k2=1:NC                            
                            @simd for k1=1:NC
                                c[k1,k2,ix,iy,iz,it] += α*a[k1,k2,ix,iy,iz,it]
                            end
                        end
                    end
                end
            end
        end
    end


    
    
    function LinearAlgebra.mul!(c::Gaugefields_4D_wing{NC},a::T1,b::T2) where {NC,T1 <: Abstractfields,T2 <: Abstractfields}
        @assert NC != 2 && NC != 3 "This function is for NC != 2,3"
        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k2=1:NC                            
                            for k1=1:NC
                                c[k1,k2,ix,iy,iz,it] = 0

                                @simd for k3=1:NC
                                    c[k1,k2,ix,iy,iz,it] += a[k1,k3,ix,iy,iz,it]*b[k3,k2,ix,iy,iz,it]
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    function LinearAlgebra.mul!(c::Gaugefields_4D_wing{NC},a::T1,b::T2,iseven::Bool) where {NC,T1 <: Abstractfields,T2 <: Abstractfields}
        @assert NC != 2 && NC != 3 "This function is for NC != 2,3"
        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        @inbounds for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        evenodd = ifelse( (ix+iy+iz+it) % 2 ==0, true,false)
                        if evenodd == iseven
                            for k2=1:NC                            
                                for k1=1:NC
                                    c[k1,k2,ix,iy,iz,it] = 0
                                    @simd for k3=1:NC
                                        c[k1,k2,ix,iy,iz,it] += a[k1,k3,ix,iy,iz,it]*b[k3,k2,ix,iy,iz,it]
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    function LinearAlgebra.mul!(c::Gaugefields_4D_wing{NC},a::T1,b::T2) where {NC,T1 <: Number,T2 <: Abstractfields}
        @assert NC != 2 && NC != 3 "This function is for NC != 2,3"
        @inbounds for i=1:length(c)
            c.U[i] = a*b.U[i]
        end
        return


        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k2=1:NC                            
                            @simd for k1=1:NC
                                c[k1,k2,ix,iy,iz,it] = a*b[k1,k2,ix,iy,iz,it]
                            end
                        end
                    end
                end
            end
        end
    end

    function LinearAlgebra.mul!(c::Gaugefields_4D_wing{NC},a::T1,b::T2,α::Ta,β::Tb) where {NC,T1 <: Abstractfields,T2 <: Abstractfields,Ta <: Number, Tb <: Number}
        @assert NC != 2 && NC != 3 "This function is for NC != 2,3"
        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k2=1:NC                            
                            for k1=1:NC
                                c[k1,k2,ix,iy,iz,it] = β*c[k1,k2,ix,iy,iz,it] 
                                @simd for k3=1:NC
                                    c[k1,k2,ix,iy,iz,it] += α*a[k1,k3,ix,iy,iz,it]*b[k3,k2,ix,iy,iz,it] 
                                end
                            end
                        end
                    end
                end
            end
        end
    end


#=
evenodd = ifelse( (ix+iy+iz+it) % 2 ==0, true,false)
                        if evenodd == iseven
=#

    function LinearAlgebra.mul!(c::Gaugefields_4D_wing{2},a::T1,b::T2,iseven::Bool) where {T1 <: Abstractfields,T2 <: Abstractfields}
        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        #println("threads = ", Threads.nthreads())
        @inbounds  for it=1:NT#,iz=1:NZ,iy=1:NY
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX
                        evenodd = ifelse( (ix+iy+iz+it) % 2 ==0, true,false)
                        if evenodd == iseven
                            a11 = a[1,1,ix,iy,iz,it]
                            a21 = a[2,1,ix,iy,iz,it]
                            a12 = a[1,2,ix,iy,iz,it]
                            a22 = a[2,2,ix,iy,iz,it]

                            b11 = b[1,1,ix,iy,iz,it]
                            b21 = b[2,1,ix,iy,iz,it]
                            b12 = b[1,2,ix,iy,iz,it]
                            b22 = b[2,2,ix,iy,iz,it]

                            c[1,1,ix,iy,iz,it] = a11*b11+a12*b21
                            c[2,1,ix,iy,iz,it] = a21*b11+a22*b21
                            c[1,2,ix,iy,iz,it] = a11*b12+a12*b22
                            c[2,2,ix,iy,iz,it] = a21*b12+a22*b22
                        end

                    end
                end
            end
        end
    end


    
    function LinearAlgebra.mul!(c::Gaugefields_4D_wing{2},a::T1,b::T2) where {T1 <: Abstractfields,T2 <: Abstractfields}
        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        #println("threads = ", Threads.nthreads())
        @inbounds  for it=1:NT#,iz=1:NZ,iy=1:NY
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX
                        a11 = a[1,1,ix,iy,iz,it]
                        a21 = a[2,1,ix,iy,iz,it]
                        a12 = a[1,2,ix,iy,iz,it]
                        a22 = a[2,2,ix,iy,iz,it]

                        b11 = b[1,1,ix,iy,iz,it]
                        b21 = b[2,1,ix,iy,iz,it]
                        b12 = b[1,2,ix,iy,iz,it]
                        b22 = b[2,2,ix,iy,iz,it]

                        c[1,1,ix,iy,iz,it] = a11*b11+a12*b21
                        c[2,1,ix,iy,iz,it] = a21*b11+a22*b21
                        c[1,2,ix,iy,iz,it] = a11*b12+a12*b22
                        c[2,2,ix,iy,iz,it] = a21*b12+a22*b22

                    end
                end
            end
        end
    end
    

    
    function LinearAlgebra.mul!(c::Gaugefields_4D_wing{2},a::T1,b::T2,α::Ta,β::Tb) where {NC,T1 <: Abstractfields,T2 <: Abstractfields,Ta <: Number, Tb <: Number}
        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        if β == zero(β)
            if α == one(α)
                mul!(c,a,b)
                return
            end
        end

        if α != 0
            @inbounds for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        @simd for ix=1:NX
                            a11 = a[1,1,ix,iy,iz,it]
                            a21 = a[2,1,ix,iy,iz,it]
                            a12 = a[1,2,ix,iy,iz,it]
                            a22 = a[2,2,ix,iy,iz,it]

                            b11 = b[1,1,ix,iy,iz,it]
                            b21 = b[2,1,ix,iy,iz,it]
                            b12 = b[1,2,ix,iy,iz,it]
                            b22 = b[2,2,ix,iy,iz,it]

                            c[1,1,ix,iy,iz,it] = (a11*b11+a12*b21)*α + β*c[1,1,ix,iy,iz,it]
                            c[2,1,ix,iy,iz,it] = (a21*b11+a22*b21)*α + β*c[2,1,ix,iy,iz,it]
                            c[1,2,ix,iy,iz,it] = (a11*b12+a12*b22)*α + β*c[1,2,ix,iy,iz,it]
                            c[2,2,ix,iy,iz,it] = (a21*b12+a22*b22)*α + β*c[2,2,ix,iy,iz,it]



                        end
                    end
                end
            end
        end
    end

    
    function LinearAlgebra.mul!(c::Gaugefields_4D_wing{3},a::T1,b::T2) where {T1 <: Abstractfields,T2 <: Abstractfields}
        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        #println("threads = ", Threads.nthreads())
        @inbounds  for it=1:NT#,iz=1:NZ,iy=1:NY
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX
                        a11 = a[1,1,ix,iy,iz,it]
                        a21 = a[2,1,ix,iy,iz,it]
                        a31 = a[3,1,ix,iy,iz,it]
                        a12 = a[1,2,ix,iy,iz,it]
                        a22 = a[2,2,ix,iy,iz,it]
                        a32 = a[3,2,ix,iy,iz,it]
                        a13 = a[1,3,ix,iy,iz,it]
                        a23 = a[2,3,ix,iy,iz,it]
                        a33 = a[3,3,ix,iy,iz,it]
                        b11 = b[1,1,ix,iy,iz,it]
                        b21 = b[2,1,ix,iy,iz,it]
                        b31 = b[3,1,ix,iy,iz,it]
                        b12 = b[1,2,ix,iy,iz,it]
                        b22 = b[2,2,ix,iy,iz,it]
                        b32 = b[3,2,ix,iy,iz,it]
                        b13 = b[1,3,ix,iy,iz,it]
                        b23 = b[2,3,ix,iy,iz,it]
                        b33 = b[3,3,ix,iy,iz,it]
                        c[1,1,ix,iy,iz,it] = a11*b11+a12*b21+a13*b31
                        c[2,1,ix,iy,iz,it] = a21*b11+a22*b21+a23*b31
                        c[3,1,ix,iy,iz,it] = a31*b11+a32*b21+a33*b31
                        c[1,2,ix,iy,iz,it] = a11*b12+a12*b22+a13*b32
                        c[2,2,ix,iy,iz,it] = a21*b12+a22*b22+a23*b32
                        c[3,2,ix,iy,iz,it] = a31*b12+a32*b22+a33*b32
                        c[1,3,ix,iy,iz,it] = a11*b13+a12*b23+a13*b33
                        c[2,3,ix,iy,iz,it] = a21*b13+a22*b23+a23*b33
                        c[3,3,ix,iy,iz,it] = a31*b13+a32*b23+a33*b33

                    end
                end
            end
        end
    end

    function LinearAlgebra.mul!(c::Gaugefields_4D_wing{3},a::T1,b::T2,iseven::Bool) where {T1 <: Abstractfields,T2 <: Abstractfields}
        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        #println("threads = ", Threads.nthreads())
        @inbounds  for it=1:NT#,iz=1:NZ,iy=1:NY
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX
                        evenodd = ifelse( (ix+iy+iz+it) % 2 ==0, true,false)
                        if evenodd == iseven

                            a11 = a[1,1,ix,iy,iz,it]
                            a21 = a[2,1,ix,iy,iz,it]
                            a31 = a[3,1,ix,iy,iz,it]
                            a12 = a[1,2,ix,iy,iz,it]
                            a22 = a[2,2,ix,iy,iz,it]
                            a32 = a[3,2,ix,iy,iz,it]
                            a13 = a[1,3,ix,iy,iz,it]
                            a23 = a[2,3,ix,iy,iz,it]
                            a33 = a[3,3,ix,iy,iz,it]
                            b11 = b[1,1,ix,iy,iz,it]
                            b21 = b[2,1,ix,iy,iz,it]
                            b31 = b[3,1,ix,iy,iz,it]
                            b12 = b[1,2,ix,iy,iz,it]
                            b22 = b[2,2,ix,iy,iz,it]
                            b32 = b[3,2,ix,iy,iz,it]
                            b13 = b[1,3,ix,iy,iz,it]
                            b23 = b[2,3,ix,iy,iz,it]
                            b33 = b[3,3,ix,iy,iz,it]
                            c[1,1,ix,iy,iz,it] = a11*b11+a12*b21+a13*b31
                            c[2,1,ix,iy,iz,it] = a21*b11+a22*b21+a23*b31
                            c[3,1,ix,iy,iz,it] = a31*b11+a32*b21+a33*b31
                            c[1,2,ix,iy,iz,it] = a11*b12+a12*b22+a13*b32
                            c[2,2,ix,iy,iz,it] = a21*b12+a22*b22+a23*b32
                            c[3,2,ix,iy,iz,it] = a31*b12+a32*b22+a33*b32
                            c[1,3,ix,iy,iz,it] = a11*b13+a12*b23+a13*b33
                            c[2,3,ix,iy,iz,it] = a21*b13+a22*b23+a23*b33
                            c[3,3,ix,iy,iz,it] = a31*b13+a32*b23+a33*b33
                        end

                    end
                end
            end
        end
    end
    

    
    function LinearAlgebra.mul!(c::Gaugefields_4D_wing{3},a::T1,b::T2,α::Ta,β::Tb) where {NC,T1 <: Abstractfields,T2 <: Abstractfields,Ta <: Number, Tb <: Number}
        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        if β == zero(β)
            if α == one(α)
                mul!(c,a,b)
                return
            end
        end

        if α != 0
            @inbounds for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        @simd for ix=1:NX
                            a11 = a[1,1,ix,iy,iz,it]
                            a21 = a[2,1,ix,iy,iz,it]
                            a31 = a[3,1,ix,iy,iz,it]
                            a12 = a[1,2,ix,iy,iz,it]
                            a22 = a[2,2,ix,iy,iz,it]
                            a32 = a[3,2,ix,iy,iz,it]
                            a13 = a[1,3,ix,iy,iz,it]
                            a23 = a[2,3,ix,iy,iz,it]
                            a33 = a[3,3,ix,iy,iz,it]
                            b11 = b[1,1,ix,iy,iz,it]
                            b21 = b[2,1,ix,iy,iz,it]
                            b31 = b[3,1,ix,iy,iz,it]
                            b12 = b[1,2,ix,iy,iz,it]
                            b22 = b[2,2,ix,iy,iz,it]
                            b32 = b[3,2,ix,iy,iz,it]
                            b13 = b[1,3,ix,iy,iz,it]
                            b23 = b[2,3,ix,iy,iz,it]
                            b33 = b[3,3,ix,iy,iz,it]
                            c[1,1,ix,iy,iz,it] = (a11*b11+a12*b21+a13*b31)*α + β*c[1,1,ix,iy,iz,it]
                            c[2,1,ix,iy,iz,it] = (a21*b11+a22*b21+a23*b31)*α + β*c[2,1,ix,iy,iz,it]
                            c[3,1,ix,iy,iz,it] = (a31*b11+a32*b21+a33*b31)*α + β*c[3,1,ix,iy,iz,it]
                            c[1,2,ix,iy,iz,it] = (a11*b12+a12*b22+a13*b32)*α + β*c[1,2,ix,iy,iz,it]
                            c[2,2,ix,iy,iz,it] = (a21*b12+a22*b22+a23*b32)*α + β*c[2,2,ix,iy,iz,it]
                            c[3,2,ix,iy,iz,it] = (a31*b12+a32*b22+a33*b32)*α + β*c[3,2,ix,iy,iz,it]
                            c[1,3,ix,iy,iz,it] = (a11*b13+a12*b23+a13*b33)*α + β*c[1,3,ix,iy,iz,it]
                            c[2,3,ix,iy,iz,it] = (a21*b13+a22*b23+a23*b33)*α + β*c[2,3,ix,iy,iz,it]
                            c[3,3,ix,iy,iz,it] = (a31*b13+a32*b23+a33*b33)*α + β*c[3,3,ix,iy,iz,it]

                        end
                    end
                end
            end
        end
    end

    function mul_skiplastindex!(c::Gaugefields_4D_wing{NC},a::T1,b::T2) where {NC,T1 <: Abstractfields,T2 <: Abstractfields}
        #@assert NC != 2 && NC != 3 "This function is for NC != 2,3"
        NT = c.NT
        NZ = c.NZ
        NY = c.NY
        NX = c.NX
        #for it=1:NT
        it = 1
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k2=1:NC                            
                            for k1=1:NC
                                c[k1,k2,ix,iy,iz,it] = 0

                                @simd for k3=1:NC
                                    c[k1,k2,ix,iy,iz,it] += a[k1,k3,ix,iy,iz,it]*b[k3,k2,ix,iy,iz,it]
                                end
                            end
                        end
                    end
                end
            end
        #end
    end

    function staggered_phase(μ,ix,iy,iz,it,NX,NY,NZ,NT)
        t = it-1
        t += ifelse(t<0,NT,0)
        t += ifelse(t ≥ NT,-NT,0)
        #boundary_factor_t = ifelse(t == NT -1,BoundaryCondition[4],1)
        z = iz-1
        z += ifelse(z<0,NZ,0)
        z += ifelse(z ≥ NZ,-NZ,0)
        #boundary_factor_z = ifelse(z == NZ -1,BoundaryCondition[3],1)
        y = iy-1
        y += ifelse(y<0,NY,0)
        y += ifelse(y ≥ NY,-NY,0)
        #boundary_factor_y = ifelse(y == NY -1,BoundaryCondition[2],1)
        x = ix-1
        x += ifelse(x<0,NX,0)
        x += ifelse(x ≥ NX,-NX,0)
        #boundary_factor_x = ifelse(x == NX -1,BoundaryCondition[1],1)
        if μ ==1
            η = 1
        elseif μ ==2
            #η = (-1.0)^(x)
            η = ifelse(x%2 == 0,1,-1)
        elseif μ ==3
            #η = (-1.0)^(x+y)
            η = ifelse((x+y)%2 == 0,1,-1)
        elseif μ ==4
            #η = (-1.0)^(x+y+z)
            η = ifelse((x+y+z)%2 == 0,1,-1)
        else
            error("η should be positive but η = $η")
        end
        return η#*boundary_factor_x*boundary_factor_y*boundary_factor_z*boundary_factor_t
    end

    function gramschmidt!(v)
        n = size(v)[1]
        for i=1:n
            @simd for j=1:i-1
                v[:,i] = v[:,i] - v[:,j]'*v[:,i]*v[:,j]
            end
            v[:,i] = v[:,i]/norm(v[:,i])
        end
    end

    function normalize_U!(u::Gaugefields_4D_wing{NC}) where NC
        NX = u.NX
        NY = u.NY
        NZ = u.NZ
        NT = u.NT

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX
                        A = u[:,:,ix,iy,iz,it]
                        gramschmidt!(A)
                        u[:,:,ix,iy,iz,it] = A[:,:]
                    end
                end
            end
        end

    end

    function normalize_U!(u::Gaugefields_4D_wing{2}) 
        NX = u.NX
        NY = u.NY
        NZ = u.NZ
        NT = u.NT

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX

                        α = u[1,1,ix,iy,iz,it]
                        β = u[2,1,ix,iy,iz,it]
                        detU = abs(α)^2 + abs(β)^2
                        u[1,1,ix,iy,iz,it] = α/detU
                        u[2,1,ix,iy,iz,it] = β/detU
                        u[1,2,ix,iy,iz,it] = -conj(β)/detU
                        u[2,2,ix,iy,iz,it] = conj(α)/detU

                    end
                end
            end
        end
    end

    function normalize_U!(u::Gaugefields_4D_wing{3}) 
        NX = u.NX
        NY = u.NY
        NZ = u.NZ
        NT = u.NT

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX

                        w1 = 0
                        w2 = 0
                        @simd for ic=1:3
                            w1 += u[2,ic,ix,iy,iz,it]*conj(u[1,ic,ix,iy,iz,it])
                            w2 += u[1,ic,ix,iy,iz,it]*conj(u[1,ic,ix,iy,iz,it])
                        end
                        zerock2 = w2
                        if zerock2 == 0 
                            println("w2 is zero  !!  (in normlz)")
                            println("u[1,1),u[1,2),u[1,3) : ",u[1,1,ix,iy,iz,it], "\t",u[1,2,ix,iy,iz,it],"\t", u[1,3,ix,iy,iz,it])
                        end

                        w1 = -w1/w2

                        x4 = (u[2,1,ix,iy,iz,it]) + w1*u[1,1,ix,iy,iz,it]
                        x5 = (u[2,2,ix,iy,iz,it]) + w1*u[1,2,ix,iy,iz,it]
                        x6 = (u[2,3,ix,iy,iz,it]) + w1*u[1,3,ix,iy,iz,it]

                        w3 = x4*conj(x4) + x5*conj(x5) + x6*conj(x6)

                        zerock3 = w3
                        if zerock3 == 0
                            println("w3 is zero  !!  (in normlz)")
                            println("x4, x5, x6 : $x4, $x5, $x6")
                            exit()
                        end

                        u[2,1,ix,iy,iz,it] = x4
                        u[2,2,ix,iy,iz,it] = x5
                        u[2,3,ix,iy,iz,it] = x6

                        w3 = 1/sqrt(w3)
                        w2 = 1/sqrt(w2)

                        u[1,1,ix,iy,iz,it] = u[1,1,ix,iy,iz,it]*w2
                        u[1,2,ix,iy,iz,it] = u[1,2,ix,iy,iz,it]*w2
                        u[1,3,ix,iy,iz,it] = u[1,3,ix,iy,iz,it]*w2
                        u[2,1,ix,iy,iz,it] = u[2,1,ix,iy,iz,it]*w3
                        u[2,2,ix,iy,iz,it] = u[2,2,ix,iy,iz,it]*w3
                        u[2,3,ix,iy,iz,it] = u[2,3,ix,iy,iz,it]*w3

                        if zerock2*zerock3 == 0 
                            println("!! devided by zero !! (in normalize)")
                            println("w2 or w3 in normlz is zero !!")
                            println("w2, w3 : $w2, $w3   ")
                            exit()
                        end
                        #println(u[:,:,ix,iy,iz,it]'*u[:,:,ix,iy,iz,it])
                    end
                end
            end
        end
        m3complv!(u)
    end

    
    function m3complv!(a::Gaugefields_4D_wing{3})
        aa = zeros(Float64,18)
        NX = a.NX
        NY = a.NY
        NZ = a.NZ
        NT = a.NT
    
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX
    
                        aa[ 1] = real( a[1,1,ix,iy,iz,it])
                        aa[ 2] = imag(a[1,1,ix,iy,iz,it])
                        aa[ 3] = real( a[1,2,ix,iy,iz,it])
                        aa[ 4] = imag(a[1,2,ix,iy,iz,it])
                        aa[ 5] = real( a[1,3,ix,iy,iz,it])
                        aa[ 6] = imag(a[1,3,ix,iy,iz,it])
                        aa[ 7] = real( a[2,1,ix,iy,iz,it])
                        aa[ 8] = imag(a[2,1,ix,iy,iz,it])
                        aa[ 9] = real( a[2,2,ix,iy,iz,it])
                        aa[10] = imag(a[2,2,ix,iy,iz,it])
                        aa[11] = real( a[2,3,ix,iy,iz,it])
                        aa[12] = imag(a[2,3,ix,iy,iz,it])
    
                        aa[13] = aa[ 3]*aa[11] - aa[ 4]*aa[12] -
                                    aa[ 5]*aa[ 9] + aa[ 6]*aa[10]
                        aa[14] = aa[ 5]*aa[10] + aa[ 6]*aa[ 9] -
                                    aa[ 3]*aa[12] - aa[ 4]*aa[11]
                        aa[15] = aa[ 5]*aa[ 7] - aa[ 6]*aa[ 8] -
                                    aa[ 1]*aa[11] + aa[ 2]*aa[12]
                        aa[16] = aa[ 1]*aa[12] + aa[ 2]*aa[11] -
                                    aa[ 5]*aa[ 8] - aa[ 6]*aa[ 7]
                        aa[17] = aa[ 1]*aa[ 9] - aa[ 2]*aa[10] -
                                    aa[ 3]*aa[ 7] + aa[ 4]*aa[ 8]
                        aa[18] = aa[ 3]*aa[ 8] + aa[ 4]*aa[ 7] -
                                    aa[ 1]*aa[10] - aa[ 2]*aa[ 9]
    
                        a[3,1,ix,iy,iz,it] = aa[13]+im*aa[14]
                        a[3,2,ix,iy,iz,it] = aa[15]+im*aa[16]
                        a[3,3,ix,iy,iz,it] = aa[17]+im*aa[18]

                        #println(a[:,:,ix,iy,iz,it]'*a[:,:,ix,iy,iz,it] )
                    end
                end
            end
        end
    end


# end