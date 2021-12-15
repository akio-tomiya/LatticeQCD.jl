module WilsonFermion_module
    using LinearAlgebra
    #using ..Fermion:FermionFields

    import ..Actions:FermiActionParam,FermiActionParam_Wilson,
                FermiActionParam_WilsonClover,FermiActionParam_Staggered
    import ..Gaugefields:GaugeFields,GaugeFields_1d,SU3GaugeFields,SU2GaugeFields,SU3GaugeFields_1d,SU2GaugeFields_1d,
                staggered_phase,SUn,SU2,SU3,SUNGaugeFields,SUNGaugeFields_1d,SU


    #import ..Fermion:FermionFields,
    #            clear!
    #using ..Fermion
    import ..AbstractFermion:FermionFields,
        Wx!,Wdagx!,clear!,substitute_fermion!,Dx!,fermion_shift!,fermion_shiftB!,add!,set_wing_fermi!,WdagWx!,apply_periodicity,Ddagx!


    """
    Struct for WilsonFermion
    """
    struct WilsonFermion <: FermionFields
        NC::Int64
        NX::Int64
        NY::Int64
        NZ::Int64
        NT::Int64
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
    end

    function Base.length(x::WilsonFermion)
        return x.NC*x.NX*x.NY*x.NZ*x.NT*4
    end

    function Base.size(x::WilsonFermion)
        return (x.NC,x.NX,x.NY,x.NZ,x.NT,4)
    end

    function Base.iterate(x::WilsonFermion, i::Int = 1)
        i == length(x.f)+1 && return nothing
        #println("i = $i ",x.f[i])
        return (x.f[i], i+1)
    end

    function WilsonFermion(NC,NX,NY,NZ,NT,fparam::FermiActionParam,BoundaryCondition) 
        r = fparam.r
        hop = fparam.hop
        eps = fparam.eps
        MaxCGstep = fparam.MaxCGstep
        return WilsonFermion(NC,NX,NY,NZ,NT,r,hop,eps,MaxCGstep,BoundaryCondition)
    end
    
    function WilsonFermion(NC,NX,NY,NZ,NT,r,hop,eps,MaxCGstep,BoundaryCondition)#r,hop,eps,MaxCGstep)
        γ,rplusγ,rminusγ = mk_gamma(r)
        hopp = zeros(ComplexF64,4)
        hopm = zeros(ComplexF64,4)
        hopp .= hop
        hopm .= hop
        Dirac_operator = "Wilson"
        return WilsonFermion(NC,NX,NY,NZ,NT,zeros(ComplexF64,NC,NX+2,NY+2,NZ+2,NT+2,4),
            γ,rplusγ,rminusγ,hop,r,hopp,hopm,eps,Dirac_operator,MaxCGstep,BoundaryCondition)
    end

    function WilsonFermion(NC,NX,NY,NZ,NT,γ,rplusγ,rminusγ,hop,r,hopp,hopm,eps,fermion,MaxCGstep,BoundaryCondition)
        return WilsonFermion(NC,NX,NY,NZ,NT,zeros(ComplexF64,NC,NX+2,NY+2,NZ+2,NT+2,4),
                    γ,rplusγ,rminusγ,hop,r,hopp,hopm,eps,fermion,MaxCGstep,BoundaryCondition)
    end

    function Base.setindex!(x::WilsonFermion,v,i1,i2,i3,i4,i5,i6) 
        x.f[i1,i2 + 1,i3 + 1,i4 + 1,i5 + 1,i6] = v
    end

    function Base.setindex!(x::WilsonFermion,v,ii)
        #ii = (((((ialpha -1)*NT+it-1)*NZ+iz-1)*NY+iy-1)*NX+ix-1)*NC+ic
        ic = (ii - 1) % x.NC + 1
        iii = (ii - ic) ÷ x.NC 
        ix = iii % x.NX + 1
        iii = (iii - ix + 1) ÷ x.NX
        iy = iii % x.NY + 1
        iii = (iii - iy + 1) ÷ x.NY
        iz = iii % x.NZ + 1
        iii = (iii - iz + 1) ÷ x.NZ
        it = iii % x.NT + 1
        iii = (iii - it + 1) ÷ x.NT
        ialpha = iii + 1
        x[ic,ix,iy,iz,it,ialpha] = v
    end

    function Base.getindex(x::WilsonFermion,i1,i2,i3,i4,i5,i6)
        return x.f[i1,i2 .+ 1,i3 .+ 1,i4 .+ 1,i5 .+ 1,i6]
    end

    function Base.getindex(x::WilsonFermion,ii)
        #ii = (((((ialpha -1)*NT+it-1)*NZ+iz-1)*NY+iy-1)*NX+ix-1)*NC+ic
        ic = (ii - 1) % x.NC + 1
        iii = (ii - ic) ÷ x.NC 
        ix = iii % x.NX + 1
        iii = (iii - ix + 1) ÷ x.NX
        iy = iii % x.NY + 1
        iii = (iii - iy + 1) ÷ x.NY
        iz = iii % x.NZ + 1
        iii = (iii - iz + 1) ÷ x.NZ
        it = iii % x.NT + 1
        iii = (iii - it + 1) ÷ x.NT
        ialpha = iii + 1
        return x[ic,ix,iy,iz,it,ialpha]
    end


    function Base.:*(a::WilsonFermion,b::WilsonFermion)
        c = 0.0im
        for α=1:4
            for it=1:a.NT
                for iz=1:a.NZ
                    for iy=1:a.NY
                        for ix=1:a.NX
                            @simd for ic=1:a.NC
                                c+= conj(a[ic,ix,iy,iz,it,α])*b[ic,ix,iy,iz,it,α]
                            end
                        end
                    end
                end
            end            
        end
        return c
    end

    function Base.:*(a::T,b::WilsonFermion) where T <: Number
        c = similar(b)
        for α=1:4
            for it=1:b.NT
                for iz=1:b.NZ
                    for iy=1:b.NY
                        for ix=1:b.NX
                            @simd for ic=1:b.NC
                                c[ic,ix,iy,iz,it,α] = a*b[ic,ix,iy,iz,it,α]
                            end
                        end
                    end
                end
            end            
        end
        return c
    end

    function LinearAlgebra.axpy!(a::T,X::WilsonFermion,Y::WilsonFermion) where T <: Number #Y = a*X+Y
        for α=1:4
            for it=1:X.NT
                for iz=1:X.NZ
                    for iy=1:X.NY
                        for ix=1:X.NX
                            @simd for ic=1:X.NC
                                Y[ic,ix,iy,iz,it,α] += a*X[ic,ix,iy,iz,it,α]
                            end
                        end
                    end
                end
            end            
        end
        return X
    end


    function LinearAlgebra.rmul!(a::WilsonFermion,b::T) where T <: Number
        for α=1:4
            for it=1:a.NT
                for iz=1:a.NZ
                    for iy=1:a.NY
                        for ix=1:a.NX
                            @simd for ic=1:a.NC
                                a[ic,ix,iy,iz,it,α] = b*a[ic,ix,iy,iz,it,α]
                            end
                        end
                    end
                end
            end            
        end
        return a
    end

    function Base.similar(x::WilsonFermion)
        return WilsonFermion(x.NC,x.NX,x.NY,x.NZ,x.NT,
                    x.γ,x.rplusγ,x.rminusγ,x.hop,x.r,x.hopp,x.hopm,x.eps,x.Dirac_operator,x.MaxCGstep,x.BoundaryCondition)
    end

    function substitute_fermion!(H,j,x::WilsonFermion)
        i = 0
        for ialpha = 1:4
            for it=1:x.NT
                for iz=1:x.NZ
                    for iy=1:x.NY
                        for ix=1:x.NX
                            @simd for ic=1:x.NC
                                i += 1
                                H[i,j] = x[ic,ix,iy,iz,it,ialpha]
                            end
                        end
                    end
                end
            end
        end
    end

    function Wx!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1}) where  {T <: FermionFields,G <: GaugeFields}
        temp = temps[4]
        temp1 = temps[1]
        temp2 = temps[2]

        clear!(temp)
        set_wing_fermi!(x)
        for ν=1:4
            fermion_shift!(temp1,U,ν,x)

            #... Dirac multiplication
            mul!(temp1,view(x.rminusγ,:,:,ν),temp1)
            
            #
            fermion_shift!(temp2,U,-ν,x)
            mul!(temp2,view(x.rplusγ,:,:,ν),temp2)

            add!(temp,x.hopp[ν],temp1,x.hopm[ν],temp2)
            
        end

        clear!(xout)
        add!(xout,1,x,-1,temp)

        #display(xout)
        #    exit()
        return
    end

    function Wx!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1},fparam::FermiActionParam_Wilson) where  {T <: FermionFields,G <: GaugeFields}
        Wx!(xout,U,x,temps)
        return
    end

    

    
    function Wx!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1},fparam::FermiActionParam_WilsonClover) where  {T <: FermionFields,G <: GaugeFields}
        Wx!(xout,U,x,temps,fparam.CloverFμν)
        return
    end

    function Wx!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1},CloverFμν::AbstractArray) where  {T <: FermionFields,G <: GaugeFields}
        temp = temps[4]
        temp1 = temps[1]
        temp2 = temps[2]

        clear!(temp)
        set_wing_fermi!(x)
        for ν=1:4
            fermion_shift!(temp1,U,ν,x)

            #... Dirac multiplication
            mul!(temp1,view(x.rminusγ,:,:,ν),temp1)
            
            #
            fermion_shift!(temp2,U,-ν,x)
            mul!(temp2,view(x.rplusγ,:,:,ν),temp2)

            add!(temp,x.hopp[ν],temp1,x.hopm[ν],temp2)
            
        end

        clear!(xout)
        add!(xout,1,x,-1,temp)

        cloverterm!(xout,CloverFμν,x)
        #println( "xout ",xout*xout)



        #display(xout)
        #    exit()
        return
    end

    function Wdagx!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1},CloverFμν::AbstractArray) where {T <: FermionFields,G <: GaugeFields}
        temp = temps[4]
        temp1 = temps[1]
        temp2 = temps[2]

        clear!(temp)
        x5 = temps[3]

        mul_γ5x!(x5,x)
        set_wing_fermi!(x5)


        for ν=1:4
            fermion_shift!(temp1,U,ν,x5)

            #... Dirac multiplication
            mul!(temp1,view(x.rminusγ,:,:,ν),temp1)
            
            #
            fermion_shift!(temp2,U,-ν,x5)
            
            mul!(temp2,view(x.rplusγ,:,:,ν),temp2)

            add!(temp,x.hopp[ν],temp1,x.hopm[ν],temp2)
        end
        clear!(temp1)
        add!(temp1,1,x5,-1,temp)

        cloverterm!(temp1,CloverFμν,x5)

        mul_γ5x!(xout,temp1)
        return
    end

    function Wdagx!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1},fparam::FermiActionParam_WilsonClover) where {T <: FermionFields,G <: GaugeFields}
        Wdagx!(xout,U,x,temps,fparam.CloverFμν)
        return

    end

    
    function Wdagx!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1}) where  {T <: FermionFields,G <: GaugeFields}
        temp = temps[4]
        temp1 = temps[1]
        temp2 = temps[2]

        clear!(temp)
        set_wing_fermi!(x)
        for ν=1:4
            fermion_shift!(temp1,U,ν,x)

            #... Dirac multiplication
            #mul!(temp1,view(x.rminusγ,:,:,ν),temp1)
            mul!(temp1,view(x.rplusγ,:,:,ν),temp1)
            
            #
            fermion_shift!(temp2,U,-ν,x)
            #mul!(temp2,view(x.rplusγ,:,:,ν),temp2)
            mul!(temp2,view(x.rminusγ,:,:,ν),temp2)

            add!(temp,x.hopp[ν],temp1,x.hopm[ν],temp2)
            
            
        end

        clear!(xout)
        add!(xout,1,x,-1,temp)

        #display(xout)
        #    exit()
        return
    end

    function Wdagx!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1},fparam::FermiActionParam_Wilson) where {T <: FermionFields,G <: GaugeFields}
        Wdagx!(xout,U,x,temps)
        return
    end

    function Dx!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1}) where  {T <: FermionFields,G <: GaugeFields}
        temp = temps[4]
        temp1 = temps[1]
        temp2 = temps[2]

        clear!(temp)
        set_wing_fermi!(x)
        for ν=1:4
            fermion_shift!(temp1,U,ν,x)

            #... Dirac multiplication
            mul!(temp1,view(x.rminusγ,:,:,ν),temp1)
            
            #
            fermion_shift!(temp2,U,-ν,x)
            mul!(temp2,view(x.rplusγ,:,:,ν),temp2)

            add!(temp,0.5,temp1,0.5,temp2)
            
        end

        clear!(xout)
        add!(xout,1/(2*x.hop),x,-1,temp)

        #display(xout)
        #    exit()
        return
    end

    function Ddagx!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1}) where  {T <: FermionFields,G <: GaugeFields}
        temp = temps[4]
        temp1 = temps[1]
        temp2 = temps[2]

        clear!(temp)
        set_wing_fermi!(x)
        for ν=1:4
            fermion_shift!(temp1,U,ν,x)

            #... Dirac multiplication
            #mul!(temp1,view(x.rminusγ,:,:,ν),temp1)
            mul!(temp1,view(x.rplusγ,:,:,ν),temp1)
            
            #
            fermion_shift!(temp2,U,-ν,x)
            #mul!(temp2,view(x.rplusγ,:,:,ν),temp2)
            mul!(temp2,view(x.rminusγ,:,:,ν),temp2)

            add!(temp,0.5,temp1,0.5,temp2)
            
            
        end

        clear!(xout)
        add!(xout,1/(2*x.hop),x,-1,temp)

        #display(xout)
        #    exit()
        return
    end

    function LinearAlgebra.mul!(xout::WilsonFermion,A::AbstractMatrix,x::WilsonFermion)
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        NC = x.NC

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

                                xout[ic,ix,iy,iz,it,1] = A[1,1]*e1+A[1,2]*e2+A[1,3]*e3+A[1,4]*e4
                                xout[ic,ix,iy,iz,it,2] = A[2,1]*e1+A[2,2]*e2+A[2,3]*e3+A[2,4]*e4
                                xout[ic,ix,iy,iz,it,3] = A[3,1]*e1+A[3,2]*e2+A[3,3]*e3+A[3,4]*e4
                                xout[ic,ix,iy,iz,it,4] = A[4,1]*e1+A[4,2]*e2+A[4,3]*e3+A[4,4]*e4



                                #=
                            for k1=1:4
                                e[k1] = x[ic,ix,iy,iz,it,k1]
                            end

                            for k1=1:4
                                f[k1] = 0
                                for k2=1:4
                                    f[k1] += A[k1,k2]*e[k2]
                                end
                                xout[ic,ix,iy,iz,it,k1] = f[k1]
                            end
                            =#

                        end
                    end
                end
            end
        end
    end

    function LinearAlgebra.mul!(xout::WilsonFermion,x::WilsonFermion,A::AbstractMatrix)
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        NC = x.NC


        #n6 = size(x.f)[6]
        #f = zeros(ComplexF64,n6)
        #e = zeros(ComplexF64,n6)
        
        for ic=1:NC
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        @simd for ix=1:NX
                            e1 = x[ic,ix,iy,iz,it,1]
                            e2 = x[ic,ix,iy,iz,it,2]
                            e3 = x[ic,ix,iy,iz,it,3]
                            e4 = x[ic,ix,iy,iz,it,4]

                            xout[ic,ix,iy,iz,it,1] = A[1,1]*e1+A[2,1]*e2+A[3,1]*e3+A[4,1]*e4
                            xout[ic,ix,iy,iz,it,2] = A[1,2]*e1+A[2,2]*e2+A[3,2]*e3+A[4,2]*e4
                            xout[ic,ix,iy,iz,it,3] = A[1,3]*e1+A[2,3]*e2+A[3,3]*e3+A[4,3]*e4
                            xout[ic,ix,iy,iz,it,4] = A[1,4]*e1+A[2,4]*e2+A[3,4]*e3+A[4,4]*e4


                            #=
                            for k1=1:n6
                                e[k1] = x[ic,ix,iy,iz,it,k1]
                            end

                            for k1=1:n6
                                f[k1] = 0
                                for k2=1:n6
                                    f[k1] += e[k2]*A[k2,k1]
                                end
                                xout[ic,ix,iy,iz,it,k1] = f[k1]
                            end
                            =#

                        end
                    end
                end
            end
        end
    end

    
    function fermion_shift!(b::WilsonFermion,u::Array{T,1},μ::Int,a::WilsonFermion) where T <: SU3GaugeFields
        if μ == 0
            substitute!(b,a)
            return
        end

        NX = a.NX
        NY = a.NY
        NZ = a.NZ
        NT = a.NT
        NC = 3#a.NC

        #NTrange = get_looprange(NT)
        #println(NTrange)
        if μ > 0
            #idel = zeros(Int64,4)
            #idel[μ] = 1

            n6 = size(a.f)[6]
            for ialpha=1:4
                #for it=NTrange
                for it=1:NT
                    it1 = it + ifelse(μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz + ifelse(μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(μ ==2,1,0) #idel[2]
                            @simd for ix=1:NX
                                ix1 = ix + ifelse(μ ==1,1,0) #idel[1]

                                b[1,ix,iy,iz,it,ialpha] = u[μ][1,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][1,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][1,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha]

                                b[2,ix,iy,iz,it,ialpha] = u[μ][2,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][2,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][2,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha]

                                b[3,ix,iy,iz,it,ialpha] = u[μ][3,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][3,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][3,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha]

                                #=
                                for k1=1:3
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:3
                                        b[k1,ix,iy,iz,it,ialpha] += u[μ][k1,k2,ix,iy,iz,it]*a[k2,ix1,iy1,iz1,it1,ialpha]
                                    end
                                end
                                =#
                            end
                        end
                    end
                end
            end
            
        elseif μ < 0
            #idel = zeros(Int64,4)
            #idel[-μ] = 1
            #n6 = size(b.f)[6]
            for ialpha =1:4
                for it=1:NT
                    it1 = it - ifelse(-μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz - ifelse(-μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy - ifelse(-μ ==2,1,0)  #idel[2]
                            @simd for ix=1:NX
                                ix1 = ix - ifelse(-μ ==1,1,0) #idel[1]

                                b[1,ix,iy,iz,it,ialpha] = conj(u[-μ][1,1,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,1,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][3,1,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha]

                                b[2,ix,iy,iz,it,ialpha] = conj(u[-μ][1,2,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,2,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][3,2,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha]

                                b[3,ix,iy,iz,it,ialpha] = conj(u[-μ][1,3,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,3,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][3,3,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha]
                                #=
                                for k1=1:3
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:3
                                        b[k1,ix,iy,iz,it,ialpha] += conj(u[-μ][k2,k1,ix1,iy1,iz1,it1])*a[k2,ix1,iy1,iz1,it1,ialpha]
                                    end
                                end
                                =#
                            end
                        end
                    end
                end
            end
        end

    end

    function fermion_shift_gamma!(b::WilsonFermion,u::Array{T,1},μ::Int,a::WilsonFermion) where T <: SU3GaugeFields
        if μ == 0
            substitute!(b,a)
            return
        end

        NX = a.NX
        NY = a.NY
        NZ = a.NZ
        NT = a.NT
        NC = 3#a.NC

        #NTrange = get_looprange(NT)
        #println(NTrange)
        if μ > 0
            #idel = zeros(Int64,4)
            #idel[μ] = 1

            n6 = size(a.f)[6]
            for ialpha=1:4
                #for it=NTrange
                for it=1:NT
                    it1 = it + ifelse(μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz + ifelse(μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(μ ==2,1,0) #idel[2]
                            @simd for ix=1:NX
                                ix1 = ix + ifelse(μ ==1,1,0) #idel[1]

                                b[1,ix,iy,iz,it,ialpha] = u[μ][1,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][1,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][1,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha]

                                b[2,ix,iy,iz,it,ialpha] = u[μ][2,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][2,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][2,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha]

                                b[3,ix,iy,iz,it,ialpha] = u[μ][3,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][3,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][3,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha]

                                #=
                                for k1=1:3
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:3
                                        b[k1,ix,iy,iz,it,ialpha] += u[μ][k1,k2,ix,iy,iz,it]*a[k2,ix1,iy1,iz1,it1,ialpha]
                                    end
                                end
                                =#
                            end
                        end
                    end
                end
            end
            
        elseif μ < 0
            #idel = zeros(Int64,4)
            #idel[-μ] = 1
            #n6 = size(b.f)[6]
            for ialpha =1:4
                for it=1:NT
                    it1 = it - ifelse(-μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz - ifelse(-μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy - ifelse(-μ ==2,1,0)  #idel[2]
                            @simd for ix=1:NX
                                ix1 = ix - ifelse(-μ ==1,1,0) #idel[1]

                                b[1,ix,iy,iz,it,ialpha] = conj(u[-μ][1,1,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,1,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][3,1,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha]

                                b[2,ix,iy,iz,it,ialpha] = conj(u[-μ][1,2,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,2,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][3,2,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha]

                                b[3,ix,iy,iz,it,ialpha] = conj(u[-μ][1,3,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,3,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][3,3,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha]
                                #=
                                for k1=1:3
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:3
                                        b[k1,ix,iy,iz,it,ialpha] += conj(u[-μ][k2,k1,ix1,iy1,iz1,it1])*a[k2,ix1,iy1,iz1,it1,ialpha]
                                    end
                                end
                                =#
                            end
                        end
                    end
                end
            end
        end

    end

    function fermion_shift!(b::WilsonFermion,u::Array{T,1},μ::Int,a::WilsonFermion) where T <: SU2GaugeFields
        if μ == 0
            substitute!(b,a)
            return
        end

        NX = a.NX
        NY = a.NY
        NZ = a.NZ
        NT = a.NT
        NC = a.NC

        #NTrange = get_looprange(NT)
        #println(NTrange)
        if μ > 0
            #idel = zeros(Int64,4)
            #idel[μ] = 1

            n6 = size(a.f)[6]
            for ialpha=1:4
                #for it=NTrange
                for it=1:NT
                    it1 = it + ifelse(μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz + ifelse(μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(μ ==2,1,0) #idel[2]
                            for ix=1:NX
                                ix1 = ix + ifelse(μ ==1,1,0) #idel[1]

                                b[1,ix,iy,iz,it,ialpha] = u[μ][1,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][1,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] 

                                b[2,ix,iy,iz,it,ialpha] = u[μ][2,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][2,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] 


                                #=
                                for k1=1:3
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:3
                                        b[k1,ix,iy,iz,it,ialpha] += u[μ][k1,k2,ix,iy,iz,it]*a[k2,ix1,iy1,iz1,it1,ialpha]
                                    end
                                end
                                =#
                            end
                        end
                    end
                end
            end
            
        elseif μ < 0
            #idel = zeros(Int64,4)
            #idel[-μ] = 1
            #n6 = size(b.f)[6]
            for ialpha =1:4
                for it=1:NT
                    it1 = it - ifelse(-μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz - ifelse(-μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy - ifelse(-μ ==2,1,0)  #idel[2]
                            for ix=1:NX
                                ix1 = ix - ifelse(-μ ==1,1,0) #idel[1]

                                b[1,ix,iy,iz,it,ialpha] = conj(u[-μ][1,1,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,1,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] 

                                b[2,ix,iy,iz,it,ialpha] = conj(u[-μ][1,2,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,2,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] 

                            end
                        end
                    end
                end
            end
        end

    end

    function fermion_shiftB!(b::WilsonFermion,u::Array{GaugeFields{SU{NC}},1},μ,a::WilsonFermion) where NC
        if μ == 0
            substitute!(b,a)
            return
        end
        

        NX = a.NX
        NY = a.NY
        NZ = a.NZ
        NT = a.NT
        #NC = a.NC


        if μ > 0
            error("""
            Sorry this case if not yet ready
            mu = $mu
            """)
            
        elseif μ < 0
            #idel = zeros(Int64,4)
            #idel[-μ] = 1
            #n6 = size(b.f)[6]
            for ialpha =1:4
                for it=1:NT
                    it1 = it + ifelse(-μ == 4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz +ifelse(-μ == 3,1,0)  #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(-μ == 2,1,0) #idel[2]
                            for ix=1:NX
                                ix1 = ix + ifelse(-μ == 1,1,0)  #idel[1]
                                                                
                                for k1=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:NC
                                        b[k1,ix,iy,iz,it,ialpha] += conj(a[k2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][k1,k2,ix,iy,iz,it])
                                    end
                                end
                                
                                
                            end
                        end
                    end
                end
            end
        end

    end

    function fermion_shiftB!(b::WilsonFermion,evensite,u::Array{GaugeFields{SU{NC}},1},μ,a::WilsonFermion)  where NC
        if μ == 0
            substitute!(b,a)
            return
        end

        ibush = ifelse(evensite,0,1)
        

        NX = a.NX
        NY = a.NY
        NZ = a.NZ
        NT = a.NT
        #NC = a.NC


        if μ > 0
            error("""
            Sorry this case if not yet ready
            mu = $mu
            """)
            
        elseif μ < 0
            #idel = zeros(Int64,4)
            #idel[-μ] = 1
            #n6 = size(b.f)[6]
            for ialpha =1:4
                for it=1:NT
                    it1 = it + ifelse(-μ == 4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz +ifelse(-μ == 3,1,0)  #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(-μ == 2,1,0) #idel[2]
                            xran =1+(1+ibush+iy+iz+it)%2:2:NX
                            for ix in xran
                            #for ix=1:NX
                                ix1 = ix + ifelse(-μ == 1,1,0)  #idel[1]
                                
                                
                                for k1=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:NC
                                        b[k1,ix,iy,iz,it,ialpha] += conj(a[k2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][k1,k2,ix,iy,iz,it])
                                    end
                                end
                                
                                
                            end
                        end
                    end
                end
            end
        end

    end

    function fermion_shiftB!(b::WilsonFermion,u::Array{SU3GaugeFields,1},μ,a::WilsonFermion) 
        if μ == 0
            substitute!(b,a)
            return
        end
        

        NX = a.NX
        NY = a.NY
        NZ = a.NZ
        NT = a.NT
        NC = a.NC


        if μ > 0
            error("""
            Sorry this case if not yet ready
            mu = $mu
            """)
            
        elseif μ < 0
            #idel = zeros(Int64,4)
            #idel[-μ] = 1
            #n6 = size(b.f)[6]
            for ialpha =1:4
                for it=1:NT
                    it1 = it + ifelse(-μ == 4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz +ifelse(-μ == 3,1,0)  #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(-μ == 2,1,0) #idel[2]
                            for ix=1:NX
                                ix1 = ix + ifelse(-μ == 1,1,0)  #idel[1]
                                

                                b[1,ix,iy,iz,it,ialpha] = conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,1,ix,iy,iz,it]) + 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,2,ix,iy,iz,it]) + 
                                                            conj(a[3,ix1,iy1,iz1,it1,ialpha])* conj(u[-μ][1,3,ix,iy,iz,it])

                                b[2,ix,iy,iz,it,ialpha] = conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,1,ix,iy,iz,it])+ 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,2,ix,iy,iz,it]) + 
                                                            conj(a[3,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,3,ix,iy,iz,it])

                                b[3,ix,iy,iz,it,ialpha] = conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][3,1,ix,iy,iz,it]) + 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][3,2,ix,iy,iz,it]) + 
                                                            conj(a[3,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][3,3,ix,iy,iz,it])

                        
                                #=
                                for k1=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:NC
                                        b[k1,ix,iy,iz,it,ialpha] += conj(a[k2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][k1,k2,ix,iy,iz,it])
                                    end
                                end
                                =#
                                
                            end
                        end
                    end
                end
            end
        end

    end

    function fermion_shiftB!(b::WilsonFermion,evensite,u::Array{SU3GaugeFields,1},μ,a::WilsonFermion) 
        if μ == 0
            substitute!(b,a)
            return
        end

        ibush = ifelse(evensite,0,1)
        

        NX = a.NX
        NY = a.NY
        NZ = a.NZ
        NT = a.NT
        NC = a.NC


        if μ > 0
            error("""
            Sorry this case if not yet ready
            mu = $mu
            """)
            
        elseif μ < 0
            #idel = zeros(Int64,4)
            #idel[-μ] = 1
            #n6 = size(b.f)[6]
            for ialpha =1:4
                for it=1:NT
                    it1 = it + ifelse(-μ == 4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz +ifelse(-μ == 3,1,0)  #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(-μ == 2,1,0) #idel[2]
                            xran =1+(1+ibush+iy+iz+it)%2:2:NX
                            for ix in xran
                            #for ix=1:NX
                                ix1 = ix + ifelse(-μ == 1,1,0)  #idel[1]
                                

                                b[1,ix,iy,iz,it,ialpha] = conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,1,ix,iy,iz,it]) + 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,2,ix,iy,iz,it]) + 
                                                            conj(a[3,ix1,iy1,iz1,it1,ialpha])* conj(u[-μ][1,3,ix,iy,iz,it])

                                b[2,ix,iy,iz,it,ialpha] = conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,1,ix,iy,iz,it])+ 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,2,ix,iy,iz,it]) + 
                                                            conj(a[3,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,3,ix,iy,iz,it])

                                b[3,ix,iy,iz,it,ialpha] = conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][3,1,ix,iy,iz,it]) + 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][3,2,ix,iy,iz,it]) + 
                                                            conj(a[3,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][3,3,ix,iy,iz,it])

                        
                                #=
                                for k1=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:NC
                                        b[k1,ix,iy,iz,it,ialpha] += conj(a[k2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][k1,k2,ix,iy,iz,it])
                                    end
                                end
                                =#
                                
                            end
                        end
                    end
                end
            end
        end

    end

    function fermion_shiftB!(b::WilsonFermion,u::Array{SU2GaugeFields,1},μ,a::WilsonFermion) 
        if μ == 0
            substitute!(b,a)
            return
        end
        

        NX = a.NX
        NY = a.NY
        NZ = a.NZ
        NT = a.NT
        NC = a.NC


        if μ > 0
            error("""
            Sorry this case if not yet ready
            mu = $mu
            """)
            
        elseif μ < 0
            #idel = zeros(Int64,4)
            #idel[-μ] = 1
            #n6 = size(b.f)[6]
            for ialpha =1:4
                for it=1:NT
                    it1 = it + ifelse(-μ == 4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz +ifelse(-μ == 3,1,0)  #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(-μ == 2,1,0) #idel[2]
                            for ix=1:NX
                                ix1 = ix + ifelse(-μ == 1,1,0)  #idel[1]
                                

                                b[1,ix,iy,iz,it,ialpha] = conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,1,ix,iy,iz,it]) + 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,2,ix,iy,iz,it])

                                b[2,ix,iy,iz,it,ialpha] = conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,1,ix,iy,iz,it])+ 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,2,ix,iy,iz,it]) 

                                
                            end
                        end
                    end
                end
            end
        end

    end

    function fermion_shiftB!(b::WilsonFermion,evensite,u::Array{SU2GaugeFields,1},μ,a::WilsonFermion) 
        if μ == 0
            substitute!(b,a)
            return
        end

        ibush = ifelse(evensite,0,1)
        

        NX = a.NX
        NY = a.NY
        NZ = a.NZ
        NT = a.NT
        NC = a.NC


        if μ > 0
            error("""
            Sorry this case if not yet ready
            mu = $mu
            """)
            
        elseif μ < 0
            #idel = zeros(Int64,4)
            #idel[-μ] = 1
            #n6 = size(b.f)[6]
            for ialpha =1:4
                for it=1:NT
                    it1 = it + ifelse(-μ == 4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz +ifelse(-μ == 3,1,0)  #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(-μ == 2,1,0) #idel[2]
                            xran =1+(1+ibush+iy+iz+it)%2:2:NX
                            for ix in xran
                                ix1 = ix + ifelse(-μ == 1,1,0)  #idel[1]
                                

                                b[1,ix,iy,iz,it,ialpha] = conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,1,ix,iy,iz,it]) + 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,2,ix,iy,iz,it])

                                b[2,ix,iy,iz,it,ialpha] = conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,1,ix,iy,iz,it])+ 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,2,ix,iy,iz,it]) 

                                
                            end
                        end
                    end
                end
            end
        end

    end




    """
    c--------------------------------------------------------------------------c
    c     y = gamma_5 * x
    c     here
    c                  ( -1       )
    c        GAMMA5 =  (   -1     )
    c                  (     +1   )
    c                  (       +1 )
    c--------------------------------------------------------------------------c
        """
        function mul_γ5x!(y::WilsonFermion,x::WilsonFermion)
            NX = x.NX
            NY = x.NY
            NZ = x.NZ
            NT = x.NT
            NC = x.NC
            for ic=1:NC
                for it=1:NT
                    for iz=1:NZ
                        for iy=1:NY
                            @simd for ix=1:NX
                                y[ic,ix,iy,iz,it,1] = -1*x[ic,ix,iy,iz,it,1]
                                y[ic,ix,iy,iz,it,2] = -1*x[ic,ix,iy,iz,it,2]
                                y[ic,ix,iy,iz,it,3] = x[ic,ix,iy,iz,it,3]
                                y[ic,ix,iy,iz,it,4] = x[ic,ix,iy,iz,it,4]
                            end
                        end
                    end
                end
            end
        end


        function mul_1plusγ5x!(y::WilsonFermion,x::WilsonFermion) #(1+gamma_5)/2
            NX = x.NX
            NY = x.NY
            NZ = x.NZ
            NT = x.NT
            NC = x.NC
            for ic=1:NC
                for it=1:NT
                    for iz=1:NZ
                        for iy=1:NY
                            @simd for ix=1:NX
                                y[ic,ix,iy,iz,it,1] = 0#-1*x[ic,ix,iy,iz,it,1]
                                y[ic,ix,iy,iz,it,2] = 0#-1*x[ic,ix,iy,iz,it,2]
                                y[ic,ix,iy,iz,it,3] = x[ic,ix,iy,iz,it,3]
                                y[ic,ix,iy,iz,it,4] = x[ic,ix,iy,iz,it,4]
                            end
                        end
                    end
                end
            end
        end

        function mul_1plusγ5x_add!(y::WilsonFermion,x::WilsonFermion,factor) #x = x +(1+gamma_5)/2
            NX = x.NX
            NY = x.NY
            NZ = x.NZ
            NT = x.NT
            NC = x.NC
            for ic=1:NC
                for it=1:NT
                    for iz=1:NZ
                        for iy=1:NY
                            @simd for ix=1:NX
                                y[ic,ix,iy,iz,it,3] += factor*x[ic,ix,iy,iz,it,3]
                                y[ic,ix,iy,iz,it,4] += factor*x[ic,ix,iy,iz,it,4]
                            end
                        end
                    end
                end
            end
        end

        function mul_1minusγ5x!(y::WilsonFermion,x::WilsonFermion) #(1-gamma_5)/2
            NX = x.NX
            NY = x.NY
            NZ = x.NZ
            NT = x.NT
            NC = x.NC
            for ic=1:NC
                for it=1:NT
                    for iz=1:NZ
                        for iy=1:NY
                            @simd for ix=1:NX
                                y[ic,ix,iy,iz,it,1] = x[ic,ix,iy,iz,it,1]
                                y[ic,ix,iy,iz,it,2] = x[ic,ix,iy,iz,it,2]
                                y[ic,ix,iy,iz,it,3] = 0#x[ic,ix,iy,iz,it,3]
                                y[ic,ix,iy,iz,it,4] = 0#x[ic,ix,iy,iz,it,4]
                            end
                        end
                    end
                end
            end
        end

        function mul_1minusγ5x_add!(y::WilsonFermion,x::WilsonFermion,factor) #+(1-gamma_5)/2
            NX = x.NX
            NY = x.NY
            NZ = x.NZ
            NT = x.NT
            NC = x.NC
            for ic=1:NC
                for it=1:NT
                    for iz=1:NZ
                        for iy=1:NY
                            @simd for ix=1:NX
                                y[ic,ix,iy,iz,it,1] += factor*x[ic,ix,iy,iz,it,1]
                                y[ic,ix,iy,iz,it,2] += factor*x[ic,ix,iy,iz,it,2]
                            end
                        end
                    end
                end
            end
        end
    

    
    """
    mk_gamma()
    c----------------------------------------------------------------------c
c     Make gamma matrix
c----------------------------------------------------------------------c
C     THE CONVENTION OF THE GAMMA MATRIX HERE
C     ( EUCLIDEAN CHIRAL REPRESENTATION )
C
C               (       -i )              (       -1 )
C     GAMMA1 =  (     -i   )     GAMMA2 = (     +1   )
C               (   +i     )              (   +1     )
C               ( +i       )              ( -1       )
C
C               (     -i   )              (     -1   )
C     GAMMA3 =  (       +i )     GAMMA4 = (       -1 )
C               ( +i       )              ( -1       )
C               (   -i     )              (   -1     )
C
C               ( -1       )
C     GAMMA5 =  (   -1     )
C               (     +1   )
C               (       +1 )
C
C     ( GAMMA_MU, GAMMA_NU ) = 2*DEL_MU,NU   FOR MU,NU=1,2,3,4   
c----------------------------------------------------------------------c
    """
    function mk_gamma(r)
        g0 = zeros(ComplexF64,4,4)
        g1 = zero(g0)
        g2 = zero(g1)
        g3 = zero(g1)
        g4 = zero(g1)
        g5 = zero(g1)
        gamma = zeros(ComplexF64,4,4,5)
        rpg = zero(gamma)
        rmg = zero(gamma)


        g0[1,1]=1.0; g0[1,2]=0.0; g0[1,3]=0.0; g0[1,4]=0.0
        g0[2,1]=0.0; g0[2,2]=1.0; g0[2,3]=0.0; g0[2,4]=0.0
        g0[3,1]=0.0; g0[3,2]=0.0; g0[3,3]=1.0; g0[3,4]=0.0
        g0[4,1]=0.0; g0[4,2]=0.0; g0[4,3]=0.0; g0[4,4]=1.0

        g1[1,1]=0.0; g1[1,2]=0.0; g1[1,3]=0.0; g1[1,4]=-im
        g1[2,1]=0.0; g1[2,2]=0.0; g1[2,3]=-im;  g1[2,4]=0.0
        g1[3,1]=0.0; g1[3,2]=+im;  g1[3,3]=0.0; g1[3,4]=0.0
        g1[4,1]=+im;  g1[4,2]=0.0; g1[4,3]=0.0; g1[4,4]=0.0

        g2[1,1]=0.0; g2[1,2]=0.0; g2[1,3]=0.0; g2[1,4]=-1.0
        g2[2,1]=0.0; g2[2,2]=0.0; g2[2,3]=1.0; g2[2,4]=0.0
        g2[3,1]=0.0; g2[3,2]=1.0; g2[3,3]=0.0; g2[3,4]=0.0
        g2[4,1]=-1.0;g2[4,2]=0.0; g2[4,3]=0.0; g2[4,4]=0.0

        g3[1,1]=0.0; g3[1,2]=0.0; g3[1,3]=-im;  g3[1,4]=0.0
        g3[2,1]=0.0; g3[2,2]=0.0; g3[2,3]=0.0; g3[2,4]=+im
        g3[3,1]=+im;  g3[3,2]=0.0; g3[3,3]=0.0; g3[3,4]=0.0
        g3[4,1]=0.0; g3[4,2]=-im;  g3[4,3]=0.0; g3[4,4]=0.0

        g4[1,1]=0.0; g4[1,2]=0.0; g4[1,3]=-1.0;g4[1,4]=0.0
        g4[2,1]=0.0; g4[2,2]=0.0; g4[2,3]=0.0; g4[2,4]=-1.0
        g4[3,1]=-1.0;g4[3,2]=0.0; g4[3,3]=0.0; g4[3,4]=0.0
        g4[4,1]=0.0; g4[4,2]=-1.0;g4[4,3]=0.0; g4[4,4]=0.0

        g5[1,1]=-1.0;g5[1,2]=0.0; g5[1,3]=0.0; g5[1,4]=0.0
        g5[2,1]=0.0; g5[2,2]=-1.0;g5[2,3]=0.0; g5[2,4]=0.0
        g5[3,1]=0.0; g5[3,2]=0.0; g5[3,3]=1.0; g5[3,4]=0.0
        g5[4,1]=0.0; g5[4,2]=0.0; g5[4,3]=0.0; g5[4,4]=1.0

        gamma[:,:,1] = g1[:,:]
        gamma[:,:,2] = g2[:,:]
        gamma[:,:,3] = g3[:,:]
        gamma[:,:,4] = g4[:,:]
        gamma[:,:,5] = g5[:,:]

        for mu=1:4
            for j=1:4
                for i=1:4
                    rpg[i,j,mu] = r*g0[i,j] + gamma[i,j,mu]
                    rmg[i,j,mu] = r*g0[i,j] - gamma[i,j,mu]
                end
            end
        end 

        return gamma,rpg,rmg


    end


    function cloverterm!(vec,CloverFμν,x)
        NT = x.NT
        NZ = x.NZ
        NY = x.NY
        NX = x.NX
        NC = x.NC

        
        i  =0
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        i += 1
                        for k1=1:NC
                            for k2=1:NC

                                c1 = x[k2,ix,iy,iz,it,1]
                                c2 = x[k2,ix,iy,iz,it,2]
                                c3 = x[k2,ix,iy,iz,it,3]
                                c4 = x[k2,ix,iy,iz,it,4]

                                vec[k1,ix,iy,iz,it,1] += CloverFμν[k1,k2,i,1]*(-   c1) + 
                                                            + CloverFμν[k1,k2,i,2]*(-im*c2) + 
                                                            + CloverFμν[k1,k2,i,3]*(-   c2) + 
                                                            + CloverFμν[k1,k2,i,4]*(-   c2) + 
                                                            + CloverFμν[k1,k2,i,5]*( im*c2) + 
                                                            + CloverFμν[k1,k2,i,6]*(-   c1)
                
                                

                                vec[k1,ix,iy,iz,it,2] += CloverFμν[k1,k2,i,1]*(   c2) + 
                                                            + CloverFμν[k1,k2,i,2]*(im*c1) + 
                                                            + CloverFμν[k1,k2,i,3]*(-   c1) + 
                                                            + CloverFμν[k1,k2,i,4]*(-   c1) + 
                                                            + CloverFμν[k1,k2,i,5]*(-im*c1) + 
                                                            + CloverFμν[k1,k2,i,6]*(   c2)

                                vec[k1,ix,iy,iz,it,3] += CloverFμν[k1,k2,i,1]*(   -c3) + 
                                                            + CloverFμν[k1,k2,i,2]*(-im*c4) + 
                                                            + CloverFμν[k1,k2,i,3]*(   c4) + 
                                                            + CloverFμν[k1,k2,i,4]*(-   c4) + 
                                                            + CloverFμν[k1,k2,i,5]*(-im*c4) + 
                                                            + CloverFμν[k1,k2,i,6]*(   c3)

                                vec[k1,ix,iy,iz,it,4] += CloverFμν[k1,k2,i,1]*(   c4) + 
                                                            + CloverFμν[k1,k2,i,2]*(im*c3) + 
                                                            + CloverFμν[k1,k2,i,3]*(   c3) + 
                                                            + CloverFμν[k1,k2,i,4]*(-   c3) + 
                                                            + CloverFμν[k1,k2,i,5]*(im*c3) + 
                                                            + CloverFμν[k1,k2,i,6]*( -  c4)


                            end
                        end
                    end
                end

            end
        end

        #println("vec = ",vec*vec)

    end




end