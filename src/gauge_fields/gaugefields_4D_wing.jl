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
        x.U[i1,i2,i3 + x.NDW,i4 + x.NDW,i5 + x.NDW,i6 + x.NDW] = v
    end

    function Base.getindex(x::Gaugefields_4D_wing,i1,i2,i3,i4,i5,i6) 
        return x.U[i1,i2,i3 .+ x.NDW,i4 .+ x.NDW,i5 .+ x.NDW,i6 .+ x.NDW]
    end


    function Base.getindex(U::Adjoint_Gaugefields{T},i1,i2,i3,i4,i5,i6) where T <: Gaugefields_4D_wing #U'
        return conj(U.parent[i2,i1,i3,i4,i5,i6])
    end

    function Base.setindex!(U::Adjoint_Gaugefields{T},v,i1,i2,i3,i4,i5,i6,μ)  where T <: Gaugefields_4D_wing
        error("type $(typeof(U)) has no setindex method. This type is read only.")
    end

    function substitute_U!(a::Array{T,1},b::Array{T,1}) where T <: Gaugefields_4D_wing
        for μ=1:4
            substitute_U!(a[μ],b[μ])
        end
    end

    function Base.similar(U::T) where T <: Gaugefields_4D_wing
        Uout = identityGaugefields_4D_wing(U.NC,U.NX,U.NY,U.NZ,U.NT,U.NDW)
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
        a.U[:,:,:,:,:,:] = copy(b.U)
        #error("substitute_U! is not implemented in type $(typeof(a)) ")
        return 
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





    
    function Base.setindex!(U::Shifted_Gaugefields{T,4},v,i1,i2,i3,i4,i5,i6)  where T <: Gaugefields_4D_wing
        error("type $(typeof(U)) has no setindex method. This type is read only.")
    end

    function Base.getindex(U::Shifted_Gaugefields{T,4},i1,i2,i3,i4,i5,i6) where T <: Gaugefields_4D_wing
        return U.parent.U[i1,i2,i3 .+ U.parent.NDW .+ U.shift[1],i4 .+ U.parent.NDW .+ U.shift[2],i5 .+ U.parent.NDW .+ U.shift[3],i6 .+ U.parent.NDW .+ U.shift[4]]
    end

    function Base.getindex(U::Adjoint_Gaugefields{Shifted_Gaugefields{T,4}},i1,i2,i3,i4,i5,i6) where T <: Gaugefields_4D_wing
        return conj(U.parent[i2,i1,i3,i4,i5,i6])
    end

    function Base.setindex!(U::Adjoint_Gaugefields{Shifted_Gaugefields{T,4}},v,i1,i2,i3,i4,i5,i6)  where T <: Gaugefields_4D_wing
        error("type $(typeof(U)) has no setindex method. This type is read only.")
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
                                u[k1,k2,-NDW+id,iy,iz,it] = u[k1,k2,NX+(id-NDW),iy,iz,it]
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
                                u[k1,k2,NX+id,iy,iz,it] = u[k1,k2,id,iy,iz,it]
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
                                u[k1,k2,ix,-NDW+id,iz,it] = u[k1,k2,ix,NY+(id-NDW),iz,it]
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
                                u[k1,k2,ix,NY+id,iz,it] = u[k1,k2,ix,id,iz,it]
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
                                u[k1,k2,ix,iy,id-NDW,it] = u[k1,k2,ix,iy,NZ+(id-NDW),it]
                                u[k1,k2,ix,iy,NZ+id,it] = u[k1,k2,ix,iy,id,it]
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
                                u[k1,k2,ix,iy,iz,id-NDW] = u[k1,k2,ix,iy,iz,NT+(id-NDW)]
                                u[k1,k2,ix,iy,iz,NT+id] = u[k1,k2,ix,iy,iz,id]
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


    """
-----------------------------------------------------c
     !!!!!   vin and vout should be different vectors

     Projectin of the etraceless antiermite part 
     vout = x/2 - Tr(x)/6
     wher   x = vin - Conjg(vin)      
-----------------------------------------------------c
    """
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

    
    


# end