module Gaugefields_4D_module
    using LinearAlgebra
    import ..AbstractGaugefields_module:AbstractGaugefields,Shifted_Gaugefields,shift_U,
                        Adjoint_Gaugefields,set_wing_U!,Abstractfields,construct_staple!,clear_U!,
                        calculate_Plaquet
    import Base

    abstract type Gaugefields_4D{NC} <: AbstractGaugefields{NC,4}
    end


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

    function clear_U!(U::Array{T,1}) where T <: Gaugefields_4D_wing
        for μ=1:4
            clear_U!(U[μ])
        end
        
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
                                Uμ[k1,k2,ix,iy,iz,it] = 0
                            end
                        end
                    end
                end
            end
        end
        set_wing_U!(Uμ)

    end

    function calculate_Plaquet(U::Array{T,1}) where T <: Gaugefields_4D
        error("calculate_Plaquet is not implemented in type $(typeof(U)) ")
    end

    #lattice shift
    function shift_U(U::Gaugefields_4D_wing,ν::T) where T <: Integer
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

    function shift_U(U::Gaugefields_4D_wing,shift::NTuple{Dim,T}) where {Dim,T <: Integer}
        Shifted_Gaugefields(U,shift,4)
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
        @inbounds for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX
                        c[1,1,ix,iy,iz,it] = a[1,1,ix,iy,iz,it]*b[1,1,ix,iy,iz,it] + 
                                                a[1,2,ix,iy,iz,it]*b[2,1,ix,iy,iz,it]
                        c[1,2,ix,iy,iz,it] = a[1,1,ix,iy,iz,it]*b[1,2,ix,iy,iz,it] + 
                                                a[1,2,ix,iy,iz,it]*b[2,2,ix,iy,iz,it]
                        c[2,1,ix,iy,iz,it] = a[2,1,ix,iy,iz,it]*b[1,1,ix,iy,iz,it] + 
                                                a[2,2,ix,iy,iz,it]*b[2,1,ix,iy,iz,it]         
                        c[2,2,ix,iy,iz,it] = a[2,1,ix,iy,iz,it]*b[1,2,ix,iy,iz,it] + 
                                                a[2,2,ix,iy,iz,it]*b[2,2,ix,iy,iz,it]                                                                               
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
        if β != 0
            @inbounds for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        @simd for ix=1:NX
                            c[1,1,ix,iy,iz,it] *= β
                            c[1,2,ix,iy,iz,it] *= β
                            c[2,1,ix,iy,iz,it] *= β    
                            c[2,2,ix,iy,iz,it] *= β                                                                            
                        end
                    end
                end
            end
        end

        if α != 0
            @inbounds for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        @simd for ix=1:NX
                            c[1,1,ix,iy,iz,it] += (a[1,1,ix,iy,iz,it]*b[1,1,ix,iy,iz,it] + 
                                                    a[1,2,ix,iy,iz,it]*b[2,1,ix,iy,iz,it])*α 
                            c[1,2,ix,iy,iz,it] += (a[1,1,ix,iy,iz,it]*b[1,2,ix,iy,iz,it] + 
                                                    a[1,2,ix,iy,iz,it]*b[2,2,ix,iy,iz,it])*α 
                            c[2,1,ix,iy,iz,it] += (a[2,1,ix,iy,iz,it]*b[1,1,ix,iy,iz,it] + 
                                                    a[2,2,ix,iy,iz,it]*b[2,1,ix,iy,iz,it])*α    
                            c[2,2,ix,iy,iz,it] += (a[2,1,ix,iy,iz,it]*b[1,2,ix,iy,iz,it] + 
                                                    a[2,2,ix,iy,iz,it]*b[2,2,ix,iy,iz,it])*α                                                                            
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
                        c[1,1,ix,iy,iz,it] = a[1,1,ix,iy,iz,it]*b[1,1,ix,iy,iz,it] + 
                                                a[1,2,ix,iy,iz,it]*b[2,1,ix,iy,iz,it]+a[1,3,ix,iy,iz,it]*b[3,1,ix,iy,iz,it]
                        c[1,2,ix,iy,iz,it] = a[1,1,ix,iy,iz,it]*b[1,2,ix,iy,iz,it] + 
                                                a[1,2,ix,iy,iz,it]*b[2,2,ix,iy,iz,it] + a[1,3,ix,iy,iz,it]*b[3,2,ix,iy,iz,it]
                        c[1,3,ix,iy,iz,it] = a[1,1,ix,iy,iz,it]*b[1,3,ix,iy,iz,it] + 
                                                a[1,2,ix,iy,iz,it]*b[2,3,ix,iy,iz,it] + a[1,3,ix,iy,iz,it]*b[3,3,ix,iy,iz,it]
                        c[2,1,ix,iy,iz,it] = a[2,1,ix,iy,iz,it]*b[1,1,ix,iy,iz,it] + 
                                                a[2,2,ix,iy,iz,it]*b[2,1,ix,iy,iz,it] + a[2,3,ix,iy,iz,it]*b[3,1,ix,iy,iz,it]         
                        c[2,2,ix,iy,iz,it] = a[2,1,ix,iy,iz,it]*b[1,2,ix,iy,iz,it] + 
                                                a[2,2,ix,iy,iz,it]*b[2,2,ix,iy,iz,it] + a[2,3,ix,iy,iz,it]*b[3,2,ix,iy,iz,it]                                                                            
                        c[2,3,ix,iy,iz,it] = a[2,1,ix,iy,iz,it]*b[1,3,ix,iy,iz,it] + 
                                                a[2,2,ix,iy,iz,it]*b[2,3,ix,iy,iz,it] + a[2,3,ix,iy,iz,it]*b[3,3,ix,iy,iz,it]                                                                            
                        c[3,1,ix,iy,iz,it] = a[3,1,ix,iy,iz,it]*b[1,1,ix,iy,iz,it] + 
                                                a[3,2,ix,iy,iz,it]*b[2,1,ix,iy,iz,it] + a[3,3,ix,iy,iz,it]*b[3,1,ix,iy,iz,it]         
                        c[3,2,ix,iy,iz,it] = a[3,1,ix,iy,iz,it]*b[1,2,ix,iy,iz,it] + 
                                                a[3,2,ix,iy,iz,it]*b[2,2,ix,iy,iz,it] + a[3,3,ix,iy,iz,it]*b[3,2,ix,iy,iz,it]                                                                            
                        c[3,3,ix,iy,iz,it] = a[3,1,ix,iy,iz,it]*b[1,3,ix,iy,iz,it] + 
                                                a[3,2,ix,iy,iz,it]*b[2,3,ix,iy,iz,it] + a[3,3,ix,iy,iz,it]*b[3,3,ix,iy,iz,it] 
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
                            c[1,1,ix,iy,iz,it] = (a[1,1,ix,iy,iz,it]*b[1,1,ix,iy,iz,it] + 
                                                    a[1,2,ix,iy,iz,it]*b[2,1,ix,iy,iz,it]+a[1,3,ix,iy,iz,it]*b[3,1,ix,iy,iz,it])*α + β*c[1,1,ix,iy,iz,it]
                            c[1,2,ix,iy,iz,it] = (a[1,1,ix,iy,iz,it]*b[1,2,ix,iy,iz,it] + 
                                                    a[1,2,ix,iy,iz,it]*b[2,2,ix,iy,iz,it] + a[1,3,ix,iy,iz,it]*b[3,2,ix,iy,iz,it])*α  + β*c[1,2,ix,iy,iz,it]
                            c[1,3,ix,iy,iz,it] = (a[1,1,ix,iy,iz,it]*b[1,3,ix,iy,iz,it] + 
                                                    a[1,2,ix,iy,iz,it]*b[2,3,ix,iy,iz,it] + a[1,3,ix,iy,iz,it]*b[3,3,ix,iy,iz,it])*α + β*c[1,3,ix,iy,iz,it]
                            c[2,1,ix,iy,iz,it] = (a[2,1,ix,iy,iz,it]*b[1,1,ix,iy,iz,it] + 
                                                    a[2,2,ix,iy,iz,it]*b[2,1,ix,iy,iz,it] + a[2,3,ix,iy,iz,it]*b[3,1,ix,iy,iz,it])*α   + β*c[2,1,ix,iy,iz,it]
                            c[2,2,ix,iy,iz,it] = (a[2,1,ix,iy,iz,it]*b[1,2,ix,iy,iz,it] + 
                                                    a[2,2,ix,iy,iz,it]*b[2,2,ix,iy,iz,it] + a[2,3,ix,iy,iz,it]*b[3,2,ix,iy,iz,it])*α    + β*c[2,2,ix,iy,iz,it]                                                                   
                            c[2,3,ix,iy,iz,it] = (a[2,1,ix,iy,iz,it]*b[1,3,ix,iy,iz,it] + 
                                                    a[2,2,ix,iy,iz,it]*b[2,3,ix,iy,iz,it] + a[2,3,ix,iy,iz,it]*b[3,3,ix,iy,iz,it] )*α   + β*c[2,3,ix,iy,iz,it]                                                                      
                            c[3,1,ix,iy,iz,it] = (a[3,1,ix,iy,iz,it]*b[1,1,ix,iy,iz,it] + 
                                                    a[3,2,ix,iy,iz,it]*b[2,1,ix,iy,iz,it] + a[3,3,ix,iy,iz,it]*b[3,1,ix,iy,iz,it])*α    + β*c[3,1,ix,iy,iz,it]    
                            c[3,2,ix,iy,iz,it] = (a[3,1,ix,iy,iz,it]*b[1,2,ix,iy,iz,it] + 
                                                    a[3,2,ix,iy,iz,it]*b[2,2,ix,iy,iz,it] + a[3,3,ix,iy,iz,it]*b[3,2,ix,iy,iz,it])*α   + β*c[3,2,ix,iy,iz,it]                                                                       
                            c[3,3,ix,iy,iz,it] = (a[3,1,ix,iy,iz,it]*b[1,3,ix,iy,iz,it] + 
                                                    a[3,2,ix,iy,iz,it]*b[2,3,ix,iy,iz,it] + a[3,3,ix,iy,iz,it]*b[3,3,ix,iy,iz,it])*α + β*c[3,3,ix,iy,iz,it]
                        end
                    end
                end
            end
        end
    end
    

end