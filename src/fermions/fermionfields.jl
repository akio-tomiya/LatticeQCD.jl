module Fermionfields
    using Random
    using LinearAlgebra
    using SparseArrays


    import ..Actions:FermiActionParam,FermiActionParam_Wilson,
                FermiActionParam_WilsonClover,FermiActionParam_Staggered
    import ..Gaugefields:GaugeFields,GaugeFields_1d,SU3GaugeFields,SU2GaugeFields,SU3GaugeFields_1d,SU2GaugeFields_1d,
                            staggered_phase,SUn,SU2,SU3,SUNGaugeFields,SUNGaugeFields_1d
    #import ..Parallel:get_looprange

    #include("cg.jl")


    abstract type FermionFields end



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

    struct StaggeredFermion <: FermionFields
        NC::Int64
        NX::Int64
        NY::Int64
        NZ::Int64
        NT::Int64
        f::Array{ComplexF64,6}

        mass::Float64
        eps::Float64
        Dirac_operator::String
        MaxCGstep::Int64
        BoundaryCondition::Array{Int8,1}
    end

    #include("boundary_fermi.jl")

    function FermionFields(NC,NX,NY,NZ,NT,fparam::FermiActionParam,BoundaryCondition)
        if findfirst("Wilson",fparam.Dirac_operator) != nothing #If string fparam.Dirac_operator includes "Wilson"
            #fparam.Dirac_operator == "Wilson"
            return WilsonFermion(NC,NX,NY,NZ,NT,fparam,BoundaryCondition)
        elseif fparam.Dirac_operator == "Staggered"
            return StaggeredFermion(NC,NX,NY,NZ,NT,fparam,BoundaryCondition)
        end
    end

    function WilsonFermion(NC,NX,NY,NZ,NT,fparam::FermiActionParam,BoundaryCondition) 
        r = fparam.r
        hop = fparam.hop
        eps = fparam.eps
        MaxCGstep = fparam.MaxCGstep
        return WilsonFermion(NC,NX,NY,NZ,NT,r,hop,eps,MaxCGstep,BoundaryCondition)
    end

    function StaggeredFermion(NC,NX,NY,NZ,NT,fparam::FermiActionParam,BoundaryCondition) 
        mass = fparam.mass
        eps = fparam.eps
        MaxCGstep = fparam.MaxCGstep
        return StaggeredFermion(NC,NX,NY,NZ,NT,mass,eps,MaxCGstep,BoundaryCondition)
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

    function StaggeredFermion(NC,NX,NY,NZ,NT,mass,eps,MaxCGstep,BoundaryCondition)#r,hop,eps,MaxCGstep)
        Dirac_operator = "Staggered"
        return StaggeredFermion(NC,NX,NY,NZ,NT,zeros(ComplexF64,NC,NX+2,NY+2,NZ+2,NT+2,1),
            mass,eps,Dirac_operator,MaxCGstep,BoundaryCondition)
    end

    function WilsonFermion(NC,NX,NY,NZ,NT,γ,rplusγ,rminusγ,hop,r,hopp,hopm,eps,fermion,MaxCGstep,BoundaryCondition)
        return WilsonFermion(NC,NX,NY,NZ,NT,zeros(ComplexF64,NC,NX+2,NY+2,NZ+2,NT+2,4),
                    γ,rplusγ,rminusγ,hop,r,hopp,hopm,eps,fermion,MaxCGstep,BoundaryCondition)
    end

    function StaggeredFermion(NC,NX,NY,NZ,NT,mass,eps,fermion,MaxCGstep,BoundaryCondition)
        return StaggeredFermion(NC,NX,NY,NZ,NT,zeros(ComplexF64,NC,NX+2,NY+2,NZ+2,NT+2,1),
                    mass,eps,fermion,MaxCGstep,BoundaryCondition)
    end

    function Base.setindex!(x::FermionFields,v,i1,i2,i3,i4,i5,i6) 
        x.f[i1,i2 + 1,i3 + 1,i4 + 1,i5 + 1,i6] = v
    end

    function Base.getindex(x::FermionFields,i1,i2,i3,i4,i5,i6)
        return x.f[i1,i2 .+ 1,i3 .+ 1,i4 .+ 1,i5 .+ 1,i6]
    end

    function Base.setindex!(x::WilsonFermion,v,i1,i2,i3,i4,i5,i6) 
        x.f[i1,i2 + 1,i3 + 1,i4 + 1,i5 + 1,i6] = v
    end

    function Base.getindex(x::WilsonFermion,i1,i2,i3,i4,i5,i6)
        return x.f[i1,i2 .+ 1,i3 .+ 1,i4 .+ 1,i5 .+ 1,i6]
    end

    function Base.:*(a::WilsonFermion,b::WilsonFermion)
        c = 0.0im
        for α=1:4
            for it=1:a.NT
                for iz=1:a.NZ
                    for iy=1:a.NY
                        for ix=1:a.NX
                            for ic=1:a.NC
                                c+= conj(a[ic,ix,iy,iz,it,α])*b[ic,ix,iy,iz,it,α]
                            end
                        end
                    end
                end
            end            
        end
        return c
    end

    function LinearAlgebra.dot(a::FermionFields,b::FermionFields)
        return a*b
    end

    function Base.:*(a::StaggeredFermion,b::StaggeredFermion)
        c = 0.0im
        α=1
        for it=1:a.NT
            for iz=1:a.NZ
                for iy=1:a.NY
                    for ix=1:a.NX
                        for ic=1:a.NC
                            c+= conj(a[ic,ix,iy,iz,it,α])*b[ic,ix,iy,iz,it,α]
                        end
                    end
                end
            end
        end            

        return c
    end


    function Base.similar(x::WilsonFermion)
        return WilsonFermion(x.NC,x.NX,x.NY,x.NZ,x.NT,
                    x.γ,x.rplusγ,x.rminusγ,x.hop,x.r,x.hopp,x.hopm,x.eps,x.Dirac_operator,x.MaxCGstep,x.BoundaryCondition)
    end

    function Base.similar(x::StaggeredFermion)
        return StaggeredFermion(x.NC,x.NX,x.NY,x.NZ,x.NT,
                    x.mass,x.eps,x.Dirac_operator,x.MaxCGstep,x.BoundaryCondition)
    end

    function clear!(a::FermionFields)
        n1,n2,n3,n4,n5,n6 = size(a.f)
        for i6=1:n6
            for i5=1:n5
                for i4=1:n4
                    for i3=1:n3
                        for i2=1:n2
                            for i1=1:n1
                                a.f[i1,i2,i3,i4,i5,i6]= 0
                            end
                        end
                    end
                end
            end
        end
    end

    function clear!(x::StaggeredFermion,evensite)
        ibush = ifelse(evensite,0,1)
        for it=1:x.NT
            for iz=1:x.NZ
                for iy=1:x.NY
                    xran =1+(1+ibush+iy+iz+it)%2:2:x.NX
                    for ix in xran
                        for ic=1:x.NC
                            x[ic,ix,iy,iz,it,1] = 0
                        end
                    end
                end
            end
        end
        return
    end



    
    function Base.display(x::FermionFields)
        
        for it=1:x.NT
            for iz=1:x.NZ
                for iy=1:x.NY
                    for ix=1:x.NX
                        for ic=1:x.NC
                            println("$ic $ix $iy $iz $it")
                            display(x[ic,ix,iy,iz,it,:])
                            println("\t")
                        end
                    end
                end

            end
        end
                        
    end

    function substitute_fermion!(H,j,x::WilsonFermion)
        i = 0
        for ialpha = 1:4
            for it=1:x.NT
                for iz=1:x.NZ
                    for iy=1:x.NY
                        for ix=1:x.NX
                            for ic=1:x.NC
                                i += 1
                                H[i,j] = x[ic,ix,iy,iz,it,ialpha]
                            end
                        end
                    end
                end
            end
        end
    end

    function substitute_fermion!(H,j,x::StaggeredFermion)
        i = 0
        for ialpha = 1:1
            for it=1:x.NT
                for iz=1:x.NZ
                    for iy=1:x.NY
                        for ix=1:x.NX
                            for ic=1:x.NC
                                i += 1
                                H[i,j] = x[ic,ix,iy,iz,it,ialpha]
                            end
                        end
                    end
                end
            end
        end
    end

    function add!(c::FermionFields,alpha::Number,a::FermionFields,beta::Number,b::FermionFields) #c = c + alpha*a + beta*b
        n1,n2,n3,n4,n5,n6 = size(a.f)

        for i6=1:n6
            for i5=1:n5
                for i4=1:n4
                    for i3=1:n3
                        for i2=1:n2
                            for i1=1:n1
                                #println(a.f[i1,i2,i3,i4,i5,i6],"\t",b.f[i1,i2,i3,i4,i5,i6] )
                                c.f[i1,i2,i3,i4,i5,i6] += alpha*a.f[i1,i2,i3,i4,i5,i6] +beta*b.f[i1,i2,i3,i4,i5,i6] 
                            end
                        end
                    end
                end
            end
        end
        return
    end

    function LinearAlgebra.axpby!(a::Number,X::FermionFields,Y::FermionFields)
        add!(Y,a,X)
        return
    end


    function add!(c::FermionFields,alpha::Number,a::FermionFields) #c = c + alpha*a 
        n1,n2,n3,n4,n5,n6 = size(a.f)

        for i6=1:n6
            for i5=1:n5
                for i4=1:n4
                    for i3=1:n3
                        for i2=1:n2
                            for i1=1:n1
                                #println(a.f[i1,i2,i3,i4,i5,i6],"\t",b.f[i1,i2,i3,i4,i5,i6] )
                                c.f[i1,i2,i3,i4,i5,i6] += alpha*a.f[i1,i2,i3,i4,i5,i6] 
                            end
                        end
                    end
                end
            end
        end
        return
    end

    function LinearAlgebra.axpby!(a::Number,X::FermionFields,b,Y::FermionFields)
        add!(b,Y,a,X)
        return
    end

    function add!(coeff::Number,c::FermionFields,alpha::Number,a::FermionFields) #c = coeff*c + alpha*a 
        n1,n2,n3,n4,n5,n6 = size(a.f)

        for i6=1:n6
            for i5=1:n5
                for i4=1:n4
                    for i3=1:n3
                        for i2=1:n2
                            for i1=1:n1
                                #println(a.f[i1,i2,i3,i4,i5,i6],"\t",b.f[i1,i2,i3,i4,i5,i6] )
                                c.f[i1,i2,i3,i4,i5,i6] = coeff*c.f[i1,i2,i3,i4,i5,i6] + alpha*a.f[i1,i2,i3,i4,i5,i6] 
                            end
                        end
                    end
                end
            end
        end
        return
    end

    function Wx!(xout::FermionFields,U::Array{GaugeFields,1},
                    x::FermionFields,temps::Array{FermionFields,1},fparam::T)  where T <: FermiActionParam
        println("This FermionFields type is not supported yet in function Wx!. The type is $(typeof(x))")
    end
    

    function Wdagx!(xout::FermionFields,U::Array{GaugeFields,1},
        x::FermionFields,temps::Array{FermionFields,1},fparam::T)  where T <: FermiActionParam
        println("This FermionFields type is not supported yet in function Wdagx!. The type is $(typeof(x))")
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

    """
    (D + m)*x
    """
    function Wx!(xout::StaggeredFermion,U::Array{G,1},
        x::StaggeredFermion,temps::Array{T,1},fparam::FermiActionParam_Staggered) where  {T <: FermionFields,G <: GaugeFields}
        
        temp = temps[4]
        Dx!(temp,U,x,temps,fparam)
        clear!(xout)
        add!(xout,fparam.mass,x,1,temp)
        set_wing_fermi!(xout)
        return
        
    end

    function Wx!(xout::StaggeredFermion,U::Array{G,1},
        x::StaggeredFermion,temps::Array{T,1},fparam::FermiActionParam_Staggered,indices) where  {T <: FermionFields,G <: GaugeFields}
        
        temp = temps[4]
        vec_indices1,vec_indices2 = Dx!(temp,U,x,temps,indices)
        clear!(xout)
        add!(xout,fparam.mass,x,1,temp)
        set_wing_fermi!(xout)
        return vec_indices1,vec_indices2
        
    end



    function Wx!(xout::StaggeredFermion,evensite,U::Array{G,1},
        x::StaggeredFermion,temps::Array{T,1},fparam::FermiActionParam_Staggered) where  {T <: FermionFields,G <: GaugeFields}
        
        
        temp = temps[4]
        Dx!(temp,evensite,U,x,temps,fparam)
        clear!(xout)
        add!(xout,fparam.mass,x,1,temp)
        set_wing_fermi!(xout)
        return
        
    end



    function Dx!(xout::StaggeredFermion,U::Array{G,1},
        x::StaggeredFermion,temps::Array{T,1}) where  {T <: FermionFields,G <: GaugeFields}
        #temp = temps[4]
        temp1 = temps[1]
        temp2 = temps[2]

        #clear!(temp)
        set_wing_fermi!(x)
        clear!(xout)
        for ν=1:4
            
            fermion_shift!(temp1,U,ν,x)

            fermion_shift!(temp2,U,-ν,x)

            add!(xout,0.5,temp1,-0.5,temp2)
            
        end
        set_wing_fermi!(xout)

        return
    end

    function Dx!(xout::StaggeredFermion,U::Array{G,1},
        x::StaggeredFermion,temps::Array{T,1},indices) where  {T <: FermionFields,G <: GaugeFields}
        #temp = temps[4]
        temp1 = temps[1]
        temp2 = temps[2]

        #clear!(temp)
        set_wing_fermi!(x)
        clear!(xout)
        vec_indices1 = Array{Tuple,1}(undef,4)
        vec_indices2 = Array{Tuple,1}(undef,4)
        
        for ν=1:4
            clear!(temp1)
            vec_indices1[ν] = fermion_shift!(temp1,U,ν,x,indices)

            #println(temp1*temp1)


            clear!(temp2)

            vec_indices2[ν] = fermion_shift!(temp2,U,-ν,x,indices)

            

            add!(xout,0.5,temp1,-0.5,temp2)
            
        end
        set_wing_fermi!(xout)

        return vec_indices1,vec_indices2
    end

    function DdagDx!(xout::StaggeredFermion,U::Array{G,1},
        x::StaggeredFermion,temps::Array{T,1},indices) where  {T <: FermionFields,G <: GaugeFields}
        #temp = temps[4]


        temp1 = temps[1]
        temp2 = temps[2]
        temp3 = temps[3]
        temp4 = temps[4]
        temp5 = temps[5]

        clear!(temp1)
        clear!(temp2)
        clear!(temp4)
        clear!(temp5)


        #clear!(temp)
        set_wing_fermi!(x)
        clear!(xout)
        vec_indices = Array{Tuple,1}(undef,8)
        #vec_indices2 = Array{Tuple,1}(undef,4)

        clear!(temp3)
        for ν=1:4
            vec_indices[2ν-1] = fermion_shift!(temp1,U,ν,x,[indices])[1]
            vec_indices[2ν] = fermion_shift!(temp2,U,-ν,x,[indices])[1]
            
        end
        add!(temp3,0.5,temp1,-0.5,temp2)

        for ν=1:4
            fermion_shift!(temp4,U,ν,temp3,vec_indices)
            fermion_shift!(temp5,U,-ν,temp3,vec_indices)
            
        end
        add!(xout,0.5,temp4,-0.5,temp5)

        set_wing_fermi!(xout)
        return

    end



    function Dx!(xout::StaggeredFermion,U::Array{G,1},
        x::StaggeredFermion,temps::Array{T,1},fparam::FermiActionParam_Staggered) where  {T <: FermionFields,G <: GaugeFields}
        Dx!(xout,U,x,temps)
        return
    end

    function Dx!(xout::StaggeredFermion,evensite,U::Array{G,1},
        x::StaggeredFermion,temps::Array{T,1},fparam::FermiActionParam_Staggered) where  {T <: FermionFields,G <: GaugeFields}
        #temp = temps[4]
        temp1 = temps[1]
        temp2 = temps[2]

        #clear!(temp)
        set_wing_fermi!(x)
        clear!(xout)
        for ν=1:4
            clear!(temp1)
            fermion_shift!(temp1,evensite,U,ν,x)

            #... Dirac multiplication
            #mul!(temp1,view(x.rminusγ,:,:,ν),temp1)
            
            #
            clear!(temp2)
            fermion_shift!(temp2,evensite,U,-ν,x)
            #mul!(temp2,view(x.rplusγ,:,:,ν),temp2)

            #add!(temp,1,temp1,-1,temp2)
            add!(xout,0.5,temp1,-0.5,temp2)
            
        end
        set_wing_fermi!(xout)

        #clear!(xout)
        #add!(xout,fparam.mass,x,1,temp)

        #display(xout)
        #    exit()
        return
    end

    

    function Wdagx!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1}) where {T <: FermionFields,G <: GaugeFields}
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

        mul_γ5x!(xout,temp1)
        return
    end

    function Wdagx!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1},fparam::FermiActionParam_Wilson) where {T <: FermionFields,G <: GaugeFields}
        Wdagx!(xout,U,x,temps)
        return
    end

    function Wdagx!(xout::StaggeredFermion,U::Array{G,1},
        x::StaggeredFermion,temps::Array{T,1},fparam::FermiActionParam_Staggered) where {T <: FermionFields,G <: GaugeFields}
        
        temp = temps[4]
        Dx!(temp,U,x,temps,fparam)
        clear!(xout)
        add!(xout,fparam.mass,x,-1,temp)
        return
        
    end

    function Wdagx!(xout::StaggeredFermion,U::Array{G,1},
        x::StaggeredFermion,temps::Array{T,1},fparam::FermiActionParam_Staggered,vec_indices1,vec_indices2) where {T <: FermionFields,G <: GaugeFields}
        
        temp = temps[4]
        Dx!(temp,U,x,temps,vec_indices1,vec_indices2)
        clear!(xout)
        add!(xout,fparam.mass,x,-1,temp)
        return
        
    end

    function WdagWx!(xout::T,U::Array{G,1},x::T,temps::Array{T,1},fparam) where {T <: FermionFields,G <: GaugeFields}

        temp = temps[5]
        Wx!(temp,U,x,temps,fparam) 
    
        Wdagx!(xout,U,temp,temps,fparam) 
    
        return
    end

    
    function WdagWx!(xout::T,U::Array{G,1},x::T,temps::Array{T,1},mass::Number,indices) where {T <: StaggeredFermion,G <: GaugeFields}

        temp = temps[6]
        #println("x = ", x*x)
        DdagDx!(temp,U,x,temps,indices)

        #println("temp ",temp*temp)
        clear!(xout)
        add!(xout,mass^2,x,-1,temp)
        #println("xout ",xout*xout)
        return
    end

    function WdagWx!(xout::T,U::Array{G,1},x::T,temps::Array{T,1},fparam::FermiActionParam_Staggered,indices) where {T <: StaggeredFermion,G <: GaugeFields}
        WdagWx!(xout,U,x,temps,fparam.mass,indices)
        return
    end

    function WdagWx!(xout::T,U::Array{G,1},x::T,temps::Array{T,1},fparam,indices) where {T <: FermionFields,G <: GaugeFields}
        WdagWx!(xout,U,x,temps,fparam)
        return
    end

    
    #=
    function WdagWx!(xout::StaggeredFermion,U::Array{G,1},
        x::StaggeredFermion,temps::Array{T,1},fparam::FermiActionParam_Staggered) where {T <: FermionFields,G <: GaugeFields}
        temp = temps[5]
        evensite = true
        Dx!(temp,evensite,U,x,temps,fparam)
        #Dx!(temp,U,x,temps,fparam)
        temp3 = temps[3]
        evensite = false
        Dx!(temp3,evensite,U,temp,temps,fparam)
        #Dx!(temp3,U,temp,temps,fparam)

        clear!(xout)
        add!(xout,fparam.mass^2,x,-1,temp3)

        return
    end
    =#
    

    function Ddagx!(xout::StaggeredFermion,U::Array{G,1},
        x::StaggeredFermion,temps::Array{T,1},fparam::FermiActionParam_Staggered) where  {T <: FermionFields,G <: GaugeFields}
        #temp = temps[4]
        temp1 = temps[1]
        temp2 = temps[2]

        #clear!(temp)
        clear!(xout)
        set_wing_fermi!(x)
        for ν=1:4
            fermion_shift!(temp1,U,ν,x)

            #... Dirac multiplication
            #mul!(temp1,view(x.rminusγ,:,:,ν),temp1)
            
            #
            fermion_shift!(temp2,U,-ν,x)
            #mul!(temp2,view(x.rplusγ,:,:,ν),temp2)

            #add!(temp,1,temp1,-1,temp2)

            add!(xout,-0.5,temp1,0.5,temp2)
            #add!(xout,1,temp1,-1,temp2)
            
        end
        set_wing_fermi!(xout)

        #clear!(xout)
        #add!(xout,fparam.mass,x,1,temp)

        #display(xout)
        #    exit()
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
                        for ix=1:NX
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





    function LinearAlgebra.mul!(xout::FermionFields,A::AbstractMatrix,x::FermionFields)
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        NC = x.NC

        n6 = size(x.f)[6]
        f = zeros(ComplexF64,n6)
        e = zeros(ComplexF64,n6)
        
        for ic=1:NC
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            for k1=1:n6
                                e[k1] = x[ic,ix,iy,iz,it,k1]
                            end

                            for k1=1:n6
                                f[k1] = 0
                                for k2=1:n6
                                    f[k1] += A[k1,k2]*e[k2]
                                end
                                xout[ic,ix,iy,iz,it,k1] = f[k1]
                            end

                        end
                    end
                end
            end
        end
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
                        for ix=1:NX
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

    function LinearAlgebra.mul!(xout::FermionFields,x::FermionFields,A::AbstractMatrix)
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        NC = x.NC

        n6 = size(x.f)[6]
        f = zeros(ComplexF64,n6)
        e = zeros(ComplexF64,n6)


        
        for ic=1:NC
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            
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
                        for ix=1:NX
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

    function LinearAlgebra.mul!(a::FermionFields,α::Number,b::FermionFields)
        n1,n2,n3,n4,n5,n6 = size(a.f)

        for i6=1:n6
            for i5=1:n5
                for i4=1:n4
                    for i3=1:n3
                        for i2=1:n2
                            for i1=1:n1
                                a.f[i1,i2,i3,i4,i5,i6] = α*b.f[i1,i2,i3,i4,i5,i6]
                            end
                        end
                    end
                end
            end
        end
        return 
    end
"""
    c-----------------------------------------------------c
c     iflag = 1
c       vv = |v1><v2| i.e., vv(a,b) = v1(a)*v2_adj(b)
c     iflag = 2
c                           vv(a,b) = v1(a)*v2(b)
c-----------------------------------------------------c
    """
    function vvmat!(vv::GaugeFields_1d{T},v1::FermionFields,v2::FermionFields,iflag) where T <: SUn 
        NX = v1.NX
        NY = v1.NY
        NZ = v1.NZ
        NT = v1.NT

        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            NC = v1.NC
            #error("NC >3 is not supported")
        end
        #else
        #    error("NC >3 is not supported")
        #end


        if iflag == 1
            
            is = 0
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            is += 1
                                for ia=1:NC
                                    for ib=1:NC
                                    vv[ia,ib,is] = 0
                                    for k=1:4
                                        vv[ia,ib,is] += v1[ia,ix,iy,iz,it,k]*conj(v2[ib,ix,iy,iz,it,k]) 
                                    end

                                end
                            end
                        end
                    end
                end
            end
        end
        if iflag == 2

                    is = 0

                    for it=1:NT
                        for iz=1:NZ
                            for iy=1:NY
                                for ix=1:NX
                                    is += 1
                                    for ia=1:NC
                                        for ib=1:NC
                                    vv[ia,ib,is] = 0
                                    for k=1:4
                                        vv[ia,ib,is] += v1[ia,ix,iy,iz,it,k]*v2[ib,ix,iy,iz,it,k]
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end


    function vvmat!(vv::GaugeFields_1d{T},v1::StaggeredFermion,v2::StaggeredFermion,iflag) where T <: SUn 
        NX = v1.NX
        NY = v1.NY
        NZ = v1.NZ
        NT = v1.NT

        if T == SU3
            NC = 3
        elseif T == SU2
            NC = 2
        else
            NC = v1.NC
            #error("NC >3 is not supported")
        end
        #else
        #    error("NC >3 is not supported")
        #end


        if iflag == 1
            
            is = 0
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            is += 1
                                for ia=1:NC
                                    for ib=1:NC
                                    vv[ia,ib,is] = 0
                                    for k=1:1
                                        vv[ia,ib,is] += v1[ia,ix,iy,iz,it,k]*conj(v2[ib,ix,iy,iz,it,k]) 
                                    end

                                end
                            end
                        end
                    end
                end
            end
        end
        if iflag == 2

                    is = 0

                    for it=1:NT
                        for iz=1:NZ
                            for iy=1:NY
                                for ix=1:NX
                                    is += 1
                                    for ia=1:NC
                                        for ib=1:NC
                                    vv[ia,ib,is] = 0
                                    for k=1:1
                                        vv[ia,ib,is] += v1[ia,ix,iy,iz,it,k]*v2[ib,ix,iy,iz,it,k]
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end


    function fermion_shift!(b::F,u::Array{T,1},μ,a::F) where {T <: GaugeFields,F <: FermionFields}
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
            idel = zeros(Int64,4)
            idel[μ] = 1

            n6 = size(a.f)[6]
            for ialpha=1:n6
                for it=1:NT
                    it1 = it + idel[4]
                    for iz=1:NZ
                        iz1 = iz + idel[3]
                        for iy=1:NY
                            iy1 = iy + idel[2]
                            for ix=1:NX
                                ix1 = ix + idel[1]
                                for k1=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:NC
                                        b[k1,ix,iy,iz,it,ialpha] += u[μ][k1,k2,ix,iy,iz,it]*a[k2,ix1,iy1,iz1,it1,ialpha]
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
        elseif μ < 0
            idel = zeros(Int64,4)
            idel[-μ] = 1
            n6 = size(b.f)[6]
            for ialpha =1:n6
                for it=1:NT
                    it1 = it - idel[4]
                    for iz=1:NZ
                        iz1 = iz -idel[3]
                        for iy=1:NY
                            iy1 = iy -idel[2]
                            for ix=1:NX
                                ix1 = ix -idel[1]
                        
                                for k1=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:NC
                                        b[k1,ix,iy,iz,it,ialpha] += conj(u[-μ][k2,k1,ix1,iy1,iz1,it1])*a[k2,ix1,iy1,iz1,it1,ialpha]
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

    end

    
    function fermion_shift!(b::F,evensite::Bool,u::Array{T,1},μ,a::F) where {T <: GaugeFields,F <: FermionFields}
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
            idel = zeros(Int64,4)
            idel[μ] = 1

            n6 = size(a.f)[6]
            for ialpha=1:n6
                for it=1:NT
                    it1 = it + idel[4]
                    for iz=1:NZ
                        iz1 = iz + idel[3]
                        for iy=1:NY
                            iy1 = iy + idel[2]
                            xran =1+(1+ibush+iy+iz+it)%2:2:NX
                            for ix in xran
                                ix1 = ix + idel[1]
                                for k1=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:NC
                                        b[k1,ix,iy,iz,it,ialpha] += u[μ][k1,k2,ix,iy,iz,it]*a[k2,ix1,iy1,iz1,it1,ialpha]
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
        elseif μ < 0
            idel = zeros(Int64,4)
            idel[-μ] = 1
            n6 = size(b.f)[6]
            for ialpha =1:n6
                for it=1:NT
                    it1 = it - idel[4]
                    for iz=1:NZ
                        iz1 = iz -idel[3]
                        for iy=1:NY
                            iy1 = iy -idel[2]
                            xran =1+(1+ibush+iy+iz+it)%2:2:NX
                            for ix in xran
                                ix1 = ix -idel[1]
                        
                                for k1=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:NC
                                        b[k1,ix,iy,iz,it,ialpha] += conj(u[-μ][k2,k1,ix1,iy1,iz1,it1])*a[k2,ix1,iy1,iz1,it1,ialpha]
                                    end
                                end
                            end
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
                            for ix=1:NX
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

    function fermion_shift!(b::StaggeredFermion,u::Array{T,1},μ::Int,a::StaggeredFermion) where T <: SU3GaugeFields
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
            for ialpha=1:1
                #for it=NTrange
                for it=1:NT
                    it1 = it + ifelse(μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz + ifelse(μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(μ ==2,1,0) #idel[2]
                            for ix=1:NX
                                ix1 = ix + ifelse(μ ==1,1,0) #idel[1]
                                η = staggered_phase(μ,ix,iy,iz,it,NX,NY,NZ,NT)

                                b[1,ix,iy,iz,it,ialpha] = η*(u[μ][1,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][1,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][1,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha])

                                b[2,ix,iy,iz,it,ialpha] = η*(u[μ][2,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][2,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][2,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha])

                                b[3,ix,iy,iz,it,ialpha] = η*(u[μ][3,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][3,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][3,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha])

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
            for ialpha =1:1
                for it=1:NT
                    it1 = it - ifelse(-μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz - ifelse(-μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy - ifelse(-μ ==2,1,0)  #idel[2]
                            for ix=1:NX
                                ix1 = ix - ifelse(-μ ==1,1,0) #idel[1]

                                η = staggered_phase(-μ,ix1,iy1,iz1,it1,NX,NY,NZ,NT)

                                b[1,ix,iy,iz,it,ialpha] = η*(conj(u[-μ][1,1,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,1,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][3,1,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha])

                                b[2,ix,iy,iz,it,ialpha] = η*(conj(u[-μ][1,2,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,2,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][3,2,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha])

                                b[3,ix,iy,iz,it,ialpha] = η*(conj(u[-μ][1,3,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,3,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][3,3,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha])
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

    function fermion_shift!(b::StaggeredFermion,u::Array{T,1},μ::Int,a::StaggeredFermion) where T <: SUNGaugeFields
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
            for ialpha=1:1
                #for it=NTrange
                for it=1:NT
                    it1 = it + ifelse(μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz + ifelse(μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(μ ==2,1,0) #idel[2]
                            for ix=1:NX
                                ix1 = ix + ifelse(μ ==1,1,0) #idel[1]
                                η = staggered_phase(μ,ix,iy,iz,it,NX,NY,NZ,NT)
                                
                                for k1=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:NC
                                        b[k1,ix,iy,iz,it,ialpha] += η*u[μ][k1,k2,ix,iy,iz,it]*a[k2,ix1,iy1,iz1,it1,ialpha]
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
            
        elseif μ < 0
            #idel = zeros(Int64,4)
            #idel[-μ] = 1
            #n6 = size(b.f)[6]
            for ialpha =1:1
                for it=1:NT
                    it1 = it - ifelse(-μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz - ifelse(-μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy - ifelse(-μ ==2,1,0)  #idel[2]
                            for ix=1:NX
                                ix1 = ix - ifelse(-μ ==1,1,0) #idel[1]

                                η = staggered_phase(-μ,ix1,iy1,iz1,it1,NX,NY,NZ,NT)

                                
                                for k1=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:NC
                                        b[k1,ix,iy,iz,it,ialpha] += η*conj(u[-μ][k2,k1,ix1,iy1,iz1,it1])*a[k2,ix1,iy1,iz1,it1,ialpha]
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
        end

    end

    function fermion_shift!(b::StaggeredFermion,evensite,u::Array{T,1},μ::Int,a::StaggeredFermion) where T <: SU3GaugeFields
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

        #NTrange = get_looprange(NT)
        #println(NTrange)
        if μ > 0
            #idel = zeros(Int64,4)
            #idel[μ] = 1

            n6 = size(a.f)[6]
            for ialpha=1:1
                #for it=NTrange
                for it=1:NT
                    it1 = it + ifelse(μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz + ifelse(μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(μ ==2,1,0) #idel[2]
                            xran =1+(1+ibush+iy+iz+it)%2:2:NX
                            for ix in xran
                            #for ix=1:NX
                                ix1 = ix + ifelse(μ ==1,1,0) #idel[1]
                                η = staggered_phase(μ,ix,iy,iz,it,NX,NY,NZ,NT)

                                b[1,ix,iy,iz,it,ialpha] = η*(u[μ][1,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][1,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][1,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha])

                                b[2,ix,iy,iz,it,ialpha] = η*(u[μ][2,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][2,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][2,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha])

                                b[3,ix,iy,iz,it,ialpha] = η*(u[μ][3,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][3,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][3,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha])

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
            for ialpha =1:1
                for it=1:NT
                    it1 = it - ifelse(-μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz - ifelse(-μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy - ifelse(-μ ==2,1,0)  #idel[2]
                            xran =1+(1+ibush+iy+iz+it)%2:2:NX
                            for ix in xran
                                ix1 = ix - ifelse(-μ ==1,1,0) #idel[1]

                                η = staggered_phase(-μ,ix1,iy1,iz1,it1,NX,NY,NZ,NT)

                                b[1,ix,iy,iz,it,ialpha] = η*(conj(u[-μ][1,1,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,1,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][3,1,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha])

                                b[2,ix,iy,iz,it,ialpha] = η*(conj(u[-μ][1,2,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,2,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][3,2,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha])

                                b[3,ix,iy,iz,it,ialpha] = η*(conj(u[-μ][1,3,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,3,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][3,3,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha])
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

    
    function fermion_shift!(b::StaggeredFermion,u::Array{T,1},μ::Int,a::StaggeredFermion) where T <: SU2GaugeFields
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
            for ialpha=1:1
                #for it=NTrange
                for it=1:NT
                    it1 = it + ifelse(μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz + ifelse(μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(μ ==2,1,0) #idel[2]
                            for ix=1:NX
                                ix1 = ix + ifelse(μ ==1,1,0) #idel[1]
                                η = staggered_phase(μ,ix,iy,iz,it,NX,NY,NZ,NT)

                                b[1,ix,iy,iz,it,ialpha] = η*(u[μ][1,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][1,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] )

                                b[2,ix,iy,iz,it,ialpha] = η*(u[μ][2,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][2,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] )


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
            for ialpha =1:1
                for it=1:NT
                    it1 = it - ifelse(-μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz - ifelse(-μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy - ifelse(-μ ==2,1,0)  #idel[2]
                            for ix=1:NX
                                ix1 = ix - ifelse(-μ ==1,1,0) #idel[1]
                                η = staggered_phase(-μ,ix1,iy1,iz1,it1,NX,NY,NZ,NT)

                                b[1,ix,iy,iz,it,ialpha] = η*(conj(u[-μ][1,1,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,1,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] )

                                b[2,ix,iy,iz,it,ialpha] = η*(conj(u[-μ][1,2,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,2,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] )

                            end
                        end
                    end
                end
            end
        end

    end

    function apply_periodicity(ix,iy,iz,it,NX,NY,NZ,NT,BC)
        
        ixflag1 = (ix < 1)
        ixflag2 = (ix > NX)
        ix += ifelse(ixflag1,NX,0)
        ix += ifelse(ixflag2,-NX,0)

        sign =  ifelse(ixflag1 || ixflag2,BC[1],1)
        #sign =  ifelse(ixflag2,BC[1],1)



        iyflag1 = (iy < 1)
        iyflag2 = (iy > NY)
        
        iy += ifelse(iyflag1,NY,0)
        iy += ifelse(iyflag2,-NY,0)

        sign = ifelse(iyflag1 || iyflag2,BC[2],1)
        #sign = ifelse(iyflag2,BC[2],1)

        izflag1 = (iz < 1)
        izflag2 = (iz > NZ)

        iz += ifelse(izflag1,NZ,0)
        iz += ifelse(izflag2,-NZ,0)

        sign = ifelse(izflag1 || izflag2,BC[3],1)
        #sign = ifelse(izflag2,BC[3],1)

        itflag1 = (it < 1)
        itflag2 = (it > NT)

        it += ifelse(itflag1,NT,0)
        it += ifelse(itflag2,-NT,0)

        sign = ifelse(itflag1 || itflag2,BC[4],1)
        #sign = ifelse(itflag2,BC[4],1)

        return ix,iy,iz,it,sign
    end


    function fermion_shift!(b::StaggeredFermion,u::Array{T,1},μ::Int,a::StaggeredFermion,vec_indices) where T <: SU2GaugeFields
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
        outindices = Array{Tuple,1}(undef,length(vec_indices))
        #clear!(b)
        if μ > 0
            for (k,indices) in enumerate(vec_indices)

                ix1 = indices[1]
                iy1 = indices[2]
                iz1 = indices[3]
                it1 = indices[4]
                ialpha = indices[5]

                ix = ix1 - ifelse(μ ==1,1,0)
                iy = iy1 - ifelse(μ ==2,1,0)
                iz = iz1 - ifelse(μ ==3,1,0)
                it = it1 - ifelse(μ ==4,1,0)

                η = staggered_phase(μ,ix,iy,iz,it,NX,NY,NZ,NT)

                ix,iy,iz,it,sign =  apply_periodicity(ix,iy,iz,it,NX,NY,NZ,NT,b.BoundaryCondition)


                outindices[k] = (ix,iy,iz,it,ialpha)

                

                b[1,ix,iy,iz,it,ialpha] += sign*η*(u[μ][1,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                            u[μ][1,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] )

                b[2,ix,iy,iz,it,ialpha] += sign*η*(u[μ][2,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                            u[μ][2,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] )
            end
            
        elseif μ < 0
            for (k,indices) in enumerate(vec_indices)

                ix1 = indices[1]
                iy1 = indices[2]
                iz1 = indices[3]
                it1 = indices[4]
                ialpha = indices[5]

                ix = ix1 + ifelse(-μ ==1,1,0)
                iy = iy1 + ifelse(-μ ==2,1,0)
                iz = iz1 + ifelse(-μ ==3,1,0)
                it = it1 + ifelse(-μ ==4,1,0)

                η = staggered_phase(-μ,ix1,iy1,iz1,it1,NX,NY,NZ,NT)

                ix,iy,iz,it,sign =  apply_periodicity(ix,iy,iz,it,NX,NY,NZ,NT,b.BoundaryCondition)

                outindices[k] = (ix,iy,iz,it,ialpha)


                b[1,ix,iy,iz,it,ialpha] += sign*η*(conj(u[-μ][1,1,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                            conj(u[-μ][2,1,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] )

                b[2,ix,iy,iz,it,ialpha] += sign*η*(conj(u[-μ][1,2,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                            conj(u[-μ][2,2,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] )
            end

        end

        return outindices

    end

    function fermion_shift!(b::StaggeredFermion,u::Array{T,1},μ::Int,a::StaggeredFermion,vec_indices) where T <: SU3GaugeFields
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
        outindices = Array{Tuple,1}(undef,length(vec_indices))
        #clear!(b)
        if μ > 0
            for (k,indices) in enumerate(vec_indices)

                ix1 = indices[1]
                iy1 = indices[2]
                iz1 = indices[3]
                it1 = indices[4]
                ialpha = indices[5]

                ix = ix1 - ifelse(μ ==1,1,0)
                iy = iy1 - ifelse(μ ==2,1,0)
                iz = iz1 - ifelse(μ ==3,1,0)
                it = it1 - ifelse(μ ==4,1,0)

                η = staggered_phase(μ,ix,iy,iz,it,NX,NY,NZ,NT)

                ix,iy,iz,it,sign =  apply_periodicity(ix,iy,iz,it,NX,NY,NZ,NT,b.BoundaryCondition)


                outindices[k] = (ix,iy,iz,it,ialpha)

                

                b[1,ix,iy,iz,it,ialpha] += sign*η*(u[μ][1,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                            u[μ][1,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                            u[μ][1,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha] 
                                            )

                b[2,ix,iy,iz,it,ialpha] += sign*η*(u[μ][2,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                            u[μ][2,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                            u[μ][2,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha] )
                b[3,ix,iy,iz,it,ialpha] += sign*η*(u[μ][3,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                            u[μ][3,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                            u[μ][3,3,ix,iy,iz,it]*a[3,ix1,iy1,iz1,it1,ialpha] )
            end
            
        elseif μ < 0
            for (k,indices) in enumerate(vec_indices)

                ix1 = indices[1]
                iy1 = indices[2]
                iz1 = indices[3]
                it1 = indices[4]
                ialpha = indices[5]

                ix = ix1 + ifelse(-μ ==1,1,0)
                iy = iy1 + ifelse(-μ ==2,1,0)
                iz = iz1 + ifelse(-μ ==3,1,0)
                it = it1 + ifelse(-μ ==4,1,0)

                η = staggered_phase(-μ,ix1,iy1,iz1,it1,NX,NY,NZ,NT)

                ix,iy,iz,it,sign =  apply_periodicity(ix,iy,iz,it,NX,NY,NZ,NT,b.BoundaryCondition)

                outindices[k] = (ix,iy,iz,it,ialpha)


                b[1,ix,iy,iz,it,ialpha] += sign*η*(conj(u[-μ][1,1,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                            conj(u[-μ][2,1,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] + 
                                            conj(u[-μ][3,1,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha] )

                b[2,ix,iy,iz,it,ialpha] += sign*η*(conj(u[-μ][1,2,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                            conj(u[-μ][2,2,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha]+ 
                                            conj(u[-μ][3,2,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha] )
                b[3,ix,iy,iz,it,ialpha] += sign*η*(conj(u[-μ][1,3,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                            conj(u[-μ][2,3,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha]+ 
                                            conj(u[-μ][3,3,ix1,iy1,iz1,it1])*a[3,ix1,iy1,iz1,it1,ialpha] )
            end

        end

        return outindices

    end

    function fermion_shift!(b::StaggeredFermion,u::Array{T,1},μ::Int,a::StaggeredFermion,vec_indices) where T <: SUNGaugeFields
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
        outindices = Array{Tuple,1}(undef,length(vec_indices))
        #clear!(b)
        if μ > 0
            for (k,indices) in enumerate(vec_indices)

                ix1 = indices[1]
                iy1 = indices[2]
                iz1 = indices[3]
                it1 = indices[4]
                ialpha = indices[5]

                ix = ix1 - ifelse(μ ==1,1,0)
                iy = iy1 - ifelse(μ ==2,1,0)
                iz = iz1 - ifelse(μ ==3,1,0)
                it = it1 - ifelse(μ ==4,1,0)

                η = staggered_phase(μ,ix,iy,iz,it,NX,NY,NZ,NT)

                ix,iy,iz,it,sign =  apply_periodicity(ix,iy,iz,it,NX,NY,NZ,NT,b.BoundaryCondition)


                outindices[k] = (ix,iy,iz,it,ialpha)

                for k1=1:NC
                    b[k1,ix,iy,iz,it,ialpha] = 0
                    for k2=1:NC
                        b[k1,ix,iy,iz,it,ialpha] += sign*η*u[μ][k1,k2,ix,iy,iz,it]*a[k2,ix1,iy1,iz1,it1,ialpha]
                    end
                end

                
            end
            
        elseif μ < 0
            for (k,indices) in enumerate(vec_indices)

                ix1 = indices[1]
                iy1 = indices[2]
                iz1 = indices[3]
                it1 = indices[4]
                ialpha = indices[5]

                ix = ix1 + ifelse(-μ ==1,1,0)
                iy = iy1 + ifelse(-μ ==2,1,0)
                iz = iz1 + ifelse(-μ ==3,1,0)
                it = it1 + ifelse(-μ ==4,1,0)

                η = staggered_phase(-μ,ix1,iy1,iz1,it1,NX,NY,NZ,NT)

                ix,iy,iz,it,sign =  apply_periodicity(ix,iy,iz,it,NX,NY,NZ,NT,b.BoundaryCondition)

                outindices[k] = (ix,iy,iz,it,ialpha)

                for k1=1:NC
                    b[k1,ix,iy,iz,it,ialpha] = 0
                    for k2=1:NC
                        b[k1,ix,iy,iz,it,ialpha] += sign*η*conj(u[-μ][k2,k1,ix1,iy1,iz1,it1])*a[k2,ix1,iy1,iz1,it1,ialpha]
                    end
                end


            end

        end

        return outindices

    end


    function fermion_shift!(b::StaggeredFermion,evensite,u::Array{T,1},μ::Int,a::StaggeredFermion) where T <: SU2GaugeFields
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

        #NTrange = get_looprange(NT)
        #println(NTrange)
        if μ > 0
            #idel = zeros(Int64,4)
            #idel[μ] = 1

            n6 = size(a.f)[6]
            for ialpha=1:1
                #for it=NTrange
                for it=1:NT
                    it1 = it + ifelse(μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz + ifelse(μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(μ ==2,1,0) #idel[2]
                            xran =1+(1+ibush+iy+iz+it)%2:2:NX
                            for ix in xran
                                ix1 = ix + ifelse(μ ==1,1,0) #idel[1]
                                η = staggered_phase(μ,ix,iy,iz,it,NX,NY,NZ,NT)

                                b[1,ix,iy,iz,it,ialpha] = η*(u[μ][1,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][1,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] )

                                b[2,ix,iy,iz,it,ialpha] = η*(u[μ][2,1,ix,iy,iz,it]*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            u[μ][2,2,ix,iy,iz,it]*a[2,ix1,iy1,iz1,it1,ialpha] )


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
            for ialpha =1:1
                for it=1:NT
                    it1 = it - ifelse(-μ ==4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz - ifelse(-μ ==3,1,0) #idel[3]
                        for iy=1:NY
                            iy1 = iy - ifelse(-μ ==2,1,0)  #idel[2]
                            xran =1+(1+ibush+iy+iz+it)%2:2:NX
                            for ix in xran
                                ix1 = ix - ifelse(-μ ==1,1,0) #idel[1]
                                η = staggered_phase(-μ,ix1,iy1,iz1,it1,NX,NY,NZ,NT)

                                b[1,ix,iy,iz,it,ialpha] = η*(conj(u[-μ][1,1,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,1,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] )

                                b[2,ix,iy,iz,it,ialpha] = η*(conj(u[-μ][1,2,ix1,iy1,iz1,it1])*a[1,ix1,iy1,iz1,it1,ialpha] + 
                                                            conj(u[-μ][2,2,ix1,iy1,iz1,it1])*a[2,ix1,iy1,iz1,it1,ialpha] )

                            end
                        end
                    end
                end
            end
        end

    end

    function fermion_shiftB!(b::WilsonFermion,u::Array{SUNGaugeFields,1},μ,a::WilsonFermion) 
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

    function fermion_shiftB!(b::WilsonFermion,evensite,u::Array{SUNGaugeFields,1},μ,a::WilsonFermion) 
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

    function fermion_shiftB!(b::StaggeredFermion,u::Array{SUNGaugeFields,1},μ,a::StaggeredFermion) 
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
            for ialpha =1:1
                for it=1:NT
                    it1 = it + ifelse(-μ == 4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz +ifelse(-μ == 3,1,0)  #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(-μ == 2,1,0) #idel[2]
                            for ix=1:NX
                                ix1 = ix + ifelse(-μ == 1,1,0)  #idel[1]
                                η = staggered_phase(-μ,ix,iy,iz,it,NX,NY,NZ,NT)
                                
                                for k1=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:NC
                                        b[k1,ix,iy,iz,it,ialpha] += η*conj(a[k2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][k1,k2,ix,iy,iz,it])
                                    end
                                end
                                
                                
                            end
                        end
                    end
                end
            end
        end

    end

    function fermion_shiftB!(b::StaggeredFermion,evensite,u::Array{SUNGaugeFields,1},μ,a::StaggeredFermion) 
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
            for ialpha =1:1
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
                                η = staggered_phase(-μ,ix,iy,iz,it,NX,NY,NZ,NT)
                                
                                for k1=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] = 0
                                    for k2=1:NC
                                        b[k1,ix,iy,iz,it,ialpha] += η*conj(a[k2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][k1,k2,ix,iy,iz,it])
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


    function fermion_shiftB!(b::StaggeredFermion,u::Array{SU3GaugeFields,1},μ,a::StaggeredFermion) 
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
            for ialpha =1:1
                for it=1:NT
                    it1 = it + ifelse(-μ == 4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz +ifelse(-μ == 3,1,0)  #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(-μ == 2,1,0) #idel[2]
                            for ix=1:NX
                                ix1 = ix + ifelse(-μ == 1,1,0)  #idel[1]
                                η = staggered_phase(-μ,ix,iy,iz,it,NX,NY,NZ,NT)

                                b[1,ix,iy,iz,it,ialpha] = η*(conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,1,ix,iy,iz,it]) + 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,2,ix,iy,iz,it]) + 
                                                            conj(a[3,ix1,iy1,iz1,it1,ialpha])* conj(u[-μ][1,3,ix,iy,iz,it]))

                                b[2,ix,iy,iz,it,ialpha] = η*(conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,1,ix,iy,iz,it])+ 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,2,ix,iy,iz,it]) + 
                                                            conj(a[3,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,3,ix,iy,iz,it]))

                                b[3,ix,iy,iz,it,ialpha] = η*(conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][3,1,ix,iy,iz,it]) + 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][3,2,ix,iy,iz,it]) + 
                                                            conj(a[3,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][3,3,ix,iy,iz,it]))

                        
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

    function fermion_shiftB!(b::StaggeredFermion,evensite,u::Array{SU3GaugeFields,1},μ,a::StaggeredFermion) 
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
            for ialpha =1:1
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
                                η = staggered_phase(-μ,ix,iy,iz,it,NX,NY,NZ,NT)

                                b[1,ix,iy,iz,it,ialpha] = η*(conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,1,ix,iy,iz,it]) + 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,2,ix,iy,iz,it]) + 
                                                            conj(a[3,ix1,iy1,iz1,it1,ialpha])* conj(u[-μ][1,3,ix,iy,iz,it]))

                                b[2,ix,iy,iz,it,ialpha] = η*(conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,1,ix,iy,iz,it])+ 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,2,ix,iy,iz,it]) + 
                                                            conj(a[3,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,3,ix,iy,iz,it]))

                                b[3,ix,iy,iz,it,ialpha] = η*(conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][3,1,ix,iy,iz,it]) + 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][3,2,ix,iy,iz,it]) + 
                                                            conj(a[3,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][3,3,ix,iy,iz,it]))

                        
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

    
    function fermion_shiftB!(b::StaggeredFermion,u::Array{SU2GaugeFields,1},μ,a::StaggeredFermion) 
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
            for ialpha =1:1
                for it=1:NT
                    it1 = it + ifelse(-μ == 4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz +ifelse(-μ == 3,1,0)  #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(-μ == 2,1,0) #idel[2]
                            for ix=1:NX
                                ix1 = ix + ifelse(-μ == 1,1,0)  #idel[1]
                                η = staggered_phase(-μ,ix,iy,iz,it,NX,NY,NZ,NT)
                                

                                b[1,ix,iy,iz,it,ialpha] = η*(conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,1,ix,iy,iz,it]) + 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,2,ix,iy,iz,it]))

                                b[2,ix,iy,iz,it,ialpha] = η*(conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,1,ix,iy,iz,it])+ 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,2,ix,iy,iz,it]) )

                                
                            end
                        end
                    end
                end
            end
        end

    end

    function fermion_shiftB!(b::StaggeredFermion,u::Array{SU2GaugeFields,1},μ,a::StaggeredFermion,indices) 
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

            ix1 = indices[1]
            iy1 = indices[2]
            iz1 = indices[3]
            it1 = indices[4]
            ialpha = indices[5]

            ix = ix1 - ifelse(-μ ==1,1,0)
            iy = iy1 - ifelse(-μ ==2,1,0)
            iz = iz1 - ifelse(-μ ==3,1,0)
            it = it1 - ifelse(-μ ==4,1,0)

            ix,iy,iz,it,sign =  apply_periodicity(ix,iy,iz,it,NX,NY,NZ,NT,b.BoundaryCondition)

            η = staggered_phase(-μ,ix,iy,iz,it,NX,NY,NZ,NT)

            b[1,ix,iy,iz,it,ialpha] = sign*η*(conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,1,ix,iy,iz,it]) + 
                                        conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,2,ix,iy,iz,it]))

            b[2,ix,iy,iz,it,ialpha] = sign*η*(conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,1,ix,iy,iz,it])+ 
                                        conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,2,ix,iy,iz,it]) )

        end

        return (ix,iy,iz,it,ialpha)

    end

    function fermion_shiftB!(b::StaggeredFermion,evensite,u::Array{SU2GaugeFields,1},μ,a::StaggeredFermion) 
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
            for ialpha =1:1
                for it=1:NT
                    it1 = it + ifelse(-μ == 4,1,0) #idel[4]
                    for iz=1:NZ
                        iz1 = iz +ifelse(-μ == 3,1,0)  #idel[3]
                        for iy=1:NY
                            iy1 = iy + ifelse(-μ == 2,1,0) #idel[2]
                            xran =1+(1+ibush+iy+iz+it)%2:2:NX
                            for ix in xran
                                ix1 = ix + ifelse(-μ == 1,1,0)  #idel[1]
                                η = staggered_phase(-μ,ix,iy,iz,it,NX,NY,NZ,NT)
                                

                                b[1,ix,iy,iz,it,ialpha] = η*(conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,1,ix,iy,iz,it]) + 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][1,2,ix,iy,iz,it]))

                                b[2,ix,iy,iz,it,ialpha] = η*(conj(a[1,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,1,ix,iy,iz,it])+ 
                                                            conj(a[2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][2,2,ix,iy,iz,it]) )

                                
                            end
                        end
                    end
                end
            end
        end

    end

    """
c-------------------------------------------------c
c     Random number function for Gaussian  Noise
    with σ^2 = 1/2
c-------------------------------------------------c
    """
    function gauss_distribution_fermi!(x::FermionFields)
        NC = x.NC
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        n6 = size(x.f)[6]
        σ = sqrt(1/2)

        for ialpha = 1:n6
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            for ic=1:NC
                                
                                x[ic,ix,iy,iz,it,ialpha] = σ*randn()+im*σ*randn()
                            end
                        end
                    end
                end
            end
        end

        set_wing_fermi!(x)

        return
    end

    """
c-------------------------------------------------c
c     Random number function for Gaussian  Noise
    with σ^2 = 1/2
c-------------------------------------------------c
    """
    function gauss_distribution_fermi!(x::FermionFields,randomfunc)
        NC = x.NC
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        n6 = size(x.f)[6]
        σ = sqrt(1/2)

        for mu = 1:n6
            for ic=1:NC
                for it=1:NT
                    for iz=1:NZ
                        for iy=1:NY
                            for ix=1:NX
                                v1 = sqrt(-log(randomfunc()+1e-10))
                                v2 = 2pi*randomfunc()

                                xr = v1*cos(v2)
                                xi = v1 * sin(v2)

                                x[ic,ix,iy,iz,it,mu] = xr + im*xi
                            end
                        end
                    end
                end
            end
        end

        set_wing_fermi!(x)

        return
    end

    """
c-------------------------------------------------c
c     Random number function Z4  Noise
c     https://arxiv.org/pdf/1611.01193.pdf
c-------------------------------------------------c
    """
    function Z4_distribution_fermi!(x::FermionFields)
        NC = x.NC
        NX = x.NX
        NY = x.NY
        NZ = x.NZ
        NT = x.NT
        n6 = size(x.f)[6]
        θ = 0.0
        N::Int32 = 4
        Ninv = Float64(1/N)
        for ialpha = 1:n6
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

        set_wing_fermi!(x)

        return
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

    #=
    function cloverterm!(vec,fparam,x)
        NT = x.NT
        NZ = x.NZ
        NY = x.NY
        NX = x.NX
        NC = x.NC

        
        i  =0
        for it=1:NT
            for iz=1:NZ
                for iy=1:NZ
                    for ix=1:NZ
                        i += 1
                        for k1=1:NC
                            for k2=1:NC

                                c1 = x[k2,ix,iy,iz,it,1]
                                c2 = x[k2,ix,iy,iz,it,2]
                                c3 = x[k2,ix,iy,iz,it,3]
                                c4 = x[k2,ix,iy,iz,it,4]

                                vec[k1,ix,iy,iz,it,1] += fparam.CloverFμν[k1,k2,i,1]*(-   c1) + 
                                                            + fparam.CloverFμν[k1,k2,i,2]*(-im*c2) + 
                                                            + fparam.CloverFμν[k1,k2,i,3]*(-   c2) + 
                                                            + fparam.CloverFμν[k1,k2,i,4]*(-   c2) + 
                                                            + fparam.CloverFμν[k1,k2,i,5]*( im*c2) + 
                                                            + fparam.CloverFμν[k1,k2,i,6]*(-   c1)
                
                                

                                vec[k1,ix,iy,iz,it,2] += fparam.CloverFμν[k1,k2,i,1]*(   c2) + 
                                                            + fparam.CloverFμν[k1,k2,i,2]*(im*c1) + 
                                                            + fparam.CloverFμν[k1,k2,i,3]*(-   c1) + 
                                                            + fparam.CloverFμν[k1,k2,i,4]*(-   c1) + 
                                                            + fparam.CloverFμν[k1,k2,i,5]*(-im*c1) + 
                                                            + fparam.CloverFμν[k1,k2,i,6]*(   c2)

                                vec[k1,ix,iy,iz,it,3] += fparam.CloverFμν[k1,k2,i,1]*(   -c3) + 
                                                            + fparam.CloverFμν[k1,k2,i,2]*(-im*c4) + 
                                                            + fparam.CloverFμν[k1,k2,i,3]*(   c4) + 
                                                            + fparam.CloverFμν[k1,k2,i,4]*(-   c4) + 
                                                            + fparam.CloverFμν[k1,k2,i,5]*(-im*c4) + 
                                                            + fparam.CloverFμν[k1,k2,i,6]*(   c3)

                                vec[k1,ix,iy,iz,it,4] += fparam.CloverFμν[k1,k2,i,1]*(   c4) + 
                                                            + fparam.CloverFμν[k1,k2,i,2]*(im*c3) + 
                                                            + fparam.CloverFμν[k1,k2,i,3]*(   c3) + 
                                                            + fparam.CloverFμν[k1,k2,i,4]*(-   c3) + 
                                                            + fparam.CloverFμν[k1,k2,i,5]*(im*c3) + 
                                                            + fparam.CloverFμν[k1,k2,i,6]*( -  c4)


                            end
                        end
                    end
                end

            end
        end

        #println("vec = ",vec*vec)

    end
    =#

    function cloverterm!(vec,CloverFμν,x)
        NT = x.NT
        NZ = x.NZ
        NY = x.NY
        NX = x.NX
        NC = x.NC

        
        i  =0
        for it=1:NT
            for iz=1:NZ
                for iy=1:NZ
                    for ix=1:NZ
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


    

    function set_wing_fermi!(a::FermionFields)
        NT = a.NT
        NZ = a.NZ
        NY = a.NY
        NX = a.NX
        NC = a.NC

        #!  X-direction
        for ialpha=1:4
            for it=1:NT
                for iz = 1:NZ
                    for iy=1:NY
                        for k=1:NC
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
                        for k=1:NC
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
                        for k=1:NC
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
                        for k=1:NC
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
                        for k=1:NC
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
                        for k=1:NC
                            a[k,ix,iy,iz,0,ialpha] = a.BoundaryCondition[4]*a[k,ix,iy,iz,NT,ialpha]
                            a[k,ix,iy,iz,NT+1,ialpha] = a.BoundaryCondition[4]*a[k,ix,iy,iz,1,ialpha]
                        end
                    end
                end
            end

        end

        

        

    end


    function set_wing_fermi!(a::StaggeredFermion)
        NT = a.NT
        NZ = a.NZ
        NY = a.NY
        NX = a.NX
        NC = a.NC

        #!  X-direction
        for ialpha=1:1
            for it=1:NT
                for iz = 1:NZ
                    for iy=1:NY
                        for k=1:NC
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
                        for k=1:NC
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
                        for k=1:NC
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
                        for k=1:NC
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
                        for k=1:NC
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
                        for k=1:NC
                            a[k,ix,iy,iz,0,ialpha] = a.BoundaryCondition[4]*a[k,ix,iy,iz,NT,ialpha]
                            a[k,ix,iy,iz,NT+1,ialpha] = a.BoundaryCondition[4]*a[k,ix,iy,iz,1,ialpha]
                        end
                    end
                end
            end

        end

        

        

    end









end