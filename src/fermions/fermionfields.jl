module Fermionfields
    using Random
    using LinearAlgebra
    using SparseArrays


    import ..Actions:FermiActionParam,FermiActionParam_Wilson,
                FermiActionParam_WilsonClover,FermiActionParam_Staggered
    import ..Gaugefields:GaugeFields,GaugeFields_1d,SU3GaugeFields,SU2GaugeFields,SU3GaugeFields_1d,SU2GaugeFields_1d,
                            staggered_phase,SUn,SU2,SU3,SUNGaugeFields,SUNGaugeFields_1d,SU
    #import ..Parallel:get_looprange

    #include("cg.jl")

    import ..AbstractFermion:FermionFields,
        Wx!,Wdagx!,clear!,substitute_fermion!,Dx!,Dxplus!,
        fermion_shift!,fermion_shiftB!,add!,set_wing_fermi!,WdagWx!,apply_periodicity,
        gauss_distribution_fermi!,gauss_distribution_fermi_Z2!,Z4_distribution_fermi!,Ddagx!

    import ..WilsonFermion_module:WilsonFermion,mul_γ5x!

    import ..StaggeredFermion_module:StaggeredFermion

    import ..DomainwallFermion_module:DomainwallFermion,D5DWx!,D5DWdagx!


#    abstract type FermionFields end


    struct RGammamatrix{ν,pm} 
        #indices::NTuple{}
        #values::
    end


    function FermionFields(NC,NX,NY,NZ,NT,fparam::FermiActionParam,BoundaryCondition)
        if findfirst("Wilson",fparam.Dirac_operator) != nothing #If string fparam.Dirac_operator includes "Wilson"
            #fparam.Dirac_operator == "Wilson"
            return WilsonFermion(NC,NX,NY,NZ,NT,fparam,BoundaryCondition)
        elseif fparam.Dirac_operator == "Staggered"
            return StaggeredFermion(NC,NX,NY,NZ,NT,fparam,BoundaryCondition)
        elseif fparam.Dirac_operator == "Domainwall"
            return DomainwallFermion(NC,NX,NY,NZ,NT,fparam,BoundaryCondition)

        else
            error(fparam.Dirac_operator,"is not supported!")
        end
    end



    #include("boundary_fermi.jl")




    function Base.setindex!(x::FermionFields,v,i1,i2,i3,i4,i5,i6) 
        x.f[i1,i2 + 1,i3 + 1,i4 + 1,i5 + 1,i6] = v
    end

    function Base.getindex(x::FermionFields,i1,i2,i3,i4,i5,i6)
        return x.f[i1,i2 .+ 1,i3 .+ 1,i4 .+ 1,i5 .+ 1,i6]
    end




    function LinearAlgebra.dot(a::FermionFields,b::FermionFields)
        return a*b
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




    function LinearAlgebra.axpby!(a::Number,X::FermionFields,Y::FermionFields)
        add!(Y,a,X)
        return
    end



    function LinearAlgebra.axpby!(a::Number,X::FermionFields,b,Y::FermionFields)
        add!(b,Y,a,X)
        return
    end




    function Wdagx_old!(xout::WilsonFermion,U::Array{G,1},
        x::WilsonFermion,temps::Array{T,1}) where {T <: FermionFields,G <: GaugeFields}
        temp = temps[4]
        temp1 = temps[1]
        temp2 = temps[2]

        clear!(temp)
        x5 = temps[3]

        #D^dag = \gamma_5 D \gamma_5


        mul_γ5x!(x5,x)
        set_wing_fermi!(x5)


        for ν=1:4
            fermion_shift!(temp1,U,ν,x5) #nu shift

            #... Dirac multiplication
            #mul!(temp1,view(x.rminusγ,:,:,ν),temp1) #(1 - \gamma)
            mul!(temp1,view(x.rplusγ,:,:,ν),temp1) #(1 + \gamma)
            
            #
            fermion_shift!(temp2,U,-ν,x5) #-nu shift
            
            #mul!(temp2,view(x.rplusγ,:,:,ν),temp2) #(1 + \gamma)
            mul!(temp2,view(x.rminusγ,:,:,ν),temp2) #(1 + \gamma)

            add!(temp,x.hopp[ν],temp1,x.hopm[ν],temp2)
        end
        clear!(temp1)
        add!(temp1,1,x5,-1,temp)

        mul_γ5x!(xout,temp1)
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
                                @simd for k2=1:n6
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
                            
                            @simd for k1=1:n6
                                e[k1] = x[ic,ix,iy,iz,it,k1]
                            end

                            for k1=1:n6
                                f[k1] = 0
                                @simd for k2=1:n6
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



    function LinearAlgebra.mul!(a::FermionFields,α::Number,b::FermionFields)
        n1,n2,n3,n4,n5,n6 = size(a.f)

        for i6=1:n6
            for i5=1:n5
                for i4=1:n4
                    for i3=1:n3
                        for i2=1:n2
                            @simd for i1=1:n1
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
                                    @simd for k=1:4
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
                                    @simd for k=1:4
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
                                    @simd for k=1:1
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
                                    @simd for k=1:1
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



    












end