"""
Generic Fermion module
"""
module AbstractFermion
    #export FermionFields,Wx!,Wdagx!,clear!,substitute_fermion!,Dx!

    #export FermionFields,
    #        Wx!,Wdagx!,clear!,substitute_fermion!,Dx!,fermion_shift!,
    #        fermion_shiftB!,add!,set_wing_fermi!,WdagWx!,apply_periodicity

    import ..Actions:FermiActionParam,FermiActionParam_Wilson,
                FermiActionParam_WilsonClover,FermiActionParam_Staggered
    import ..Gaugefields:GaugeFields,GaugeFields_1d,SU3GaugeFields,SU2GaugeFields,SU3GaugeFields_1d,SU2GaugeFields_1d,
                staggered_phase,SUn,SU2,SU3,SUNGaugeFields,SUNGaugeFields_1d,SU

    abstract type FermionFields end



    
    function clear!(a::FermionFields)
        n1,n2,n3,n4,n5,n6 = size(a.f)
        for i6=1:n6
            for i5=1:n5
                for i4=1:n4
                    for i3=1:n3
                        for i2=1:n2
                            @simd for i1=1:n1
                                a.f[i1,i2,i3,i4,i5,i6]= 0
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
                            @simd for i1=1:n1
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

    function add!(c::FermionFields,alpha::Number,a::FermionFields) #c = c + alpha*a 
        n1,n2,n3,n4,n5,n6 = size(a.f)

        for i6=1:n6
            for i5=1:n5
                for i4=1:n4
                    for i3=1:n3
                        for i2=1:n2
                            @simd for i1=1:n1
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

    function add!(coeff::Number,c::FermionFields,alpha::Number,a::FermionFields) #c = coeff*c + alpha*a 
        n1,n2,n3,n4,n5,n6 = size(a.f)

        for i6=1:n6
            for i5=1:n5
                for i4=1:n4
                    for i3=1:n3
                        for i2=1:n2
                            @simd for i1=1:n1
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
        error("This FermionFields type is not supported yet in function Wx!. The type is $(typeof(x))")
    end

    function Dx!(xout::T,U::Array{G,1},
        x::T,temps::Array{T,1}) where  {T <: FermionFields,G <: GaugeFields}
        error("This FermionFields type is not supported yet in function Dx!. The type is $(typeof(x))")

        return
    end

    function fermion_shiftB!(b::T,evensite,u::Array{GaugeFields{SU{NC}},1},μ,a::T)  where {NC,T <: FermionFields}
        error("This FermionFields type is not supported yet in function fermion_shiftB!. The type is $(typeof(x))")
    end


    function Wdagx!(xout::FermionFields,U::Array{GaugeFields,1},
        x::FermionFields,temps::Array{FermionFields,1},fparam::T)  where T <: FermiActionParam
        error("This FermionFields type is not supported yet in function Wdagx!. The type is $(typeof(x))")
    end

    function WdagWx!(xout::T,U::Array{G,1},x::T,temps::Array{T,1},fparam) where {T <: FermionFields,G <: GaugeFields}
        temp = temps[5]
        Wx!(temp,U,x,temps,fparam) 
    
        Wdagx!(xout,U,temp,temps,fparam) 
    
        return
    end



    function WdagWx!(xout::T,U::Array{G,1},x::T,temps::Array{T,1},fparam,indices) where {T <: FermionFields,G <: GaugeFields}
        WdagWx!(xout,U,x,temps,fparam)
        return
    end



    function substitute_fermion!(H,j,x::FermionFields)
        error("substitute_fermion! is not implemented for $(typeof(x))")
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

    
    function fermion_shift!(b::F,u::Array{GaugeFields{SU{NC}},1},μ,a::F) where {NC,F <: FermionFields}
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
                                    @simd for k2=1:NC
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
                                    @simd for k2=1:NC
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

    
    function fermion_shift!(b::F,evensite::Bool,u::Array{GaugeFields{SU{NC}},1},μ,a::F) where {NC,F <: FermionFields}
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
                                    @simd for k2=1:NC
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
                                    @simd for k2=1:NC
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



end