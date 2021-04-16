module Gaugefields
    using LinearAlgebra
    using Random
    using JLD
    import ..Actions:GaugeActionParam,GaugeActionParam_standard,GaugeActionParam_autogenerator,SmearingParam_single,SmearingParam_multi,SmearingParam
    import ..Wilsonloops:Wilson_loop_set,calc_coordinate,make_plaq_staple_prime,calc_shift,make_plaq,make_plaq_staple,
                            Tensor_wilson_lines_set,Tensor_wilson_lines,Tensor_derivative_set,
                            get_leftstartposition,get_rightstartposition,Wilson_loop
    import ..SUN_generator:Generator
    

    abstract type SUn end

    abstract type SU{N} <: SUn
    end

    const U1 = SU{1}
    const SU2 = SU{2}
    const SU3 = SU{3}


    struct GaugeFields{T <: SUn} 
        g::Array{ComplexF64,6}
        NX::Int64
        NY::Int64
        NZ::Int64
        NT::Int64
        NC::Int64
        NDW::Int64
        NV::Int64


        function GaugeFields(NC,NDW,NX,NY,NZ,NT)
            sutype = SU{NC}
                
            NV = NX*NY*NZ*NT
            g = zeros(ComplexF64,NC,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW)
            return new{sutype}(g,NX,NY,NZ,NT,NC,NDW,NV)
        end

    end



    const U1GaugeFields  = GaugeFields{U1}
    const SU2GaugeFields  = GaugeFields{SU2}
    const SU3GaugeFields  = GaugeFields{SU3}
    const SUNGaugeFields{N}  = GaugeFields{SU{N}}



    struct Adjoint_GaugeFields{T <: SUn} 
        parent::GaugeFields{T}
    end

    const Adjoint_U1GaugeFields  = Adjoint_GaugeFields{U1}
    const Adjoint_SU2GaugeFields  = Adjoint_GaugeFields{SU2}
    const Adjoint_SU3GaugeFields  = Adjoint_GaugeFields{SU3}

    
    function Base.adjoint(x::GaugeFields{T}) where T <: SUn
        Adjoint_GaugeFields{T}(x)
    end

    struct GaugeFields_1d{T <: SUn} 
        g::Array{ComplexF64,3}
        NC::Int64
        NV::Int64
    
        function GaugeFields_1d(NC,NX,NY,NZ,NT)
            sutype = SU{NC}

            NV = NX*NY*NZ*NT
            g = zeros(ComplexF64,NC,NC,NV)
            return new{sutype}(g,NC,NV)
        end
    
        function GaugeFields_1d(NC,NV) 
            sutype = SU{NC}
            g = zeros(ComplexF64,NC,NC,NV)
            return new{sutype}(g,NC,NV)
        end
    end

    const U1GaugeFields_1d  = GaugeFields_1d{U1}
    const SU2GaugeFields_1d  = GaugeFields_1d{SU2}
    const SU3GaugeFields_1d  = GaugeFields_1d{SU3}
    const SUNGaugeFields_1d{N}  = GaugeFields_1d{SU{N}}

    struct Adjoint_GaugeFields_1d{T <: SUn} 
        parent::GaugeFields_1d{T}
    end

    const Adjoint_1U1GaugeFields_1d  = Adjoint_GaugeFields_1d{U1}
    const Adjoint_SU2GaugeFields_1d  = Adjoint_GaugeFields_1d{SU2}
    const Adjoint_SU3GaugeFields_1d  = Adjoint_GaugeFields_1d{SU3}

    function Base.adjoint(x::GaugeFields_1d{T}) where T <: SUn
        Adjoint_GaugeFields_1d{T}(x)
    end

    function Base.setindex!(x::GaugeFields,v,i1,i2,i3,i4,i5,i6) 
        x.g[i1,i2,i3 .+ x.NDW,i4 + x.NDW,i5 + x.NDW,i6 + x.NDW] = v
    end

    function Base.getindex(x::GaugeFields,i1,i2,i3,i4,i5,i6)
        return x.g[i1,i2,i3 .+ x.NDW,i4 .+ x.NDW,i5 .+ x.NDW,i6 .+ x.NDW]
    end



    function Base.setindex!(x::GaugeFields_1d,v,i1,i2,i3) 
        x.g[i1,i2,i3] = v
    end

    function Base.getindex(x::GaugeFields_1d,i1,i2,i3)
        return x.g[i1,i2,i3]
    end

    struct Loops{T_1d}
        loopset::Wilson_loop_set
        temps::Array{T_1d,1}
        function Loops(U::Array{GaugeFields{S},1},loopset::Wilson_loop_set) where {S <: SUn}
            NC = U[1].NC
            NX = U[1].NX
            NY = U[1].NY
            NZ = U[1].NZ
            NT = U[1].NT
            temp1 = GaugeFields_1d(NC,NX,NY,NZ,NT)
            temp2 = GaugeFields_1d(NC,NX,NY,NZ,NT)
            temp3 = GaugeFields_1d(NC,NX,NY,NZ,NT)
            T_1d = typeof(temp1)
            return new{T_1d}(loopset,[temp1,temp2,temp3])
        end

        function Loops(U::Array{GaugeFields{S},1},loopset::Wilson_loop_set,temps) where {S <: SUn}
            T_1d = typeof(temps[1])
            return new{T_1d}(loopset,temps)
        end
    end

    function evaluate_loops(loops::Loops,U::Array{GaugeFields{S},1})  where {S <: SUn}
        NC = U[1].NC
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        xout = GaugeFields_1d(NC,NX,NY,NZ,NT)
        evaluate_loops!(xout,loops,U)
        return xout
    end

    function evaluate_loops!(xout::T_1d,loops::Loops,U::Array{GaugeFields{S},1})  where {S <: SUn,T_1d <: GaugeFields_1d}
        evaluate_wilson_loops!(xout,loops.loopset,U,loops.temps)
    end

    

    function evaluate_loops!(V,loops::Loops,U::Array{T,1},ix,iy,iz,it) where T <: GaugeFields
        evaluate_wilson_loops!(V,loops.loopset,U,ix,iy,iz,it)
    end

    function evaluate_loops(loops::Loops,U::Array{T,1},ix,iy,iz,it) where T <: GaugeFields
        NC = U[1].NC
        V = zeros(ComplexF64,NC,NC)
        evaluate_loops!(V,loops,U,ix,iy,iz,it)
        return V
    end

    function evaluate_wilson_loops!(xout::T_1d,w::Wilson_loop_set,U::Array{GaugeFields{S},1}) where {S <: SUn,T_1d <: GaugeFields_1d}
        NC = U[1].NC
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT

        temp1 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp2 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp3 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        evaluate_wilson_loops!(xout,w,U,[temp1,temp2,temp3])

    end




    function evaluate_wilson_loops!(xout::T_1d,w::Wilson_loop_set,U::Array{GaugeFields{S},1},temps::Array{T_1d,1}) where {S <: SUn,T_1d <: GaugeFields_1d}

        num = length(w)
        clear!(xout)
        temp1 = temps[1]
        temp2 = temps[2]
        temp3 = temps[3]

        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        NC = U[1].NC



        for i=1:num
            wi = w[i]
            numloops = length(wi)    
    
            shifts = calc_shift(wi)
            #println("shift ",shifts)
            
            loopk = wi[1]
            k = 1
            #println("k = $k shift: ",shifts[k])
            gauge_shift_all!(temp1,shifts[1],U[loopk[1]])

            
            loopk1_2 = loopk[2]
            for k=2:numloops
                loopk = wi[k]
                #println("k = $k shift: ",shifts[k])
                #println("gauge_shift!(temp2,$(shifts[k]),$(loopk[1]) )")
                clear!(temp2)
                gauge_shift_all!(temp2,shifts[k],U[loopk[1]])

                multiply_12!(temp3,temp1,temp2,k,loopk,loopk1_2)
                temp1,temp3 = temp3,temp1

                
            end
            add!(xout,temp1)
            
        end
    end

    function evaluate_wilson_loops!(xout::T_1d,w::Wilson_loop_set,U::Array{GaugeFields{S},1},temps::Array{T_1d,1}) where {S <: SUn,T_1d <: GaugeFields_1d}

        num = length(w)
        clear!(xout)
        temp1 = temps[1]
        temp2 = temps[2]
        temp3 = temps[3]

        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        NC = U[1].NC



        for i=1:num
            wi = w[i]
            numloops = length(wi)    
    
            shifts = calc_shift(wi)
            #println("shift ",shifts)
            
            loopk = wi[1]
            k = 1
            #println("k = $k shift: ",shifts[k])
            gauge_shift_all!(temp1,shifts[1],U[loopk[1]])

            
            loopk1_2 = loopk[2]
            for k=2:numloops
                loopk = wi[k]
                #println("k = $k shift: ",shifts[k])
                #println("gauge_shift!(temp2,$(shifts[k]),$(loopk[1]) )")
                clear!(temp2)
                gauge_shift_all!(temp2,shifts[k],U[loopk[1]])

                multiply_12!(temp3,temp1,temp2,k,loopk,loopk1_2)
                temp1,temp3 = temp3,temp1

                
            end
            add!(xout,temp1)
            
        end
    end

    function shift_xyzt(shift,ix,iy,iz,it)
        ix1 = ix + shift[1]
        iy1 = iy + shift[2]
        iz1 = iz + shift[3]
        it1 = it + shift[4]
        return ix1,iy1,iz1,it1
    end

    function multiply_12!(temp3,temp1,temp2,k,loopk,loopk1_2)
        if loopk[2] == 1
            if k==2
                if loopk1_2 == 1
                    mul!(temp3,temp1,temp2)
                else
                    mul!(temp3,temp1',temp2)
                end
            else
                mul!(temp3,temp1,temp2)
            end
        elseif loopk[2] == -1
            if k==2
                if loopk1_2 == 1
                    mul!(temp3,temp1,temp2')
                else
                    mul!(temp3,temp1',temp2')
                end
            else
                mul!(temp3,temp1,temp2')
            end
        else
            error("Second element should be 1 or -1 but now $(loopk)")
        end
        return
    end

    function evaluate_wilson_loops!(V,w::Wilson_loop_set,U::Array{T,1},ix,iy,iz,it) where T <: GaugeFields
        num = length(w)

        V .= 0 #zeros(ComplexF64,NC,NC)
        temp1 = zero(V)
        temp2 = zero(V)
        temp3 = zero(V)
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        NDW = U[1].NDW


        for i=1:num
            wi = w[i]
            numloops = length(wi)    
            coordinates = calc_coordinate(wi)
            #println("positions ",coordinates)        
            shifts = calc_shift(wi)
            #println("shift ",shifts)
            
            loopk = wi[1]
            ix1,iy1,iz1,it1 = shift_xyzt(shifts[1],ix,iy,iz,it)
            ix1,iy1,iz1,it1 =periodiccheck(ix1,iy1,iz1,it1,NX,NY,NZ,NT,NDW)

            
            temp1[:,:] = U[loopk[1]][:,:,ix1,iy1,iz1,it1]


            loopk1_2 = loopk[2]

            #gauge_shift_all!(temp1,shifts[1],U[loopk[1]])
            for k=2:numloops
                loopk = wi[k]
                
                #gauge_shift_all!(temp2,shifts[k],U[loopk[1]])
                ix1,iy1,iz1,it1 = shift_xyzt(shifts[k],ix,iy,iz,it)
                ix1,iy1,iz1,it1 =periodiccheck(ix1,iy1,iz1,it1,NX,NY,NZ,NT,NDW)
                
                temp2[:,:] = U[loopk[1]][:,:,ix1,iy1,iz1,it1]

                multiply_12!(temp3,temp1,temp2,k,loopk,loopk1_2)

                temp1,temp3 = temp3,temp1
            end
            @. V += temp1
            #add!(xout,temp1)
            
        end
    end

    function calc_coordinate_shift(loopk,coordinates)
        if loopk[2] == 1
            shifts = coordinates
        elseif loopk[2] == -1
            if loopk[1] == 1
                shifts = (coordinates[1]+loopk[2],coordinates[2],coordinates[3],coordinates[4])
            elseif loopk[1] == 2
                shifts = (coordinates[1],coordinates[2]+loopk[2],coordinates[3],coordinates[4])
            elseif loopk[1] == 3
                shifts = (coordinates[1],coordinates[2],coordinates[3]+loopk[2],coordinates[4])
            elseif loopk[1] == 4
                shifts = (coordinates[1],coordinates[2],coordinates[3],coordinates[4]+loopk[2])
            end
        end

        if loopk[1] == 1
            coordinates = (coordinates[1]+loopk[2],coordinates[2],coordinates[3],coordinates[4])
        elseif loopk[1] == 2
            coordinates = (coordinates[1],coordinates[2]+loopk[2],coordinates[3],coordinates[4])
        elseif loopk[1] == 3
            coordinates = (coordinates[1],coordinates[2],coordinates[3]+loopk[2],coordinates[4])
        elseif loopk[1] == 4
            coordinates = (coordinates[1],coordinates[2],coordinates[3],coordinates[4]+loopk[2])
        end

        return coordinates,shifts
    end


    function evaluate_wilson_loops!(V,w::Wilson_loop_set,U::Array{GaugeFields{SU{NC}},1},ix,iy,iz,it,temps) where NC
        num = length(w)

        V .= 0 #zeros(ComplexF64,NC,NC)
        temp1 = temps[1]
        temp2 = temps[2]
        temp3 = temps[3]
        temp1 .= 0
        temp2 .= 0
        temp3 .= 0

        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        NDW = U[1].NDW
        #NC = U[1].NDW



        for i=1:num
            wi = w[i]
            numloops = length(wi)    
            #println("d")
            #@time coordinates = calc_coordinate(wi)
            #println("positions ",coordinates)     




            #shifts = calc_shift(wi)
            #println("shift ",shifts)
            
            loopk = wi[1]


            
            coordinates = wi.origin

            coordinates,shifts2 = calc_coordinate_shift(loopk,coordinates)
            #println(shifts)
            #println(shifts[1],"\t",shifts2)

            
            
            
            


            #println("d3")

            ix1,iy1,iz1,it1 = shift_xyzt(shifts2,ix,iy,iz,it)

            #ix1,iy1,iz1,it1 = shift_xyzt(shifts[1],ix,iy,iz,it)
            #ix1,iy1,iz1,it1 = shift_xyzt(shifts,ix,iy,iz,it)
            #println("d4")
            ix1,iy1,iz1,it1 =periodiccheck(ix1,iy1,iz1,it1,NX,NY,NZ,NT,NDW)

            #println("d5")
            for k=1:NC
                for j=1:NC
                    temp1[j,k] = U[loopk[1]][j,k,ix1,iy1,iz1,it1]
                end
            end
            #println("t ",temp1)


            loopk1_2 = loopk[2]

            #gauge_shift_all!(temp1,shifts[1],U[loopk[1]])
            for k=2:numloops
                loopk = wi[k]
                
                #gauge_shift_all!(temp2,shifts[k],U[loopk[1]])
                #println("d6")

                coordinates,shifts2 = calc_coordinate_shift(loopk,coordinates)

                #println(shifts[k],"\t",shifts2)

                ix1,iy1,iz1,it1 = shift_xyzt(shifts2,ix,iy,iz,it)
                
                #ix1,iy1,iz1,it1 = shift_xyzt(shifts[k],ix,iy,iz,it)
                #ix1,iy1,iz1,it1 = shift_xyzt(shifts,ix,iy,iz,it)
                #println("d7")
                ix1,iy1,iz1,it1 =periodiccheck(ix1,iy1,iz1,it1,NX,NY,NZ,NT,NDW)
                
                #println("d8")
                for k=1:NC
                    for j=1:NC
                        temp2[j,k] = U[loopk[1]][j,k,ix1,iy1,iz1,it1]
                    end
                end

                #@time temp2[:,:] = U[loopk[1]][:,:,ix1,iy1,iz1,it1]

                #println("d9")
                multiply_12!(temp3,temp1,temp2,k,loopk,loopk1_2)

                #println("d10")
                temp1,temp3 = temp3,temp1
            end
            @. V += temp1
            #add!(xout,temp1)
            
        end
        #exit()
    end


    function set_wing!(u::Array{T,1}) where T <: GaugeFields
        for μ=1:4
            set_wing!(u[μ])
        end
    end



    function set_wing!(u::GaugeFields{SU{NC}},ix,iy,iz,it) where NC#where T <: SUn
        NT = u.NT
        NY = u.NY
        NZ = u.NZ
        NX = u.NX
        NDW = u.NDW


        function update!(u,ixp,iyp,izp,itp,ix,iy,iz,it)
            for k2=1:NC
                for k1=1:NC
                    u[k1,k2,ixp,iyp,izp,itp] = u[k1,k2,ix,iy,iz,it]
                end
            end
            return
        end

        isxwing = 0
        isywing = 0
        iszwing = 0
        istwing = 0

        if ix ≤ NDW
            isxwing = 1
        elseif ix ≥ NX + 1- NDW
            isxwing = -1
        end

        if iy ≤ NDW
            isywing = 1
        elseif iy ≥ NY + 1- NDW
            isywing = -1
        end

        if iz ≤ NDW
            iszwing = 1
        elseif iz ≥ NZ + 1- NDW
            iszwing = -1
        end

        if it ≤ NDW
            istwing = 1
        elseif it ≥ NT + 1- NDW
            istwing = -1
        end

        if isxwing == 0 && isywing == 0 && iszwing == 0 && istwing == 0
            return
        end


        ixp = ix
        iyp = iy
        izp = iz
        itp = it
    
        if isxwing == 1
            ixp = NX+ix
        elseif isxwing == -1
            ixp = ix-NX
        end

        if isywing == 1
            iyp = NY+iy
        elseif isywing == -1
            iyp = iy-NY
        end

        if iszwing == 1
            izp = NZ+iz
        elseif iszwing == -1
            izp = iz-NZ
        end

        if istwing == 1
            itp = NT+it
        elseif istwing == -1
            itp = it-NT
        end

        if isxwing != 0
            update!(u,ixp,iy,iz,it,ix,iy,iz,it)

            if isywing != 0
                update!(u,ixp,iyp,iz,it,ix,iy,iz,it)
                if iszwing != 0
                    update!(u,ixp,iyp,izp,it,ix,iy,iz,it)

                    if istwing != 0
                        update!(u,ixp,iyp,izp,itp,ix,iy,iz,it)
                    end
                end
                if istwing != 0
                    update!(u,ixp,iyp,iz,itp,ix,iy,iz,it)
                end
            end
            
            if iszwing != 0
                update!(u,ixp,iy,izp,it,ix,iy,iz,it)
                if istwing != 0
                    update!(u,ixp,iy,izp,itp,ix,iy,iz,it)
                end
            end
            
            if istwing != 0
                update!(u,ixp,iy,iz,itp,ix,iy,iz,it)
            end
            
        end

        if isywing != 0
            update!(u,ix,iyp,iz,it,ix,iy,iz,it)

            if iszwing != 0
                update!(u,ix,iyp,izp,it,ix,iy,iz,it)
                if istwing != 0
                    update!(u,ix,iyp,izp,itp,ix,iy,iz,it)
                end
            end

            if istwing != 0
                update!(u,ix,iyp,iz,itp,ix,iy,iz,it)
            end
        

        end

        if iszwing != 0
            update!(u,ix,iy,izp,it,ix,iy,iz,it)
            if istwing != 0
                update!(u,ix,iy,izp,itp,ix,iy,iz,it)
            end

        end

        if istwing != 0
            update!(u,ix,iy,iz,itp,ix,iy,iz,it)
        end


    end

    #function set_wing!(u::GaugeFields{T}) where T <: SUn
    function set_wing!(u::GaugeFields{SU{NC}}) where NC
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
    
    function Base.display(x::GaugeFields{SU{NC}}) where NC
        NX=x.NX
        NY=x.NY
        NZ=x.NZ
        NT=x.NT


        icum = 0
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        icum += 1
                        println("i = $icum")
                        println("ix,iy,iz,it = $ix $iy $iz $it ")
                        display(x[:,:,ix,iy,iz,it])
                        println("\t")
                    end
                end
            end
        end
    end

    function Base.display(a::GaugeFields_1d) 
        for i=1:a.NV
            println("i = $i")
            display(a[:,:,i])
            println("\t")
        end
    end

    function substitute!(a::Array{T,1},b::Array{T,1}) where T <: GaugeFields
        for μ=1:4
            substitute!(a[μ],b[μ])
        end
    end

    function substitute!(a::GaugeFields,b::GaugeFields)
        a.g[:,:,:,:,:,:] = copy(b.g)
        return 
    end


    function substitute!(a::GaugeFields{SU{NC}},b::GaugeFields_1d{SU{NC}}) where NC
        NT = a.NT
        NY = a.NY
        NZ = a.NZ
        NX = a.NX


        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX               
                        icum = (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix 
                        for k2=1:NC
                            @simd for k1=1:NC
                                a[k1,k2,ix,iy,iz,it] = b[k1,k2,icum]
                            end
                        end

                        #func!(a,b,icum,ix,iy,iz,it)

                    end
                end
            end
        end

    end

    function substitute!(a::GaugeFields_1d{SU{NC}},b::GaugeFields{SU{NC}}) where NC

        NT = b.NT
        NY = b.NY
        NZ = b.NZ
        NX = b.NX


        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX               
                        icum = (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix 


                        for j=1:NC
                            @simd for i=1:NC
                                a[i,j,icum] = b[i,j,ix,iy,iz,it]
                            end
                        end

                    end
                end
            end
        end

    end

    function substitute!(a::GaugeFields_1d,α::Number,b::GaugeFields_1d)
        n1,n2,n3 = size(a.g)

        for i3=1:n3
            for i2=1:n2
                @simd for i1=1:n1
                    a.g[i1,i2,i3]= α*b.g[i1,i2,i3]
                end
            end
        end

    end

    function clear!(a::GaugeFields)
        @. a.g = 0
        return 
    end

    function clear!(a::GaugeFields_1d)
        @. a.g = 0
        return 
    end


    function IdentityGauges(NC,NX,NY,NZ,NT,NDW)
        U = GaugeFields(NC,NDW,NX,NY,NZ,NT)
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
        return U
    end

    
    function RandomGauges(NC,NX,NY,NZ,NT,NDW)
        U = GaugeFields(NC,NDW,NX,NY,NZ,NT)
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
        normalize!(U)
        return U
    end

    function Oneinstanton(NC,NX,NY,NZ,NT,NDW)
        U = Array{SU2GaugeFields,1}(undef,4)
        for μ = 1:4
            U[μ] = GaugeFields(NC,NDW,NX,NY,NZ,NT)
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

        #normalize!(U)
        return U

    end

        



    function Base.similar(x::GaugeFields) 
        return GaugeFields(x.NC,x.NDW,x.NX,x.NY,x.NZ,x.NT)
    end

    function Base.similar(x::GaugeFields_1d) 
        return GaugeFields_1d(x.NC,x.NV) 
    end

    function Base.similar(x::Array{T,1}) where T <: GaugeFields
        xout = Array{T,1}(undef,4)
        for μ=1:4
            xout[μ] = similar(x[μ])
        end
        return xout
    end

    function Base.copyto!(Unew::GaugeFields,U::GaugeFields)
        Unew = deepcopy(U)
        return
    end

    function elementwise_tr!(s,u::GaugeFields_1d{SU{NC}},v::GaugeFields_1d{SU{NC}}) where NC

        NV=u.NV

        for i=1:NV
            for k1=1:NC
                @simd for k2=1:NC
                    s[i] += real(u[k1,k2,i]*conj(v[k1,k2,i]))
                end
            end
        end
        return

    end


    function LinearAlgebra.tr(a::GaugeFields{SU{NC}}) where NC
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
                        end
                    end
                end
            end
        end
        return s

    end

    function LinearAlgebra.tr(a::GaugeFields_1d{SU{NC}}) where NC
        NV=a.NV

        s = 0

        for i=1:NV

            @simd for k=1:NC
                s += a[k,k,i]
            end

        end
        return s

    end

    function elementwise_apply!(a::GaugeFields{SU{NC}},func!::Function) where NC
        NT = a.NT
        NY = a.NY
        NZ = a.NZ
        NX = a.NX


        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX   
                        for j=1:NC
                            @simd for i=1:NC
                                func!(a,i,j,ix,iy,iz,it,b...)
                            end
                        end
                    end
                end
            end
        end
        return
    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{SU{NC}},a::GaugeFields_1d{SU{NC}},b::GaugeFields_1d{SU{NC}}) where NC


        NV=a.NV
        
        for i=1:NV
            #mulabc!(a,b,c,i)
            
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] = 0
                    @simd for k3=1:NC
                        c[k1,k2,i] += a[k1,k3,i]*b[k3,k2,i]
                    end
                end
            end
            
        end

    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{SU{NC}},a::GaugeFields{SU{NC}},b::GaugeFields_1d{SU{NC}}) where NC


        NV=a.NV
        NT = a.NT
        NZ = a.NZ
        NY = a.NY
        NX = a.NX
        #NC = b.NC

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        i= (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix
                        for k2=1:NC                            
                            for k1=1:NC
                                c[k1,k2,i] = 0
                                @simd for k3=1:NC
                                    c[k1,k2,i] += a[k1,k3,ix,iy,iz,it]*b[k3,k2,i]
                                end
                            end
                        end
                    end
                end
            end
        end

        

    end



    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::GaugeFields_1d{T},b::GaugeFields_1d{T}) where T <: SU3
        NV=a.NV
        #NC=a.NC
        #mulabc! = NCmul(NC)

        @simd for i=1:NV
            c[1,1,i] = a[1,1,i] * b[1,1,i] +
                a[1,2,i] * b[2,1,i] +
                a[1,3,i] * b[3,1,i]
            c[1,2,i] = a[1,1,i] * b[1,2,i] + 
                    a[1,2,i] * b[2,2,i] +
                    a[1,3,i] * b[3,2,i]
            c[1,3,i] = a[1,1,i] * b[1,3,i] +
                    a[1,2,i] * b[2,3,i] +
                    a[1,3,i] * b[3,3,i]

            c[2,1,i] = a[2,1,i] * b[1,1,i] +
                    a[2,2,i] * b[2,1,i] +
                    a[2,3,i] * b[3,1,i]
            c[2,2,i] = a[2,1,i] * b[1,2,i] + 
                    a[2,2,i] * b[2,2,i] +
                    a[2,3,i] * b[3,2,i]
            c[2,3,i] = a[2,1,i] * b[1,3,i] + 
                    a[2,2,i] * b[2,3,i] +
                    a[2,3,i] * b[3,3,i]

            c[3,1,i] = a[3,1,i] * b[1,1,i] + 
                    a[3,2,i] * b[2,1,i] +
                    a[3,3,i] * b[3,1,i]
            c[3,2,i] = a[3,1,i] * b[1,2,i] + 
                    a[3,2,i] * b[2,2,i] +
                    a[3,3,i] * b[3,2,i]
            c[3,3,i] = a[3,1,i] * b[1,3,i] + 
                    a[3,2,i] * b[2,3,i] + 
                    a[3,3,i] * b[3,3,i]
        end

    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::GaugeFields_1d{T},b::GaugeFields_1d{T}) where T <: SU2
        NV=a.NV
        #NC=a.NC
        #mulabc! = NCmul(NC)

        @simd for i=1:NV
            c[1,1,i] = a[1,1,i] * b[1,1,i] +
                a[1,2,i] * b[2,1,i] 

            c[1,2,i] = a[1,1,i] * b[1,2,i] + 
                    a[1,2,i] * b[2,2,i] 


            c[2,1,i] = a[2,1,i] * b[1,1,i] +
                    a[2,2,i] * b[2,1,i] 
            c[2,2,i] = a[2,1,i] * b[1,2,i] + 
                    a[2,2,i] * b[2,2,i] 

        end

    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{SU{NC}},a::GaugeFields_1d{SU{NC}},b::Adjoint_GaugeFields_1d{SU{NC}}) where NC

        NV=a.NV


        for i=1:NV
            #mulabc!(a,b,c,i)
            
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] = 0
                    @simd for k3=1:NC
                        c[k1,k2,i] += a[k1,k3,i]*conj(b.parent[k2,k3,i])
                    end
                end
            end
            
                        
        end
        

    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{SU{NC}},a::Number) where NC
        NV=c.NV


        #NC=a.NC
        #mulabc! = NCmul_aconjb(NC)

        for i=1:NV
            for k2= 1:NC
                @simd for k1=1:NC
                    c[k1,k2,i] *= a
                end
            end                        
        end
        return
    
    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::GaugeFields_1d{T},b::Adjoint_GaugeFields_1d{T}) where T <: SU3
        
        NV=a.NV
        #NC=a.NC
        #mulabc! = NCmul_aconjb(NC)

        @simd for i=1:NV
            c[1,1,i] = a[1,1,i] * conj(b.parent[1,1,i]) +
                a[1,2,i] * conj(b.parent[1,2,i]) +
                a[1,3,i] * conj(b.parent[1,3,i])
            c[1,2,i] = a[1,1,i] * conj(b.parent[2,1,i]) + 
                    a[1,2,i] * conj(b.parent[2,2,i]) +
                    a[1,3,i] * conj(b.parent[2,3,i])
            c[1,3,i] = a[1,1,i] * conj(b.parent[3,1,i])+
                    a[1,2,i] * conj(b.parent[3,2,i]) +
                    a[1,3,i] * conj(b.parent[3,3,i])

            c[2,1,i] = a[2,1,i] * conj(b.parent[1,1,i]) +
                    a[2,2,i] * conj(b.parent[1,2,i]) +
                    a[2,3,i] * conj(b.parent[1,3,i])
            c[2,2,i] = a[2,1,i] * conj(b.parent[2,1,i]) + 
                    a[2,2,i] * conj(b.parent[2,2,i]) +
                    a[2,3,i] * conj(b.parent[2,3,i])
            c[2,3,i] = a[2,1,i] * conj(b.parent[3,1,i]) + 
                    a[2,2,i] * conj(b.parent[3,2,i]) +
                    a[2,3,i] * conj(b.parent[3,3,i])

            c[3,1,i] = a[3,1,i] * conj(b.parent[1,1,i]) + 
                    a[3,2,i] * conj(b.parent[1,2,i]) +
                    a[3,3,i] * conj(b.parent[1,3,i])
            c[3,2,i] = a[3,1,i] * conj(b.parent[2,1,i]) + 
                    a[3,2,i] * conj(b.parent[2,2,i]) +
                    a[3,3,i] * conj(b.parent[2,3,i])
            c[3,3,i] = a[3,1,i] * conj(b.parent[3,1,i]) + 
                    a[3,2,i] * conj(b.parent[3,2,i]) + 
                    a[3,3,i] * conj(b.parent[3,3,i])

                        
        end
    
    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::GaugeFields_1d{T},b::Adjoint_GaugeFields_1d{T}) where T <: SU2
        
        NV=a.NV
        #NC=a.NC
        #mulabc! = NCmul_aconjb(NC)

        @simd for i=1:NV
            c[1,1,i] = a[1,1,i] * conj(b.parent[1,1,i]) +
                a[1,2,i] * conj(b.parent[1,2,i]) 
            c[1,2,i] = a[1,1,i] * conj(b.parent[2,1,i]) + 
                    a[1,2,i] * conj(b.parent[2,2,i]) 

            c[2,1,i] = a[2,1,i] * conj(b.parent[1,1,i]) +
                    a[2,2,i] * conj(b.parent[1,2,i]) 
            c[2,2,i] = a[2,1,i] * conj(b.parent[2,1,i]) + 
                    a[2,2,i] * conj(b.parent[2,2,i]) 

                        
        end
    
    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::Adjoint_GaugeFields_1d{T},b::Adjoint_GaugeFields_1d{T}) where T <: SU2
        
        NV=a.parent.NV
        #NC=a.NC
        #mulabc! = NCmul_aconjb(NC)

        @simd for i=1:NV
            c[1,1,i] = conj(a.parent[1,1,i]) * conj(b.parent[1,1,i]) +
                conj(a.parent[2,1,i] )* conj(b.parent[1,2,i]) 
            c[1,2,i] = conj(a.parent[1,1,i]) * conj(b.parent[2,1,i]) + 
                conj(a.parent[2,1,i]) * conj(b.parent[2,2,i]) 

            c[2,1,i] = conj(a.parent[1,2,i]) * conj(b.parent[1,1,i]) +
                conj(a.parent[2,2,i]) * conj(b.parent[1,2,i]) 
            c[2,2,i] = conj(a.parent[1,2,i]) * conj(b.parent[2,1,i]) + 
                conj(a.parent[2,2,i]) * conj(b.parent[2,2,i])        
        end
    
    end


    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::Adjoint_GaugeFields_1d{T},b::Adjoint_GaugeFields_1d{T}) where T <: SU3
        
        NV=a.parent.NV
        #NC=a.NC
        #mulabc! = NCmul_aconjb(NC)

        @simd for i=1:NV
            c[1,1,i] = conj(a.parent[1,1,i]) * conj(b.parent[1,1,i]) +
                conj(a.parent[2,1,i] )* conj(b.parent[1,2,i]) +
                conj(a.parent[3,1,i] )* conj(b.parent[1,3,i]) 
            c[1,2,i] = conj(a.parent[1,1,i]) * conj(b.parent[2,1,i]) + 
                conj(a.parent[2,1,i]) * conj(b.parent[2,2,i]) +
                conj(a.parent[3,1,i]) * conj(b.parent[2,3,i]) 

            c[1,3,i] = conj(a.parent[1,1,i]) * conj(b.parent[3,1,i]) + 
                conj(a.parent[2,1,i]) * conj(b.parent[3,2,i]) +
                conj(a.parent[3,1,i]) * conj(b.parent[3,3,i]) 

            c[2,1,i] = conj(a.parent[1,2,i]) * conj(b.parent[1,1,i]) +
                conj(a.parent[2,2,i]) * conj(b.parent[1,2,i]) +
                conj(a.parent[3,2,i]) * conj(b.parent[1,3,i]) 

            c[2,2,i] = conj(a.parent[1,2,i]) * conj(b.parent[2,1,i]) + 
                conj(a.parent[2,2,i]) * conj(b.parent[2,2,i])      +
                conj(a.parent[3,2,i]) * conj(b.parent[2,3,i])   

            c[2,3,i] = conj(a.parent[1,2,i]) * conj(b.parent[3,1,i]) + 
                conj(a.parent[2,2,i]) * conj(b.parent[3,2,i]) +
                conj(a.parent[3,2,i]) * conj(b.parent[3,3,i]) 

            c[3,1,i] = conj(a.parent[1,3,i]) * conj(b.parent[1,1,i]) +
                conj(a.parent[2,3,i]) * conj(b.parent[1,2,i]) +
                conj(a.parent[3,3,i]) * conj(b.parent[1,3,i]) 

            c[3,2,i] = conj(a.parent[1,3,i]) * conj(b.parent[2,1,i]) + 
                conj(a.parent[2,3,i]) * conj(b.parent[2,2,i])      +
                conj(a.parent[3,3,i]) * conj(b.parent[2,3,i])   

            c[3,3,i] = conj(a.parent[1,3,i]) * conj(b.parent[3,1,i]) + 
                conj(a.parent[2,3,i]) * conj(b.parent[3,2,i]) +
                conj(a.parent[3,3,i]) * conj(b.parent[3,3,i]) 
        end
    
    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{SU{NC}},a::Adjoint_GaugeFields_1d{SU{NC}},b::Adjoint_GaugeFields_1d{SU{NC}}) where NC

        NV=c.NV
        #NC=c.NC
        #println(NC)


        for i=1:NV
            #mulabc!(a,b,c,i)
            
    
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] = 0
                    @simd for k3=1:NC
                        c[k1,k2,i] += conj(a.parent[k3,k1,i])*conj(b.parent[k2,k3,i])
                    end
                end
            end
            
                        
        end

    end


    function LinearAlgebra.mul!(c::GaugeFields_1d{SU{NC}},a::Adjoint_GaugeFields_1d{SU{NC}},b::GaugeFields_1d{SU{NC}}) where NC

        NV=c.NV
        #NC=c.NC


        for i=1:NV
            #mulabc!(a,b,c,i)
            
    
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] = 0
                    @simd for k3=1:NC
                        c[k1,k2,i] += conj(a.parent[k3,k1,i])*b[k3,k2,i]
                    end
                end
            end
            
                        
        end

    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{SU{NC}},Udag::Adjoint_GaugeFields{SU{NC}}) where NC #Udag*c
        NV=c.NV
        NT = Udag.parent.NT
        NZ = Udag.parent.NZ
        NY = Udag.parent.NY
        NX = Udag.parent.NX
        #NC = b.NC
        #println(NC)
        tmp = zeros(ComplexF64,NC,NC)


        
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        i= (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix
                        for k2=1:NC                            
                            for k1=1:NC
                                tmp[k1,k2] = 0
                                @simd for k3=1:NC
                                    tmp[k1,k2] += conj(Udag.parent[k3,k1,ix,iy,iz,it])*c[k3,k2,i]
                                end
                            end
                        end
                        @. c[:,:,i] = tmp[:,:]
                    end
                end
            end
        end


    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{SU{NC}},U::GaugeFields{SU{NC}}) where NC #U*c
        NV=c.NV
        NT = U.NT
        NZ = U.NZ
        NY = U.NY
        NX = U.NX
        #NC = b.NC
        #println(NC)
        tmp = zeros(ComplexF64,NC,NC)


        
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        i= (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix
                        for k2=1:NC                            
                            for k1=1:NC
                                tmp[k1,k2] = 0
                                @simd for k3=1:NC
                                    tmp[k1,k2] += U[k1,k3,ix,iy,iz,it]*c[k3,k2,i]
                                end
                            end
                        end
                        @. c[:,:,i] = tmp[:,:]
                    end
                end
            end
        end


    end


    function LinearAlgebra.mul!(c::GaugeFields_1d{SU{NC}},a::Adjoint_GaugeFields{SU{NC}},b::GaugeFields_1d{SU{NC}}) where NC

        NV=c.NV
        #NC=c.NC

        NT = a.parent.NT
        NZ = a.parent.NZ
        NY = a.parent.NY
        NX = a.parent.NX

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        i= (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix
                        for k2=1:NC                            
                            for k1=1:NC
                                c[k1,k2,i] = 0
                                @simd for k3=1:NC
                                    c[k1,k2,i] += conj(a.parent[k3,k1,ix,iy,iz,it])*b[k3,k2,i]
                                end
                            end
                        end
                    end
                end
            end
        end


    end


    function LinearAlgebra.mul!(c::GaugeFields_1d{SU{NC}},a::Adjoint_GaugeFields_1d{SU{NC}},b::GaugeFields_1d{SU{NC}}) where NC

        NV=c.NV
        #NC=c.NC


        for i=1:NV
            #mulabc!(a,b,c,i)
            
    
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] = 0
                    @simd for k3=1:NC
                        c[k1,k2,i] += conj(a.parent[k3,k1,i])*b[k3,k2,i]
                    end
                end
            end
            
                        
        end

    end


    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::Adjoint_GaugeFields_1d{T},b::GaugeFields_1d{T}) where T <: SU3
        NV=c.NV
        #NC=c.NC

        #mulabc! = NCmul_conjab(c.NC)

        @simd for i=1:NV

            c[1,1,i] = conj(a.parent[1,1,i]) * b[1,1,i] +
                conj(a.parent[2,1,i]) * b[2,1,i] +
                conj(a.parent[3,1,i]) * b[3,1,i]
            c[1,2,i] = conj(a.parent[1,1,i]) * b[1,2,i] + 
                    conj(a.parent[2,1,i]) * b[2,2,i] +
                    conj(a.parent[3,1,i]) * b[3,2,i]
            c[1,3,i] = conj(a.parent[1,1,i]) * b[1,3,i] +
                    conj(a.parent[2,1,i]) * b[2,3,i] +
                    conj(a.parent[3,1,i]) * b[3,3,i]

            c[2,1,i] = conj(a.parent[1,2,i]) * b[1,1,i] +
                    conj(a.parent[2,2,i]) * b[2,1,i] +
                    conj(a.parent[3,2,i]) * b[3,1,i]
            c[2,2,i] = conj(a.parent[1,2,i]) * b[1,2,i] + 
                    conj(a.parent[2,2,i]) * b[2,2,i] +
                    conj(a.parent[3,2,i]) * b[3,2,i]
            c[2,3,i] = conj(a.parent[1,2,i]) * b[1,3,i] + 
                    conj(a.parent[2,2,i]) * b[2,3,i] +
                    conj(a.parent[3,2,i]) * b[3,3,i]

            c[3,1,i] = conj(a.parent[1,3,i]) * b[1,1,i] + 
                    conj(a.parent[2,3,i]) * b[2,1,i] +
                    conj(a.parent[3,3,i]) * b[3,1,i]
            c[3,2,i] = conj(a.parent[1,3,i]) * b[1,2,i] + 
                    conj(a.parent[2,3,i]) * b[2,2,i] +
                    conj(a.parent[3,3,i]) * b[3,2,i]
            c[3,3,i] = conj(a.parent[1,3,i]) * b[1,3,i] + 
                    conj(a.parent[2,3,i]) * b[2,3,i] + 
                    conj(a.parent[3,3,i]) * b[3,3,i]
            #=
            for k2=1:NC                            
                for k1=1:NC
                    c[k1,k2,i] = 0
                    for k3=1:NC
                        c[k1,k2,i] += conj(a.parent[k3,k1,i])*b[k3,k2,i]
                    end
                end
            end
            =#
                        
        end

    end

    function LinearAlgebra.mul!(c::GaugeFields_1d{T},a::Adjoint_GaugeFields_1d{T},b::GaugeFields_1d{T}) where T <: SU2
        NV=c.NV
        #NC=c.NC

        #mulabc! = NCmul_conjab(c.NC)

        @simd for i=1:NV

            c[1,1,i] = conj(a.parent[1,1,i]) * b[1,1,i] +
                conj(a.parent[2,1,i]) * b[2,1,i] 
            c[1,2,i] = conj(a.parent[1,1,i]) * b[1,2,i] + 
                    conj(a.parent[2,1,i]) * b[2,2,i] 


            c[2,1,i] = conj(a.parent[1,2,i]) * b[1,1,i] +
                    conj(a.parent[2,2,i]) * b[2,1,i] 
            c[2,2,i] = conj(a.parent[1,2,i]) * b[1,2,i] + 
                    conj(a.parent[2,2,i]) * b[2,2,i] 


                        
        end

    end

    function add!(c::GaugeFields_1d{SU{NC}},a::GaugeFields_1d{SU{NC}}) where NC

        NV=c.NV

        for i=1:NV
            #ncadd!(a,c,i)
            
            for k2=1:NC                            
                @simd for k1=1:NC
                    c[k1,k2,i] += a[k1,k2,i] 
                end
            end
            
        end


    end

    function muladd!(c::GaugeFields_1d{SU{NC}},α::Number,a::GaugeFields_1d{SU{NC}}) where NC

        NV=c.NV

        for i=1:NV
            #ncadd!(a,c,i)
            
            for k2=1:NC                            
                @simd for k1=1:NC
                    c[k1,k2,i] += α*a[k1,k2,i] 
                end
            end
            
        end


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

    function gauge_shift!(a::GaugeFields_1d{SU{NC}},ν::N,b::GaugeFields{SU{NC}}) where {N <: Int,NC}


        if ν == 0
            substitute!(a,b)
            return
        end

        idel = zeros(Int8,4)
        if ν > 0
            idel[ν] = +1
        else
            idel[-ν] = -1
        end

        
        NT = b.NT
        NZ = b.NZ
        NY = b.NY
        NX = b.NX
        #NC = b.NC
        #println(NC)


        
        for it=1:NT
            it1 = it + idel[4]
            for iz=1:NZ
                iz1 = iz + idel[3]
                for iy=1:NY
                    iy1 = iy+idel[2]
                    for ix=1:NX
                        ix1 = ix+idel[1]
                        icum = (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix
                        #func!(a,b,icum,ix1,iy1,iz1,it1)

                        
                        
                        for k2=1:NC
                            @simd for k1=1:NC
                                a[k1,k2,icum] = b[k1,k2,ix1,iy1,iz1,it1]
                            end
                        end
                        
                        
                        
                    end
                end
            end
        end
        return
    end

    function gauge_shift!(a::GaugeFields_1d{SU{NC}},dir::Tuple,b::GaugeFields{SU{NC}}) where NC 

        #....  Stop for the exceptional case   ...
        if dir[1] != 0 && dir[2] != 0
            if dir[1] == dir[2]
                println("Sorry this case is not yet considered!")
                println("dir : ",dir[1],"\t",dir[2])
                exit()
            end
        end

        if dir[1] == 0 && dir[2] == 0
            subst!(ibush,b,a)
            return
        end

        idel1 = zeros(Int64,4)
        idel2 = zeros(Int64,4)
        if dir[1] > 0
            idel1[+dir[1]] = +1
        elseif dir[1] < 0
            idel1[-dir[1]] = -1
        end

        if dir[2] > 0
            idel2[+dir[2]] = +1
        elseif dir[2] < 0
            idel2[-dir[2]] = -1
        end


        
        NT = b.NT
        NZ = b.NZ
        NY = b.NY
        NX = b.NX


        #func! = NCsubstitute_63(b.NC)
        #        NC = b.NC
        

        
        for it=1:NT
            it1 = it + idel1[4] + idel2[4]
            for iz=1:NZ
                iz1 = iz + idel1[3]+idel2[3]
                for iy=1:NY
                    iy1 = iy+idel1[2]+idel2[2]
                    for ix=1:NX
                        ix1 = ix+idel1[1] + idel2[1]
                        icum = (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix
                        #func!(a,b,icum,ix1,iy1,iz1,it1)

                        
                        for k2=1:NC
                            @simd for k1=1:NC
                                a[k1,k2,icum] = b[k1,k2,ix1,iy1,iz1,it1]
                            end
                        end
                        
                    end
                end
            end
        end
        return
    end

    struct Shift_set
        shift::NTuple{4,Int8}
    end

    function gauge_shift_all!(a::GaugeFields_1d{T},dir,b::GaugeFields{T}) where T <: SUn 
        shift = Shift_set(dir)
        gauge_shift!(a,shift,b) 
        return
    end

    function periodiccheck(ix,iy,iz,it,NX,NY,NZ,NT,NDW)
        return periodiccheck(ix,NX,NDW),periodiccheck(iy,NY,NDW),periodiccheck(iz,NZ,NDW),periodiccheck(it,NT,NDW)
    end

    function periodiccheck(ix,NX,NDW)
        ix1 = ix + ifelse(ix > NX+NDW,-NX,0)
        ix1 += ifelse(ix < 1-NDW,NX,0)
        return ix1
    end

    function gauge_shift!(a::GaugeFields_1d{SU{NC}},shift::Shift_set,b::GaugeFields{SU{NC}}) where NC 

        
        NT = b.NT
        NZ = b.NZ
        NY = b.NY
        NX = b.NX
        NDW= b.NDW


        #func! = NCsubstitute_63(b.NC)
        #        NC = b.NC
        

        
        for it=1:NT
            it1 = it + shift.shift[4]
            it1 = periodiccheck(it1,NT,NDW)
            for iz=1:NZ
                iz1 = iz + shift.shift[3]
                iz1 = periodiccheck(iz1,NZ,NDW)
                for iy=1:NY
                    iy1 = iy+ shift.shift[2]
                    iy1 = periodiccheck(iy1,NY,NDW)
                    for ix=1:NX
                        ix1 = ix+ shift.shift[1]
                        ix1 = periodiccheck(ix1,NX,NDW)
                        
                        icum = (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix
                        #func!(a,b,icum,ix1,iy1,iz1,it1)

                        
                        for k2=1:NC
                            @simd for k1=1:NC
                                a[k1,k2,icum] = b[k1,k2,ix1,iy1,iz1,it1]
                            end
                        end
                        
                    end
                end
            end
        end
        return
    end

    function make_staple!(staple::GaugeFields_1d,U,μ,temp1,temp2,temp3)
        clear!(staple)

        idir2 = zeros(Int64,2)

        for ν=1:4
            if ν == μ
                continue
            end

                        #=
            c       x+nu temp2
            c        .---------.
            c        I         I
            c  temp1 I         I
            c        I         I
            c        .         .
            c        x        x+mu
            =#
            substitute!(temp1,U[ν])
            gauge_shift!(temp2,ν,U[μ])
            #mul!(temp3,U[ν],temp2)
            mul!(temp3,temp1,temp2)

            gauge_shift!(temp1,μ,U[ν])
            mul!(temp2,temp3,temp1')
            add!(staple,temp2)



        end
    end



    function make_staple_double!(staple::GaugeFields_1d,U,μ,temp1,temp2,temp3)
        clear!(staple)
        #println("mu = ",μ)

        #loops = make_plaq_staple(μ)
        loops = make_plaq_staple_prime(μ)
        evaluate_wilson_loops!(staple,loops,U,[temp1,temp2,temp3])
        return


        clear!(staple)
        
        for ν=1:4
            if ν == μ
                continue
            end

                        #=
            c       x+nu temp2
            c        .---------.
            c        I         I
            c  temp1 I         I
            c        I         I
            c        .         .
            c        x        x+mu
            =#
            substitute!(temp1,U[ν])
            #println("[]")
            gauge_shift!(temp2,ν,U[μ])

            #println(ν)
            
            mul!(temp3,temp1,temp2)

            gauge_shift!(temp1,μ,U[ν])

            #println(μ)
            mul!(temp2,temp3,temp1')



            add!(staple,temp2)

            #=
            c        x
            c        .         .
            c        I         I
            c  temp1 I         I
            c        I         I
            c        .---------.
            c      x-nu        x+mu-nu
            =#
            gauge_shift!(temp1,-ν,U[ν])
            #println(-ν)
            gauge_shift!(temp2,-ν,U[μ])
            #println(-ν)

            mul!(temp3,temp1',temp2)




            #idir2 = Ivec2(μ,-ν)
            
            #println((μ,-ν))
            gauge_shift!(temp1,(μ,-ν),U[ν])

            mul!(temp2,temp3,temp1)

            add!(staple,temp2)
            #println("mu")

            #if ν == 4
            #    display(staple)
            #    exit()

            #end



        end
        
    end

    


    function TA(vin::T) where T <: GaugeFields_1d
        vout = deepcopy(vin)
        projlink!(vout,vin)
        return vout
    end

"""
    c-----------------------------------------------------c
c     !!!!!   vin and vout should be different vectors
c
c     Projectin of the etraceless antiermite part 
c     vout = x/2 - Tr(x)/6
c     wher   x = vin - Conjg(vin)      
c-----------------------------------------------------c
    """
    function projlink!(vout::SU3GaugeFields,vin::SU3GaugeFields)
        
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

    const sr3 = sqrt(3)
    const sr3i = 1/sr3

    function projlink!(vout::SU3GaugeFields_1d,vin::SU3GaugeFields_1d)
        fac13 = 1/3
        nv = vin.NV

        @simd for i=1:nv
            v11 = vin[1,1,i]
            v22 = vin[2,2,i]
            v33 = vin[3,3,i]

            tri = fac13*(imag(v11)+imag(v22)+imag(v33))

            vout[1,1,i] = (imag(v11)-tri)*im
            vout[2,2,i] = (imag(v22)-tri)*im
            vout[3,3,i] = (imag(v33)-tri)*im

        end

        @simd for i=1:nv
            v12 = vin[1,2,i]
            v13 = vin[1,3,i]
            v21 = vin[2,1,i]
            v23 = vin[2,3,i]
            v31 = vin[3,1,i]
            v32 = vin[3,2,i]

            x12 = v12 - conj(v21)
            x13 = v13 - conj(v31)
            x23 = v23 - conj(v32)
        
            x21 = - conj(x12)
            x31 = - conj(x13)
            x32 = - conj(x23)

            vout[1,2,i] = 0.5  * x12
            vout[1,3,i] = 0.5  * x13
            vout[2,1,i] = 0.5  * x21
            vout[2,3,i] = 0.5  * x23
            vout[3,1,i] = 0.5  * x31
            vout[3,2,i] = 0.5  * x32
            
        end
    end

    function TA!(vout::Array{T,2},vin::Array{T,2}) where T <: Number
        NC,_ = size(vout)
        fac1N = 1/NC
        tri = 0.0
        @simd for k=1:NC
            tri += imag(vin[k,k])
        end
        tri *= fac1N
        @simd for k=1:NC
            vout[k,k] = (imag(vin[k,k])-tri)*im
        end
        for k1=1:NC
            @simd for k2=k1+1:NC
                vv = 0.5*(vin[k1,k2] - conj(vin[k2,k1]))
                vout[k1,k2] = vv
                vout[k2,k1] = -conj(vv)
            end
        end

    end

    function projlink!(vout::SUNGaugeFields_1d{NC},vin::SUNGaugeFields_1d{NC}) where NC
        #NC = vout.NC
        fac1N = 1/NC
        nv = vin.NV

        for i=1:nv
            tri = 0.0
            @simd for k=1:NC
                tri += imag(vin[k,k,i])
            end
            tri *= fac1N
            @simd for k=1:NC
                vout[k,k,i] = (imag(vin[k,k,i])-tri)*im
            end

        end

        for i=1:nv
            for k1=1:NC
                @simd for k2=k1+1:NC
                    vv = 0.5*(vin[k1,k2,i] - conj(vin[k2,k1,i]))
                    vout[k1,k2,i] = vv
                    vout[k2,k1,i] = -conj(vv)
                end
            end
            
        end
    end


    function projlink!(vout::U1GaugeFields_1d,vin::U1GaugeFields_1d) 
        for i=1:vin.NV
            vout[1,1,i] = 2*imag(vin[1,1,i])
            #imag(vin[1,1,i])
        end
    end
    """
    c-----------------------------------------------------c
    c     !!!!!   vin and vout should be different vectors
    c
    c     Projection of the traceless antiermite part 
    c     vout = x/2 - Tr(x)/6
    c     where   x = vin - Conjg(vin)      
    c-----------------------------------------------------c
    """
    function projlink!(vout::SU2GaugeFields,vin::SU2GaugeFields)
        
        fac12 = 1/2
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

                        tri = fac12*(imag(v11)+imag(v22))

                        vout[1,1,ix,iy,iz,it] = (imag(v11)-tri)*im
                        vout[2,2,ix,iy,iz,it] = (imag(v22)-tri)*im
                    end
                end
            end
        end

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    @simd for ix=1:NX

                        v12 = vin[1,2,ix,iy,iz,it]
                        #v13 = vin[1,3,ix,iy,iz,it]
                        v21 = vin[2,1,ix,iy,iz,it]

                        x12 = v12 - conj(v21)

                        x21 = - conj(x12)

                        vout[1,2,ix,iy,iz,it] = 0.5  * x12

                        vout[2,1,ix,iy,iz,it] = 0.5  * x21


                    end
                end
            end
        end

    end



    function projlink!(vout::SU2GaugeFields_1d,vin::SU2GaugeFields_1d)
        fac12 = 1/2
        nv = vin.NV

        @simd for i=1:nv
            v11 = vin[1,1,i]
            v22 = vin[2,2,i]

            tri = fac12*(imag(v11)+imag(v22))

            vout[1,1,i] = (imag(v11)-tri)*im
            vout[2,2,i] = (imag(v22)-tri)*im

        end

        @simd for i=1:nv
            v12 = vin[1,2,i]
            v21 = vin[2,1,i]


            x12 = v12 - conj(v21)

        
            x21 = - conj(x12)

            vout[1,2,i] = 0.5  * x12
            vout[2,1,i] = 0.5  * x21

            
        end

    end

    const sr3ih = 0.5*sr3i
    const sqr3inv = sr3i


    function lambdamul(b::SU3GaugeFields_1d,a::SU3GaugeFields_1d,k,generator::Nothing)
        #=
        c----------------------------------------------------------------------c
        c     b = (lambda_k/2)*a
        C             lambda_k : GellMann matrices. k=1, 8 
        c----------------------------------------------------------------------c
        =#
        NV = a.NV

        if k==1
            @simd for i=1:NV
                b[1,1,i] = 0.5 * a[2,1,i] 
                b[1,2,i] = 0.5 * a[2,2,i]
                b[1,3,i] = 0.5 * a[2,3,i]
                b[2,1,i] = 0.5 * a[1,1,i]
                b[2,2,i] = 0.5 * a[1,2,i]
                b[2,3,i] = 0.5 * a[1,3,i]
                b[3,1,i] = 0
                b[3,2,i] = 0
                b[3,3,i] = 0
            end
        elseif k==2
            @simd for i=1:NV
                b[1,1,i] = -0.5*im * a[2,1,i] 
                b[1,2,i] = -0.5*im * a[2,2,i]
                b[1,3,i] = -0.5*im * a[2,3,i]
                b[2,1,i] =  0.5*im * a[1,1,i]
                b[2,2,i] =  0.5*im * a[1,2,i]
                b[2,3,i] =  0.5*im * a[1,3,i]
                b[3,1,i] = 0
                b[3,2,i] = 0
                b[3,3,i] = 0
            end
        elseif k==3
            @simd for i=1:NV
                b[1,1,i] =  0.5 * a[1,1,i] 
                b[1,2,i] =  0.5 * a[1,2,i]
                b[1,3,i] =  0.5 * a[1,3,i]
                b[2,1,i] = -0.5 * a[2,1,i]
                b[2,2,i] = -0.5 * a[2,2,i]
                b[2,3,i] = -0.5 * a[2,3,i]
                b[3,1,i] = 0
                b[3,2,i] = 0
                b[3,3,i] = 0
            end
        elseif k==4
            @simd for i=1:NV
                b[1,1,i] = 0.5 * a[3,1,i] 
                b[1,2,i] = 0.5 * a[3,2,i]
                b[1,3,i] = 0.5 * a[3,3,i]
                b[2,1,i] = 0
                b[2,2,i] = 0
                b[2,3,i] = 0
                b[3,1,i] = 0.5 * a[1,1,i]
                b[3,2,i] = 0.5 * a[1,2,i]
                b[3,3,i] = 0.5 * a[1,3,i]
            end
        elseif k==5
            @simd for i=1:NV
                b[1,1,i] = -0.5*im * a[3,1,i] 
                b[1,2,i] = -0.5*im * a[3,2,i]
                b[1,3,i] = -0.5*im * a[3,3,i]
                b[2,1,i] = 0
                b[2,2,i] = 0
                b[2,3,i] = 0
                b[3,1,i] =  0.5*im * a[1,1,i]
                b[3,2,i] =  0.5*im * a[1,2,i]
                b[3,3,i] =  0.5*im * a[1,3,i]
            end
        elseif k==6
            @simd for i=1:NV
                b[1,1,i] = 0
                b[1,2,i] = 0
                b[1,3,i] = 0
                b[2,1,i] = 0.5 * a[3,1,i] 
                b[2,2,i] = 0.5 * a[3,2,i]
                b[2,3,i] = 0.5 * a[3,3,i]
                b[3,1,i] = 0.5 * a[2,1,i]
                b[3,2,i] = 0.5 * a[2,2,i]
                b[3,3,i] = 0.5 * a[2,3,i]
            end
        elseif k==7
            @simd for i=1:NV
                b[1,1,i] = 0
                b[1,2,i] = 0
                b[1,3,i] = 0
                b[2,1,i] = -0.5*im * a[3,1,i] 
                b[2,2,i] = -0.5*im * a[3,2,i]
                b[2,3,i] = -0.5*im * a[3,3,i]
                b[3,1,i] =  0.5*im * a[2,1,i]
                b[3,2,i] =  0.5*im * a[2,2,i]
                b[3,3,i] =  0.5*im * a[2,3,i]
            end
        elseif k==8
            @simd for i=1:NV
                b[1,1,i] =  sr3ih * a[1,1,i] 
                b[1,2,i] =  sr3ih * a[1,2,i]
                b[1,3,i] =  sr3ih * a[1,3,i]
                b[2,1,i] =  sr3ih * a[2,1,i] 
                b[2,2,i] =  sr3ih * a[2,2,i]
                b[2,3,i] =  sr3ih * a[2,3,i]
                b[3,1,i] = -sqr3inv * a[3,1,i]
                b[3,2,i] = -sqr3inv * a[3,2,i]
                b[3,3,i] = -sqr3inv * a[3,3,i]
            end
        else
            error("k should be k <= 8 but k = $k")
        end

        return
    end

    function lambdamul(b::SU2GaugeFields_1d,a::SU2GaugeFields_1d,k,generator::Nothing)
        #=
        c----------------------------------------------------------------------c
        c     b = (lambda_k/2)*a
        C             lambda_k : GellMann matrices. k=1, 8 
        c----------------------------------------------------------------------c
        =#
        NV = a.NV
    
    
        if k==1
            @simd for i=1:NV
                b[1,1,i] = -0.5*im* a[2,1,i]*im
                b[1,2,i] = -0.5*im * a[2,2,i]*im
    
                b[2,1,i] = -0.5*im * a[1,1,i]*im
                b[2,2,i] = -0.5*im * a[1,2,i]*im
    
    
            end
        elseif k==2
            @simd for i=1:NV
                b[1,1,i] = -0.5 * a[2,1,i] *im
                b[1,2,i] = -0.5 * a[2,2,i]*im
    
                b[2,1,i] =  0.5 * a[1,1,i]*im
                b[2,2,i] =  0.5 * a[1,2,i]*im
    
            end
        elseif k==3
            @simd for i=1:NV
                b[1,1,i] =  -0.5*im * a[1,1,i] *im
                b[1,2,i] =  -0.5*im * a[1,2,i]*im
    
                b[2,1,i] = 0.5*im * a[2,1,i]*im
                b[2,2,i] = 0.5*im * a[2,2,i]*im
    
            end
        else
            error("k should be k <= 3 but k = $k")
        end
    
        return
    end

    function lambdamul(b::SUNGaugeFields_1d{NC},a::SUNGaugeFields_1d{NC},k,generator) where NC
        #=
        c----------------------------------------------------------------------c
        c     b = (lambda_k/2)*a
        C             lambda_k : GellMann matrices. k=1, 8 
        c----------------------------------------------------------------------c
        =#
        NV = a.NV
        #NC = generator.NC
        matrix = generator.generator[k]
        for i=1:NV
            for k2=1:NC
                for k1=1:NC
                    b[k1,k2,i] = 0
                    @simd for l=1:NC
                        b[k1,k2,i] += matrix[k1,l]*a[l,k2,i]/2
                    end
                end
            end
        end
    

        return
    end
    

    #function calc_Polyakov(u::Array{T,1}) where T <: GaugeFields
    function calc_Polyakov(u::Array{GaugeFields{SU{NC}},1}) where NC
        NX = u[1].NX
        NY = u[1].NY
        NZ = u[1].NZ
        NT = u[1].NT
        #NC = u[1].NC

        Pol = zeros(ComplexF64,NC,NC,NX,NY,NZ)
        
        tmp1= zeros(ComplexF64,NC,NC)
        tmp2= zero(tmp1)

        #set_wing!(u)

        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    for k2=1:NC
                        @simd for k1=1:NC
                            Pol[k1,k2,ix,iy,iz] = u[4][k1,k2,ix,iy,iz,1]
                        end
                    end
                end
            end
        end

        for it=2:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        for k2=1:NC
                            @simd for k1=1:NC
                                tmp1[k1,k2] = Pol[k1,k2,ix,iy,iz]
                                tmp2[k1,k2] = u[4][k1,k2,ix,iy,iz,it]
                            end
                        end
                        Pol[:,:,ix,iy,iz] = tmp1*tmp2
                    end
                end
            end
        end

        avePol = 0.0
        for iz=1:NX
            for iy=1:NY
                for ix=1:NX
                    @simd for k1=1:NC
                        avePol += Pol[k1,k1,ix,iy,iz]
                    end
                end
            end
        end
        avePol /= NX*NY*NZ


        return avePol
    end

    function calc_Plaq_notrace_1d(U::Array{T,1},μ,ν,origin) where T <: GaugeFields

        NC = U[1].NC
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT

        temp1 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp2 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp3 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        plaq = GaugeFields_1d(NC,NX,NY,NZ,NT)

        loop = make_plaq(μ,ν,origin)

        evaluate_wilson_loops!(plaq,loop ,U,[temp1,temp2,temp3])

        return plaq
    end

    function calc_Plaq_notrace_1d(U::Array{T,1},μ,ν) where T <: GaugeFields
        origin = [0,0,0,0]
        return calc_Plaq_notrace_1d(U,μ,ν,origin)
    end

    function calc_Plaq_notrace_1d(U::Array{T,1}) where T <: GaugeFields

        NC = U[1].NC
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT

        temp1 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp2 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp3 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        plaq = GaugeFields_1d(NC,NX,NY,NZ,NT)

        loop = make_plaq()

        evaluate_wilson_loops!(plaq,loop ,U,[temp1,temp2,temp3])

        return plaq
    end


    function calc_Plaq(U::Array{T,1}) where T <: GaugeFields
        NC = U[1].NC
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT

        temp1 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp2 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        temp3 = GaugeFields_1d(NC,NX,NY,NZ,NT)
        staple = GaugeFields_1d(NC,NX,NY,NZ,NT)

        return calc_Plaq!(U::Array{T,1},temp1,temp2,temp3,staple)
    end

    function calc_Plaq!(U::Array{T,1},temp1,temp2,temp3,staple) where T <: GaugeFields
        plaq = 0
        for μ=1:4
            clear!(staple)
            make_staple!(staple,U,μ,temp1,temp2,temp3)
            substitute!(temp1,U[μ])

            mul!(temp2,temp1,staple')
            #mul!(temp2,U[μ],staple')
            plaq += tr(temp2)
            #println(plaq)
            #exit()
        end
        #plaq /= 
        return plaq*0.5
    end

    function calc_Plaq!(U::Array{T,1},temps::Array{T_1d,1}) where {T <: GaugeFields,T_1d <: GaugeFields_1d}
        plaq = 0
        temp1 = temps[1]
        temp2 = temps[2]
        temp3 = temps[3]
        staple = temps[4]

        return calc_Plaq!(U,temp1,temp2,temp3,staple)
    end

    function calc_Plaq!(U::Array{T,1},temps::Array{T_1d,1}) where {T <: GaugeFields,T_1d <: GaugeFields}
        plaq = 0
        temp1 = temps[1]
        temp2 = temps[2]
        temp3 = temps[3]
        staple = temps[4]

        for μ=1:4
            clear!(staple)
            make_staple!(staple,U,μ,temp1,temp2,temp3)
            substitute!(temp1,U[μ])

            mul!(temp2,temp1,staple')
            #mul!(temp2,U[μ],staple')
            plaq += tr(temp2)
            #println(plaq)
            #exit()
        end
        #plaq /= 
        return plaq*0.5
    end



    function calc_GaugeAction(U::Array{T,1},gparam::GaugeActionParam_standard,temps::Array{T_1d,1}) where {T <: GaugeFields,T_1d <: GaugeFields_1d}
        plaq = calc_Plaq!(U,temps)

        Sg = -plaq*gparam.β/gparam.NTRACE
        return Sg,plaq
    end

    function calc_GaugeAction(U::Array{T,1},gparam::GaugeActionParam_standard,temps::Array{T_1d,1}) where {T <: GaugeFields,T_1d <: GaugeFields}
        plaq = calc_Plaq!(U,temps)

        Sg = -plaq*gparam.β/gparam.NTRACE
        return Sg,plaq
    end


    function calc_GaugeAction(U::Array{T,1},gparam::GaugeActionParam_autogenerator,temps::Array{T_1d,1}) where {T <: GaugeFields,T_1d <: GaugeFields_1d}
        Sg = 0
        loopaction = temps[4]        
        trs = zeros(ComplexF64,gparam.numactions)
        
        for i = 1:gparam.numactions
            #if gparam.βs[i] != 0
                evaluate_wilson_loops!(loopaction,gparam.loops[i],U,temps[1:3])
                trs[i] = tr(loopaction)
                #println("i = $i ",trs[i])
                sg = (-gparam.βs[i]/gparam.NTRACE)*trs[i] #/2
                #println("$i-th actions ",sg)
                Sg += sg#(-gparam.βs[i]/gparam.NTRACE)*tr(loopaction)/2       
            #end     
        end

        return Sg,trs
        
    end


    function normalize!(u::Array{ComplexF64,2})
        NC,_ = size(u)
        if NC == 2
            normalize2!(u)
        elseif NC == 3
            normalize3!(u)
        else
            normalizeN!(u)
        end
        return 
    end

    function normalize2!(u)
        α = u[1,1]
        β = u[2,1]
        detU = abs(α)^2 + abs(β)^2
        u[1,1] = α/detU
        u[2,1] = β/detU
        u[1,2] = -conj(β)/detU
        u[2,2] = conj(α)/detU
    end

    function normalizeN!(u)
        gramschmidt!(u)
    end

    function normalize3!(u)
        w1 = 0
        w2 = 0
        for ic=1:3
            w1 += u[2,ic]*conj(u[1,ic])
            w2 += u[1,ic]*conj(u[1,ic])
        end
        zerock2 = w2
        if zerock2 == 0 
            println("w2 is zero  !!  (in normlz)")
            println("u[1,1),u[1,2),u[1,3) : ",u[1,1], "\t",u[1,2],"\t", u[1,3])
        end

        w1 = -w1/w2

        x4 = (u[2,1]) + w1*u[1,1]
        x5 = (u[2,2]) + w1*u[1,2]
        x6 = (u[2,3]) + w1*u[1,3]

        w3 = x4*conj(x4) + x5*conj(x5) + x6*conj(x6)

        zerock3 = w3
        if zerock3 == 0
            println("w3 is zero  !!  (in normlz)")
            println("x4, x5, x6 : $x4, $x5, $x6")
            exit()
        end

        u[2,1] = x4
        u[2,2] = x5
        u[2,3] = x6

        w3 = 1/sqrt(w3)
        w2 = 1/sqrt(w2)

        u[1,1] = u[1,1]*w2
        u[1,2] = u[1,2]*w2
        u[1,3] = u[1,3]*w2
        u[2,1] = u[2,1]*w3
        u[2,2] = u[2,2]*w3
        u[2,3] = u[2,3]*w3

        if zerock2*zerock3 == 0 
            println("!! devided by zero !! (in normalize)")
            println("w2 or w3 in normlz is zero !!")
            println("w2, w3 : $w2, $w3   ")
            exit()
        end

        m3complv3!(u)
    end

    """
    normalize!(u)
    ----------------------------------------------------------------------c
         normalizes links                                                 c
         input  x : arbitrary 3*3 complex matrix                          c
         output x : SU(3) matrix 
    """
    function normalize!(u::GaugeFields{T}) where T <: SU3
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

    """
    normalize!(u)
    ----------------------------------------------------------------------c
         normalizes links                                                 c
         input  x : arbitrary 2*2 complex matrix                          c
         output x : SU(2) matrix 
    """
    function normalize!(u::GaugeFields{T}) where T <: SU2
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

    function gramschmidt!(v)
        n = size(v)[1]
        for i=1:n
            @simd for j=1:i-1
                v[:,i] = v[:,i] - v[:,j]'*v[:,i]*v[:,j]
            end
            v[:,i] = v[:,i]/norm(v[:,i])
        end
    end

    function normalize!(u::GaugeFields{T}) where T <: SUn
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


    function m3complv3!(a)
        aa = zeros(Float64,18)
        aa[ 1] = real( a[1,1])
        aa[ 2] = imag(a[1,1])
        aa[ 3] = real( a[1,2])
        aa[ 4] = imag(a[1,2])
        aa[ 5] = real( a[1,3])
        aa[ 6] = imag(a[1,3])
        aa[ 7] = real( a[2,1])
        aa[ 8] = imag(a[2,1])
        aa[ 9] = real( a[2,2])
        aa[10] = imag(a[2,2])
        aa[11] = real( a[2,3])
        aa[12] = imag(a[2,3])

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

        a[3,1] = aa[13]+im*aa[14]
        a[3,2] = aa[15]+im*aa[16]
        a[3,3] = aa[17]+im*aa[18]
        return
    end


    function m3complv!(a::SU3GaugeFields)
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

    function get_insertposition(tensor,ix,iy,iz,it,NX,NY,NZ,NT,NDW)
        ix1 = ix + tensor.position[1]
        iy1 = iy + tensor.position[2]
        iz1 = iz + tensor.position[3]
        it1 = it + tensor.position[4]
        ix1,iy1,iz1,it1 =periodiccheck(ix1,iy1,iz1,it1,NX,NY,NZ,NT,NDW)
        return ix1,iy1,iz1,it1
    end

    function calc_coordinate_tensor(w,origin)
        numloops = length(w)
        coordinates = Array{NTuple{4,Int8},1}(undef,numloops)
        if numloops == 1
            coordinates[1] = origin
        end
        
        for k=1:numloops-1
            if k == 1
                coordinates[k] = origin
            end
            loopk = w[k]
            if loopk[1] == 1
                coordinates[k+1] = (coordinates[k][1]+loopk[2],coordinates[k][2],coordinates[k][3],coordinates[k][4])
            elseif loopk[1] == 2
                coordinates[k+1] = (coordinates[k][1],coordinates[k][2]+loopk[2],coordinates[k][3],coordinates[k][4])
            elseif loopk[1] == 3
                coordinates[k+1] = (coordinates[k][1],coordinates[k][2],coordinates[k][3]+loopk[2],coordinates[k][4])
            elseif loopk[1] == 4
                coordinates[k+1] = (coordinates[k][1],coordinates[k][2],coordinates[k][3],coordinates[k][4]+loopk[2])
            end
        end
        return coordinates
    end

    function evaluate_tensor_matrix!(temp2,temp3,U,tensorlist,origin,ix,iy,iz,it,NX,NY,NZ,NT,NDW)

        # #=
        wi = tensorlist
        
        numloops = length(wi)    
        #println(wi)
        coordinates = calc_coordinate_tensor(wi,origin)
        shifts = calc_shift_tensor(wi,origin)

        #println("coordinates ",coordinates)
        #println("shift ",shifts)
       
        #exit()
        loopk = wi[1]
        ix1,iy1,iz1,it1 = shift_xyzt(shifts[1],ix,iy,iz,it)
        ix1,iy1,iz1,it1 =periodiccheck(ix1,iy1,iz1,it1,NX,NY,NZ,NT,NDW)

        
        matrixout= U[loopk[1]][:,:,ix1,iy1,iz1,it1]

        loopk1_2 = loopk[2]

        

        if numloops == 1
            if loopk1_2 == -1
                matrixout[:,:] = U[loopk[1]][:,:,ix1,iy1,iz1,it1]'
            end
            #println("matrixout 1 ",matrixout)
        else
            #println("matrixout 2 ",matrixout)
            #gauge_shift_all!(temp1,shifts[1],U[loopk[1]])
            for k=2:numloops
                loopk = wi[k]
                
                #gauge_shift_all!(temp2,shifts[k],U[loopk[1]])
                ix1,iy1,iz1,it1 = shift_xyzt(shifts[k],ix,iy,iz,it)
                ix1,iy1,iz1,it1 =periodiccheck(ix1,iy1,iz1,it1,NX,NY,NZ,NT,NDW)
                #println(loopk[1])
                
                temp2[:,:] = U[loopk[1]][:,:,ix1,iy1,iz1,it1]
                #println("temp2 ",temp2)

                multiply_12!(temp3,matrixout,temp2,k,loopk,loopk1_2)

                matrixout,temp3 = temp3,matrixout
            end
        end

        return  matrixout

        
    end

    function calc_shift_tensor(w,origin)
        numloops = length(w)
        coordinates = calc_coordinate_tensor(w,origin)
        shifts = Array{NTuple{4,Int8},1}(undef,numloops)
        for k=1:numloops
            loopk = w[k]
            if loopk[2] == 1
                shifts[k] = coordinates[k]
            elseif loopk[2] == -1
                if loopk[1] == 1
                    shifts[k] = (coordinates[k][1]+loopk[2],coordinates[k][2],coordinates[k][3],coordinates[k][4])
                elseif loopk[1] == 2
                    shifts[k] = (coordinates[k][1],coordinates[k][2]+loopk[2],coordinates[k][3],coordinates[k][4])
                elseif loopk[1] == 3
                    shifts[k] = (coordinates[k][1],coordinates[k][2],coordinates[k][3]+loopk[2],coordinates[k][4])
                elseif loopk[1] == 4
                    shifts[k] = (coordinates[k][1],coordinates[k][2],coordinates[k][3],coordinates[k][4]+loopk[2])
                end
            end
        end
        return shifts
    end

    function evaluate_tensor_matrix!(matrixout,temp2,temp3,U::Array{GaugeFields{SU{NC}}},tensorlist,origin,ix,iy,iz,it,NX,NY,NZ,NT,NDW) where NC

        # #=
        wi = tensorlist
        
        numloops = length(wi)    
        #println(wi)
        coordinates = calc_coordinate_tensor(wi,origin)
        shifts = calc_shift_tensor(wi,origin)

        #println("coordinates ",coordinates)
        #println("shift ",shifts)
       
        #exit()
        loopk = wi[1]
        ix1,iy1,iz1,it1 = shift_xyzt(shifts[1],ix,iy,iz,it)
        ix1,iy1,iz1,it1 =periodiccheck(ix1,iy1,iz1,it1,NX,NY,NZ,NT,NDW)

        for i=1:NC
            for j=1:NC
                matrixout[j,i]= U[loopk[1]][j,i,ix1,iy1,iz1,it1]
            end
        end

        loopk1_2 = loopk[2]

        

        if numloops == 1
            if loopk1_2 == -1
                for i=1:NC
                    for j=1:NC
                        matrixout[j,i]= conj(U[loopk[1]][i,j,ix1,iy1,iz1,it1])
                    end
                end

                #matrixout[:,:] = U[loopk[1]][:,:,ix1,iy1,iz1,it1]'
            end
            #println("matrixout 1 ",matrixout)
        else
            #println("matrixout 2 ",matrixout)
            #gauge_shift_all!(temp1,shifts[1],U[loopk[1]])
            for k=2:numloops
                loopk = wi[k]
                
                #gauge_shift_all!(temp2,shifts[k],U[loopk[1]])
                ix1,iy1,iz1,it1 = shift_xyzt(shifts[k],ix,iy,iz,it)
                ix1,iy1,iz1,it1 =periodiccheck(ix1,iy1,iz1,it1,NX,NY,NZ,NT,NDW)
                #println(loopk[1])

                for i=1:NC
                    for j=1:NC
                        temp2[j,i]= U[loopk[1]][j,i,ix1,iy1,iz1,it1]
                    end
                end
                
                #temp2[:,:] = U[loopk[1]][:,:,ix1,iy1,iz1,it1]
                #println("temp2 ",temp2)

                multiply_12!(temp3,matrixout,temp2,k,loopk,loopk1_2)

                #matrixout,temp3 = temp3,matrixout
                for i=1:NC
                    for j=1:NC
                        matrixout[j,i]= temp3[j,i]
                    end
                end

            end
        end

        #return  matrixout

        
    end

    function evaluate_tensor!(leftmatrix,rightmatrix,tensor::Tensor_wilson_lines,U::Array{GaugeFields{SU{NC}},1},ix,iy,iz,it,mat_tmps) where NC
        #println(tensor.left)
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        NDW = U[1].NDW
        #NC = U[1].NC
        #leftmatrix .= 0 
        #rightmatrix .= 0
        
        temp1 = mat_tmps[1] #zeros(ComplexF64,NC,NC)
        temp2 = mat_tmps[2] #zero(temp1)
        temp3 = mat_tmps[3] #zero(temp1)
        
        loopk1_2 = (0,0)

        if tensor.left[1] == (0,0)
            #leftmatrix = zero(temp1)
            leftmatrix .=  0 #zero(temp1)
            for i=1:NC 
                leftmatrix[i,i] = 1
            end
            #println(leftmatrix)
        else
            leftstart = get_leftstartposition(tensor)
            #println(tensor.left)
            evaluate_tensor_matrix!(leftmatrix,temp2,temp3,U,tensor.left,leftstart,ix,iy,iz,it,NX,NY,NZ,NT,NDW)
            #leftmatrix = evaluate_tensor_matrix!(temp2,temp3,U,tensor.left,leftstart,ix,iy,iz,it,NX,NY,NZ,NT,NDW)
        end
        #println("left")
        #println(leftmatrix)

        if tensor.right[1] == (0,0)
            rightmatrix .= 0#zero(temp1)
            #rightmatrix = zero(temp1)
            for i=1:NC 
                rightmatrix[i,i] = 1
            end
            #println(rightmatrix)
        else
            #println("right")
            #println(tensor.right)
            rightstart = get_rightstartposition(tensor)
            evaluate_tensor_matrix!(rightmatrix,temp2,temp3,U,tensor.right,rightstart,ix,iy,iz,it,NX,NY,NZ,NT,NDW)
            #rightmatrix = evaluate_tensor_matrix!(temp2,temp3,U,tensor.right,rightstart,ix,iy,iz,it,NX,NY,NZ,NT,NDW)
        end
        #println("right")
        #println(rightmatrix)
        return #leftmatrix,rightmatrix
        
    
    end

    function calc_Qmatrix(staple_nu_set,nu,ρs,U,NC,ix,iy,iz,it)
        C = zeros(ComplexF64,NC,NC)
        Ctmp = zeros(ComplexF64,NC,NC)
        
        #evaluate_wilson_loops!(C,loops,U,ix,iy,iz,it)
        for i=1:length(ρs)
            Ctmp .= 0
            evaluate_wilson_loops!(Ctmp,staple_nu_set[i][nu],U,ix,iy,iz,it)
            C[:,:] .+= ρs[i]*Ctmp[:,:] 
            

        end
        #Ω = ρ*C[:,:]*U[nu][:,:,ix,iy,iz,it]'
        #if (ix,iy,iz,it) == (2,2,2,2) && nu == 1
        #    println("ρs = ",ρs)
        #    println("C")
        #    display(C)
        #    println("\t")
        #end

        Ω = C[:,:]*U[nu][:,:,ix,iy,iz,it]'
        if NC == 1
            #Q .= (-1/2)*(Ω' .- Ω)
            #Q .= (1/2)*(Ω' .- Ω)
            #Q = (Ω' .- Ω)
            Q = -(Ω' .- Ω)
        else
            trterm =-(1/(2NC))*tr(Ω' - Ω)
            Q = -(1/2)*(Ω' - Ω) 
            for i=1:NC
                Q[i,i] += trterm
            end
            #+ (1/(2NC))*tr(Ω' - Ω)*I0_2
            #Q = (1/2)*(Ω' - Ω) - (1/(2NC))*tr(Ω' - Ω)*I0_2
        end
        return Q
    end
    
    function calc_Qmatrix!(Q,staple_nu_set,nu,ρs,U::Array{GaugeFields{SU{NC}},1},ix,iy,iz,it,tmp_matrices) where NC
        C = tmp_matrices[1] #zeros(ComplexF64,NC,NC)
        Ctmp = tmp_matrices[2] #zeros(ComplexF64,NC,NC)
        Ω = tmp_matrices[3]
        
        
        #evaluate_wilson_loops!(C,loops,U,ix,iy,iz,it)
        C .= 0
        for i=1:length(ρs)
            Ctmp .= 0
            #evaluate_wilson_loops!(Ctmp,staple_nu_set[i][nu],U,ix,iy,iz,it)
            evaluate_wilson_loops!(Ctmp,staple_nu_set[i][nu],U,ix,iy,iz,it,tmp_matrices[4:6])

            C[:,:] .+= ρs[i]*Ctmp #[:,:] 
        end

        for i=1:NC
            for j=1:NC
                tmp_matrices[4][j,i] = conj(U[nu][i,j,ix,iy,iz,it])
            end
        end

        #mul!(Ω,C,U[nu][:,:,ix,iy,iz,it]')
        mul!(Ω,C,tmp_matrices[4])
        #Ω = C[:,:]*U[nu][:,:,ix,iy,iz,it]'
        if NC == 1
            #Q .= (-1/2)*(Ω' .- Ω)
            #Q .= (1/2)*(Ω' .- Ω)
            #Q = (Ω' .- Ω)
            Q[:,:] = -(Ω' .- Ω)
        #elseif NC == 2
        else

            @. Q[:,:] = -(1/2)*(Ω' - Ω)
            #tr1 = (1/(2NC))*tr(Ω')
            tr2 = (1/(2NC))*tr(Ω)
            for i=1:NC
                Q[i,i] += conj(tr2)-tr2
                #Q[i,i] += tr1-tr2
            end
            #Q = -(1/2)*(Ω' - Ω) + (1/(2NC))*tr(Ω' - Ω)*I0_2
            #Q = (1/2)*(Ω' - Ω) - (1/(2NC))*tr(Ω' - Ω)*I0_2
        end
        return 
        #return Q
    end

    function calc_Bmatrix(q,Q,NC)

        #q = sqrt((-1/2)*tr(Q^2))
        B = -(-sin(q)*I0_2 +(cos(q)/q -sin(q)/q^2 )*Q)*(1/2q)
    end

    function calc_Bmatrix!(B,q,Q,NC)
        @assert NC == 2 "NC should be 2! now $NC"
        mul!(B,cos(q)/q -sin(q)/q^2,Q)
        for i=1:NC
            B[i,i] += -sin(q)
        end
        B .*= -1/2q
        #B[:,:] .= (cos(q)/q -sin(q)/q^2 )*Q

        #q = sqrt((-1/2)*tr(Q^2))
        #B = -(-sin(q)*I0_2 +(cos(q)/q -sin(q)/q^2 )*Q)*(1/2q)
    end

    const eps_Q = 1e-18
    const I0_1 = zeros(1,1) .+ 1
    const I0_2 = [1 0;0 1]


    function calc_Mmatrix!(M,dSdUnu,staple_nu_set,nu,ρs,U::Array{GaugeFields{SU{NC}},1},ix,iy,iz,it,tmp_matrices) where NC
        #M = zeros(ComplexF64,NC,NC)
        Q = tmp_matrices[1]
        #println("Q in M")
        #@time Q = calc_Qmatrix(staple_nu_set,nu,ρs,U,NC,ix,iy,iz,it)
        calc_Qmatrix!(Q,staple_nu_set,nu,ρs,U,ix,iy,iz,it,tmp_matrices[2:end])
        B = tmp_matrices[2]
        UdSdU = tmp_matrices[3]
        trQ2 = tr(Q^2)
        if abs(trQ2) > eps_Q
            #println("Q ")
            #display(Q)
            #println("\t")
            if NC == 2
                q = sqrt((-1/2)*trQ2)
                calc_Bmatrix!(B,q,Q,NC)
                #println(B)

                #B = calc_Bmatrix(q,Q,NC)
                #println(B)
                #exit()
                #println("B = ",B)
                #Unu = U[nu][:,:,ix,iy,iz,it]
                

                for i=1:NC
                    for j=1:NC
                        tmp_matrices[4][j,i] = U[nu][j,i,ix,iy,iz,it]
                    end
                end

                mul!(UdSdU,tmp_matrices[4],dSdUnu)
                #println("Q in M ",Q)
                #println("q in M ",q)
                #M[:,:] = tr(Unu*dSdUnu*B)*Q + (sin(q)/q)*Unu*dSdUnu
                #M[:,:] = (sin(q)/q)*UdSdU
                trsum = 0.0im
                for i=1:NC
                    for j=1:NC
                        trsum += UdSdU[i,j]*B[j,i]
                    end
                end
                #M[:,:] .+= trsum*Q

                for i=1:NC
                    for j=1:NC
                        M[j,i] = (sin(q)/q)*UdSdU[j,i] + trsum*Q[j,i]
                    end
                end
                
                

                #M[:,:] = tr(Unu*dSdUnu*B)*Q + (sin(q)/q)*Unu*dSdUnu
                #Q = calc_Qmatrix(staple_nu,nu,ρ,U,NC,ix,iy,iz,it)
            elseif NC == 1
                M[:,:] = U[nu][:,:,ix,iy,iz,it]*dSdUnu*exp.(Q)
            else
                #=
                @time begin
                M .= 0
                A = U[nu][:,:,ix,iy,iz,it]*dSdUnu
                nmax = 3
                for n=1:nmax
                    fn = factorial(n)
                    for m=0:n-1
                        M += Q^(n-1-m)*A*Q^m/fn
                    end
                end
                println(M)
                end
                =#
                e,v = eigen(Q) #

                #=

                A = tmp_matrices[4]
                A .= 0
                for i=1:NC
                    for j=1:NC
                        for l=1:NC
                            A[j,i] += U[nu][j,l,ix,iy,iz,it]*dSdUnu[l,i]
                        end
                    end
                end
                B .= 0
                for i=1:NC
                    for j=1:NC
                        B[j,i] += conj(v[j])*A[j,i]*v[i]
                    end
                end
                =#

                for i=1:NC
                    for j=1:NC
                        tmp_matrices[4][j,i] = U[nu][j,i,ix,iy,iz,it]
                    end
                end
                mul!(UdSdU,tmp_matrices[4],dSdUnu)


                #=
                A = tmp_matrices[4]
                A .= 0
                for i=1:NC
                    for j=1:NC
                        for l=1:NC
                            A[j,i] += U[nu][j,l,ix,iy,iz,it]*dSdUnu[l,i]
                        end
                    end
                end
                =#
                mul!(tmp_matrices[4],UdSdU,v)
                mul!(B,v',tmp_matrices[4])
                
                #=
                B .= 0
                for i=1:NC
                    for j=1:NC
                        for k=1:NC
                            for l=1:NC
                                B[j,i] += conj(v[l,j])*UdSdU[l,k]*v[k,i]
                            end
                        end
                    end
                end
                =#
                
                #println(B)
                #println(A)

                #A = U[nu][:,:,ix,iy,iz,it]*dSdUnu
                #println(A)
                #B = v'*A*v

                #println(B)
                #exit()
                #println(e)
                #M0 = zero(M)
                tmp_matrices[4] .= 0
                #M0 .= 0

                nmax = 3
                for n=1:nmax
                    fn = factorial(n)
                    for l=1:NC
                        for k=1:NC
                            factor = B[k,l]*e[k]^(n-1)/fn
                            elek = e[l]/e[k]
                            for m=0:n-1
                                tmp_matrices[4][k,l] += elek^m*factor
                            end
                        end
                    end
                    #println("n = $n")
                    #display(v*M0*v')
                    #println("\t")
                end
                #M = deepcopy(M0)
                #=
                M .= 0
                for i=1:NC
                    for j=1:NC
                        M[j,i] += v[j]*M0[j,i]*conj(v[i])
                    end
                end
                =#
                #M[:,:] = v*M0*v'
                #println(M)

                mul!(UdSdU,tmp_matrices[4],v')
                mul!(M,v,UdSdU)

                #=
                M .= 0
                for i=1:NC
                    for j=1:NC
                        for k=1:NC
                            for l=1:NC
                                M[j,i] += v[j,l]*tmp_matrices[4][l,k]*conj(v[i,k])
                            end
                        end
                    end
                end
                =#
                #println(M)
                #exit()
                

                


                #error("not supoorted yet")
            end
        else
            mul!(M,U[nu][:,:,ix,iy,iz,it],dSdUnu)
            #M[:,:] =  U[nu][:,:,ix,iy,iz,it]*dSdUnu
        end

        #return

        #println("M")
        #display(M)
        #println("\t")
        

        
        


        #exit()



        return 
        #return M
    end


    function calc_matrix_tensors!(leftmatrix,rightmatrix,tensor,U,ix,iy,iz,it,mat_tmps)
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        NDW = U[1].NDW
        evaluate_tensor!(leftmatrix,rightmatrix,tensor,U,ix,iy,iz,it,mat_tmps) 
        ix1,iy1,iz1,it1 = get_insertposition(tensor,ix,iy,iz,it,NX,NY,NZ,NT,NDW)
        #display(tensor)
        #println("insert ", ix1,iy1,iz1,it1)
        return  ix1,iy1,iz1,it1
    end

    function get_1dindex(ix,iy,iz,it,NX,NY,NZ,NT)
        #println("before ",(ix,iy,iz,it))

        ix1 = ix + ifelse(ix < 1,NX,0) + ifelse(ix > NX,-NX,0)
        iy1 = iy + ifelse(iy < 1,NY,0) + ifelse(iy > NY,-NY,0)
        iz1 = iz + ifelse(iz < 1,NZ,0) + ifelse(iz > NZ,-NZ,0)
        it1 = it + ifelse(it < 1,NT,0) + ifelse(it > NT,-NT,0)

        #println("after ",(ix1,iy1,iz1,it1))

        icum = (((it1-1)*NZ+iz1-1)*NY+iy1-1)*NX+ix1 

        #println("1d index ", icum)

        return icum
    end

    function calc_Λmatrix!(Λ,M,NC)
        #println("M= ", M)
        if NC == 1
            #Λ = -(M - M')
            @. Λ[:,:] = (M - M')
        elseif NC == 2
            #Λ = (1/2)*(M - M') - (1/(2NC))*tr(M - M')*I0_2
            @. Λ[:,:] = (1/2)*(M - M')
            trM = (1/(2NC))*(M[1,1]-conj(M[1,1]) + M[2,2] - conj(M[2,2]))#  tr(M - M')
            #trM = (1/(2NC))*tr(M - M')
            for i=1:NC
                Λ[i,i] += - trM
            end
            #Λ = 2*Λ
        else
            @. Λ[:,:] = (1/2)*(M - M')
            trM = (1/(2NC))*tr(M - M')
            for i=1:NC
                Λ[i,i] += - trM
            end
        end
        #display(Λ)
        #println("\t")
        #exit()
        return 
#        return Λ
    end

    function evaluate_tensor_lines(mu,dSdU::Array{T_1d,1},smearing,U::Array{GaugeFields{SU{NC}},1},umu::Tuple{I,I},ix,iy,iz,it,ρs) where {NC,I <: Int,T_1d <: GaugeFields_1d}
        
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        NDW = U[1].NDW
        #NC = U[1].NC
        V = zeros(ComplexF64,NC,NC)
        Λ = zero(V)
        M = zero(V)
        tmpmat = zero(V)
        leftmatrix = zero(V)
        rightmatrix = zero(V)
        dSdUnu = zero(V)
        n = 7
        tmp_matrices = Array{Array{ComplexF64,2},1}(undef,n)
        for i=1:n
            tmp_matrices[i] = zero(V)
        end

        t = smearing.tensor_derivative
        num = length(t)
        #ρs = smearing.ρs


        for i=1:num
            Ti = t[i]
            for nu=1:4
                #staple_nu = gparam.smearing.staples_for_stout[nu]
                staple_nu_set = smearing.staples_for_stout
                timu = Ti[nu]

                #println((ix,iy,iz,it))

                for (k,timu_k) in enumerate(timu)
                    if haskey(timu_k,umu)
                        timuumu = timu_k[umu]
                        for m =1:length(timuumu)
                            tensor = timuumu[m]
                            #println("tensors")
                            #leftmatrix,rightmatrix,ix1,iy1,iz1,it1 = calc_matrix_tensors(tensor,U,ix,iy,iz,it)
                            ix1,iy1,iz1,it1 = calc_matrix_tensors!(leftmatrix,rightmatrix,tensor,U,ix,iy,iz,it,tmp_matrices)
                            icum = get_1dindex(ix1,iy1,iz1,it1,NX,NY,NZ,NT)
                            for i=1:NC
                                for j=1:NC
                                    dSdUnu[j,i] = dSdU[nu][j,i,icum]
                                end
                            end

                            #dSdUnu = dSdU[nu][:,:,icum]
                            #println("M")
                            #@time M = calc_Mmatrix(dSdUnu,staple_nu_set,nu,ρs,NC,U,ix1,iy1,iz1,it1)
                            calc_Mmatrix!(M,dSdUnu,staple_nu_set,nu,ρs,U,ix1,iy1,iz1,it1,tmp_matrices)
            
                            #println("Λ")
                            calc_Λmatrix!(Λ,M,NC)
                            #@time Λ = calc_Λmatrix(M,NC)
                            #@time Umu = U[nu][:,:,ix1,iy1,iz1,it1]
                            #@time Vtmp = ρs[i]*rightmatrix*Umu'*Λ*leftmatrix

                            mul!(tmp_matrices[1],Λ,leftmatrix)
                            for i=1:NC
                                for j=1:NC
                                    tmp_matrices[3][j,i] = conj(U[nu][i,j,ix1,iy1,iz1,it1]) #U[nu][:,:,ix1,iy1,iz1,it1]'
                                end
                            end
                            mul!(tmp_matrices[2],tmp_matrices[3],tmp_matrices[1])
                            mul!(tmp_matrices[1],rightmatrix,tmp_matrices[2])
                            tmp_matrices[1] .*= ρs[i]

                            #println("Vtmp $i,$nu")
                            #display(Vtmp/ρ)
                            #println("\t")


                            @. V += tmp_matrices[1]
                            #V += Vtmp

                            #=
                            println("$i $nu $k")
                            display(leftmatrix)
                            println("\t")
                            =#
                            #println("Vtmp, ", Vtmp)
                        end
                    end
                end
            end
            #exit()
            
        end

        t = smearing.tensor_derivative_dag
        num = length(t)


        for i=1:num
            Ti = t[i]
            for nu=1:4
                staple_nu_set = smearing.staples_for_stout
                #staple_nu = gparam.smearing.staples_for_stout[nu]
                timu = Ti[nu]

                #println((ix,iy,iz,it))

                for (k,timu_k) in enumerate(timu)
                    if haskey(timu_k,umu)
                        timuumu = timu_k[umu]
                        for m =1:length(timuumu)
                            tensor = timuumu[m]
                            #leftmatrix,rightmatrix,ix1,iy1,iz1,it1 = calc_matrix_tensors(tensor,U,ix,iy,iz,it)
                            ix1,iy1,iz1,it1 = calc_matrix_tensors!(leftmatrix,rightmatrix,tensor,U,ix,iy,iz,it,tmp_matrices)
                            icum = get_1dindex(ix1,iy1,iz1,it1,NX,NY,NZ,NT)
                            for i=1:NC
                                for j=1:NC
                                    dSdUnu[j,i] = dSdU[nu][j,i,icum]
                                end
                            end
                            #dSdUnu = dSdU[nu][:,:,icum]
                            calc_Mmatrix!(M,dSdUnu,staple_nu_set,nu,ρs,U,ix1,iy1,iz1,it1,tmp_matrices)
                            #M = calc_Mmatrix(dSdUnu,staple_nu_set,nu,ρs,NC,U,ix1,iy1,iz1,it1)
            
                            calc_Λmatrix!(Λ,M,NC)
                            #Λ = calc_Λmatrix(M,NC)
                            #Umu = U[nu][:,:,ix1,iy1,iz1,it1]
                            #Vtmp = -ρs[i]*rightmatrix*Λ*Umu*leftmatrix

                            for i=1:NC
                                for j=1:NC
                                    tmp_matrices[3][j,i] = U[nu][j,i,ix1,iy1,iz1,it1] #U[nu][:,:,ix1,iy1,iz1,it1]
                                end
                            end

                            mul!(tmp_matrices[1],tmp_matrices[3],leftmatrix)
                            mul!(tmp_matrices[2],Λ,tmp_matrices[1])
                            mul!(tmp_matrices[1],rightmatrix,tmp_matrices[2])
                            tmp_matrices[1] .*= -ρs[i]

                            @. V += tmp_matrices[1]



                            #println("Vtmp $i,$nu")
                            #display(Vtmp/ρ)
                            #println("\t")

                            #V += Vtmp
                            #println("Vtmp, ", Vtmp)
                        end
                    end
                end
            end
            #exit()
            
        end

        #println("Theta old = ", V)
        

        y = ix,iy,iz,it
        ix1,iy1,iz1,it1 = y
        #ix1,iy1,iz1,it1 = periodiccheck(ix1,iy1,iz1,it1,NX,NY,NZ,NT,NDW)
        icum = get_1dindex(ix1,iy1,iz1,it1,NX,NY,NZ,NT)
        #icum = (((it1-1)*NZ+iz1-1)*NY+iy1-1)*NX+ix1

        staple_mu_set = smearing.staples_for_stout
        
        
        dSdUnu = dSdU[mu][:,:,icum]
        #M = calc_Mmatrix(dSdUnu,staple_mu_set,mu,ρs,NC,U,y...)
        calc_Mmatrix!(M,dSdUnu,staple_mu_set,mu,ρs,U,y...,tmp_matrices)
        calc_Λmatrix!(Λ,M,NC)
        #Λ = calc_Λmatrix(M,NC)
        

        dSdρs = zero(ρs)
        Umu = U[mu][:,:,ix,iy,iz,it]

        C = zeros(ComplexF64,NC,NC)
        Ctmp = zero(C)
        for i=1:num
            Ctmp .= 0
            staple_mu = staple_mu_set[i][mu]
            #evaluate_wilson_loops!(Ctmp,staple_nu_set[i][nu],U,ix,iy,iz,it,tmp_matrices[4:6])
            evaluate_wilson_loops!(Ctmp,staple_mu,U,y...,tmp_matrices[4:6])
            C[:,:] .+= ρs[i]*Ctmp[:,:]
            dSdρs[i] += 2*real(tr(Umu'*Λ*Ctmp))
            
        end


        
        #Vtmp = -ρ*C'*M
        Vtmp = -C'*Λ


        V += Vtmp
        #println("Vtmp, ", Vtmp)

        Q = calc_Qmatrix(staple_mu_set,mu,ρs,U,NC,ix,iy,iz,it)
        #println("Q ", Q)
        Vtmp = dSdUnu*exp(Q)

        #=
        if y == (2,2,2,2) && mu == 1
            println("U,M,Λ,Q")
            display(Umu)
            println("\t")
            display(M)
            println("\t")
            display(Λ)
            println("\t")
            display(Q)
            println("\t")
        end
        =#

        #println("dSdUnu ",dSdUnu)
        #println("exp(Q) = ",exp(Q))
        V+= Vtmp 
        #println("Vtmp, ", Vtmp)
        

        return V,dSdρs
    end


    function apply_smearing(U::Array{GaugeFields{SU{NC}},1},smearing::T) where {NC,T <: SmearingParam_single}
        Uout = similar(U)
        calc_stout_multi!(Uout,U,smearing.ρs,smearing.staples_for_stout) 
        return Uout
    end


    function apply_smearing(U::Array{GaugeFields{SU{NC}},1},smearing::T) where {NC,T <: SmearingParam_multi}
        Uout_multi = calc_stout_multi(U,smearing.ρs,smearing.staples_for_stout) 
        return Uout_multi
    end

    function calc_stout(U::Array{GaugeFields{SU{NC}},1},ρ) where NC
        Uout = similar(U)
        calc_stout!(Uout,U,ρ)
        return Uout
    end

    function calc_stout!(Uout::Array{GaugeFields{SU{NC}},1},ρ) where NC
        Uin = deepcopy(Uout)
        calc_stout!(Uout,Uin,ρ)
    end

    function calc_stout!(Uout::Array{GaugeFields{SU{NC}},1},U::Array{GaugeFields{SU{NC}},1},ρ) where NC
        @assert Uout != U "input U and output U should not be same!"

        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        Q = zeros(ComplexF64,NC,NC)
        C = zeros(ComplexF64,NC,NC)

        I0N = zeros(ComplexF64,NC,NC)
        for i=1:NC
            I0N[i,i] = 1
        end

        for μ=1:4
            loops = make_plaq_staple_prime(μ)
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            C .= 0
                            evaluate_wilson_loops!(C,loops,U,ix,iy,iz,it)
                            #display(C)
                            #println("\t")

                            Ω = ρ*C[:,:]*U[μ][:,:,ix,iy,iz,it]'
                            Q = (im/2)*(Ω' .- Ω) .- (im/(2NC))*tr(Ω' .- Ω)*I0N
                            Uout[μ][:,:,ix,iy,iz,it] = exp(im*Q)*U[μ][:,:,ix,iy,iz,it]
                            set_wing!(Uout[μ],ix,iy,iz,it)

                            
                        end
                    end
                end
            end
        end
        #display(Uout[1][:,:,1,1,1,1])
        return
    end

    const eps_Q = 1e-18

    function calc_stout_multi(Uin::Array{GaugeFields{SU{NC}},1},ρs::Array{Array{T,1},1},staples)  where {NC,T <: Number}
        numlayer = length(ρs)
        #println("numlayer = ",numlayer,"\t",ρs)
        Uout_multi = Array{Array{GaugeFields{SU{NC}},1}}(undef,numlayer)
        for i=1:numlayer
            Uout_multi[i] = similar(Uin)
        end
        calc_stout_multi!(Uout_multi,Uin,ρs,staples)

        return Uout_multi
    end

    function calc_stout_multi!(Uout_multi::Array{Array{GaugeFields{SU{NC}},1},1},Uin::Array{GaugeFields{SU{NC}},1},ρs::Array{Array{T,1},1},staples)  where {NC,T <: Number}
        numlayer = length(ρs)
        Utmp = similar(Uin)
        #Uout_multi = Array{Array{GaugeFields{SU{NC}},1}}(undef,numlayer)
        U = deepcopy(Uin)
        for i = 1:numlayer
            if i != numlayer
                calc_stout_multi!(Utmp,U,ρs[i],staples)
                set_wing!(Utmp)
                Uout_multi[i] = deepcopy(Utmp)
                Utmp,U = U,Utmp            
            else
                calc_stout_multi!(Uout_multi[i],U,ρs[i],staples)
                set_wing!(Uout_multi[i])
            end
        end
    end

    function calc_stout_multi!(Uout::Array{GaugeFields{SU{NC}},1},Uin::Array{GaugeFields{SU{NC}},1},ρs::Array{Array{T,1},1},staples)  where {NC,T <: Number}
        numlayer = length(ρs)
        Utmp = similar(Uin)
        #Uout_multi = Array{Array{GaugeFields{SU{NC}},1}}(undef,numlayer)
        U = deepcopy(Uin)
        for i = 1:numlayer
            if i != numlayer
                calc_stout_multi!(Utmp,U,ρs[i],staples)
                set_wing!(Utmp)
                Utmp,U = U,Utmp            
            else
                calc_stout_multi!(Uout,U,ρs[i],staples)
                set_wing!(Uout)
            end
        end
        
    end


    function calc_stout_multi!(Uout::Array{GaugeFields{SU{NC}},1},
                U::Array{GaugeFields{SU{NC}},1},ρs::Array{T,1},staples) where {NC,T <: Number}
        @assert Uout != U "input U and output U should not be same!"

        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        Q = zeros(ComplexF64,NC,NC)
        C = zeros(ComplexF64,NC,NC)
        Ctmp = zeros(ComplexF64,NC,NC)
        Ω = zeros(ComplexF64,NC,NC)
        I0N = zeros(ComplexF64,NC,NC)
        for i=1:NC
            I0N[i,i] = 1
        end
        Umu = zeros(ComplexF64,NC,NC)
        Umutmp = zeros(ComplexF64,NC,NC)
        tmp_matrices = Array{Array{ComplexF64,2},1}(undef,3)
        for k=1:3
            tmp_matrices[k] = zero(C)
        end


        num = length(ρs)
        

        for μ=1:4
            #loops = make_plaq_staple_prime(μ)
            #display(loops)
            #exit()
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            Ω .= 0
                            C .= 0
                            for i=1:num
                                Ctmp .= 0
                                evaluate_wilson_loops!(Ctmp,staples[i][μ],U,ix,iy,iz,it,tmp_matrices)
                                #C[:,:] .+= ρs[i]*Ctmp[:,:]
                                for j=1:NC
                                    for k=1:NC
                                        C[k,j] += ρs[i]*Ctmp[k,j]
                                    end
                                end
                                #@. C[:,:] .+= ρs[i]*Ctmp
                            end
                            for i=1:NC
                                for j=1:NC
                                    Umu[j,i] =   U[μ][j,i,ix,iy,iz,it]
                                end
                            end

                            for i=1:NC
                                for j=1:NC
                                    Ω[j,i] = 0
                                    for k=1:NC
                                        Ω[j,i] += C[j,k]*conj(Umu[i,k])
                                    end
                                end
                            end
                            #Ω[:,:] = C[:,:]*U[μ][:,:,ix,iy,iz,it]'

                            #display(Ω)
                            #println("\t")
                            if NC == 1
                                #Q = (im/2)*(Ω' .- Ω) 
                                #Q = (1/2)*(Ω' .- Ω) 
                                Q = -(Ω' .- Ω)
                                #Q = (Ω' .- Ω)
                                #Q = im*(Ω' .- Ω) 
                            else

                                #Q = (im/2)*(Ω' .- Ω) .- (im/(2NC))*tr(Ω' .- Ω)
                                #Q = (1/2)*(Ω' .- Ω) .- (1/(2NC))*tr(Ω' .- Ω)
                                #
                                for i=1:NC
                                    for j=1:NC
                                        Q[j,i] = -(1/2)*(conj(Ω[i,j])-Ω[j,i])
                                    end
                                end
                                trΩ = tr(Ω)
                                for i=1:NC
                                    Q[i,i] += (1/(2NC))*(conj(trΩ)-trΩ)
                                end
                                #Q = -(1/2)*(Ω' - Ω) + (1/(2NC))*tr(Ω' - Ω)*I0N
                                #Q = (1/2)*(Ω' - Ω) - (1/(2NC))*tr(Ω' - Ω)*I0N
                            end
                            #println(Q)
                            trQ2 = 0im
                            for i=1:NC
                                for k=1:NC
                                    trQ2 += Q[i,k]*Q[k,i]
                                end
                            end
                            #trQ2 = tr(Q^2)
                            #Uout[μ][:,:,ix,iy,iz,it] = exp(im*Q)*U[μ][:,:,ix,iy,iz,it]
                            if abs(trQ2) > eps_Q
                                if NC == 2
                                    #exp(Q) = cosq I + sinq/q Q
                                    q = sqrt((-1/2)*trQ2)
                                    for i=1:NC
                                        for j=1:NC
                                            tmp_matrices[1][j,i] = (sin(q)/q)*Q[j,i]
                                        end
                                    end
                                    for i=1:NC
                                        tmp_matrices[1][i,i] += cos(q)
                                    end
                                    mul!(Umutmp,tmp_matrices[1],Umu)
                                    
                                else
                                    mul!(Umutmp,exp(Q),Umu)
                                end
                                for i=1:NC
                                    for j=1:NC
                                        Uout[μ][j,i,ix,iy,iz,it] = Umutmp[j,i]
                                    end
                                end
                                #Uout[μ][:,:,ix,iy,iz,it] = exp(Q)*U[μ][:,:,ix,iy,iz,it]
                            else
                                Uout[μ][:,:,ix,iy,iz,it] = U[μ][:,:,ix,iy,iz,it]
                            end
                            #display(Uout[μ][:,:,ix,iy,iz,it] )
                            #println("\t")
                            set_wing!(Uout[μ],ix,iy,iz,it)

                            
                        end
                    end
                end
            end
        end


        
    end
end