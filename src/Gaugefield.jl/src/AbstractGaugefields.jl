module AbstractGaugefields_module
    using LinearAlgebra
    import ..Wilsonloops_module:Wilson_loop_set,calc_coordinate,make_plaq_staple_prime,calc_shift,make_plaq,make_plaq_staple,
                        Tensor_wilson_lines_set,Tensor_wilson_lines,Tensor_derivative_set,
                        get_leftstartposition,get_rightstartposition,Wilson_loop,calc_loopset_μν_name   
    import Wilsonloop:loops_staple_prime,Wilsonline,get_position,get_direction,Adjoint_GLink,GLink
    

    using MPI
    using InteractiveUtils
    
    abstract type Abstractfields end

    abstract type AbstractGaugefields{NC,Dim} <: Abstractfields
    end


    abstract type Adjoint_fields{T} <: Abstractfields
    end

    struct Adjoint_Gaugefields{T} <: Adjoint_fields{T}
        parent::T
    end

    function Base.adjoint(U::T) where T <: Abstractfields
        Adjoint_Gaugefields{T}(U)
    end

    abstract type Unit_Gaugefield{NC,Dim} <: Abstractfields
    end

    function Base.adjoint(U::Unit_Gaugefield) where T <: Unit_Gaugefield
        return U
    end

    abstract type Shifted_Gaugefields{NC,Dim} <: Abstractfields
    end

    struct Staggered_Gaugefields{T,μ} <: Abstractfields
        parent::T

        function Staggered_Gaugefields(u::T,μ) where T <: Abstractfields
            return new{T,μ}(u)
        end
    end

    function staggered_U(u::T,μ) where T <: Abstractfields
        return Staggered_Gaugefields(u,μ)
    end

    function Base.setindex!(U::T,v...)  where T <: Staggered_Gaugefields
        error("type $(typeof(U)) has no setindex method. This type is read only.")
    end

    function Base.setindex!(U::T,v...)  where T <:  Unit_Gaugefield
        error("type $(typeof(U)) has no setindex method. This type is read only.")
    end

    struct Gaugefield_latticeindices{Dim,NC,T}
        NN::NTuple{Dim,Int64}
        NC::Int8

        function Gaugefield_latticeindices(u::AbstractGaugefields{Dim,NC}) where {Dim,NC}
            _,_,NN... = size(u)
            return new{Dim,NC,typeof{u}}(NN,NC)
        end
    end

    function Gaugefield_latticeindices(U::Array{<: AbstractGaugefields{Dim,NC},1}) where {Dim,NC}
        return Gaugefield_latticeindices(U[1])
    end

    function Base.iterate(g::Gaugefield_latticeindices{Dim,NC,T}) where {Dim,NC,T}
        N<1 && return nothing
        return (1, 2)
    end

    include("./4D/gaugefields_4D.jl")
    include("TA_Gaugefields.jl")

    function LinearAlgebra.mul!(C,A::T,B) where T <: Unit_Gaugefield
        substitute_U!(C,B)
    end

    function LinearAlgebra.mul!(C,A,B::T) where T <: Unit_Gaugefield
        substitute_U!(C,A)
    end   
    

    function Staggered_Gaugefields(u::AbstractGaugefields{NC,Dim}) where {NC,Dim}
        if Dim == 4
            return Staggered_Gaugefields_4D(u)
        else
            error("Dim = $Dim is not supported")
        end
    end
    
    
    #include("./gaugefields_4D_wing_mpi.jl")


    function Base.similar(U::T) where T <: AbstractGaugefields
        error("similar! is not implemented in type $(typeof(U)) ")
    end

    function substitute_U!(a::Array{T1,1},b::Array{T2,1}) where {T1 <: AbstractGaugefields,T2 <: AbstractGaugefields}
        error("substitute_U! is not implemented in type $(typeof(a)) and $(typeof(b))")
    end

    function substitute_U!(a::T1,b::T2) where {T1 <: AbstractGaugefields,T2 <: AbstractGaugefields}
        error("substitute_U! is not implemented in type $(typeof(a)) and $(typeof(b))")
        return 
    end

    function substitute_U!(a::Array{T1,1},b::Array{T2,1},iseven::Bool) where {T1 <: AbstractGaugefields,T2 <: AbstractGaugefields}
        error("substitute_U! is not implemented in type $(typeof(a)) and $(typeof(b))")
    end

    function substitute_U!(a::T1,b::T2,iseven::Bool) where {T1 <: AbstractGaugefields,T2 <: Abstractfields}
        error("substitute_U! is not implemented in type $(typeof(a)) and $(typeof(b))")
        return 
    end


    function RandomGauges(NC,NDW,NN...;mpi = false,PEs=nothing,mpiinit = nothing)
        dim = length(NN)
        if mpi
            if PEs == nothing || mpiinit == nothing
                error("not implemented yet!")
            else
                if dim==4
                    U = randomGaugefields_4D_wing_mpi(NC,NN[1],NN[2],NN[3],NN[4],NDW,PEs,mpiinit = mpiinit)
                else
                    error("not implemented yet!")
                end
                
            end
        else
            if dim == 4
                U = randomGaugefields_4D_wing(NC,NN[1],NN[2],NN[3],NN[4],NDW)
            else
                error("not implemented yet!")
            end
        end
        return U
    end

    function IdentityGauges(NC,NDW,NN...;mpi = false,PEs=nothing,mpiinit = nothing)
        dim = length(NN)
        if mpi
            if PEs == nothing || mpiinit == nothing
                error("not implemented yet!")
            else
                if dim==4
                    U = identityGaugefields_4D_wing_mpi(NC,NN[1],NN[2],NN[3],NN[4],NDW,PEs,mpiinit = mpiinit)
                else
                    error("not implemented yet!")
                end
                
            end
        else
            if dim == 4
                U = identityGaugefields_4D_wing(NC,NN[1],NN[2],NN[3],NN[4],NDW)
            else
                error("not implemented yet!")
            end
        end
        return U
    end

    function Oneinstanton(NC,NDW,NN...;mpi = false,PEs=nothing,mpiinit = nothing)
        dim = length(NN)
        if mpi
            if PEs == nothing || mpiinit == nothing
                error("not implemented yet!")
            else
                if dim==4
                    U = Oneinstanton_4D_wing_mpi(NC,NN[1],NN[2],NN[3],NN[4],NDW,PEs,mpiinit = mpiinit)
                else
                    error("not implemented yet!")
                end
                
            end
        else
            if dim == 4
                U = Oneinstanton_4D_wing(NC,NN[1],NN[2],NN[3],NN[4],NDW)
            else
                error("not implemented yet!")
            end
        end
        return U
    end

    function construct_gauges(NC,NDW,NN...;mpi = false,PEs=nothing,mpiinit = nothing)
        dim = length(NN)
        if mpi
            if PEs == nothing || mpiinit == nothing
                error("not implemented yet!")
            else
                if dim==4
                    U = identityGaugefields_4D_wing_mpi(NC,NN[1],NN[2],NN[3],NN[4],NDW,PEs,mpiinit = mpiinit)
                else
                    error("not implemented yet!")
                end
                
            end
        else
            if dim == 4
                U = identityGaugefields_4D_wing(NC,NN[1],NN[2],NN[3],NN[4],NDW)
            else
                error("not implemented yet!")
            end
        end
        return U
    end



    function clear_U!(U::T) where T <: AbstractGaugefields
        error("clear_U! is not implemented in type $(typeof(U)) ")
    end

    function clear_U!(U::Array{<: AbstractGaugefields{NC,Dim},1}) where {NC,Dim}
        for μ=1:Dim
            clear_U!(U[μ])
        end
    end

    function clear_U!(U::T,iseven::Bool) where T <: AbstractGaugefields
        error("clear_U! is not implemented in type $(typeof(U)) ")
    end

    function clear_U!(U::Array{<: AbstractGaugefields{NC,Dim},1},iseven::Bool) where {NC,Dim}
        for μ=1:Dim
            clear_U!(U[μ],iseven)
        end
    end

    function unit_U!(U::T) where T <: AbstractGaugefields
        error("unit_U! is not implemented in type $(typeof(U)) ")
    end

    function unit_U!(U::Array{<: AbstractGaugefields{NC,Dim},1}) where {NC,Dim}
        for μ=1:Dim
            unit_U!(U[μ])
        end
    end

    
    function shift_U(U::AbstractGaugefields{NC,Dim},ν) where {NC,Dim}
        error("shift_U is not implemented in type $(typeof(U)) ")
        return nothing
    end

    function map_U!(U::AbstractGaugefields{NC,Dim},f::Function,V::AbstractGaugefields{NC,Dim},iseven::Bool) where {NC,Dim, T<: Abstractfields}
        error("map_U! is not implemented in type $(typeof(U)) ")
        return nothing
    end
    

    function identitymatrix(U::T) where T <: AbstractGaugefields
        error("identitymatrix is not implemented in type $(typeof(U)) ")
    end

    function set_wing_U!(U::Array{<: AbstractGaugefields{NC,Dim},1}) where {NC,Dim}
        for μ=1:Dim
            set_wing_U!(U[μ])
        end
    end

    function set_wing_U!(U::T) where T <: AbstractGaugefields
        error("set_wing_U! is not implemented in type $(typeof(U)) ")
    end

    function evaluate_gaugelinks_evenodd!(uout::T,w::Wilsonline{Dim},U::Array{T,1},temps::Array{T,1},iseven) where {T<: AbstractGaugefields,Dim}

        #Uold = temps[1]
        Unew = temps[1]
        #Utemp2 = temps[2]
        #clear_U!(uout)

        glinks = w
        numlinks = length(glinks)
        if numlinks == 0
            unit_U!(uout)
            return
        end

        j = 1    
        U1link = glinks[1]
        direction = get_direction(U1link)
        position = get_position(U1link)
        isU1dag = ifelse(typeof(U1link) <: Adjoint_GLink,true,false)


        if numlinks == 1
            substitute_U!(Unew,U[direction])
            Ushift1 = shift_U(Unew,position)
            if isU1dag 
                #println("Ushift1 ",Ushift1'[1,1,1,1,1,1])
                substitute_U!(uout,Ushift1',iseven)
            else
                substitute_U!(uout,Ushift1,iseven)
            end
            return
        end

        substitute_U!(Unew,U[direction])
        Ushift1 = shift_U(Unew,position)

        for j=2:numlinks
            Ujlink = glinks[j]
            isUkdag = ifelse(typeof(Ujlink) <: Adjoint_GLink,true,false)
            position = get_position(Ujlink)
            direction = get_direction(Ujlink)

            Ushift2 = shift_U(U[direction],position)
            multiply_12!(uout,Ushift1,Ushift2,j,isUkdag,isU1dag,iseven)


            substitute_U!(Unew,uout)

            Ushift1 = shift_U(Unew,(0,0,0,0))
        end


    end

    function evaluate_gaugelinks!(uout::T,w::Wilsonline{Dim},U::Array{T,1},temps::Array{T,1}) where {T<: AbstractGaugefields,Dim}

        #Uold = temps[1]
        Unew = temps[1]
        #Utemp2 = temps[2]
        #clear_U!(uout)

        glinks = w
        numlinks = length(glinks)
        if numlinks == 0
            unit_U!(uout)
            return
        end

        j = 1    
        U1link = glinks[1]
        direction = get_direction(U1link)
        position = get_position(U1link)
        isU1dag = ifelse(typeof(U1link) <: Adjoint_GLink,true,false)

        #show(glinks)   
        #println("in evaluate_gaugelinks!")
        #show(w)
        #println("numlinks = $numlinks")
        if numlinks == 1
            substitute_U!(Unew,U[direction])
            Ushift1 = shift_U(Unew,position)
            if isU1dag 
                #println("Ushift1 ",Ushift1'[1,1,1,1,1,1])
                substitute_U!(uout,Ushift1')
            else
                substitute_U!(uout,Ushift1)
            end
            return
        end

        #j = 1    
        #U1link = glinks[1]
        #direction = get_direction(U1link)
        #position = get_position(U1link)
        #println("i = $i j = $j position = $position")
        substitute_U!(Unew,U[direction])
        Ushift1 = shift_U(Unew,position)

        #ix,iy,iz,it=(2,2,2,2)
        #println("posotion = $position")
        #pos = Tuple([ix,iy,iz,it] .+ collect(position))
        #U1 = Unew[:,:,pos...]
        #println("U1, ",Unew[:,:,pos...])
        #isU1dag = ifelse(typeof(U1link) <: Adjoint_GLink,true,false)

        


        for j=2:numlinks
            Ujlink = glinks[j]
            isUkdag = ifelse(typeof(Ujlink) <: Adjoint_GLink,true,false)
            position = get_position(Ujlink)
            direction = get_direction(Ujlink)
            #println("j = $j position = $position")
            #println("a,b, $isUkdag , $isU1dag")
            Ushift2 = shift_U(U[direction],position)
            multiply_12!(uout,Ushift1,Ushift2,j,isUkdag,isU1dag)

            #pos = Tuple([ix,iy,iz,it] .+ collect(position))
            #U2 = U[direction][:,:,pos...]
            #println("U1U2dag ", U1*U2')
            substitute_U!(Unew,uout)
            
            #println("Unew ", Unew[:,:,ix,iy,iz,it])

            Ushift1 = shift_U(Unew,(0,0,0,0))
            #println("uout ", uout[:,:,ix,iy,iz,it])
        end


        #println("uout2 ", uout[:,:,ix,iy,iz,it])
        
        
    end

    function evaluate_gaugelinks_evenodd!(xout::T,w::Array{<:Wilsonline{Dim},1},U::Array{T,1},temps::Array{T,1},iseven) where {T<: AbstractGaugefields,Dim}
        num = length(w)
        temp1 = temps[1]
        temp2 = temps[2]

        #ix,iy,iz,it=(2,2,2,2)
        
        clear_U!(xout,iseven)
        for i=1:num
            glinks = w[i]
            evaluate_gaugelinks_evenodd!(temp2,glinks,U,[temp1],iseven)
            #println("uout2 ", temp2[:,:,ix,iy,iz,it])
            add_U!(xout,temp2,iseven)
            #println("xout ", xout[:,:,ix,iy,iz,it])
        end

        #println("xout2 ", xout[:,:,ix,iy,iz,it])
        return


    end

    function evaluate_gaugelinks!(xout::T,w::Array{<:Wilsonline{Dim},1},U::Array{T,1},temps::Array{T,1}) where {T<: AbstractGaugefields,Dim}
        num = length(w)
        temp1 = temps[1]
        temp2 = temps[2]

        #ix,iy,iz,it=(2,2,2,2)
        
        clear_U!(xout)
        for i=1:num
            glinks = w[i]
            evaluate_gaugelinks!(temp2,glinks,U,[temp1])
            #println("uout2 ", temp2[:,:,ix,iy,iz,it])
            add_U!(xout,temp2)
            #println("xout ", xout[:,:,ix,iy,iz,it])
        end

        #println("xout2 ", xout[:,:,ix,iy,iz,it])
        return


    end

    #=
    function evaluate_gaugelinks_inside!(U,glinks,Ushift1,Uold,Unew,numlinks)
        for k=2:numlinks
            position =get_position(glinks[k])
            direction = get_direction(glinks[k])
            Ushift2 = shift_U(U[direction],position)
            multiply_12!(Unew,Ushift1,Ushift2,k,loopk,loopk1_2)

            Unew,Uold = Uold,Unew
            Ushift1 = shift_U(Uold,(0,0,0,0))
        end
    end
    =#
    
    function multiply_12!(temp3,temp1,temp2,k,isUkdag::Bool,isU1dag::Bool)
        if k==2
            if isUkdag
                if isU1dag
                    mul!(temp3,temp1',temp2')
                else
                    mul!(temp3,temp1,temp2')
                end
            else
                if isU1dag
                    mul!(temp3,temp1',temp2)
                else
                    mul!(temp3,temp1,temp2)
                end
            end
        else
            if isUkdag
                mul!(temp3,temp1,temp2')
            else
                mul!(temp3,temp1,temp2)
            end
        end
        return
    end

        
    function multiply_12!(temp3,temp1,temp2,k,isUkdag::Bool,isU1dag::Bool,iseven)
        if k==2
            if isUkdag
                if isU1dag
                    mul!(temp3,temp1',temp2',iseven)
                else
                    mul!(temp3,temp1,temp2',iseven)
                end
            else
                if isU1dag
                    mul!(temp3,temp1',temp2,iseven)
                else
                    mul!(temp3,temp1,temp2,iseven)
                end
            end
        else
            if isUkdag
                mul!(temp3,temp1,temp2',iseven)
            else
                mul!(temp3,temp1,temp2,iseven)
            end
        end
        return
    end




    function evaluate_wilson_loops!(xout::T,w::Wilson_loop_set,U::Array{T,1},temps::Array{T,1}) where T<: AbstractGaugefields
        num = length(w)
        clear_U!(xout)
        Uold = temps[1]
        Unew = temps[2]

        for i=1:num
            wi = w[i]
            numloops = length(wi)    

            shifts = calc_shift(wi)
            #println("shift ",shifts)
            
            loopk = wi[1]
            k = 1
            #println("k = $k shift: ",shifts[k])
            substitute_U!(Uold,U[loopk[1]])
            Ushift1 = shift_U(Uold,shifts[1])

            #gauge_shift_all!(temp1,shifts[1],U[loopk[1]])

            
            loopk1_2 = loopk[2]
            evaluate_wilson_loops_inside!(U,shifts,wi,Ushift1,Uold,Unew,numloops,loopk,loopk1_2)
            
            #=
            for k=2:numloops
                
                loopk = wi[k]
                #println("k = $k shift: ",shifts[k])
                #println("gauge_shift!(temp2,$(shifts[k]),$(loopk[1]) )")
                #clear!(temp2)
                Ushift2 = shift_U(U[loopk[1]],shifts[k])
                #gauge_shift_all!(temp2,shifts[k],U[loopk[1]])

                #multiply_12!(temp3,temp1,temp2,k,loopk,loopk1_2)
                multiply_12!(Unew,Ushift1,Ushift2,k,loopk,loopk1_2)
                Unew,Uold = Uold,Unew
                #Ushift1 = shift_U(Uold,(0,0,0,0))
                Ushift1 = Uold
                #temp1,temp3 = temp3,temp1
            end
            =#
            add_U!(xout,Uold)
            #println("i = $i ",Uold[1,1,1,1,1,1])
            #add_U!(xout,Ushift1)
            #add!(xout,temp1)
            
        end
    end

    function evaluate_wilson_loops_inside!(U,shifts,wi,Ushift1,Uold,Unew,numloops,loopk,loopk1_2)
        for k=2:numloops
                
            loopk = wi[k]
            #println("k = $k shift: ",shifts[k])
            #println("gauge_shift!(temp2,$(shifts[k]),$(loopk[1]) )")
            #clear!(temp2)
            Ushift2 = shift_U(U[loopk[1]],shifts[k])
            #gauge_shift_all!(temp2,shifts[k],U[loopk[1]])

            #multiply_12!(temp3,temp1,temp2,k,loopk,loopk1_2)
            multiply_12!(Unew,Ushift1,Ushift2,k,loopk,loopk1_2)
            
            Unew,Uold = Uold,Unew
            Ushift1 = shift_U(Uold,(0,0,0,0))
            
            #Ushift1 = Uold
            #temp1,temp3 = temp3,temp1
        end
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


    function calculate_Plaquette(U::Array{T,1}) where T <: AbstractGaugefields
        error("calculate_Plaquette is not implemented in type $(typeof(U)) ")
    end

    function calculate_Plaquette(U::Array{T,1},temps::Array{T1,1}) where {T <: AbstractGaugefields,T1 <: AbstractGaugefields}
        return calculate_Plaquette(U,temps[1],temps[2])
    end

    function calculate_Plaquette(U::Array{T,1},temp::AbstractGaugefields{NC,Dim},staple::AbstractGaugefields{NC,Dim}) where {NC,Dim,T <: AbstractGaugefields}
        plaq = 0
        V = staple
        for μ=1:Dim
            construct_staple!(V,U,μ,temp)
            mul!(temp,U[μ],V')
            plaq += tr(temp)
            
        end



        return real(plaq*0.5)
    end

    function construct_staple!(staple::AbstractGaugefields,U,μ) where T <: AbstractGaugefields
        error("construct_staple! is not implemented in type $(typeof(U)) ")
    end

    function add_force!(F::Array{T1,1},U::Array{T2,1},temps::Array{<: AbstractGaugefields{NC,Dim},1};
        plaqonly = false,staplefactors::Union{Array{<: Number,1},Nothing} = nothing,factor = 1) where {NC,Dim,T1 <: AbstractGaugefields,T2 <: AbstractGaugefields}
        error("add_force! is not implemented in type $(typeof(F)) ")
    end
    
    function add_force!(F::Array{T1,1},U::Array{T2,1},temps::Array{<: AbstractGaugefields{NC,Dim},1};
                plaqonly = false,staplefactors::Union{Array{<: Number,1},Nothing} = nothing,factor = 1) where {NC,Dim,T1 <: TA_Gaugefields,T2 <: AbstractGaugefields}
        @assert length(temps) >= 3 "length(temps) should be >= 3. But $(length(temps))"
        #println("add force, plaqonly = $plaqonly")
        
        V = temps[3]  
        temp1 = temps[1]
        temp2 = temps[2]    

        for μ=1:Dim
            if plaqonly

                construct_double_staple!(V,U,μ,temps[1:2])

                mul!(temp1,U[μ],V') #U U*V
            else
                clear_U!(V)
                for i=1:gparam.numactions
                    loops = gparam.staples[i][μ]
                    evaluate_wilson_loops!(temp3,loops,U,[temp1,temp2])
                    add_U!(V,staplefactors[i],temp3)
                    #add_U!(V,gparam.βs[i]/gparam.β,temp3)
                end
                mul!(temp1,U[μ],V) #U U*V
            end

            Traceless_antihermitian_add!(F[μ],factor,temp1)
            #add_U!(F[μ],factor,temp2)
        end

    end

    #=
    function add_force!(F::Array{T,1},U::Array{T,1},temps::Array{<: AbstractGaugefields{NC,Dim},1},factor = 1) where {NC,Dim,T <: AbstractGaugefields,GP}
        @assert length(temps) >= 3 "length(temps) should be >= 3. But $(length(temps))"
        clear_U!(F)
        V = temps[3]  
        temp1 = temps[1]
        temp2 = temps[2]    

        for μ=1:Dim
            construct_double_staple!(V,U,μ,temps[1:2])

            mul!(temp1,U[μ],V') #U U*V

            a = temp1[:,:,1,1,1,1]
            println(a'*a)

            Traceless_antihermitian!(temp2,temp1)
            #println(temp2[1,1,1,1,1,1])
            a = temp2[:,:,1,1,1,1]
            println(a'*a)
            error("a")
            add_U!(F[μ],factor,temp2)
        end
    end
    =#
    

    function exptU!(uout::T,t::N,f::T1,temps::Array{T,1}) where {N <: Number, T <: AbstractGaugefields,T1 <: AbstractGaugefields} #uout = exp(t*u)
        error("expUt! is not implemented in type $(typeof(f)) ")
    end

    function exptU!(uout::T,f::T1,temps::Array{T,1}) where {T <: AbstractGaugefields,T1 <: AbstractGaugefields}
        expU!(uout,1,f,temps)
    end

    function exp_aF_U!(W::Array{<: AbstractGaugefields{NC,Dim},1},a::N,F::Array{T1,1},U::Array{T,1},temps::Array{T,1}) where {NC,Dim,N <: Number, T <: AbstractGaugefields,T1 <: AbstractGaugefields} #exp(a*F)*U
        @assert a != 0 "Δτ should not be zero in expF_U! function!"
        expU = temps[1]
        temp1 = temps[2]
        temp2 = temps[3]
        #clear_U!(temp1)
        #clear_U!(temp2)
        #clear_U!(expU)
        #clear_U!(W)

        for μ=1:Dim
            exptU!(expU,a,F[μ],[temp1,temp2])
            mul!(W[μ],expU,U[μ])
        end

        set_wing_U!(W)
    end


    function staple_prime()
        loops_staple_prime = Array{Wilson_loop_set,2}(undef,4,4)
        for Dim=1:4
            for μ=1:Dim
                loops_staple_prime[Dim,μ] = make_plaq_staple_prime(μ,Dim)
            end
        end
        return loops_staple_prime
    end
    const loops_staple_prime_old = staple_prime()


    function construct_double_staple!(staple::AbstractGaugefields{NC,Dim},U::Array{T,1},μ,temps::Array{<: AbstractGaugefields{NC,Dim},1}) where {NC,Dim,T <: AbstractGaugefields}
        #println("mu = ",μ)
        #loops = loops_staple_prime_old[Dim,μ] #make_plaq_staple_prime(μ,Dim)
        #println("staple")
        #@time evaluate_wilson_loops!(staple,loops,U,temps)
        #println(staple[1,1,1,1,1,1])
        loops = loops_staple_prime[(Dim,μ)]
        evaluate_gaugelinks!(staple,loops,U,temps)
        #println(staple[1,1,1,1,1,1])
        #error("construct!!")
    end


    function construct_staple!(staple::AbstractGaugefields{NC,Dim},U::Array{T,1},μ,temp::AbstractGaugefields{NC,Dim}) where {NC,Dim,T <: AbstractGaugefields}
        U1U2 = temp
        firstterm = true

        for ν=1:Dim
            if ν == μ
                continue
            end
            
            #=
                    x+nu temp2
                    .---------.
                    I         I
            temp1 I         I
                    I         I
                    .         .
                    x        x+mu
            =#
            U1 = U[ν]    
            U2 = shift_U(U[μ],ν)
            #println(typeof(U1))
            mul!(U1U2,U1,U2)
            #error("test")
            
            U3 = shift_U(U[ν],μ)
            #mul!(staple,temp,Uμ')
            #  mul!(C, A, B, α, β) -> C, A B α + C β
            if firstterm
                β = 0
                firstterm = false
            else
                β = 1
            end
            mul!(staple,U1U2,U3',1,β) #C = alpha*A*B + beta*C

            #println("staple ",staple[1,1,1,1,1,1])
            

            #mul!(staple,U0,Uν,Uμ')
        end
        set_wing_U!(staple)
    end

    function Base.size(U::T) where T <: AbstractGaugefields
        error("Base.size is not implemented in type $(typeof(U)) ")
    end

    function calculate_Polyakov_loop(U::Array{T,1},temp1::AbstractGaugefields{NC,Dim},temp2::AbstractGaugefields{NC,Dim}) where {NC,Dim,T <: AbstractGaugefields}
        Uold = temp1
        Unew = temp2
        shift = zeros(Int64,Dim)
        
        μ = Dim
        _,_,NN... = size(U[1]) #NC,NC,NX,NY,NZ,NT 4D case
        lastaxis = NN[end]
        #println(lastaxis)

        substitute_U!(Uold,U[μ])
        for i=2:lastaxis
            shift[μ] = i-1
            U1 = shift_U(U[μ],Tuple(shift))
            mul_skiplastindex!(Unew,Uold,U1)
            Uold,Unew = Unew,Uold
        end

        set_wing_U!(Uold)
        poly = partial_tr(Uold,μ)/prod(NN[1:Dim-1])
        return poly

    end





    function mul_skiplastindex!(c::T,a::T1,b::T2) where {T <: AbstractGaugefields, T1 <: Abstractfields,T2 <: Abstractfields}
        error("mul_skiplastindex! is not implemented in type $(typeof(c)) ")
    end


    function Traceless_antihermitian(vin::T) where T <: AbstractGaugefields
        #error("Traceless_antihermitian is not implemented in type $(typeof(vin)) ")
        vout = deepcopy(vin)
        Traceless_antihermitian!(vout,vin)
        return vout
    end

    function Traceless_antihermitian_add!(U::T,factor,temp1) where T <: AbstractGaugefields
        error("Traceless_antihermitian_add! is not implemented in type $(typeof(U)) ")
    end

    function Traceless_antihermitian!(vout::T,vin::T) where T <: AbstractGaugefields
        error("Traceless_antihermitian! is not implemented in type $(typeof(vout)) ")
    end

    function add_U!(c::T,a::T1) where {T<: AbstractGaugefields,T1 <: Abstractfields}
        error("add_U! is not implemented in type $(typeof(c)) ")
    end

    function add_U!(c::Array{<: AbstractGaugefields{NC,Dim},1},α::N,a::Array{T1,1}) where {NC,Dim,T1 <: Abstractfields, N<:Number}
        for μ=1:Dim
            add_U!(c[μ],α,a[μ])
        end
    end

    function add_U!(c::T,a::T1,iseven) where {T<: AbstractGaugefields,T1 <: Abstractfields}
        error("add_U! is not implemented in type $(typeof(c)) ")
    end

    function add_U!(c::Array{<: AbstractGaugefields{NC,Dim},1},α::N,a::Array{T1,1},iseven) where {NC,Dim,T1 <: Abstractfields, N<:Number}
        for μ=1:Dim
            add_U!(c[μ],α,a[μ],iseven)
        end
    end

    function add_U!(c::T,α::N,a::T1) where {T<: AbstractGaugefields,T1 <: Abstractfields, N<:Number}
        error("add_U! is not implemented in type $(typeof(c)) ")
    end

    function LinearAlgebra.mul!(c::T,a::T1,b::T2,α::Ta,β::Tb) where {T<: AbstractGaugefields,T1 <: Abstractfields,T2 <: Abstractfields,Ta <: Number, Tb <: Number}
        error("LinearAlgebra.mul! is not implemented in type $(typeof(c)) ")
    end

    function LinearAlgebra.mul!(c::T,a::N,b::T2) where {T<: AbstractGaugefields,N <: Number ,T2 <: Abstractfields}
        error("LinearAlgebra.mul! is not implemented in type $(typeof(c)) ")
    end

    function partial_tr(a::T,μ) where T<: Abstractfields
        error("partial_tr is not implemented in type $(typeof(a)) ")
    end

    function LinearAlgebra.tr(a::T) where T<: Abstractfields
        error("LinearAlgebra.tr! is not implemented in type $(typeof(a)) ")
    end

    """
        Tr(A*B)
    """
    function LinearAlgebra.tr(a::T,b::T) where T<: Abstractfields
        error("LinearAlgebra.tr! is not implemented in type $(typeof(a)) ")
    end

    """
    M = (U*δ_prev) star (dexp(Q)/dQ)
    Λ = TA(M)
    """
    function construct_Λmatrix_forSTOUT!(Λ,δ_prev::T,Q,u::T) where T <: AbstractGaugefields
        error("construct_Λmatrix_forSTOUT! is not implemented in type $(typeof(u)) ")
    end

    const eps_Q = 1e-18

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
    

    function calc_Mmatrix!(Mn,δn_prev,Qn,Un,u::AbstractGaugefields{2,Dim},tempmatrices) where {Dim}
        Unδn = tempmatrices[1]
        B = tempmatrices[2]
        tmp_matrix1 = tempmatrices[3]

        trQ2 = 0.0
        for i=1:2
            for j=1:2
                trQ2 += Qn[i,j]*Qn[j,i]
            end
        end

        if abs(trQ2) > eps_Q
            q = sqrt((-1/2)*trQ2)
            calc_Bmatrix!(B,q,Qn,NC)
            for i=1:2
                for j=1:2
                    tmp_matrix1[j,i] = Un[j,i]
                end
            end

            mul!(Unδn,tmp_matrix1,δn_prev)
            trsum = 0.0im
            for i=1:2
                for j=1:2
                    trsum += Unδn[i,j]*B[j,i]
                end
            end
            for i=1:2
                for j=1:2
                    Mn[j,i] = (sin(q)/q)*Unδn[j,i] + trsum*Qn[j,i]
                end
            end
        end
    end

    function calc_Mmatrix!(Mn,δn_prev,Qn,Un,u::AbstractGaugefields{3,Dim},tempmatrices) where {Dim}
        Unδn = tempmatrices[1]
        tmp_matrix1 = tempmatrices[2]
        tmp_matrix2 = tempmatrices[3]
        #println("Qn ", Qn)
        trQ2 = 0.0
        for i=1:3
            for j=1:3
                trQ2 += Qn[i,j]*Qn[j,i]
            end
        end
        #println("tr", trQ2)

        if abs(trQ2) > eps_Q
            Qn ./= im
            #println("Qn b ",Qn)
            f0,f1,f2,b10,b11,b12,b20,b21,b22 = calc_coefficients_Q(Qn)
            for i=1:3
                for j=1:3
                    tmp_matrix1[j,i] = Un[j,i]
                end
            end
            mul!(Unδn,tmp_matrix1,δn_prev)
            
            B1 = tmp_matrix1
            B1 .= 0
            B2 = tmp_matrix2
            B2 .= 0
            for i=1:3
                B1[i,i] = b10
                B2[i,i] = b20
            end
            for j=1:3
                for i=1:3
                    B1[i,j] += b11*Qn[i,j]
                    B2[i,j] += b21*Qn[i,j]
                    for k=1:3  
                        B1[i,j] += b12*Qn[i,k]*Qn[k,j]
                        B2[i,j] += b22*Qn[i,k]*Qn[k,j]
                    end
                end
            end
            #println("coeff, ",(f0,f1,f2,b10,b11,b12,b20,b21,b22))
            #println("B1 ",B1)
            #println("B2 ",B2)

            trB1 = 0.0
            trB2 = 0.0
            for i=1:3
                for j=1:3
                    trB1 += Unδn[i,j]*B1[j,i]
                    trB2 += Unδn[i,j]*B2[j,i]
                end
            end

            for j=1:3
                for i=1:3
                    Mn[i,j] = trB1*Qn[i,j] + f1*Unδn[i,j]
                    for k=1:3  
                        Mn[i,j] += trB2*Qn[i,k]*Qn[k,j]+f2*(Qn[i,k]*Unδn[k,j]+Unδn[i,k]*Qn[k,j])
                    end
                end
            end
            Mn ./= im
        end
    end

    function calc_Mmatrix!(Mn,δn_prev,Qn,Un,u::AbstractGaugefields{NC,Dim},tempmatrices) where {NC,Dim}
        error("not supported yet")

        @assert NC > 3 "NC > 3 not NC = $NC"
        Unδn = tempmatrices[1]
        B = tempmatrices[2]
        tempmatrix = tempmatrices[3]



        trQ2 = 0.0
        for i=1:NC
            for j=1:NC
                trQ2 += Qn[i,j]*Qn[j,i]
            end
        end

        if abs(trQ2) > eps_Q
            e,v = eigen(Qn) 
            mul!(Unδn,Un,δn_prev)
            #=
                    A star dexp(Q)/dQ = \sum_{n=0}^{infty} \frac{1}{n!} 
                                            \sum_{k=0}^{n-1} i^{n-1-k}P^+ D^{n-1-k} P A P^+ D^k P i^k
                                    = P^+ (\sum_{n=0}^{infty} (i^{n-1}/n!) sum_{k=0}^{n-1} D^{n-1-k} B D^k)  P
                                    B = P A P+
                =#

            mul!(tempmatrix,Unδn,v)
            mul!(B,v',tempmatrix)

            


        end

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

    function calc_coefficients_Q(Q)
        @assert size(Q) == (3,3)
        c0 = Q[1,1]*Q[2,2]*Q[3,3]+Q[1,2]*Q[2,3]*Q[3,1]+Q[1,3]*Q[2,1]*Q[3,2]-Q[1,3]*Q[2,2]*Q[3,1]-Q[1,2]*Q[2,1]*Q[3,3]-Q[1,1]*Q[2,3]*Q[3,2]
        #@time cdet = det(Q)
        ##println(c0,"\t",cdet)
        #exit() 
        
        c1 = 0.0
        for i=1:3
            for j=1:3
                c1 += Q[i,j]*Q[j,i]
            end
        end
        c1 /= 2
        c0max = 2*(c1/3)^(3/2)
        θ = acos(c0/c0max)
        u = sqrt(c1/3)*cos(θ/3)
        w = sqrt(c1)*sin(θ/3)
        ξ0 = sin(w)/w
        ξ1 = cos(w)/w^2 - sin(w)/w^3
    
        emiu = exp(-im*u)
        e2iu = exp(2*im*u)
    
        h0 = (u^2-w^2)*e2iu + emiu*(
            8u^2*cos(w)+2*im*u*(3u^2+w^2)* ξ0
        )
        h1 = 2u*e2iu-emiu*(
            2u*cos(w)-im*(3u^2-w^2)* ξ0
        )
        h2 = e2iu - emiu*(cos(w)+3*im*u*ξ0)
    
        denom = 9u^2-w^2
        
        f0 = h0/denom
        f1 = h1/denom
        f2 = h2/denom
    
        r10 = 2*(u+im*(u^2-w^2))*e2iu + 
                2*emiu*(
                    4u*(2-im*u)*cos(w) + 
                    im*(9u^2+w^2-im*u*(3u^2+w^2))*ξ0
                )
        r11 = 2*(1+2*im*u)*e2iu+ 
                emiu*(
                    -2*(1-im*u)*cos(w)+
                    im*(6u+im*(w^2-3u^2))*ξ0
                )
        r12 = 2*im*e2iu + im*emiu*(
            cos(w) -3*(1-im*u)*ξ0
        )
        r20 = -2*e2iu+2*im*u*emiu*(
            cos(w)+(1+4*im*u)*ξ0+3u^2*ξ1
        )
        r21 = -im*emiu*(
            cos(w)+(1+2*im*u)*ξ0 - 
            3*u^2*ξ1
        )
        r22 = emiu*(
            ξ0-3*im*u*ξ1
        )
        b10 = (
            2*u*r10+(3u^2-w^2)*r20 - 2*(15u^2+w^2)*f0
            )/(
                2*denom^2
            )
        
        b11 = (
                2*u*r11+(3u^2-w^2)*r21 - 2*(15u^2+w^2)*f1
                )/(
                    2*denom^2
                )
        b12 = (
            2*u*r12+(3u^2-w^2)*r22 - 2*(15u^2+w^2)*f2
            )/(
                2*denom^2
            )
        b20 = (
            r10 - 3*u*r20 - 24*u*f0
            )/(2*denom^2)
        b21 = (
                r11 - 3*u*r21 - 24*u*f1
                )/(2*denom^2)
        b22 = (
            r12 - 3*u*r22 - 24*u*f2
            )/(2*denom^2)
    
        return f0,f1,f2,b10,b11,b12,b20,b21,b22
    end
    

    function staggered_phase(μ,iii...)
        error("staggered_phase is not implemented")
    end



end