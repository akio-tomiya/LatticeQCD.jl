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



    include("./4D/gaugefields_4D.jl")
    include("TA_Gaugefields.jl")

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

    function substitute_U!(a::Array{T,1},b::Array{T,1}) where T <: AbstractGaugefields
        error("substitute_U! is not implemented in type $(typeof(a)) ")
    end

    function substitute_U!(a::T,b::T) where T <: AbstractGaugefields
        error("substitute_U! is not implemented in type $(typeof(a)) ")
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

    
    function shift_U(U::AbstractGaugefields{NC,Dim},ν) where {NC,Dim}
        error("shift_U is not implemented in type $(typeof(U)) ")
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


    function evaluate_gaugelinks!(xout::T,w::Array{<:Wilsonline{Dim}},U::Array{T,1},temps::Array{T,1}) where {T<: AbstractGaugefields,Dim}
        num = length(w)
        clear_U!(xout)
        Uold = temps[1]
        Unew = temps[2]
        for i=1:num
            glinks = w[i]
            numlinks = length(glinks)
            #show(glinks)        
            j = 1    
            U1link = glinks[1]
            direction = get_direction(U1link)
            position = get_position(U1link)
            #println("i = $i j = $j position = $position")
            substitute_U!(Uold,U[direction])
            Ushift1 = shift_U(Uold,position)
            isU1dag = ifelse(typeof(U1link) <: Adjoint_GLink,true,false)

            for j=2:numlinks
                Ujlink = glinks[j]
                isUkdag = ifelse(typeof(Ujlink) <: Adjoint_GLink,true,false)
                position = get_position(Ujlink)
                direction = get_direction(Ujlink)
                #println("i = $i j = $j position = $position")
                Ushift2 = shift_U(U[direction],position)
                multiply_12!(Unew,Ushift1,Ushift2,j,isUkdag,isU1dag)

                Unew,Uold = Uold,Unew
                Ushift1 = shift_U(Uold,(0,0,0,0))
            end
            add_U!(xout,Uold)
            #println("i = $i ",Uold[1,1,1,1,1,1])
        end
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


    function staggered_phase(μ,iii...)
        error("staggered_phase is not implemented")
    end



end