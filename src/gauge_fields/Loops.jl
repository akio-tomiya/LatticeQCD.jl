module Loops_module
    using LinearAlgebra
    import ..AbstractGaugefields_module:AbstractGaugefields,clear_U!,shift_U,add_U!,substitute_U!,construct_gauges
    import ..Wilsonloops:Wilson_loop_set,calc_coordinate,make_plaq_staple_prime,calc_shift,make_plaq,make_plaq_staple,
                            Tensor_wilson_lines_set,Tensor_wilson_lines,Tensor_derivative_set,
                            get_leftstartposition,get_rightstartposition,Wilson_loop    

    struct Loops{T}
        loopset::Wilson_loop_set
        temps::Array{T,1}
        function Loops(U::Array{T,1},loopset::Wilson_loop_set;mpi = false,PEs=nothing,mpiinit = nothing) where T<: AbstractGaugefields
            NC,NC,NN... = size(U[1]) #NC,NC,NX,NY,NZ,NT 4D case
            NDW = U[1].NDW
            #NC = U[1].NDW

            temp1 = construct_gauges(NC,NDW,NN...;mpi = mpi,PEs=PEs,mpiinit = mpiinit)
            temp2 = construct_gauges(NC,NDW,NN...;mpi = mpi,PEs=PEs,mpiinit = mpiinit)
            temp3 = construct_gauges(NC,NDW,NN...;mpi = mpi,PEs=PEs,mpiinit = mpiinit)

            return new{T}(loopset,[temp1,temp2,temp3])
        end

        function Loops(U::Array{T,1},loopset::Wilson_loop_set,temps) where T<: AbstractGaugefields
            return new{T}(loopset,temps)
        end
    end

    function evaluate_loops(loops::Loops,U::Array{T,1})  where T<: AbstractGaugefields
        xout = deepcopy(loops.temps[1])
        evaluate_loops!(xout,loops,U)
        return xout
    end

    function evaluate_loops!(xout::T,loops::Loops,U::Array{T,1})  where {T<: AbstractGaugefields}
        evaluate_wilson_loops!(xout,loops.loopset,U,loops.temps)
    end

    

    function evaluate_loops!(V,loops::Loops,U::Array{T,1},ix,iy,iz,it) where T<: AbstractGaugefields
        evaluate_wilson_loops!(V,loops.loopset,U,ix,iy,iz,it)
    end

    function evaluate_loops(loops::Loops,U::Array{T,1},ix,iy,iz,it) where T<: AbstractGaugefields
        NC = U[1].NC
        V = zeros(ComplexF64,NC,NC)
        evaluate_loops!(V,loops,U,ix,iy,iz,it)
        return V
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
                Ushift1 = Uold
                #temp1,temp3 = temp3,temp1

                
            end
            add_U!(xout,Ushift1)
            #add!(xout,temp1)
            
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

    function calc_large_wiloson_loop!(temp_Wmat::Array{<: AbstractGaugefields{NC,Dim},2},W_operator,U::Array{T,1}) where {T <: AbstractGaugefields,NC,Dim}
        W = temp_Wmat
        for μ=1:Dim
            for ν=1:Dim
                if μ == ν
                    continue
                end
                #println(typeof(μ)," ",typeof(ν))
                #exit()
                loopset = Loops(U,W_operator[μ,ν])
                W[μ,ν] = evaluate_loops(loopset,U)
            end
        end
        return 
    end

end