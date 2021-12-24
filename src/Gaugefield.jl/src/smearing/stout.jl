mutable struct Stoutsmearing <: Abstractsmearing
    staples_for_stout::Array{Array{Wilson_loop_set,1},1}
    tensor_derivative::Array{Tensor_derivative_set,1}
    staples_for_stout_dag::Array{Array{Wilson_loop_set,1},1}
    tensor_derivative_dag::Array{Tensor_derivative_set,1}
    ρs::Array{Array{Float64,1},1}
end

function Stoutsmearing(loops_smearing,ρs)
    numloops = length(loops_smearing)
    staplesforsmear_set, tensor_derivative, staplesforsmear_dag_set, tensor_derivative_dag = make_loops_Stoutsmearing(loops_smearing,rand(numloops))
    return Stoutsmearing(
        staplesforsmear_set, 
        tensor_derivative,
        staplesforsmear_dag_set,
        tensor_derivative_dag,
        ρs
    )
end

function make_loops_Stoutsmearing(loops_smearing,ρs)
    num = length(loops_smearing)
    @assert num == length(ρs) "number of ρ elements in stout smearing scheme should be $num. Now $(length(ρs)). "
    staplesforsmear_set = Array{Wilson_loop_set,1}[]
    staplesforsmear_dag_set = Array{Wilson_loop_set,1}[]
    println("staple for stout smearing")

    tensor_derivative = Array{Tensor_derivative_set,1}(undef,num)
    tensor_derivative_dag = Array{Tensor_derivative_set,1}(undef,num)
    


    for i=1:num
        loop_smearing = loops_smearing[i]

        staplesforsmear = Array{Wilson_loop_set,1}(undef,4)
        staplesforsmear_dag = Array{Wilson_loop_set,1}(undef,4)
        


        staple = make_staples(loop_smearing)

        for μ=1:4
            staplesforsmear_dag[μ] = staple[μ]##make_plaq_staple(μ)
            staplesforsmear[μ] = staplesforsmear_dag[μ]'
            #println("$μ -direction")
            #display(staplesforsmear[μ])
            #staplesforsmear[μ] = make_plaq_staple(μ)
            #staplesforsmear_dag[μ] = staplesforsmear[μ]'
            #println("dagger: $μ -direction")
            #display(staplesforsmear_dag[μ])
        end
        tensor_derivative[i] = Tensor_derivative_set(staplesforsmear)
        tensor_derivative_dag[i] = Tensor_derivative_set(staplesforsmear_dag)

        push!(staplesforsmear_set,staplesforsmear )
        push!(staplesforsmear_dag_set,staplesforsmear_dag)

    end

    return staplesforsmear_set, 
        tensor_derivative, 
        staplesforsmear_dag_set, 
        tensor_derivative_dag 
        
end
