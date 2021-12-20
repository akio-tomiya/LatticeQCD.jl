module Abstractsmearing_module
    import ..Wilsonloops:Wilson_loop_set,make_staples,Wilson_loop_set,
            make_cloverloops,Tensor_derivative_set, make_loops
    import ..AbstractGaugefields_module:AbstractGaugefields,Abstractfields,initialize_TA_Gaugefields,add_force!,
                        exp_aF_U!,clear_U!,add_U!

    abstract type Abstractsmearing end

    struct Nosmearing <: Abstractsmearing 
    end

    include("./stout.jl")
    include("./gradientflow.jl")

    function construct_smearing(smearingparameters,loops_list,L,coefficients,numlayers)
        if smearingparameters == "nothing"
            smearing = Nosmearing()
        elseif smearingparameters == "stout"
            @assert loops_list != nothing "loops should be put if you want to use smearing schemes"
            loops = make_loops(loops_list,L)

            @assert coefficients != nothing "coefficients should be put if you want to use smearing schemes"
            println("stout smearing will be used")
            if numlayers == 1
                input_coefficients = [coefficients]
            else
                input_coefficients = coefficients
            end

            smearing = Stoutsmearing(loops_smearing,input_coefficients)
        else
            
        end
        return smearing
    end

    function calc_smearedU(Uin::Array{T,1},smearing;calcdSdU = false,temps = nothing) where T<: AbstractGaugefields
        if smearing != nothing && typeof(smearing) != Nosmearing
            println(smearing)
            if typeof(smearing) <: SmearingParam_single
                Uout_multi = nothing
                U = apply_smearing_U(Uin,smearing)
            elseif typeof(smearing) <: SmearingParam_multi
                Uout_multi = apply_smearing_U(Uin,smearing)
                U = Uout_multi[end]
            else
                error("something is wrong in calc_smearingU")
            end
            set_wing!(U)  #we want to remove this.
            if calcdSdU 
                dSdU = [temps[end-3],temps[end-2],temps[end-1],temps[end]]    
            else
                dSdU = nothing
            end
        else
            dSdU = nothing
            Uout_multi = nothing
            U = Uin
        end
        return U,Uout_multi,dSdU
    end

    function apply_smearing_U(Uin::Array{T,1},smearing) where T<: Abstractfields
        error("apply_smearing_U is not implemented in type $(typeof(Uin)) ")
    end

    

    #=


    abstract type SmearingParam end

    struct SmearingParam_nosmearing <: SmearingParam 
    end

    abstract type SmearingParam_single <: SmearingParam
    end

    abstract type SmearingParam_multi <: SmearingParam
    end

    mutable struct SmearingParam_stout <: SmearingParam_single
        staples_for_stout::Array{Array{Wilson_loop_set,1},1}
        tensor_derivative::Array{Tensor_derivative_set,1}
        staples_for_stout_dag::Array{Array{Wilson_loop_set,1},1}
        tensor_derivative_dag::Array{Tensor_derivative_set,1}
        ρs::Array{Float64,1}
        #ρs::Array{Float64,1}
    end

    mutable struct SmearingParam_stout_multi <: SmearingParam_multi
        staples_for_stout::Array{Array{Wilson_loop_set,1},1}
        tensor_derivative::Array{Tensor_derivative_set,1}
        staples_for_stout_dag::Array{Array{Wilson_loop_set,1},1}
        tensor_derivative_dag::Array{Tensor_derivative_set,1}
        ρs::Array{Array{Float64,1},1}
        #ρs::Array{Float64,1}
    end

    const Nosmearing = SmearingParam_nosmearing
    const Stout = SmearingParam_stout

    =#
end