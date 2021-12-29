module Abstractsmearing_module
    using LinearAlgebra
    import ..Wilsonloops_module:Wilson_loop_set,make_staples,Wilson_loop_set,
            make_cloverloops,Tensor_derivative_set, make_loops
    import ..AbstractGaugefields_module:AbstractGaugefields,Abstractfields,initialize_TA_Gaugefields,add_force!,
                        exp_aF_U!,clear_U!,add_U!,evaluate_wilson_loops!,exptU!,
                        Traceless_antihermitian_add!,set_wing_U!,Traceless_antihermitian,evaluate_gaugelinks!,
                        construct_Λmatrix_forSTOUT!,Traceless_antihermitian!
    import Wilsonloop:Wilsonline,DwDU,make_loopforactions,make_Cμ,derive_U,derive_Udag

    abstract type Abstractsmearing end

    struct Nosmearing <: Abstractsmearing 
    end

    abstract type CovLayer{Dim} end

    struct CovNeuralnet{Dim} <: Abstractsmearing
        numlayers::Int64
        layers::Array{CovLayer{Dim},1}
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

            smearing = Stoutsmearing(loops,input_coefficients)
        elseif smearingparameters == "covnet_stout"
            if numlayers == 1
                input_coefficients = [coefficients]
            else
                input_coefficients = coefficients
            end
            println("covnet verion of the stout smearing will be used")
            smearing = CovNeuralnet_STOUT(loops_list,input_coefficients,L;Dim=length(L))
        else
            error("smearing = $smearing is not supported")
        end
        return smearing
    end

    function calc_smearedU(Uin::Array{T,1},smearing;calcdSdU = false,temps = nothing) where T<: AbstractGaugefields
        if smearing != nothing && typeof(smearing) != Nosmearing
            #println(smearing)
            println(typeof(smearing))
            Uout_multi = apply_smearing_U(Uin,smearing)
            U = Uout_multi[end]

            #=
            if typeof(smearing) <: SmearingParam_single
                Uout_multi = nothing
                U = apply_smearing_U(Uin,smearing)
            elseif typeof(smearing) <: SmearingParam_multi
                Uout_multi = apply_smearing_U(Uin,smearing)
                U = Uout_multi[end]
            else
                error("something is wrong in calc_smearingU")
            end
            =#
            set_wing_U!(U)  #we want to remove this.
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
        error("apply_smearing_U is not implemented in type $(typeof(Uin)) and smearing type $(typeof(smearing))")
    end

    function apply_smearing_U(Uin::Array{T,1},smearing::CovNeuralnet{Dim}) where {T<: Abstractfields,Dim}
        temp1  = similar(Uin[1])
        temp2  = similar(Uin[1])
        temp3  = similar(Uin[1])
        temp4  = similar(Uin[1])
        F0 = initialize_TA_Gaugefields(Uin[1])
        numlayers = smearing.numlayers
        Uout_multi = Array{typeof(Uin),1}(undef,numlayers)
        for i=1:numlayers
            Uout_multi[i] = similar(Uin)
        end
        apply_neuralnet!(Uout_multi,smearing,Uin,[temp1,temp2,temp3,temp4],[F0])
        return Uout_multi

        error("apply_smearing_U is not implemented in type $(typeof(Uin)) and smearing type $(typeof(smearing))")
    end

    function apply_neuralnet!(Uout_multi,net::CovNeuralnet{Dim},Uin,temps,temps_F) where {Dim}
        layer = net.layers[1]
        apply_layer!(Uout_multi[1],layer,Uin,temps,temps_F)
        for i=2:net.numlayers
            layer = net.layers[i]
            apply_layer!(Uout_multi[i],layer,Uout_multi[i-1],temps,temps_F)
        end
    end

    function apply_layer!(Uout,layer::T,Uin,temps,temps_F) where T <: CovLayer
        error("apply_layer!(Uout,layer,Uin,,temps,temps_F) is not implemented with type $(typeof(layer)) of layer.")
    end

    """
    layer_pullback!(δ_prev,δ_next,layer::T,Uprev,temps,tempf) 
    This is a function for a back propagation
        δ_next,Uprev -> δ_prev
    """
    function layer_pullback!(δ_prev,δ_next,layer::T,Uprev,temps,tempf) where T <: CovLayer
        error("layer_pullback!(δ_prev,δ_next,layer,Uprev,temps,tmeps_F) is not implemented with type $(typeof(layer)) of layer.")
    end
    

    function apply_smearing_U(U::Array{T,1},smearing::Stoutsmearing) where T<: Abstractfields
        numlayer = length(smearing.ρs)
        Uout_multi = Array{typeof(U),1}(undef,numlayer)
        for i=1:numlayer
            Uout_multi[i] = similar(U)
        end
        #println("smearing.ρs ", smearing.ρs)
        #println("type U ",typeof(Uout_multi))
        #Uout = similar(U)
        calc_stout_multi!(Uout_multi,U,smearing.ρs,smearing.staples_for_stout) 
        #Uout_multi = calc_stout_multi!(U,smearing.ρs,smearing.staples_for_stout) 
        return Uout_multi
    end

    #=
    function calc_stout_multi(Uin::Array{<: AbstractGaugefields{NC,Dim},1},ρs::Array{Array{T,1},1},staples)  where {NC,Dim,T <: Number}
        numlayer = length(ρs)
        #println("numlayer = ",numlayer,"\t",ρs)
        Utype = eltype(Uin)
        Uout_multi = Array{Array{Utype,1}}(undef,numlayer)
        for i=1:numlayer
            Uout_multi[i] = similar(Uin)
        end
        calc_stout_multi!(Uout_multi,Uin,ρs,staples)

        return Uout_multi
    end
    =#

    function calc_stout_multi!(Uout_multi::Vector{<: Vector{<: AbstractGaugefields{NC,Dim}}},Uin::Array{<: AbstractGaugefields{NC,Dim},1},ρs::Array{Array{T,1},1},staples)  where {NC,Dim,T <: Number}
        numlayer = length(ρs)
        Utmp = similar(Uin)
        #Uout_multi = Array{Array{GaugeFields{SU{NC}},1}}(undef,numlayer)
        U = deepcopy(Uin)
        for i = 1:numlayer
            if i != numlayer
                apply_stout_smearing!(Utmp,U,ρs[i],staples)
                Uout_multi[i] = deepcopy(Utmp)
                Utmp,U = U,Utmp            
            else
                apply_stout_smearing!(Uout_multi[i],U,ρs[i],staples)
            end
        end
    end

    function calc_stout_multi!(Uout::Array{<: AbstractGaugefields{NC,Dim},1},Uin::Array{<: AbstractGaugefields{NC,Dim},1},ρs::Array{Array{T,1},1},staples)  where {NC,Dim,T <: Number}
        numlayer = length(ρs)
        Utmp = similar(Uin)
        #Uout_multi = Array{Array{GaugeFields{SU{NC}},1}}(undef,numlayer)
        U = deepcopy(Uin)
        for i = 1:numlayer
            if i != numlayer
                apply_stout_smearing!(Utmp,U,ρs[i],staples)
                Utmp,U = U,Utmp            
            else
                apply_stout_smearing!(Uout,U,ρs[i],staples)
            end
        end
        
    end

    function apply_stout_smearing!(Uout::Array{<: AbstractGaugefields{NC,Dim},1},
        U::Array{<: AbstractGaugefields{NC,Dim},1},ρs::Array{T,1},staples) where {NC,Dim,T <: Number}
        @assert Uout != U "input U and output U should not be same!"
        V  = similar(U[1])
        temp1  = similar(U[1])
        temp2  = similar(U[1])
        temp3  = similar(U[1])
        F0 = initialize_TA_Gaugefields(U[1])

        num = length(ρs)

        for μ=1:Dim
            clear_U!(V)
            for i=1:num
                loops = staples[i][μ]
                evaluate_wilson_loops!(temp3,loops,U,[temp1,temp2])
                add_U!(V,ρs[i],temp3)
            end
            mul!(temp1,V,U[μ]') #U U*V
            clear_U!(F0)
            Traceless_antihermitian_add!(F0,1,temp1)
            
            exptU!(temp3,1,F0,[temp1,temp2])
            
            
            mul!(Uout[μ],temp3,U[μ])        
        end
        set_wing_U!(Uout)


        #error("ee")

        #error("not implemented")
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