module Measure_topological_charge_module
    using LinearAlgebra
    import ..AbstractMeasurement_module:AbstractMeasurement,measure
    import ..Gaugefield:AbstractGaugefields,
                                        Traceless_antihermitian!,Traceless_antihermitian,
                                        initialize_TA_Gaugefields
    import ..Gaugefields:calculate_Plaquette,calculate_Polyakov_loop,substitute_U!,calculate_energy_density
    import ..Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3,println_verbose2,println_verbose1,
            print_verbose1,print_verbose2,print_verbose3
    import ..Gaugefield:Wilson_loop,Wilson_loop_set,calc_loopset_μν_name            
    import ..Gaugefield:Loops,evaluate_loops,TA_Gaugefields
    import ..Smearing:gradientflow!,calc_stout!,calc_fatlink_APE!,calc_stout,calc_fatlink_APE,calc_multihit!
    import ..Actions:GaugeActionParam_autogenerator,GaugeActionParam


    mutable struct Measure_topological_charge{T} <: AbstractMeasurement
        filename::String
        fp::IOStream
        tempU::Array{T,1}
        printvalues::Bool
        Nflowsteps::Int64
        eps_flow::Float64
        Usmr::Array{T,1}
        temp_UμνTA::Array{T,2}
        smearing_type::String
        gparam::GaugeActionParam
        numflow::Int64

        function Measure_topological_charge(filename,U::Array{T,1},params,gparam;printvalues = true) where T
            fp = open(filename,"w")
            
            tempU = Array{T,1}(undef,3)
            for i=1:3
                tempU[i] = similar(U[1])
            end

            Usmr = deepcopy(U)

            Nflowsteps = params["Nflowsteps"]#50
            eps_flow =  params["eps_flow"] #0.01

            temp_UμνTA = Array{T,2}(undef,4,4)
            if haskey(params,"smearing_type")
                smearing_type = params["smearing_type"]
            else
                smearing_type = "gradient_flow"
            end

            numflow = params["numflow"]

            

            m = new{T}(filename,fp,tempU,printvalues,
                    Nflowsteps,
                    eps_flow,
                    Usmr,
                    temp_UμνTA,
                    smearing_type,
                    gparam,
                    numflow
            )
            finalizer(m) do m
                close(m.fp)
            end
            return m
        end

    end


    function measure(m::M,itrj,U::Array{<: AbstractGaugefields{NC,Dim},1};verbose = Verbose_2()) where {M <: Measure_topological_charge,NC,Dim}
        #
        temp1 = m.tempU[1]
        temp2 = m.tempU[2]
        temp3 = m.tempU[3]

        println_verbose2(verbose,"# epsilon for the Wilson flow is $(m.eps_flow)")
        substitute_U!(m.Usmr,U)       

        τ = 0.0
        
        plaq = calculate_Plaquette(m.Usmr,temp1,temp2)
        Qplaq = calculate_topological_charge_plaq(m.Usmr,m.temp_UμνTA)
        #println("Qplaq = ",Qplaq)
        Qclover= calculate_topological_charge_clover(m.Usmr,m.temp_UμνTA)
        
        #println("Qclover = ",Qclover)
        Qimproved= calculate_topological_charge_improved(m.Usmr,m.temp_UμνTA,Qclover)
        #println("Qimproved = ",Qimproved)
        clov = calculate_energy_density(m.Usmr)

        if m.printvalues
            println_verbose1(verbose,"$itrj $τ $plaq $clov $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #flow itrj flowtime plaq E Qplaq Qclov Qimproved")
            println(m.fp,"$itrj $τ $plaq $clov $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #flow itrj flowtime plaq E Qplaq Qclov Qimproved")
            flush(stdout)
        end

        F0 = initialize_TA_Gaugefields(U)

        Ftemps = Array{typeof(F0),1}(undef,4)
        Ftemps[1] = F0
        for i=2:4
            Ftemps[i] = initialize_TA_Gaugefields(U)
        end

        temps = Array{typeof(U),1}(undef,2)
        for i=1:2
            temps[i] = similar(U)
        end


        if m.smearing_type == "gradient_flow"
            for iflow = 1:m.numflow#5000 # eps=0.01: t_flow = 50
                @time gradientflow!(m.Usmr,m.gparam,
                    Ftemps,temps,[temp1,temp2,temp3],
                    m.Nflowsteps,m.eps_flow)

                
                #println(typeof(m.Usmr[1]))
                #println(m.Usmr[1][1,1,1,1,1,1])
                #error("m.Usmr")

                plaq = calculate_Plaquette(m.Usmr,temp1,temp2)
                Qplaq = calculate_topological_charge_plaq(m.Usmr,m.temp_UμνTA)
                Qclover= calculate_topological_charge_clover(m.Usmr,m.temp_UμνTA)
                Qimproved= calculate_topological_charge_improved(m.Usmr,m.temp_UμνTA,Qclover)
                clov = calculate_energy_density(m.Usmr)
                #@time Q = calc_topological_charge(Usmr)
                τ = iflow*m.eps_flow*m.Nflowsteps
                if m.printvalues
                    println_verbose1(verbose,"$itrj $(round(τ, digits=3)) $plaq $clov $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #flow itrj flowtime plaq E Qplaq Qclov Qimproved")
                    println(m.fp,"$itrj $(round(τ, digits=3)) $plaq $clov $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #flow itrj flowtime plaq E Qplaq Qclov Qimproved")
                    #if iflow%10 == 0
                    flush(stdout)
                end
            end
        else
            error("$(m.smearig_type) is not suppoorted")
        end

        return 
    end



    function calculate_topological_charge_plaq(U::Array{T,1},temp_UμνTA) where T <: AbstractGaugefields
        UμνTA = temp_UμνTA
        numofloops = calc_UμνTA!(UμνTA,"plaq",U)
        Q = calc_Q(UμνTA,numofloops,U)
        return Q
    end

    function calculate_topological_charge_clover(U::Array{T,1},temp_UμνTA) where T <: AbstractGaugefields
        UμνTA = temp_UμνTA
        numofloops = calc_UμνTA!(UμνTA,"clover",U)
        Q = calc_Q(UμνTA,numofloops,U)
        return Q
    end

    function calculate_topological_charge_improved(U::Array{T,1},temp_UμνTA,Qclover) where T <: AbstractGaugefields
        UμνTA = temp_UμνTA
        #numofloops = calc_UμνTA!(UμνTA,"clover",U)
        #Qclover = calc_Q(UμνTA,numofloops,U)

        numofloops = calc_UμνTA!(UμνTA,"rect",U)
        Qrect = 2*calc_Q(UμνTA,numofloops,U)
        c1 = -1/12
        c0 = 5/3
        Q = c0*Qclover + c1*Qrect
        return Q
    end

    function calc_UμνTA!(temp_UμνTA,name::String,U)
        loops_μν,numofloops = calc_loopset_μν_name(name)
        calc_UμνTA!(temp_UμνTA,loops_μν,U)
        return numofloops
    end

    function calc_UμνTA!(temp_UμνTA,loops_μν,U)
        UμνTA = temp_UμνTA
        for μ=1:4
            for ν=1:4
                if ν == μ
                    continue
                end
                loopset = Loops(U,loops_μν[μ,ν])
                UμνTA[μ,ν] = evaluate_loops(loopset,U)
                UμνTA[μ,ν] = Traceless_antihermitian(UμνTA[μ,ν])
            end
        end
        return 
    end



        #=
    implementation of topological charge is based on
    https://arxiv.org/abs/1509.04259
    =#
    function calc_Q(UμνTA,numofloops,U::Array{<: AbstractGaugefields{NC,Dim},1}) where {NC,Dim}
        Q = 0.0
        if Dim == 4
            ε(μ,ν,ρ,σ) = epsilon_tensor(μ,ν,ρ,σ)  
        else
            error("Dimension $Dim is not supported")
        end
        for μ=1:Dim
            for ν=1:Dim
                if ν == μ
                    continue
                end
                Uμν = UμνTA[μ,ν]                 
                for ρ =1:Dim
                    for σ=1:Dim
                        if ρ == σ
                            continue
                        end
                        Uρσ = UμνTA[ρ,σ]
                        s = tr(Uμν,Uρσ)
                        Q += ε(μ,ν,ρ,σ)*s/numofloops^2
                    end
                end
            end
        end

        return -Q/(32*(π^2))
    end




    #topological charge
    function epsilon_tensor(mu::Int,nu::Int,rho::Int,sigma::Int) 
        sign=1 # (3) 1710.09474 extended epsilon tensor
        if mu < 0
            sign*=-1
            mu=-mu
        end
        if nu < 0
            sign*=-1
            nu=-nu
        end
        if rho < 0
            sign*=-1
            rh=-rho
        end
        if sigma < 0
            sign*=-1
            sigma=-sigma
        end
        epsilon = zeros(Int,4,4,4,4)
        epsilon[ 1, 2, 3, 4 ] = 1
        epsilon[ 1, 2, 4, 3 ] = -1
        epsilon[ 1, 3, 2, 4 ] = -1
        epsilon[ 1, 3, 4, 2 ] = 1
        epsilon[ 1, 4, 2, 3 ] = 1
        epsilon[ 1, 4, 3, 2 ] = -1
        epsilon[ 2, 1, 3, 4 ] = -1
        epsilon[ 2, 1, 4, 3 ] = 1
        epsilon[ 2, 3, 1, 4 ] = 1
        epsilon[ 2, 3, 4, 1 ] = -1
        epsilon[ 2, 4, 1, 3 ] = -1
        epsilon[ 2, 4, 3, 1 ] = 1
        epsilon[ 3, 1, 2, 4 ] = 1
        epsilon[ 3, 1, 4, 2 ] = -1
        epsilon[ 3, 2, 1, 4 ] = -1
        epsilon[ 3, 2, 4, 1 ] = 1
        epsilon[ 3, 4, 1, 2 ] = 1
        epsilon[ 3, 4, 2, 1 ] = -1
        epsilon[ 4, 1, 2, 3 ] = -1
        epsilon[ 4, 1, 3, 2 ] = 1
        epsilon[ 4, 2, 1, 3 ] = 1
        epsilon[ 4, 2, 3, 1 ] = -1
        epsilon[ 4, 3, 1, 2 ] = -1
        epsilon[ 4, 3, 2, 1 ] = 1
        return epsilon[mu,nu,rho,sigma]*sign
    end
end