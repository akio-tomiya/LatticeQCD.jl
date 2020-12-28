module Mainrun
    using Dates
    import ..LTK_universe:Universe,show_parameters,make_WdagWmatrix,calc_Action,set_β!,set_βs!,get_β
    import ..Actions:Setup_Gauge_action,Setup_Fermi_action,GaugeActionParam_autogenerator
    import ..Measurements:calc_plaquette,measure_correlator,Measurement,calc_polyakovloop,measure_chiral_cond,calc_topological_charge,
                measurements,Measurement_set
    import  ..MD:md_initialize!,MD_parameters_standard,md!,metropolis_update!,construct_MD_parameters
    import ..System_parameters:Params,print_parameters,parameterloading,Params_set#,parameterloading2
    import ..Print_config:write_config
    import ..Smearing:gradientflow!
    import ..ILDG_format:ILDG,load_gaugefield,load_gaugefield!
    import ..Heatbath:heatbath!
    import ..Wilsonloops:make_plaq
    import ..IOmodule:saveU,loadU,loadU!
    import ..SLMC:SLMC_data,show_effbeta,update_slmcdata!
    import ..Gaugefields:calc_GaugeAction

    import ..Actions:GaugeActionParam_standard,
                    GaugeActionParam,
                    GaugeActionParam_autogenerator


    import ..System_parameters:system,actions,md,cg,wilson,staggered,measurement

    function run_LQCD(filename::String)



        include(pwd()*"/"*filename)
        params_set = Params_set(system,actions,md,cg,wilson,staggered,measurement)



        
        run_LQCD(params_set)
    end


    function run_LQCD()
        parameters = parameterloading()

        load_gaugefield()

        univ = Universe(parameters)

        run_LQCD!(univ,parameters)
        return 
    end

    function run_LQCD(parameters::Params)
        univ = Universe(parameters)
        run_LQCD!(univ,parameters)

        return 
    end

    function run_LQCD(params_set::Params_set)
        parameters = parameterloading(params_set)
        univ = Universe(parameters)
        run_LQCD!(univ,parameters)

        return 
    end

    function run_LQCD!(univ::Universe,parameters::Params)
        println("# ",pwd())
        println("# ",Dates.now())

        show_parameters(univ)


        mdparams = construct_MD_parameters(parameters)

        measset = Measurement_set(univ,measurement_methods=parameters.measurement_methods)

        run_core!(parameters,univ,mdparams,measset)


    end

    function run_core!(parameters,univ,mdparams,meas)

        # If an algorithm uses fermion integration (trlog(D+m)),
        # it is flagged
        if parameters.update_method == "IntegratedHMC" || 
                    parameters.update_method == "SLHMC" || 
                    parameters.update_method == "IntegratedHB" ||
                    parameters.update_method == "SLMC" 
            isIntegratedFermion = true
        else
            isIntegratedFermion = false
        end



        Nsteps = parameters.Nsteps
        if parameters.update_method == "Fileloading"
            println("load U from ",parameters.loadU_dir)
            if parameters.loadU_format == "JLD"
                filename_load =  filter(f -> contains(f,".jld"),readdir("./$(parameters.loadU_dir)"))
            elseif parameters.loadU_format == "ILDG"
                filename_load =  filter(f -> contains(f,"ildg"),readdir("./$(parameters.loadU_dir)"))
            end
            #filename = filter(f -> isfile(f), readdir("./$(parameters.loadU_dir)"))
            #println(filename)
            numfiles = length(filename_load)
            println("Num of files = $numfiles")
            Nsteps = numfiles-1
            filename_i = filename_load[1]
            if parameters.loadU_format == "JLD"
                loadU!(parameters.loadU_dir*"/"*filename_i,univ.U)
            elseif parameters.loadU_format == "ILDG"
                ildg = ILDG(parameters.loadU_dir*"/"*filename_i)
                i = 1
                load_gaugefield!(univ.U,i,ildg,parameters.L,parameters.NC)
                measurements(0,univ.U,univ,meas) 
                exit()
            end

        elseif parameters.update_method == "SLHMC" || parameters.update_method == "SLMC"
            #slmc_data = SLMC_data(1,univ.NC)
            if typeof(univ.gparam) == GaugeActionParam_autogenerator
                slmc_data = SLMC_data(length(univ.gparam.couplinglist),univ.NC)
            else
                slmc_data = SLMC_data(1,univ.NC)
            end
            
        end


        numaccepts = 0

        measurements(0,univ.U,univ,meas) # check consistency of preparation.
        if parameters.saveU_format != nothing && parameters.update_method != "Fileloading"
            itrj = 0
            itrjstring = lpad(itrj,8,"0")
            itrjsavecount = 0

            println("save gaugefields U every $(parameters.saveU_every) trajectory")
        end


        if isIntegratedFermion
            Sfold = nothing
        end
        
        for itrj=1:Nsteps
            # Update for different updaters
            # HMC: Hybrid Monte-Carlo
            if parameters.update_method == "HMC"
                Hold = md_initialize!(univ)

                @time Hnew = md!(univ,mdparams)
                accept = metropolis_update!(univ,Hold,Hnew)
                numaccepts += ifelse(accept,1,0)
                println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")

            # Heatbath
            elseif parameters.update_method == "Heatbath"
                @time heatbath!(univ)

            # Integrated HMC
            # HMC with S = -tr(log(D+m)), instead of the pseudo-fermins
            elseif parameters.update_method == "IntegratedHMC"
                Sgold = md_initialize!(univ)
            
                @time Sgnew,Sfnew,Sgold,Sfold = md!(univ,Sfold,Sgold,mdparams)
                Sold = Sgold + Sfold
                Snew = Sgnew + Sfnew
                accept = metropolis_update!(univ,Sold,Snew)
                numaccepts += ifelse(accept,1,0)
                println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                
                Sfold = ifelse(accept,Sfnew,Sfold)

            # File loading
            # This actually loads files, and performs measurements
            elseif parameters.update_method == "Fileloading"
                filename_i = filename_load[itrj+1]
                loadU!(parameters.loadU_dir*"/"*filename_i,univ.U)

            # SLHMC: Self-learing Hybrid Monte-Carlo
            # This uses gluonic effective action for MD
            elseif parameters.update_method == "SLHMC"

                Sgold = md_initialize!(univ)
            
                @time Sgnew,Sfnew,Sgold,Sfold,plaq = md!(univ,Sfold,Sgold,mdparams)
                Sold = Sgold + Sfold
                Snew = Sgnew + Sfnew

            elseif parameters.update_method == "SLMC"
                Sgold = md_initialize!(univ)
            
                @time Sgnew,Sfnew,Sgeffnew,Sgold,Sfold,Sgeffold = md!(univ,Sfold,Sgold,mdparams)
                Sold = Sgold + Sfold -Sgeffold
                Snew = Sgnew + Sfnew -Sgeffnew

            elseif parameters.update_method == "IntegratedHB"
                Sgold = md_initialize!(univ)
            
                @time Sgnew,Sfnew,Sgeffnew,Sgold,Sfold,Sgeffold = md!(univ,Sfold,Sgold,mdparams)
                Sold = Sgold + Sfold -Sgeffold
                Snew = Sgnew + Sfnew -Sgeffnew
                accept = metropolis_update!(univ,Sold,Snew)
                numaccepts += ifelse(accept,1,0)
                println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                
                Sfold = ifelse(accept,Sfnew,Sfold)
            end

            if parameters.update_method == "SLHMC" || parameters.update_method == "SLMC"
                S2,plaqetc = calc_GaugeAction(univ)
                plaq = plaqetc
                #S2,plaq = calc_Action(univ)

                if typeof(plaqetc) == Float64
                    plaqetc = [plaq]
                end
                #println(plaq)
                #exit()
                Sf = Sfnew

                outputdata = univ.gparam.β*plaqetc[1]*slmc_data.factor + Sf
                update_slmcdata!(slmc_data,plaqetc,outputdata)
                #update_slmcdata!(slmc_data,[plaq],outputdata)
                βeffs,Econst,IsSucs = show_effbeta(slmc_data,univ.gparam)
                #βeffs,Econst,IsSucs = show_effbeta(slmc_data)

                println("#S = ",outputdata)
                println("#Estimated Seff = ",Econst + sum(βeffs[:].*plaqetc[:])*slmc_data.factor)
                if IsSucs && itrj ≥ parameters.firstlearn 
                    #mdparams.βeff = βeffs[1]
                    mdparams.βeff = βeffs[:]
                end


                accept = metropolis_update!(univ,Sold,Snew)
                numaccepts += ifelse(accept,1,0)
                println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                Sfold = ifelse(accept,Sfnew,Sfold)
            end# update end

            measurements(itrj,univ.U,univ,meas)

            if itrj % parameters.saveU_every == 0 && parameters.saveU_format != nothing && parameters.update_method != "Fileloading"
                itrjsavecount += 1
                itrjstring = lpad(itrjsavecount,8,"0")
                filename = parameters.saveU_dir*"/conf_$(itrjstring).jld"
                saveU(filename,univ.U)
            end

            
            println("-------------------------------------")
            flush(stdout)
        end
        return
    end

end