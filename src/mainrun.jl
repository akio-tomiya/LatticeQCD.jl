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
    import ..ILDG_format:ILDG,load_gaugefield
    import ..Heatbath:heatbath!
    import ..Wilsonloops:make_plaq
    import ..IOmodule:saveU,loadU,loadU!
    import ..SLMC:SLMC_data,show_effbeta,update_slmcdata!
    import ..Gaugefields:calc_GaugeAction

    import ..System_parameters:system,actions,md,cg,wilson,staggered,measurement

    function run_LQCD(filename::String)

        #=
        filename= "/Users/yuki/ILDG/ckpoint_lat.ildg.1000"
        ildg = ILDG(filename)
        i = 1
        U = load_gaugefield(i,ildg)
        plaq = calc_plaquette(U)
            println("-------------------------------------")
            println("$i-th plaq = ",plaq)
            println("-------------------------------------")
        exit()
        =#

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
        if parameters.upgrade_method == "IntegratedHMC" || 
                    parameters.upgrade_method == "SLHMC" || 
                    parameters.upgrade_method == "IntegratedHB" ||
                    parameters.upgrade_method == "SLMC" 
            isIntegratedFermion = true
        else
            isIntegratedFermion = false
        end

        Nsteps = parameters.Nsteps
        if parameters.upgrade_method == "Fileloading"
            println("load U from ",parameters.loadU_dir)
            filename_load =  filter(f -> contains(f,".jld"),readdir("./$(parameters.loadU_dir)"))
            #filename = filter(f -> isfile(f), readdir("./$(parameters.loadU_dir)"))
            #println(filename)
            numfiles = length(filename_load)
            println("Num of files = $numfiles")
            Nsteps = numfiles-1
            filename_i = filename_load[1]
            loadU!(parameters.loadU_dir*"/"*filename_i,univ.U)
        elseif parameters.upgrade_method == "SLHMC" || parameters.upgrade_method == "SLMC"
            slmc_data = SLMC_data(1,univ.NC)
        end


        numaccepts = 0

        measurements(0,univ.U,univ,meas) # check consistency
        if parameters.saveU_format != nothing && parameters.upgrade_method != "Fileloading"
            itrj = 0
            itrjstring = lpad(itrj,8,"0")
            itrjsavecount = 0

            println("save gaugefields U every $(parameters.saveU_every) trajectory")
        end


        if isIntegratedFermion
            Sfold = nothing
        end
        
        for itrj=1:Nsteps
            
            if parameters.upgrade_method == "HMC"
                Hold = md_initialize!(univ)

                @time Hnew = md!(univ,mdparams)
                accept = metropolis_update!(univ,Hold,Hnew)
                numaccepts += ifelse(accept,1,0)
                println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")

            elseif parameters.upgrade_method == "Heatbath"
                @time heatbath!(univ)
            elseif parameters.upgrade_method == "IntegratedHMC"
                Sgold = md_initialize!(univ)
            
                @time Sgnew,Sfnew,Sgold,Sfold = md!(univ,Sfold,Sgold,mdparams)
                Sold = Sgold + Sfold
                Snew = Sgnew + Sfnew
                accept = metropolis_update!(univ,Sold,Snew)
                numaccepts += ifelse(accept,1,0)
                println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                
                Sfold = ifelse(accept,Sfnew,Sfold)
            elseif parameters.upgrade_method == "Fileloading"
                filename_i = filename_load[itrj+1]
                loadU!(parameters.loadU_dir*"/"*filename_i,univ.U)
            elseif parameters.upgrade_method == "SLHMC"


                Sgold = md_initialize!(univ)
            
                @time Sgnew,Sfnew,Sgold,Sfold,plaq = md!(univ,Sfold,Sgold,mdparams)
                Sold = Sgold + Sfold
                Snew = Sgnew + Sfnew

                S2,plaq2 = calc_Action(univ)
                Sf = Sfnew
                outputdata = univ.gparam.β*plaq2*slmc_data.factor + Sf
                update_slmcdata!(slmc_data,[plaq],outputdata)
                βeffs,Econst,IsSucs = show_effbeta(slmc_data)


                accept = metropolis_update!(univ,Sold,Snew)
                numaccepts += ifelse(accept,1,0)
                
                
                
                
                #println("Sg = ",Sg,"\t",univ.gparam.β*plaq/univ.gparam.NTRACE,
                #"\t",univ.gparam.β*plaq2/univ.gparam.NTRACE,"Sg2 ",S2)
                println("#S = ",outputdata)
                println("#Estimated Seff = ",Econst + βeffs[1]*plaq*slmc_data.factor)
                if IsSucs && itrj ≥ parameters.firstlearn 
                    mdparams.βeff = βeffs[1]
                end
                Sfold = ifelse(accept,Sfnew,Sfold)
                println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
            elseif parameters.upgrade_method == "IntegratedHB"
                Sgold = md_initialize!(univ)
            
                @time Sgnew,Sfnew,Sgeffnew,Sgold,Sfold,Sgeffold = md!(univ,Sfold,Sgold,mdparams)
                Sold = Sgold + Sfold -Sgeffold
                Snew = Sgnew + Sfnew -Sgeffnew
                accept = metropolis_update!(univ,Sold,Snew)
                numaccepts += ifelse(accept,1,0)
                println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                
                Sfold = ifelse(accept,Sfnew,Sfold)
            elseif parameters.upgrade_method == "SLMC"
                Sgold = md_initialize!(univ)
            
                @time Sgnew,Sfnew,Sgeffnew,Sgold,Sfold,Sgeffold = md!(univ,Sfold,Sgold,mdparams)
                Sold = Sgold + Sfold -Sgeffold
                Snew = Sgnew + Sfnew -Sgeffnew

                S2,plaq = calc_GaugeAction(univ)
                Sf = Sfnew

                outputdata = univ.gparam.β*plaq*slmc_data.factor + Sf
                update_slmcdata!(slmc_data,[plaq],outputdata)
                βeffs,Econst,IsSucs = show_effbeta(slmc_data)

                println("#S = ",outputdata)
                println("#Estimated Seff = ",Econst + βeffs[1]*plaq*slmc_data.factor)
                if IsSucs && itrj ≥ parameters.firstlearn 
                    mdparams.βeff = βeffs[1]
                end

                accept = metropolis_update!(univ,Sold,Snew)
                numaccepts += ifelse(accept,1,0)
                println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                Sfold = ifelse(accept,Sfnew,Sfold)
            end

            measurements(itrj,univ.U,univ,meas)


            #dH = Sold -Snew

            

            if itrj % parameters.saveU_every == 0 && parameters.saveU_format != nothing && parameters.upgrade_method != "Fileloading"
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