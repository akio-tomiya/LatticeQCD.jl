module Mainrun
    using Dates
    using InteractiveUtils
    import ..LTK_universe:Universe,show_parameters,make_WdagWmatrix,calc_Action,set_β!,set_βs!,get_β,
                            Wilsonloops_actions,calc_looptrvalues,calc_trainingdata,calc_looptrvalues_site
    import ..Actions:Setup_Gauge_action,Setup_Fermi_action,GaugeActionParam_autogenerator
    import ..Measurements:calc_plaquette,measure_correlator,Measurement,calc_polyakovloop,measure_chiral_cond,calc_topological_charge,
                measurements,Measurement_set
    import  ..MD:md_initialize!,MD_parameters_standard,md!,metropolis_update!,construct_MD_parameters
    import ..System_parameters:Params,print_parameters,parameterloading,Params_set#,parameterloading2
    import ..Print_config:write_config
    import ..Smearing:gradientflow!
    import ..ILDG_format:ILDG,load_gaugefield,load_gaugefield!,save_binarydata
    import ..Heatbath:heatbath!
    import ..Wilsonloops:make_plaq,make_loopforactions,make_plaqloops,make_rectloops,make_polyakovloops
    import ..IOmodule:saveU,loadU,loadU!
    import ..SLMC:SLMC_data,show_effbeta,update_slmcdata!
    import ..Gaugefields:calc_GaugeAction

    import ..Actions:GaugeActionParam_standard,
                    GaugeActionParam,
                    GaugeActionParam_autogenerator
    import ..Verbose_print:println_verbose1,println_verbose2,Verbose_1


    import ..System_parameters:system,actions,md,cg,wilson,staggered,measurement

    


    function run_LQCD(filenamein::String)
        #if isdemo
        #    params_set = Params_set(Demo.system,Demo.actions,Demo.md,Demo.cg,Demo.wilson,Demo.staggered,Demo.measurement)
        #else
        filename = filenamein
        include(pwd()*"/"*filename)
        params_set = Params_set(system,actions,md,cg,wilson,staggered,measurement)
        #end

        plaq = run_LQCD(params_set)
    end



    #=
    function run_LQCD()
        parameters = parameterloading()

        load_gaugefield()

        univ = Universe(parameters)

        run_LQCD!(univ,parameters)
        return 
    end
    =#

    function run_LQCD(parameters::Params)
        univ = Universe(parameters)
        plaq = run_LQCD!(univ,parameters)

        return plaq
    end

    function run_LQCD(params_set::Params_set)
        parameters = parameterloading(params_set)
        univ = Universe(parameters)
        plaq = run_LQCD!(univ,parameters)

        return  plaq
    end

    function run_LQCD!(univ::Universe,parameters::Params)
        verbose = univ.kind_of_verboselevel
        #println("# ",pwd())
        #println("# ",Dates.now())
        println_verbose1(verbose,"# ",pwd())
        println_verbose1(verbose,"# ",Dates.now())
        versioninfo(verbose)

        #show_parameters(univ)


        mdparams = construct_MD_parameters(parameters)

        measset = Measurement_set(univ,parameters.measuredir,measurement_methods=parameters.measurement_methods)

        #if isdemo
        #run_demo!(parameters,univ,measset)
        #else
        plaq = run_core!(parameters,univ,mdparams,measset)
        #end


    end

    function run_init_Fileloading!(parameters,univ,mdparams,meas)
        ildg = nothing
        println_verbose1(verbose,"load U from ",parameters.loadU_dir)
        if parameters.loadU_format == "JLD"
            filename_load =  filter(f -> contains(f,".jld"),readdir("./$(parameters.loadU_dir)"))
        elseif parameters.loadU_format == "ILDG"
            filename_load =  filter(f -> contains(f,"ildg"),readdir("./$(parameters.loadU_dir)"))
        end
        #filename = filter(f -> isfile(f), readdir("./$(parameters.loadU_dir)"))
        #println(filename)
        numfiles = length(filename_load)
        println_verbose1(verbose,"Num of files = $numfiles")
        #println("Num of files = $numfiles")
        Nsteps = numfiles-1
        filename_i = filename_load[1]
        if parameters.loadU_format == "JLD"
            loadU!(parameters.loadU_dir*"/"*filename_i,univ.U)
        elseif parameters.loadU_format == "ILDG"
            ildg = ILDG(parameters.loadU_dir*"/"*filename_i)
            i = 1
            load_gaugefield!(univ.U,i,ildg,parameters.L,parameters.NC)

        end
        return Nsteps,numfiles,filename_load,ildg

    end

    
    function print_trainingdata(trs,Sg,Sf)
        for i=1:length(trs)
            print(real(trs[i]),"\t",imag(trs[i]),"\t")
        end
        println(Sg,"\t",Sf,"\t #trainingdata")
    end

    function print_trainingdata(fp,trs,Sg,Sf)
        for i=1:length(trs)
            print(fp,real(trs[i]),"\t",imag(trs[i]),"\t")
        end
        println(fp,Sg,"\t",Sf,"\t #trainingdata")
    end

    function print_trainingdata_site(fp,trs,Sg,Sf)
        return 
        
        numloops,NX,NY,NZ,NT = size(trs)
        println(fp,Sg,"\t",Sf," $numloops $NX $NY $NZ $NT")
        for it = 1:NT
            for iz = 1:NZ
                for iy = 1:NY
                    for ix=1:NX
                        print(fp,"$ix $iy $iz $it ")
                        for i=1:numloops
                            print(fp,real(trs[i,ix,iy,iz,it]),"\t",imag(trs[i,ix,iy,iz,it]),"\t")
                        end
                        println(fp,"\t")
                        
                    end
                end
            end
        end
        #println(fp,Sg,"\t",Sf,"\t #trainingdata_site")

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
        verbose = univ.kind_of_verboselevel


        if parameters.integratedFermionAction
            loopactions = Wilsonloops_actions(univ)
            trainingfp = open(parameters.training_data_name,"w")
            trainingfp2 = open("trainset.txt","w")
            couplinglist = loopactions.couplinglist
            print(trainingfp,"# ")
            for i=1:loopactions.numloops
                print(trainingfp,"Re($(loopactions.couplinglist[i])) Im($(loopactions.couplinglist[i])) ")
            end
            println(trainingfp,"Sg Sf")
            #println(trainingfp,"#Re(plaq) Im(plaq) Re(rect) Im(rect) Re(polyx) Im(polyy) Re(polyz) Im(polyz) Re(polyt) Im(polyt) Sg Sf")
        end


        Nsteps = parameters.Nsteps
        if parameters.update_method == "Fileloading"
            Nsteps,numfiles,filename_load,ildg = run_init_Fileloading!(parameters,univ,mdparams,meas)
        elseif parameters.update_method == "SLHMC" || parameters.update_method == "SLMC"
            #slmc_data = SLMC_data(1,univ.NC)
            if typeof(univ.gparam) == GaugeActionParam_autogenerator
                slmc_data = SLMC_data(length(univ.gparam.couplinglist),univ.NC)
            else
                slmc_data = SLMC_data(1,univ.NC)
            end
        end


        numaccepts = 0

        if parameters.Nthermalization ≤ 0
            plaq,poly = measurements(0,univ.U,univ,meas;verbose = univ.kind_of_verboselevel) # check consistency of preparation.
        end

        if parameters.integratedFermionAction
            trs,Sg,Sf = calc_trainingdata(loopactions,univ)
            print_trainingdata(trs,Sg,Sf)
            print_trainingdata(trainingfp,trs,Sg,Sf)
            trs_site = calc_looptrvalues_site(loopactions,univ)
            print_trainingdata_site(trainingfp2,trs_site,Sg,Sf)
        end


        if parameters.saveU_format != nothing && parameters.update_method != "Fileloading"
            itrj = 0
            itrjstring = lpad(itrj,8,"0")
            itrjsavecount = 0
            println_verbose1(verbose,"save gaugefields U every $(parameters.saveU_every) trajectory")
            #println("save gaugefields U every $(parameters.saveU_every) trajectory")
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

                if parameters.integratedFermionAction
                    trs,Sg,Sf = calc_trainingdata(loopactions,univ)
                    print_trainingdata(trs,Sg,Sf)
                    print_trainingdata(trainingfp,trs,Sg,Sf)
                    trs_site = calc_looptrvalues_site(loopactions,univ)
                    print_trainingdata_site(trainingfp2,trs_site,Sg,Sf)
                end

                accept = metropolis_update!(univ,Hold,Hnew)
                numaccepts += ifelse(accept,1,0)
                #println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                println_verbose1(verbose,"Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")

            # Heatbath
            elseif parameters.update_method == "Heatbath"
                @time heatbath!(univ)

                if parameters.integratedFermionAction
                    trs,Sg,Sf = calc_trainingdata(loopactions,univ)
                    print_trainingdata(trs,Sg,Sf)
                    print_trainingdata(trainingfp,trs,Sg,Sf)
                    trs_site = calc_looptrvalues_site(loopactions,univ)
                    print_trainingdata_site(trainingfp2,trs_site,Sg,Sf)
                end

            # Integrated HMC
            # HMC with S = -tr(log(D+m)), instead of the pseudo-fermins
            elseif parameters.update_method == "IntegratedHMC"
                Sgold = md_initialize!(univ)
            
                @time Sgnew,Sfnew,Sgold,Sfold = md!(univ,Sfold,Sgold,mdparams)
                Sold = Sgold + Sfold
                Snew = Sgnew + Sfnew

                if parameters.integratedFermionAction
                    trs,Sg,Sf = calc_trainingdata(loopactions,univ)
                    print_trainingdata(trs,Sg,Sf)
                    print_trainingdata(trainingfp,trs,Sg,Sf)
                    trs_site = calc_looptrvalues_site(loopactions,univ)
                    print_trainingdata_site(trainingfp2,trs_site,Sg,Sf)
                end

                accept = metropolis_update!(univ,Sold,Snew)
                numaccepts += ifelse(accept,1,0)
                #println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                println_verbose1(verbose,"Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                
                
                Sfold = ifelse(accept,Sfnew,Sfold)

            # File loading
            # This actually loads files, and performs measurements
            elseif parameters.update_method == "Fileloading"
                filename_i = filename_load[itrj+1]
                if parameters.loadU_format == "JLD"
                    loadU!(parameters.loadU_dir*"/"*filename_i,univ.U)
                elseif parameters.loadU_format == "ILDG"
                    ildg = ILDG(parameters.loadU_dir*"/"*filename_i)
                    i = 1
                    load_gaugefield!(univ.U,i,ildg,parameters.L,parameters.NC)
                end

                if parameters.integratedFermionAction
                    trs,Sg,Sf = calc_trainingdata(loopactions,univ)
                    print_trainingdata(trs,Sg,Sf)
                    print_trainingdata(trainingfp,trs,Sg,Sf)
                    trs_site = calc_looptrvalues_site(loopactions,univ)
                    print_trainingdata_site(trainingfp2,trs_site,Sg,Sf)
                end

                #loadU!(parameters.loadU_dir*"/"*filename_i,univ.U)

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

                if parameters.integratedFermionAction
                    trs,Sg,Sf = calc_trainingdata(loopactions,univ)
                    print_trainingdata(trs,Sg,Sf)
                    print_trainingdata(trainingfp,trs,Sg,Sf)
                    trs_site = calc_looptrvalues_site(loopactions,univ)
                    print_trainingdata_site(trainingfp2,trs_site,Sg,Sf)
                end


                accept = metropolis_update!(univ,Sold,Snew)
                numaccepts += ifelse(accept,1,0)
                #println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                println_verbose1(verbose,"Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                
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
                βeffs,Econst,IsSucs = show_effbeta(slmc_data,univ.gparam)
                println("βeffs = ",βeffs)

                println_verbose1(verbose,"#S = ",outputdata)
                println_verbose1(verbose,"#Estimated Seff = ",Econst + sum(βeffs[:].*plaqetc[:])*slmc_data.factor)
                if IsSucs && itrj ≥ parameters.firstlearn 
                    mdparams.βeff = βeffs[:]
                end


                if parameters.integratedFermionAction
                    trs,Sg,Sf = calc_trainingdata(loopactions,univ)
                    print_trainingdata(trs,Sg,Sf)
                    print_trainingdata(trainingfp,trs,Sg,Sf)
                    trs_site = calc_looptrvalues_site(loopactions,univ)
                    print_trainingdata_site(trainingfp2,trs_site,Sg,Sf)
                end

                accept = metropolis_update!(univ,Sold,Snew)
                numaccepts += ifelse(accept,1,0)
                println_verbose1(verbose,"Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                #println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
                Sfold = ifelse(accept,Sfnew,Sfold)
            end# update end

            if itrj ≥ parameters.Nthermalization
                plaq,poly = measurements(itrj,univ.U,univ,meas;verbose = univ.kind_of_verboselevel)
            end


            if itrj % parameters.saveU_every == 0 && parameters.saveU_format != nothing && parameters.update_method != "Fileloading"
                itrjsavecount += 1
                itrjstring = lpad(itrjsavecount,8,"0")
                if parameters.saveU_format == "JLD"
                    filename = parameters.saveU_dir*"/conf_$(itrjstring).jld"
                    saveU(filename,univ.U)
                elseif parameters.saveU_format == "ILDG"
                    filename = parameters.saveU_dir*"/conf_$(itrjstring).ildg"
                    save_binarydata(univ.U,filename)
                end
            end

            println_verbose1(verbose,"-------------------------------------")
            #println("-------------------------------------")
            flush(stdout)
            flush(verbose)
            if parameters.integratedFermionAction
                flush(trainingfp)
                flush(trainingfp2)
            end
        end
        if parameters.integratedFermionAction
            close(trainingfp)
        end
        return plaq
    end

end