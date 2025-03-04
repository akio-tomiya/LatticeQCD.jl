module Mainrun
using Dates
using InteractiveUtils
import Gaugefields: Gradientflow, println_verbose_level1, get_myrank

import ..Universe_module: Univ, get_gauge_action, is_quenched
import ..AbstractMD_module: MD, runMD!
import ..AbstractUpdate_module: Updatemethod, update!
import ..AbstractMeasurement_module:
    Plaquette_measurement,
    measure,
    Polyakov_measurement,
    Topological_charge_measurement,
    Energy_density_measurement,
    Measurements_set,
    measure_withflow

import ..LTK_universe:
    Universe,
    show_parameters,
    make_WdagWmatrix,
    calc_Action,
    set_β!,
    set_βs!,
    get_β,
    Wilsonloops_actions,
    calc_looptrvalues,
    calc_trainingdata,
    calc_looptrvalues_site
import ..Actions: Setup_Gauge_action, Setup_Fermi_action, GaugeActionParam_autogenerator
import ..Measurements:
    calc_plaquette,
    measure_correlator,
    Measurement,
    calc_polyakovloop,
    measure_chiral_cond,
    calc_topological_charge,
    measurements,
    Measurement_set
import ..MD:
    md_initialize!, MD_parameters_standard, md!, metropolis_update!, construct_MD_parameters
import ..System_parameters: Params
import ..Print_config: write_config
import ..Smearing: gradientflow!, calc_fatlink_APE, calc_stout
import ..ILDG_format: ILDG, load_gaugefield, load_gaugefield!, save_binarydata
import ..Bridge_format: load_BridgeText!, save_textdata
import ..Heatbath: heatbath!, overrelaxation!
import ..Wilsonloops:
    make_plaq, make_loopforactions, make_plaqloops, make_rectloops, make_polyakovloops
import ..IOmodule: saveU, loadU, loadU!
import ..SLMC: SLMC_data, show_effbeta, update_slmcdata!
import ..Gaugefields: calc_GaugeAction

import ..Actions:
    GaugeActionParam_standard, GaugeActionParam, GaugeActionParam_autogenerator
import ..Verbose_print: println_verbose1, println_verbose2, Verbose_1


import ..System_parameters: system, actions, md, cg, wilson, staggered, measurement
import ..Transform_oldinputfile: transform_to_toml
import ..Parameters_TOML: construct_Params_from_TOML
import ..LQCD: run_LQCD_file




function run_LQCD(filenamein::String)
    plaq = run_LQCD_file(filenamein)
    return plaq
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
    plaq = run_LQCD!(univ, parameters)

    return plaq
end



function run_LQCD_new!(univ::Univ, parameters::Params)
    println_verbose_level1(univ.U[1], "# " * pwd())
    println_verbose_level1(univ.U[1], "# " * string(Dates.now()))

    if get_myrank(univ.U) == 0
        InteractiveUtils.versioninfo()
    end


    #md = MD(parameters,univ)
    #update_Uold!(univ)
    #substitute_U!(univ.Uold,univ.U)
    gauge_action = get_gauge_action(univ)
    quench = is_quenched(univ)
    updatemethod = Updatemethod(
        univ.U,
        gauge_action,
        parameters.update_method,
        quench,
        parameters.Δτ,
        parameters.MDsteps,
        fermi_action=univ.fermi_action,
        SextonWeingargten=parameters.SextonWeingargten,
        loadU_dir=parameters.loadU_dir,
        loadU_format=parameters.loadU_format,
        isevenodd=parameters.isevenodd,
        β=parameters.β,
        ITERATION_MAX=parameters.ITERATION_MAX,
        numOR=parameters.numOR,
        useOR=parameters.useOR,
    )
    #runMD!(univ.U,md)
    eps_flow = 0.01
    numflow = 500
    Nflow = 1

    meas = Measurements_set(univ.U, parameters.measuredir, parameters.measurement_methods)
    gradientflow = Gradientflow(univ.U, Nflow=1, eps=eps_flow)
    #=
    plaq_m = Plaquette_measurement(univ.U,filename="plaq.txt")
    poly_m = Polyakov_measurement(univ.U,filename="poly.txt")
    topo_m = Topological_charge_measurement(univ.U, filename="topo.txt",
            TC_methods = ["plaquette","clover"])
    energy_m = Energy_density_measurement(univ.U, filename="energydensity.txt")
    =#

    #=
    plaq = measure(plaq_m,0,univ.U)
    poly = measure(poly_m,0,univ.U)
    topo = measure(topo_m,0,univ.U)
    energy = measure(energy_m,0,univ.U)
    =#

    dτ = Nflow * eps_flow
    measure(meas, 0, univ.U)
    measure_withflow(meas, 0, univ.U, gradientflow, numflow, dτ)


    for itrj = parameters.initialtrj:parameters.Nsteps
        println_verbose_level1(univ.U[1], "# itrj = $itrj")
        @time update!(updatemethod, univ.U)
        measure(meas, itrj, univ.U)
        measure_withflow(meas, itrj, univ.U, gradientflow, numflow, dτ)
        #=
        plaq = measure(plaq_m,itrj,univ.U)
        poly = measure(poly_m,itrj,univ.U)
        topo = measure(topo_m,itrj,univ.U)
        energy = measure(energy_m,itrj,univ.U)
        =#

    end

    error("error in run_LQCD_new!")
end

function run_LQCD!(univ::Universe, parameters::Params)
    verbose = univ.kind_of_verboselevel
    #println("# ",pwd())
    #println("# ",Dates.now())
    println_verbose1(verbose, "# ", pwd())
    println_verbose1(verbose, "# ", Dates.now())
    if univ.U[1].myrank == 0
        InteractiveUtils.versioninfo()
    end

    #show_parameters(univ)


    mdparams = construct_MD_parameters(parameters)

    measset = Measurement_set(
        univ,
        parameters.measuredir,
        measurement_methods=parameters.measurement_methods,
    )

    #if isdemo
    #run_demo!(parameters,univ,measset)
    #else
    plaq = run_core!(parameters, univ, mdparams, measset)
    #end


end

function checkfile(filename, list)
    for name in list

    end
end



function run_init_Fileloading!(parameters, univ, mdparams, meas, verbose)
    ildg = nothing
    println_verbose1(verbose, "load U from ", parameters.loadU_dir)
    if parameters.loadU_format == "JLD"
        datatype = "JLD"
        filename_load =
            filter(f -> contains(f, ".jld2"), readdir("./$(parameters.loadU_dir)"))
        #filename_load =  filter(f -> contains(f,".jld"),readdir("./$(parameters.loadU_dir)"))
    elseif parameters.loadU_format == "ILDG"
        datatype = "ILDG"
        filename_load =
            filter(f -> contains(f, "ildg"), readdir("./$(parameters.loadU_dir)"))
    elseif parameters.loadU_format == "BridgeText"
        datatype = "BridgeText"
        filename_load =
            filter(f -> contains(f, "txt"), readdir("./$(parameters.loadU_dir)"))
    else
        error("loadU_format should be JLD or ILDG")
    end

    if parameters.loadU_fromfile
        data = readlines("./$(parameters.loadU_dir)/$(parameters.loadU_filename)")
        datafromlist = String[]

        for (i, datai) in enumerate(data)
            havesharp = findfirst('#', datai)
            #println("havesharp: $havesharp")
            if havesharp == 1
                datau = ""
            elseif havesharp != nothing
                datau = datai[begin:havesharp-1]
            else
                datau = datai
            end

            datau = strip(datau)

            if datatype == "JLD"
                doescontains = contains(datau, ".jld2")
                #doescontains = contains(datau,".jld")
            elseif datatype == "ILDG"
                doescontains = contains(datau, ".ildg")
            elseif datatype == "BridgeText"
                doescontains = contains(datau, ".txt")
            end



            #println(doescontains)
            if doescontains
                if isfile("./$(parameters.loadU_dir)/$(datau)")
                    push!(datafromlist, datau)
                else
                    println(
                        "Warning! $(datau) does not exist in $(parameters.loadU_dir)! We ignore this.",
                    )
                end
                #println(datau)

            end
        end
        filename_load = datafromlist
        #println(datafromlist)

    end

    #exit()


    #filename = filter(f -> isfile(f), readdir("./$(parameters.loadU_dir)"))
    #println(filename)
    numfiles = length(filename_load)
    println_verbose1(verbose, "Num of files = $numfiles")
    for file in filename_load
        println_verbose1(verbose, "$file")
    end
    #println("Num of files = $numfiles")
    Nsteps = numfiles - 1
    filename_i = filename_load[1]
    if parameters.loadU_format == "JLD"
        loadU!(parameters.loadU_dir * "/" * filename_i, univ.U)
    elseif parameters.loadU_format == "ILDG"
        ildg = ILDG(parameters.loadU_dir * "/" * filename_i)
        i = 1
        load_gaugefield!(univ.U, i, ildg, parameters.L, parameters.NC)
    elseif parameters.loadU_format == "BridgeText"
        load_BridgeText!(
            parameters.loadU_dir * "/" * filename_i,
            univ.U,
            parameters.L,
            parameters.NC,
        )
    end
    return Nsteps, numfiles, filename_load, ildg

end


function print_trainingdata(trs, Sg, Sf)
    for i = 1:length(trs)
        print(real(trs[i]), "\t", imag(trs[i]), "\t")
    end
    println(Sg, "\t", Sf, "\t #trainingdata")
end

function print_trainingdata(fp, trs, Sg, Sf)
    for i = 1:length(trs)
        print(fp, real(trs[i]), "\t", imag(trs[i]), "\t")
    end
    println(fp, Sg, "\t", Sf, "\t #trainingdata")
end

function print_trainingdata_site(fp, trs, Sg, Sf)
    return

    numloops, NX, NY, NZ, NT = size(trs)
    println(fp, Sg, "\t", Sf, " $numloops $NX $NY $NZ $NT")
    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    print(fp, "$ix $iy $iz $it ")
                    for i = 1:numloops
                        print(
                            fp,
                            real(trs[i, ix, iy, iz, it]),
                            "\t",
                            imag(trs[i, ix, iy, iz, it]),
                            "\t",
                        )
                    end
                    println(fp, "\t")

                end
            end
        end
    end
    #println(fp,Sg,"\t",Sf,"\t #trainingdata_site")

end



function run_core!(parameters, univ, mdparams, meas)

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
        trainingfp = open(parameters.training_data_name, "w")
        trainingfp2 = open("trainset.txt", "w")
        couplinglist = loopactions.couplinglist
        print(trainingfp, "# ")
        for i = 1:loopactions.numloops
            print(
                trainingfp,
                "Re($(loopactions.couplinglist[i])) Im($(loopactions.couplinglist[i])) ",
            )
        end
        println(trainingfp, "Sg Sf")
        #println(trainingfp,"#Re(plaq) Im(plaq) Re(rect) Im(rect) Re(polyx) Im(polyy) Re(polyz) Im(polyz) Re(polyt) Im(polyt) Sg Sf")
    end


    Nsteps = parameters.Nsteps
    if parameters.update_method == "Fileloading"
        Nsteps, numfiles, filename_load, ildg =
            run_init_Fileloading!(parameters, univ, mdparams, meas, verbose)
    elseif parameters.update_method == "SLHMC" || parameters.update_method == "SLMC"
        #slmc_data = SLMC_data(1,univ.NC)
        if typeof(univ.gparam) == GaugeActionParam_autogenerator
            slmc_data = SLMC_data(length(univ.gparam.couplinglist), univ.NC)
        else
            slmc_data = SLMC_data(1, univ.NC)
        end
    end


    numaccepts = 0

    if (parameters.Nthermalization ≤ 0 && parameters.initialtrj == 1) ||
       parameters.update_method == "Fileloading"
        plaq, poly =
            measurements(0, univ.U, univ, meas; verbose=univ.kind_of_verboselevel) # check consistency of preparation.
    end

    if parameters.integratedFermionAction
        trs, Sg, Sf = calc_trainingdata(loopactions, univ)
        print_trainingdata(trs, Sg, Sf)
        print_trainingdata(trainingfp, trs, Sg, Sf)
        trs_site = calc_looptrvalues_site(loopactions, univ)
        print_trainingdata_site(trainingfp2, trs_site, Sg, Sf)
    end


    if parameters.saveU_format != nothing && parameters.update_method != "Fileloading"
        itrj = 0
        itrjstring = lpad(itrj, 8, "0")
        itrjsavecount = 0
        println_verbose1(
            verbose,
            "save gaugefields U every $(parameters.saveU_every) trajectory",
        )
        #println("save gaugefields U every $(parameters.saveU_every) trajectory")
    end


    if isIntegratedFermion
        Sfold = nothing
    end

    for itrj = parameters.initialtrj:Nsteps
        println("# itrj = $itrj")
        # Update for different updaters
        # HMC: Hybrid Monte-Carlo
        if parameters.update_method == "HMC"
            Hold = md_initialize!(univ)

            @time Hnew = md!(univ, mdparams)

            if parameters.integratedFermionAction
                trs, Sg, Sf = calc_trainingdata(loopactions, univ)
                print_trainingdata(trs, Sg, Sf)
                print_trainingdata(trainingfp, trs, Sg, Sf)
                trs_site = calc_looptrvalues_site(loopactions, univ)
                print_trainingdata_site(trainingfp2, trs_site, Sg, Sf)
            end


            accept = metropolis_update!(univ, Hold, Hnew)
            numaccepts += ifelse(accept, 1, 0)
            #println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
            println_verbose1(
                verbose,
                "Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %",
            )

            # Heatbath
        elseif parameters.update_method == "Heatbath"
            @time heatbath!(univ)
            if parameters.useOR
                for ior = 1:parameters.numOR
                    @time overrelaxation!(univ)
                end
            end

            if parameters.integratedFermionAction
                trs, Sg, Sf = calc_trainingdata(loopactions, univ)
                print_trainingdata(trs, Sg, Sf)
                print_trainingdata(trainingfp, trs, Sg, Sf)
                trs_site = calc_looptrvalues_site(loopactions, univ)
                print_trainingdata_site(trainingfp2, trs_site, Sg, Sf)
            end

            # Integrated HMC
            # HMC with S = -tr(log(D+m)), instead of the pseudo-fermins
        elseif parameters.update_method == "IntegratedHMC"
            Sgold = md_initialize!(univ)

            @time Sgnew, Sfnew, Sgold, Sfold = md!(univ, Sfold, Sgold, mdparams)
            Sold = Sgold + Sfold
            Snew = Sgnew + Sfnew

            if parameters.integratedFermionAction
                trs, Sg, Sf = calc_trainingdata(loopactions, univ)
                print_trainingdata(trs, Sg, Sf)
                print_trainingdata(trainingfp, trs, Sg, Sf)
                trs_site = calc_looptrvalues_site(loopactions, univ)
                print_trainingdata_site(trainingfp2, trs_site, Sg, Sf)
            end

            accept = metropolis_update!(univ, Sold, Snew)
            numaccepts += ifelse(accept, 1, 0)
            #println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
            println_verbose1(
                verbose,
                "Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %",
            )


            Sfold = ifelse(accept, Sfnew, Sfold)

            # File loading
            # This actually loads files, and performs measurements
        elseif parameters.update_method == "Fileloading"
            filename_i = filename_load[itrj+1]
            if parameters.loadU_format == "JLD"
                loadU!(parameters.loadU_dir * "/" * filename_i, univ.U)
            elseif parameters.loadU_format == "ILDG"
                ildg = ILDG(parameters.loadU_dir * "/" * filename_i)
                i = 1
                load_gaugefield!(univ.U, i, ildg, parameters.L, parameters.NC)
            elseif parameters.loadU_format == "BridgeText"
                load_BridgeText!(
                    parameters.loadU_dir * "/" * filename_i,
                    univ.U,
                    parameters.L,
                    parameters.NC,
                )
            end

            if parameters.integratedFermionAction
                trs, Sg, Sf = calc_trainingdata(loopactions, univ)
                print_trainingdata(trs, Sg, Sf)
                print_trainingdata(trainingfp, trs, Sg, Sf)
                trs_site = calc_looptrvalues_site(loopactions, univ)
                print_trainingdata_site(trainingfp2, trs_site, Sg, Sf)
            end

            #loadU!(parameters.loadU_dir*"/"*filename_i,univ.U)

            # SLHMC: Self-learing Hybrid Monte-Carlo
            # This uses gluonic effective action for MD
        elseif parameters.update_method == "SLHMC"

            Sgold = md_initialize!(univ)

            @time Sgnew, Sfnew, Sgold, Sfold, plaq = md!(univ, Sfold, Sgold, mdparams)
            Sold = Sgold + Sfold
            Snew = Sgnew + Sfnew

        elseif parameters.update_method == "SLMC"
            Sgold = md_initialize!(univ)

            @time Sgnew, Sfnew, Sgeffnew, Sgold, Sfold, Sgeffold =
                md!(univ, Sfold, Sgold, mdparams)
            Sold = Sgold + Sfold - Sgeffold
            Snew = Sgnew + Sfnew - Sgeffnew

        elseif parameters.update_method == "IntegratedHB"
            Sgold = md_initialize!(univ)

            @time Sgnew, Sfnew, Sgeffnew, Sgold, Sfold, Sgeffold =
                md!(univ, Sfold, Sgold, mdparams)
            Sold = Sgold + Sfold - Sgeffold
            Snew = Sgnew + Sfnew - Sgeffnew

            if parameters.integratedFermionAction
                trs, Sg, Sf = calc_trainingdata(loopactions, univ)
                print_trainingdata(trs, Sg, Sf)
                print_trainingdata(trainingfp, trs, Sg, Sf)
                trs_site = calc_looptrvalues_site(loopactions, univ)
                print_trainingdata_site(trainingfp2, trs_site, Sg, Sf)
            end


            accept = metropolis_update!(univ, Sold, Snew)
            numaccepts += ifelse(accept, 1, 0)
            #println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
            println_verbose1(
                verbose,
                "Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %",
            )

            Sfold = ifelse(accept, Sfnew, Sfold)
        end

        if parameters.update_method == "SLHMC" || parameters.update_method == "SLMC"
            S2, plaqetc = calc_GaugeAction(univ)
            plaq = plaqetc
            #S2,plaq = calc_Action(univ)

            if typeof(plaqetc) == Float64
                plaqetc = [plaq]
            end
            #println(plaq)
            #exit()
            Sf = Sfnew

            outputdata = univ.gparam.β * plaqetc[1] * slmc_data.factor + Sf

            update_slmcdata!(slmc_data, plaqetc, outputdata)
            βeffs, Econst, IsSucs = show_effbeta(slmc_data, univ.gparam)
            println("βeffs = ", βeffs)

            println_verbose1(verbose, "#S = ", outputdata)
            println_verbose1(
                verbose,
                "#Estimated Seff = ",
                Econst + sum(βeffs[:] .* plaqetc[:]) * slmc_data.factor,
            )
            if IsSucs && itrj ≥ parameters.firstlearn
                mdparams.βeff = βeffs[:]
            end


            if parameters.integratedFermionAction
                trs, Sg, Sf = calc_trainingdata(loopactions, univ)
                print_trainingdata(trs, Sg, Sf)
                print_trainingdata(trainingfp, trs, Sg, Sf)
                trs_site = calc_looptrvalues_site(loopactions, univ)
                print_trainingdata_site(trainingfp2, trs_site, Sg, Sf)
            end

            accept = metropolis_update!(univ, Sold, Snew)
            numaccepts += ifelse(accept, 1, 0)
            println_verbose1(
                verbose,
                "Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %",
            )
            #println("Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %")
            Sfold = ifelse(accept, Sfnew, Sfold)
        end# update end

        if (itrj ≥ parameters.Nthermalization) || parameters.update_method == "Fileloading"
            plaq, poly =
                measurements(itrj, univ.U, univ, meas; verbose=univ.kind_of_verboselevel)
        end


        if itrj % parameters.saveU_every == 0 &&
           parameters.saveU_format != nothing &&
           parameters.update_method != "Fileloading"
            itrjsavecount += 1

            itrjstring = lpad(itrj, 8, "0")
            #itrjstring = lpad(itrjsavecount,8,"0")
            if parameters.saveU_format == "JLD"
                #filename = parameters.saveU_dir*"/conf_$(itrjstring).jld"
                filename = parameters.saveU_dir * "/conf_$(itrjstring).jld2"
                saveU(filename, univ.U)
            elseif parameters.saveU_format == "ILDG"
                filename = parameters.saveU_dir * "/conf_$(itrjstring).ildg"
                save_binarydata(univ.U, filename)
            elseif parameters.saveU_format == "BridgeText"
                filename = parameters.saveU_dir * "/conf_$(itrjstring).txt"
                save_textdata(univ.U, filename)
            else
                error("$(parameters.saveU_format) is not supported")
            end
        end

        println_verbose1(verbose, "-------------------------------------")
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
