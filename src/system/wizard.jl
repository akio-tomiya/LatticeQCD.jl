module Wizard
    using REPL.TerminalMenus
    import ..System_parameters:Params_set,print_parameters,Params

    function print_wizard_logo(outs)
        blue    = "\033[34m"
        red     = "\033[31m"
        green   = "\033[32m"
        magenta = "\033[35m"
        normal  = "\033[0m\033[0m"

        logo = raw"""
    --------------------------------------------------------------------------------  
    run_wizard       
    　　　　　格　　　　　　　格　　　　　　　
    　　　　　色　　　　　　　格　　　　
    　　　　色色色　　　　　　格　　　
    　子子色色色色色子子子子子格子子子子
    　　　　色色色　　　　　　格　　　　
    　　　　　色　　　　　　　格　　　
    　　　　　格　　　　　　　格　　　
    　　　　　力　　　　　　　学　　　　　　LatticeQCD.jl
    　　　　力力力　　　　　学学学　　　
    　子子力力力力力子子子学学学学学子子　　
    　　　　力力力　　　　　学学学　　　　　
    　　　　　力　　　　　　　学　　　　　　
    　　　　　格　　　　　　　格　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　
        """


        logo = replace(logo, "Q" => "$(red)Q$(normal)")
        logo = replace(logo, "C" => "$(blue)C$(normal)")
        logo = replace(logo, "D" => "$(green)D$(normal)")
        logo = replace(logo, "色" => "$(red)色$(normal)")
        logo = replace(logo, "力" => "$(blue)力$(normal)")
        logo = replace(logo, "学" => "$(green)学$(normal)")

        println(outs, logo)
        println(outs,
            "Welcome to a wizard for Lattice QCD.\n",
            "We'll get you set up simulation parameters in no time."
        )
        println("--------------------------------------------------------------------------------")
        println("If you leave the prompt empty, a default value will be used.")
        println("To exit, press Ctrl + c.")
    end

    
    function run_wizard()        
        print_wizard_logo(stdout)
        simple = request("Choose wizard mode",RadioMenu([
            "simple",
            "expert",
        ]))
        if simple == 1
            isexpert = false
            #return run_wizard_simple()
        else
            isexpert = true
        end

        system,md,actions = initialize()
        filename = Base.prompt("put the name of the parameter file you make", default="my_parameters.jl")
        set_verboselevel!(system,isexpert)
        set_randomseed!(system,isexpert)
        set_Lattice_size!(system,isexpert)
        set_gaugegroup!(system,isexpert)
        set_beta!(system,isexpert)

        system["initialtrj"] = 1

        if check_isfileloading()
            wilson,cg,staggered = set_loadingformat!(system)
        else
            set_initialconfs!(system)

            wilson,cg,staggered,system = set_dynamicalfermion!(system,isexpert)
            set_update_method!(system,isexpert)
            set_nthermalization!(system,isexpert)
            set_totaltrajectoy!(system,isexpert)

            if system["update_method"] == "HMC" || system["update_method"] == "IntegratedHMC" || system["update_method"] == "SLHMC"
                set_MDparams!(system,md,isexpert)
            end
        end

        measurement = Dict()
        options = ["Plaquette","Polyakov_loop","Topological_charge","Chiral_condensate","Pion_correlator"]

        set_measurementmethods!(system,measurement,staggered,wilson,options,isexpert)
        
        headername = make_headername(system,staggered,wilson)

        set_measurementdir!(measurement,headername)
        set_logdir!(system,headername)
        set_saveU!(system,headername,isexpert)

        params_set = Params_set(system,actions,md,cg,wilson,staggered,measurement)

        
        print_parameters(filename,params_set)

        #p = Params(params_set)

        println("""
        --------------------------------------------------------------------------------  
        run_wizard is done. 
        
        The returned value in this run_wizard() is params_set.
        If you want to run a simulation in REPL or other Julia codes,  just do

        run_LQCD(params_set)

        or 

        run_LQCD("$filename")

        The output parameter file is $filename. 
        If you want to run a simulation, just do

        julia run.jl $filename

        --------------------------------------------------------------------------------  
        """)
        
        return params_set
    end

    function wilson_wizard_simple!(system)
        wilson = Dict()
        cg = Dict()
        staggered = Dict()
        wtype = 1

        if wtype == 1
            println("Standard Wilson fermion action will be used")
            system["Dirac_operator"] = "Wilson"
        else
            println("Wilson+Clover fermion action will be used")
            system["Dirac_operator"] = "WilsonClover"
        end

        #println("Inut the hopping parameter κ: typical value = 0.141139")
        hop = parse(Float64,Base.prompt("Input the hopping parameter κ", default="0.141139"))
        #hop = parse(Float64,readline(stdin))
        if hop <= 0
            error("Invalid value for κ=$hop. This has to be strictly positive.")
        end
        wilson["hop"] = hop
        println("κ = $hop")
        wilson["r"] = 1
        wilson["Clover_coefficient"] = 1.5612

        eps = 1e-19
        MaxCGstep = 3000

        cg["eps"] = eps
        cg["MaxCGstep"] =MaxCGstep 
        return wilson,cg,staggered,system
    end

    function wilson_wizard!(system)
        wilson = Dict()
        cg = Dict()
        staggered = Dict()
        wtype = request("Choose Wilson fermion type",RadioMenu([
                    "Standard Wilson fermion action",
                    "Wilson+Clover fermion action",
                ]))
        if wtype == 1
            println("Standard Wilson fermion action will be used")
            system["Dirac_operator"] = "Wilson"
        else
            println("Wilson+Clover fermion action will be used")
            system["Dirac_operator"] = "WilsonClover"
        end

        #println("Inut the hopping parameter κ: typical value = 0.141139")
        hop = parse(Float64,Base.prompt("Input the hopping parameter κ", default="0.141139"))
        #hop = parse(Float64,readline(stdin))
        if hop <= 0
            error("Invalid value for κ=$hop. This has to be strictly positive.")
        end
        wilson["hop"] = hop
        println("κ = $hop")
        wilson["r"] = 1
        wilson["Clover_coefficient"] = 1.5612

        eps = parse(Float64,Base.prompt("relative error in CG loops", default="1e-19"))
        MaxCGstep = parse(Int64,Base.prompt("Maximum iteration steps in CG loops", default="3000"))
        if eps<=0
            error("Invalid value for eps=$eps. This has to be strictly positive.")
        end
        if MaxCGstep<=0
            error("Invalid value for MaxCGstep=$MaxCGstep. This has to be strictly positive.")
        end
        cg["eps"] = eps
        cg["MaxCGstep"] =MaxCGstep 
        return wilson,cg,staggered,system
    end

    function staggered_wizard!(system)
        wilson = Dict()
        cg = Dict()
        staggered = Dict()
        system["Dirac_operator"] = "Staggered"

        mass = parse(Float64,Base.prompt("Input mass", default="0.5"))
        if mass<=0
            error("Invalid value for mass=$mass. This has to be strictly positive.")
        end
        staggered["mass"] = mass
        Nftype = request("Choose the number of flavors(tastes)",RadioMenu([
                    "2 (RHMC will be used)",
                    "3 (RHMC will be used)",
                    "4 (HMC will be used)",
                    "8 (HMC will be used)",
                    "1 (RHMC will be used)"
                ]))

        if Nftype == 1
            staggered["Nf"] = 2
        elseif Nftype == 2
            staggered["Nf"] = 3
        elseif Nftype == 3
            staggered["Nf"] = 4
        elseif Nftype == 4
            staggered["Nf"] = 8
        elseif Nftype == 5
            staggered["Nf"] = 1
        end

        eps = parse(Float64,Base.prompt("relative error in CG loops", default="1e-19"))
        MaxCGstep = parse(Int64,Base.prompt("Maximum iteration steps in CG loops", default="3000"))
        if eps<=0
            error("Invalid value for eps=$eps. This has to be strictly positive.")
        end
        if MaxCGstep<=0
            error("Invalid value for MaxCGstep=$MaxCGstep. This has to be strictly positive.")
        end
        cg["eps"] = eps
        cg["MaxCGstep"] = MaxCGstep 
        return wilson,cg,staggered,system
    end

    function simple_wizard()
        system,md,actions = initialize()

        filename = Base.prompt("put the name of the parameter file you make", default="my_parameters.jl")

        system["verboselevel"] = 1



    end

    function initialize()
        system = Dict()
        md = Dict()
        actions = Dict()
        actions["use_autogeneratedstaples"] = false
        actions["couplinglist"] = []
        actions["couplingcoeff"] = []
        system["BoundaryCondition"] = [1,1,1,-1]
        system["Nwing"] = 1

        return system,md,actions
    end

    function set_verboselevel!(system,isexpert)
        if isexpert 
            verboselevel = parse(Int64,Base.prompt("verbose level ?", default="2"))
        else
            verboselevel = 1
        end

        if 1 ≤ verboselevel ≤ 3
            println("verbose level = ",verboselevel)
        else
            error("verbose level should be 1 ≤ verboselevel ≤ 3")
        end
        system["verboselevel"] = verboselevel
    end

    function set_randomseed!(system,isexpert)
        if isexpert 
            system["randomseed"] = parse(Int64,Base.prompt("Input random seed.", default="111"))
        else
            system["randomseed"] = 111
        end
    end

    function set_Lattice_size!(system,isexpert)
        if isexpert 
            println("Input Lattice size, L=(Nx,Ny,Nz,Nt)")
            NX = parse(Int64,Base.prompt("Nx ?", default="4"))
            NY = parse(Int64,Base.prompt("Ny ?", default="4"))
            NZ = parse(Int64,Base.prompt("Nz ?", default="4"))
            NT = parse(Int64,Base.prompt("Nt ?", default="4"))
            #NT = parse(Int64,readline(stdin))
            L = (NX,NY,NZ,NT)
            system["L"] = L
            if (NX<= 0)|(NY<= 0)|(NZ<= 0)|(NT<= 0)
                error("Invalid parameter L=$L, elements must be positive integers")
            end
        else
            NX = parse(Int64,Base.prompt("Input spatial lattice size ", default="4"))
            NT = parse(Int64,Base.prompt("Input temporal lattice size ", default="4"))
            L = (NX,NX,NX,NT)
            system["L"] = L
            if (NX<= 0)|(NT<= 0)
                error("Invalid parameter L=$L, elements must be positive integers")
            end
        end
        println("Lattice is $L")
    end

    function set_gaugegroup!(system,isexpert)
        if isexpert 
            SNC = request("Choose a gauge group",RadioMenu([
                        "SU(3)",
                        "SU(2)",
                    ]))
            NC = ifelse(SNC == 1,3,2)
        else
            NC = 3
        end
        system["NC"] = NC
        println("SU($NC) will be used")
    end

    function set_beta!(system,isexpert)
        NC = system["NC"]
        if NC == 3
            β = parse(Float64,Base.prompt("β ?", default="5.7"))
        elseif NC == 2
            β = parse(Float64,Base.prompt("β ?", default="2.7"))
        end
        system["β"] = β
        if β<0
            error("Invalid value for β=$β. This has to be positive or zero")
        end
    end

    function check_isfileloading()
        fileloading = request("Do you perform only measurements on configurations in a directory? (no update)",RadioMenu([
            "No",
            "Yes",
        ]))
        if fileloading == 2
            isfileloading = true
        else
            isfileloading = false
        end
        return isfileloading
    end

    function set_loadingformat!(system)
        loadtype = request("Choose a configuration format for loading",RadioMenu([
                    "JLD",
                    "ILDG",
            ]))
        system["update_method"] = "Fileloading"

        if loadtype == 1
            system["loadU_format"] = "JLD"
            
        elseif loadtype == 2
            system["loadU_format"] = "ILDG"
        end

        if system["loadU_format"] ≠ nothing
            system["loadU_dir"] = String(Base.prompt("Loading directory", default="./confs"))
        end


        filelist = request("Which configurations do you use?",RadioMenu([
                    "All configurations in the directory",
                    "Configurations written in a list",
                ]))
        
        if filelist==1
            system["loadU_fromfile"] = false
        else
            system["loadU_fromfile"] = true

            system["loadU_filename"] = String(Base.prompt("name of the list in $(system["loadU_dir"])", default="filelist.txt"))
        end

        system["initial"] = "cold"
        system["Dirac_operator"] = nothing
        system["quench"] = true

        cg = Dict()
        wilson = Dict()
        staggered = Dict()

        return wilson,cg,staggered

    end

    function set_initialconfs!(system)
        NC = system["NC"]

        if NC == 3
            initialconf = request("Choose initial configurations",RadioMenu([
                        "cold start",
                        "hot start",
                        "start from a file",
                    ]))
        elseif NC == 2
            initialconf = request("Choose initial configurations",RadioMenu([
                    "cold start",
                    "hot start",
                    "start from a file",
                    "start from one instanton (Radius is half of Nx)",
                ]))
        end
        if initialconf == 1
            system["initial"] = "cold"
        elseif initialconf == 2
            system["initial"] = "hot"
        elseif initialconf == 3
            loadtype = request("Choose a configuration format for loading",RadioMenu([
                "JLD",
                "ILDG",
            ]))

            if loadtype == 1
                system["loadU_format"] = "JLD"
                system["initial"] = String(Base.prompt("Input the file name that you want to use",default="./confs/conf_00000001.jld"))
                
            elseif loadtype == 2
                system["loadU_format"] = "ILDG"
                system["initial"] = String(Base.prompt("Input the file name that you want to use",default="./confs/conf_00000001.ildg"))
            end
            system["initialtrj"] = parse(Int64,Base.prompt("Start trj number?", default="1"))

        elseif initialconf == 4
            system["initial"] = "Start from one instanton"
        end
    end

    function set_smearing!(system)
        smearing = request("Choose smearing scheme for fermions",RadioMenu([
                        "nothing",
                        "stout",
                    ]))
        if smearing == 1
            system["smearing_for_fermion"] = "nothing"
        elseif smearing == 2
            system["smearing_for_fermion"] = "stout"
            system["stout_numlayers"] = 1
            #system["stout_numlayers"] = parse(Int64,Base.prompt("How many stout layers do you consider?", default="1"))
            if system["stout_numlayers"] == 1
                kindsof_loops = ["plaquette","rectangular","chair","polyakov_x","polyakov_y","polyakov_z","polyakov_t"]
                stout_menu = MultiSelectMenu(kindsof_loops)
                choices = request("Select the kinds of loops you want to add in stout smearing:", stout_menu)
                count = 0
                ρs = Float64[]
                loops = String[]
                for  i in choices
                    count += 1
                    ρ = parse(Float64,Base.prompt("coefficient ρ for $(kindsof_loops[i]) loop?", default="0.1"))
                    push!(ρs,ρ)
                    push!(loops,kindsof_loops[i])
                end
                #println(ρs)
                system["stout_ρ"] = ρs
                system["stout_loops"] = loops
                #system["stout_ρ"] = [parse(Float64,Base.prompt("stout parameter ρ ?", default="0.1"))]
            else
                error("system[\"stout_numlayers\"] = $(system["stout_numlayers"]) is not supported yet!")
            end
        end
    end

    function set_dynamicalfermion!(system,isexpert)
        if isexpert 
            ftype = request("Choose a dynamical fermion",RadioMenu([
                        "Nothing (quenched approximation)",
                        "Wilson Fermion (2-flavor)",
                        "Staggered Fermion",
                    ]))
            if ftype == 1
                cg = Dict()
                wilson = Dict()
                staggered = Dict()
                system["Dirac_operator"] = nothing
                system["quench"] = true

                
            elseif ftype == 2
                wilson,cg,staggered,system = wilson_wizard!(system)
                system["quench"] = false
                set_smearing!(system)
            elseif ftype == 3
                wilson,cg,staggered,system = staggered_wizard!(system)
                system["quench"] = false
                set_smearing!(system)
            end
            
        else
            system["Dirac_operator"] = "Wilson"

            wilson,cg,staggered,system = wilson_wizard_simple!(system)
            system["quench"] = false
        end

        return wilson,cg,staggered,system 

    end

    function set_update_method!(system,isexpert)
        if isexpert
            if system["quench"] == true
                methodtype = request("Choose an update method",RadioMenu([
                    "Heatbath",
                    "Hybrid Monte Carlo",
                ]))
                if methodtype == 2
                    system["update_method"] = "HMC"
                else
                    system["update_method"] = "Heatbath"
                    or = request("Use overrelazation method?",RadioMenu([
                            "true",
                            "false",
                        ]))
                    system["useOR"] = ifelse(or==1,true,false)
                    if system["useOR"]
                        system["numOR"] = parse(Int64,Base.prompt("How many times do you want to do the OR?", default="3"))
                    end

                end
            else
                methodtype = request("Choose an update method",RadioMenu([
                    "Hybrid Monte Carlo",
                    "Integrated HMC",
                    "Self-learning Hybrid Monte Carlo (SLHMC)",
                    "Self-learning Monte Carlo (SLMC)",
                ]))
                if methodtype == 1
                    system["update_method"] = "HMC"
                elseif methodtype == 2
                    system["update_method"] = "IntegratedHMC"
                elseif methodtype == 3
                    system["update_method"] = "SLHMC"
                    system["βeff"] = parse(Float64,Base.prompt("Input initial effective β", default="$β"))
                    system["firstlearn"] = parse(Int64,Base.prompt("When do you want to start updating the effective action?", default="10"))
                    system["quench"] = true
                elseif methodtype == 4
                    system["update_method"] = "SLMC"
                    system["βeff"] = parse(Float64,Base.prompt("Input initial effective β", default="$β"))
                    system["firstlearn"] = parse(Int64,Base.prompt("When do you want to start updating the effective action?", default="10"))
                    system["quench"] = true
                end
            end
        else
            system["update_method"] = "HMC"

        end
    end

    function set_nthermalization!(system,isexpert)
        Nthermalization = 0
        if isexpert
            if system["update_method"] == "HMC" || system["update_method"] == "IntegratedHMC" || system["update_method"] == "SLHMC"|| system["update_method"] == "Heatbath" || system["update_method"] == "SLMC"
                Nthermalization = parse(Int64,Base.prompt("Input the number of thermalization steps (no mearurement)", default="0"))
                if Nthermalization<0
                    error("Invalid value for Nthermalization=$Nthermalization. This has to be positive/zero.")
                end
            end
        end
        system["Nthermalization"] = Nthermalization
    end

    function set_totaltrajectoy!(system,isexpert)
        if isexpert
            if system["update_method"] == "HMC" || system["update_method"] == "IntegratedHMC" || system["update_method"] == "SLHMC"|| system["update_method"] == "Heatbath" || system["update_method"] == "SLMC"
                Nsteps = parse(Int64,Base.prompt("Input the number of total trajectories afterthermalization", default="$(100+system["initialtrj"])"))
                if Nsteps<=0
                    error("Invalid value for Nsteps=$Nsteps. This has to be strictly positive.")
                end
            end
        else
            Nsteps = parse(Int64,Base.prompt("Input number of total trajectories", default="100"))
            if Nsteps<=0
                error("Invalid value for Nsteps=$Nsteps. This has to be strictly positive.")
            end
        end
        system["Nsteps"] = Nsteps
    end

    function set_MDparams!(system,md,isexpert)
        if isexpert
            println("Choose parameters for MD")
            MDsteps = parse(Int64,Base.prompt("Input MD steps", default="20"))
            Δτ = parse(Float64,Base.prompt("Input Δτ", default="$(1/MDsteps)"))
            
            #SextonWeingargten = parse(Bool,Base.prompt("Use SextonWeingargten method? true or false", default="false"))

            if system["Dirac_operator"] != nothing
                SW = request("Use SextonWeingargten method? multi-time scale",RadioMenu([
                        "false",
                        "true",
                    ]))
                SextonWeingargten = ifelse(SW==1,false,true)
                
                if SextonWeingargten
                    N_SextonWeingargten = parse(Int64,Base.prompt("Input number of SextonWeingargten steps", default="2"))
                else
                    N_SextonWeingargten = 2
                end
                md["SextonWeingargten"] = SextonWeingargten
                md["N_SextonWeingargten"] = N_SextonWeingargten
            end

            if MDsteps<=0
                error("Invalid value for MDsteps=$MDsteps. This has to be strictly positive.")
            end
            if Δτ<=0
                error("Invalid value for Δτ=$Δτ. This has to be strictly positive.")
            end

            md["MDsteps"] = MDsteps
            md["Δτ"] = Δτ
        else
            MDsteps = 20
            Δτ = 1/MDsteps

            SextonWeingargten = false
            N_SextonWeingargten = 2

            md["MDsteps"] = MDsteps
            md["Δτ"] = Δτ
            md["SextonWeingargten"] = SextonWeingargten
            md["N_SextonWeingargten"] = N_SextonWeingargten

        end
    end

    function set_measurementmethods!(system,measurement,staggered,wilson,options,isexpert)
        if isexpert
            L = system["L"]
            measurementmenu = MultiSelectMenu(options)
            choices = request("Select the measurement methods you want to do:", measurementmenu)
            nummeasurements  = length(choices)
            #println(choices)
            measurement_methods = Array{Dict,1}(undef,nummeasurements)
            count = 0
            
            for  i in choices
                count += 1
                if i == 1
                    measurement_methods[count] = plaq_wizard()
                elseif i == 2
                    measurement_methods[count] = poly_wizard()
                elseif i == 3
                    measurement_methods[count] = topo_wizard(L,system)
                elseif i == 4
                    measurement_methods[count] = chiral_wizard(staggered,system)
                elseif i == 5
                    #measurement_methods[count] = pion_wizard(wilson,system)
                    measurement_methods[count] = pion_wizard(wilson,staggered,system)
                end
            end
        else
            choices =  [1,2,5]
            nummeasurements  = length(choices)
            #println(choices)
            measurement_methods = Array{Dict,1}(undef,nummeasurements)
            count = 0
            
            for  i in choices
                count += 1
                if i == 1
                    measurement_methods[count] =Dict()
                    println("You measure plaquette")
                    measurement_methods[count]["methodname"] = "Plaquette"
                    measurement_methods[count]["measure_every"] = 1#parse(Int64,Base.prompt("How often measure Plaquette loops?", default="1"))
                    measurement_methods[count]["fermiontype"] = nothing

                elseif i == 2
                    measurement_methods[count] =Dict()
                    println("You measure Polyakov loop")
                    measurement_methods[count]["methodname"] = "Polyakov_loop"
                    measurement_methods[count]["measure_every"] = 1#parse(Int64,Base.prompt("How often measure Plaquette loops?", default="1"))
                    measurement_methods[count]["fermiontype"] = nothing
                elseif i == 5
                    measurement_methods[count] =Dict()
                    measurement_methods[count] = pion_wizard_simple(wilson)
                end
            end
        end

        measurement["measurement_methods"] = measurement_methods
    end

    function make_headername(system,staggered,wilson)
        L = system["L"]
        headername = system["update_method"]*"_L"*string(L[1],pad=2)*string(L[2],pad=2)*string(L[3],pad=2)*string(L[4],pad=2)*"_beta"*string(system["β"])

        if system["update_method"] == "HMC"
            if system["quench"] == true
                headername *= "_quenched"
            else
                headername *= "_"*system["Dirac_operator"]
                if system["Dirac_operator"] == "Staggered"
                    headername *= "_mass"*string(staggered["mass"])*"_Nf"*string(staggered["Nf"])
                elseif system["Dirac_operator"] == "Wilson" || system["Dirac_operator"] == "WilsonClover"
                    headername *= "_kappa"*string(wilson["hop"])
                end
            end
        elseif system["update_method"] == "Heatbath"
            headername *= "_quenched"
        else
            if system["Dirac_operator"] != nothing
                headername *= "_"*system["Dirac_operator"]
                if system["Dirac_operator"] == "Staggered"
                    headername *= "_mass"*string(staggered["mass"])*"_Nf"*string(staggered["Nf"])
                elseif system["Dirac_operator"] == "Wilson" || system["Dirac_operator"] == "WilsonClover"
                    headername *= "_kappa"*string(wilson["hop"])
                end
            else
            end
        end
        return headername
    end

    function set_measurementdir!(measurement,headername)
        nummeasurements  =length(measurement["measurement_methods"])
        if nummeasurements != 0
            measurement["measurement_basedir"] = String(Base.prompt("base directory for measurements", default="./measurements"))
            measurement["measurement_dir"] = String(Base.prompt("directory for measurements in $(measurement["measurement_basedir"])/", default=headername))
        end
    end

    function set_logdir!(system,headername)
        system["log_dir"] = String(Base.prompt("log directory", default="./logs"))
        system["logfile"] = String(Base.prompt("logfile name", default=headername*".txt"))
    end

    function set_saveU!(system,headername,isexpert)
        if isexpert
            savetype = request("Choose a configuration format for saving",RadioMenu([
                        "no save",
                        "JLD",
                        "ILDG",
                ]))
            if savetype == 2
                system["saveU_format"] = "JLD"
                
            elseif savetype == 1
                system["saveU_format"] = nothing
                system["saveU_dir"] = ""
            elseif savetype == 3
                system["saveU_format"] = "ILDG"
            end

            if system["saveU_format"] ≠ nothing
                
                system["saveU_every"] = parse(Int64,Base.prompt("How often do you save a configuration in file (Save every)?", default="10"))
                #system["saveU_basedir"] = String(Base.prompt("base directory for saving configuration", default="./confs"))
                system["saveU_dir"] = String(Base.prompt("Saving directory", default="./confs_$(headername)"))
                #system["saveU_dir"] = system["saveU_basedir"]*"/"*system["saveU_dir"]
            end
        else
            system["saveU_format"] = nothing
            system["saveU_dir"] = ""
        end
    end



    function  plaq_wizard()
        method = Dict()

        println("You measure Plaquette loops")
        method["methodname"] = "Plaquette"
        method["measure_every"] = parse(Int64,Base.prompt("How often measure Plaquette loops?", default="1"))
        method["fermiontype"] = nothing

        return method
    end

    function  poly_wizard()
        method = Dict()

        println("You measure Polyakov loops")
        method["methodname"] = "Polyakov_loop"
        method["measure_every"] = parse(Int64,Base.prompt("How often measure Polyakov loops?", default="1"))
        method["fermiontype"] = nothing

        return method
    end

    function  topo_wizard(L,system)
        method = Dict()

        if system["update_method"] == "Fileloading"
            defaultvalue = 1
        else
            defaultvalue = 10
        end


        println("You measure a topological charge")
        method["methodname"] = "Topological_charge"
        method["measure_every"] = parse(Int64,Base.prompt("How often measure a topological charge?", default="$defaultvalue"))
        method["fermiontype"] = nothing
        method["numflow"]  = parse(Int64,Base.prompt("How many times do you want to flow gauge fields to measure the topological charge?", default="10"))

        Nflowsteps = L[1]
        eps_flow = 0.01

        method["Nflowsteps"]  = parse(Int64,Base.prompt("Nflowsteps?", default="$Nflowsteps"))
        method["eps_flow"]  = parse(Float64,Base.prompt("eps_flow?", default="$eps_flow"))



        return method
    end

    function  chiral_wizard(staggered,system)
        if haskey(staggered,"mass")
            mass_default = staggered["mass"]
        else
            mass_default = 0.5
        end


        method = Dict()
        println("You measure chiral condensates with the statteggred fermion")

        method["methodname"] = "Chiral_condensate" 

        if system["update_method"] == "Fileloading"
            defaultvalue = 1
        else
            defaultvalue = 10
        end

        method["measure_every"] = parse(Int64,Base.prompt("How often measure chiral condensates?", default="$defaultvalue"))
        method["fermiontype"] = "Staggered"
        method["mass"] = parse(Float64,Base.prompt("Input mass for the measurement of chiral condensates", default="$mass_default"))

        #Nfsystem = system["Nf"]
        #Nf = parse(Int64,Base.prompt("Number of flavors (tastes) for the measurement of chiral condensates", default="$Nfsystem"))
        #method["Nf"] = Nf
        method["Nf"] = 4
        println("Number of flavors (tastes) for the measurement of chiral condensates is $(method["Nf"])")
        
        

        eps = parse(Float64,Base.prompt("relative error in CG loops", default="1e-19"))
        MaxCGstep = parse(Int64,Base.prompt("Maximum iteration steps in CG loops", default="3000"))
        if eps<=0
            error("Invalid value for eps=$eps. This has to be strictly positive.")
        end
        if MaxCGstep<=0
            error("Invalid value for MaxCGstep=$MaxCGstep. This has to be strictly positive.")
        end
        method["eps"] = eps
        method["MaxCGstep"] = MaxCGstep

        set_smearing!(method)
    

        return method
    end


    function  pion_wizard(wilson,staggered,system)
        println("You measure Pion_correlator")
        method = Dict()
        method["methodname"] = "Pion_correlator" 

        if system["update_method"] == "Fileloading"
            defaultvalue = 1
        else
            defaultvalue = 10
        end


        method["measure_every"] = parse(Int64,Base.prompt("How often measure Pion_correlator?", default="$defaultvalue"))


        wtype = request("Choose fermion type for the measurement of Pion_correlator",RadioMenu([
            "Standard Wilson fermion action",
            "Wilson+Clover fermion action",
            "Staggered fermion action",
        ]))

        
        if wtype == 1
            println("Standard Wilson fermion action will be used for the measurement")
            method["fermiontype"] = "Wilson"
            isWilson = true
        elseif wtype == 2
            println("Wilson+Clover fermion action will be used for the measurement")
            method["fermiontype"] = "WilsonClover"
            isWilson = true
        elseif wtype == 3
            println("Staggered fermion action will be used for the measurement")
            method["fermiontype"] = "Staggered"
            isWilson = false
        end

        if isWilson
            if haskey(wilson,"hop")
                hop_default = wilson["hop"]
            else
                hop_default = 0.141139
            end

            hop = parse(Float64,Base.prompt("Input the hopping parameter κ for the measurement", default="$hop_default"))
            if hop <= 0
                error("Invalid parameter κ=$hop")
            end
            method["hop"] = hop
            method["r"] = 1
        else
            if haskey(staggered,"mass")
                mass_default = staggered["mass"]
            else
                mass_default = 0.5
            end
            method["mass"] = parse(Float64,Base.prompt("Input mass for the measurement of Pion_correlator", default="$mass_default"))
            method["Nf"] = 4
        end

        eps = parse(Float64,Base.prompt("relative error in CG loops", default="1e-19"))
        MaxCGstep = parse(Int64,Base.prompt("Maximum iteration steps in CG loops", default="3000"))
        if eps<=0
            error("Invalid value for eps=$eps. This has to be strictly positive.")
        end
        if MaxCGstep<=0
            error("Invalid value for MaxCGstep=$MaxCGstep. This has to be strictly positive.")
        end
        method["eps"] = eps
        method["MaxCGstep"] = MaxCGstep

        set_smearing!(method)

        return method
    end

    function  pion_wizard(wilson,system)
        if haskey(wilson,"hop")
            hop_default = wilson["hop"]
        else
            hop_default = 0.141139
        end

        if system["update_method"] == "Fileloading"
            defaultvalue = 1
        else
            defaultvalue = 10
        end


        method = Dict()

        println("You measure Pion_correlator with the Wilson quark operator")

        method["methodname"] = "Pion_correlator" 
        method["measure_every"] = parse(Int64,Base.prompt("How often measure Pion_correlator?", default="$defaultvalue"))
        wtype = request("Choose Wilson fermion type for the measurement of Pion_correlator",RadioMenu([
                    "Standard Wilson fermion action",
                    "Wilson+Clover fermion action",
                ]))
        if wtype == 1
            println("Standard Wilson fermion action will be used for the measurement")
            method["fermiontype"] = "Wilson"
        else
            println("Wilson+Clover fermion action will be used for the measurement")
            method["fermiontype"] = "WilsonClover"
        end

        hop = parse(Float64,Base.prompt("Input the hopping parameter κ for the measurement", default="$hop_default"))
        if hop <= 0
            error("Invalid parameter κ=$hop")
        end
        method["hop"] = hop
        method["r"] = 1

        eps = parse(Float64,Base.prompt("relative error in CG loops", default="1e-19"))
        MaxCGstep = parse(Int64,Base.prompt("Maximum iteration steps in CG loops", default="3000"))
        if eps<=0
            error("Invalid value for eps=$eps. This has to be strictly positive.")
        end
        if MaxCGstep<=0
            error("Invalid value for MaxCGstep=$MaxCGstep. This has to be strictly positive.")
        end
        method["eps"] = eps
        method["MaxCGstep"] = MaxCGstep

        set_smearing!(method)

        return method
    end

    function  pion_wizard_simple(wilson)
        if haskey(wilson,"hop")
            hop_default = wilson["hop"]
        else
            hop_default = 0.141139
        end


        method = Dict()

        println("You measure Pion_correlator with the Wilson quark operator")

        method["methodname"] = "Pion_correlator" 
        method["measure_every"] = 10
        println("Standard Wilson fermion action will be used for the measurement")
        method["fermiontype"] = "Wilson"

        hop = hop_default
        if hop <= 0
            error("Invalid parameter κ=$hop")
        end
        method["hop"] = hop
        method["r"] = 1

        return method
    end


end

#using .Wizard
#Wizard.run_wizard()