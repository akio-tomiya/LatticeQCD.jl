module System_parameters
    using Random
    export Params



    defaultmeasures = Array{Dict,1}(undef,2)
    for i=1:length(defaultmeasures)
        defaultmeasures[i] = Dict()
    end
    defaultmeasures[1]["methodname"] = "Plaquette"
    defaultmeasures[1]["measure_every"] = 1
    defaultmeasures[1]["fermiontype"] = nothing
    defaultmeasures[2]["methodname"] = "Polyakov_loop"
    defaultmeasures[2]["measure_every"] = 1
    defaultmeasures[2]["fermiontype"] = nothing

    system = Dict()
    actions = Dict()
    md = Dict()
    cg = Dict()
    wilson = Dict()
    staggered = Dict()
    measurement = Dict()

    system["Dirac_operator"]  = "Wilson"
    system["quench"] = false
    system["initial"] = "cold"
    system["BoundaryCondition"]=[1,1,1,-1]
    system["Nwing"]= 1
    system["randomseed"] = 111
    #system["IntegratedHMC"] = false
    #system["Heatbath"] = false
    system["verboselevel"] = 2

    system["saveU_format"] = "JLD"
    system["saveU_every"] = 1
    system["saveU_dir"] = "./confs"

    #system["loadU_format"] = "JLD"
    #system["loadU_dir"] = "./confs"

    system["update_method"] = "HMC"

    actions["use_autogeneratedstaples"] = false
    actions["couplinglist"] = []
    actions["couplingcoeff"] = []
    actions["coupling_loops"] = nothing

    md["Δτ"] = 0.05
    md["MDsteps"] = 20
    md["SextonWeingargten"] = false
    system["Nsteps"] = 100
    system["Nthermalization"] = 0

    #md["N_SextonWeingargten"] 

    cg["eps"] = 1e-19
    cg["MaxCGstep"] = 3000

    #wilson["hop"] 
    wilson["r"] = 1
    wilson["Clover_coefficient"] =  1.5612
    

    #staggered["mass"]
    staggered["Nf"] = 4

     

    measurement["measurement_methods"] = defaultmeasures

    mutable struct Params_set
        system::Dict
        actions::Dict
        md::Dict
        cg::Dict
        wilson::Dict
        staggered::Dict
        measurement::Dict

        function Params_set()
            return new(system,actions,md,cg,wilson,staggered,measurement)
        end

        function Params_set(system,actions,md,cg,wilson,staggered,measurement)
            return new(system,actions,md,cg,wilson,staggered,measurement)
        end
    end


    struct Params
        L::Tuple  # = (2,2,2,2) # Mandatory
        β::Float64 # = 6 # Mandatory
        NC::Int64 # = 3 
        Dirac_operator::Union{Nothing,String} # = "Wilson"
        quench::Bool # = false
        initial::String # = "cold"
        BoundaryCondition::Array{Int8,1} # =[1,1,1,-1]
        Nwing::Int64  # = 1 
        randomseed::Int64 # = 111

        #For actions
        use_autogeneratedstaples::Bool # = false
        couplinglist::Array{String,1} # = []
        couplingcoeff::Array{Float64,1} # = []
        coupling_loops::Union{Nothing,Array{Array{Array{Tuple{Int,Int},1},1},1}}  
        
        #For MD
        Δτ::Float64 # = 0.05
        MDsteps::Int64 # = 20
        SextonWeingargten::Bool # = false
        Nsteps::Int64 # = 100
        Nthermalization::Int64 # = 10
        
        #For Sexton-Weingargten
        N_SextonWeingargten::Int64 # = 10 # Mandatory

        #For IntegratedHMC
        #IntegratedHMC::Bool # = false

        #For Heatbath
        #Heatbath::Bool # = false
        
        #For CG
        eps::Float64 # = 1e-19
        MaxCGstep::Int64 # = 3000

        #For Wilson

        hop::Float64 # = 0.141139#Hopping parameter # Mandatory
        r::Float64 # = 1 #Wilson term
        #For WilsonClover
        Clover_coefficient::Float64  # = 1.5612

        #For Staggered
        mass::Float64 # = 0.5 # Mandatory
        Nf::Int8 # = 4

        #verbose
        verboselevel::Int8 # = 2

        #For measurement
        measurement_methods::Array{Dict,1} # = defaultmeasures

        saveU_format::Union{String,Nothing}
        saveU_every::Int64
        saveU_dir::String

        loadU_format::Union{String,Nothing}
        loadU_dir::String

        update_method::String
        βeff::Union{Float64,Array{Float64,1}}
        firstlearn::Int64
        log_dir::String
        logfile::String
        load_fp::IOStream
        measuredir::String
        integratedFermionAction::Bool
        training_data_name::String




        
        function Params(system,actions,md,cg,wilson,staggered,measurement)
            #pnames = fieldnames(Params)
            #display(pnames)

            if haskey(system,"L")
                L = system["L"] #::Tuple  = (2,2,2,2) # Mandatory
            else
                error("system[\"L\"] should be set. eg. L = (4,4,4,4)")
            end

            if haskey(system,"β")
                β = system["β"] #::Float64 # = 6 # Mandatory
            else
                error("system[\"β\"] should be set. eg. β = 2.1")
            end

            
            NC = system["NC"] #::Int64 = 3 
            Dirac_operator = system["Dirac_operator"] #::Union{Nothing,String} = "Wilson"
            quench = system["quench"] #::Bool = false
            initial = system["initial"] #::String = "cold"
            BoundaryCondition = system["BoundaryCondition"] #::Array{Int8,1}=[1,1,1,-1]
            Nwing = system["Nwing"] #::Int64 = 1
            randomseed = system["randomseed"] #::Int64 = 111


            use_autogeneratedstaples = actions["use_autogeneratedstaples"] #::Bool = false
            couplinglist = actions["couplinglist"] #::Array{String,1} = []
            couplingcoeff = actions["couplingcoeff"] #::Array{Float64,1} = []

            coupling_loops = actions["coupling_loops"]

            Δτ = md["Δτ"] #::Float64 = 0.05
            MDsteps = md["MDsteps"]#::Int64 = 20
            SextonWeingargten = md["SextonWeingargten"] #::Bool = false
            Nsteps = system["Nsteps"] #::Int64 = 100
            Nthermalization = system["Nthermalization"] #::Int64 = 10
            if SextonWeingargten
                if haskey(md,"N_SextonWeingargten")
                    N_SextonWeingargten = md["N_SextonWeingargten"]
                else
                    error("md[\"N_SextonWeingargten\"] should be set.")
                end
            else
                N_SextonWeingargten = 0
            end

            
            #For Heatbath
            #Heatbath = system["Heatbath"] #::Bool = false
            
            #For CG
            if quench == false
                eps = cg["eps"] #::Float64= 1e-19
                MaxCGstep = cg["MaxCGstep"] #::Int64 = 3000
            else
                eps = 1e-19
                MaxCGstep = 3000
            end

            #For Wilson
            if Dirac_operator == "Wilson" || Dirac_operator == "WilsonClover"
                println("$Dirac_operator fermion is used")
                if haskey(wilson,"hop")
                    hop = wilson["hop"] #::Float64 # = 0.141139#Hopping parameter # Mandatory
                else
                    error("wilson[\"hop\"] should be set.")
                end           
                r = wilson["r"] #::Float64 = 1 #Wilson term     
            else
                wilson["hop"]  = 0
                hop = 0
                r = 0
            end
            
            #For WilsonClover
            if Dirac_operator == "WilsonClover"
                Clover_coefficient = wilson["Clover_coefficient"] #::Float64  = 1.5612
            else
                wilson["Clover_coefficient"] = 0
                Clover_coefficient = 0
            end

            #For Staggered
            if Dirac_operator == "Staggered"
                println("Staggered fermion is used")
                if haskey(staggered,"mass")
                    mass = staggered["mass"] #::Float64 # = 0.5 # Mandatory
                else
                    error("staggered[\"mass\"] should be set.")
                end
                Nf = staggered["Nf"] #::Int8 = 4
            else
                staggered["mass"] = 0
                mass = 0
                staggered["Nf"] = 0
                Nf = staggered["Nf"] #::Int8 = 4
            end
            

            verboselevel = system["verboselevel"]

            measurement_methods = measurement["measurement_methods"]

            saveU_format = system["saveU_format"]
            seveU_every = system["saveU_every"]
            saveU_dir = system["saveU_dir"]

            if saveU_format ≠ nothing
                if isdir(saveU_dir) == false
                    mkdir(saveU_dir)
                end
                println("$saveU_dir is used for saving configurations")
            end

            βeff = 0.0
            firstlearn = 0

            update_method = system["update_method"]
            if update_method == "HMC"
                println("HMC will be used")
            elseif update_method == "Heatbath"
                println("Heatbath will be used")

                if system["quench"] == false
                    error("system[\"quench\"] = false. The Heatbath method needs the quench update. Put the other system[\"update_method\"] != \"Heatbath\" or system[\"quench\"] = true")
                end
            elseif update_method == "IntegratedHMC"
                println("IntegratedHMC will be used")

                #IntegratedHMC = system["IntegratedHMC"]
                if system["quench"] == false
                    println("system[\"quench\"] = true is set")
                    system["quench"] = true
                    #error("system[\"quench\"] = false. The IntegratedHMC needs the quench update. Put the other system[\"update_method\"] != \"IntegratedHMC\" or system[\"quench\"] = true")
                end
    
            elseif update_method == "Fileloading"
                println("No update will be used (read-measure mode)")
                system["quench"] = true
            elseif update_method == "SLHMC"
                println("SLHMC will be used")
                if system["quench"] == false
                    println("system[\"quench\"] = true is set")
                    system["quench"] = true
                    #error("system[\"quench\"] = false. The SLHMC needs the quench update. Put the other system[\"update_method\"] != \"SLHMC\" or system[\"quench\"] = true")
                end

                if haskey(system,"βeff")
                    βeff = system["βeff"]
                else
                    error("system[\"βeff\"] should be set when you want to do SLHMC.")
                end

                if haskey(system,"firstlearn")
                    firstlearn = system["firstlearn"]
                else
                    error("system[\"firstlearn\"] should be set when you want to do SLHMC.")
                end
            elseif update_method == "SLMC"
                println("SLMC will be used")
                if system["quench"] == false
                    println("system[\"quench\"] = true is set")
                    system["quench"] = true
                    #error("system[\"quench\"] = false. The SLMC needs the quench update. Put the other system[\"update_method\"] != \"SLMC\" or system[\"quench\"] = true")
                end

                if haskey(system,"βeff")
                    βeff = system["βeff"]
                else
                    error("system[\"βeff\"] should be set when you want to do SLMC.")
                end

                if haskey(system,"firstlearn")
                    firstlearn = system["firstlearn"]
                else
                    error("system[\"firstlearn\"] should be set when you want to do SLMC.")
                end
            elseif update_method == "IntegratedHB"
                println("IntegratedHB will be used")

                #IntegratedHMC = system["IntegratedHMC"]
                if system["quench"] == false
                    println("system[\"quench\"] = true is set")
                    system["quench"] = true
                    #error("system[\"quench\"] = false. The IntegratedHB needs the quench update. Put the other system[\"update_method\"] != \"IntegratedHB\" or system[\"quench\"] = true")
                end
            else
                error("""
                system[\"update_method\"] = $update_method is not supported.
                Supported methods are 
                HMC
                Heatbath
                IntegratedHMC
                Fileloading
                """)
            end

            if update_method == "Fileloading"
                if haskey(system,"loadU_format")
                    loadU_format = system["loadU_format"]
                else
                    error("""system[\"loadU_format\"] should be set.
                    Supported formats are 
                    JLD
                    """)
                end

                if haskey(system,"loadU_dir")
                    loadU_dir = system["loadU_dir"]
                else
                    error("system[\"loadU_dir\"] should be set")
                end
                
            else
                loadU_format = nothing
                loadU_dir = "."
            end

            if haskey(system,"log_dir")
                log_dir = system["log_dir"]
                if isdir(log_dir) == false
                    mkdir(log_dir)
                end
            else
                error("system[\"log_dir\"] should be set")
            end

            if haskey(system,"logfile")
                logfile = pwd()*"/"*log_dir*"/"*system["logfile"]
            else
                error("system[\"logfile\"] should be set")
            end
            load_fp = open(logfile,"w")


            if haskey(measurement,"measurement_basedir")
                measurement_basedir = measurement["measurement_basedir"]
                if isdir(measurement_basedir) == false
                    mkdir(measurement_basedir)
                end
            else
                measurement_basedir = "nothing"
            end

            if haskey(measurement,"measurement_dir")
                measurement_dir = measurement["measurement_dir"]
                if isdir(pwd()*"/"*measurement_basedir*"/"*measurement_dir) == false
                    mkdir(pwd()*"/"*measurement_basedir*"/"*measurement_dir)
                end
            else
                measurement_dir = "nothing"
            end
            measuredir = pwd()*"/"*measurement_basedir*"/"*measurement_dir

            if haskey(actions,"IntegratedFermionAction")
                integratedFermionAction = actions["IntegratedFermionAction"]
            else
                integratedFermionAction = false
            end

            if haskey(actions,"training_data_name")
                training_data_name = pwd()*"/"*actions["training_data_name"]
            else
                training_data_name = pwd()*"/"*"trainingdata.txt"
            end



            
            return new(
                L,# = system["L"] #::Tuple  = (2,2,2,2) # Mandatory
                β,# = system["β"] #::Float64 # = 6 # Mandatory
                NC,# = system["NC"] #::Int64 = 3 
                Dirac_operator,# = system["Dirac_operator"] #::Union{Nothing,String} = "Wilson"
                quench,# = system["quench"] #::Bool = false
                initial,# = system["initial"] #::String = "cold"
                BoundaryCondition,# = system["BoundaryCondition"] #::Array{Int8,1}=[1,1,1,-1]
                Nwing,# = system["Nwing"] #::Int64 = 1
                randomseed,# = system["randomseed"] #::Int64 = 111
                use_autogeneratedstaples,# = actions["use_autogeneratedstaples"] #::Bool = false
                couplinglist,# = actions["couplinglist"] #::Array{String,1} = []
                couplingcoeff,# = actions["couplingcoeff"] #::Array{Float64,1} = []
                coupling_loops,
                Δτ,# = md["Δτ"] #::Float64 = 0.05
                MDsteps,# = md["MDsteps"]#::Int64 = 20
                SextonWeingargten,# = md["SextonWeingargten"] #::Bool = false
                Nsteps,# = system["Nsteps"] #::Int64 = 100
                Nthermalization,# = system["Nthermalization"] #::Int64 = 10
                N_SextonWeingargten,# = md["N_SextonWeingargten"]
                #IntegratedHMC,# = system["IntegratedHMC"]
                #Heatbath,# = system["Heatbath"] #::Bool = false
                eps,# = cg["cg"] #::Float64= 1e-19
                MaxCGstep,# = cg["MaxCGstep"] #::Int64 = 3000
                hop,# = wilson["hop"] #::Float64 # = 0.141139#Hopping parameter # Mandatory
                r,# = wilson["r"] #::Float64 = 1 #Wilson term
                Clover_coefficient,# = wilson["Clover_coefficient"] #::Float64  = 1.5612
                mass,# = staggered["mass"] #::Float64 # = 0.5 # Mandatory
                Nf,# = staggered["Nf"] #::Int8 = 4
                verboselevel,
                measurement_methods,
                saveU_format,
                seveU_every,
                saveU_dir,
                loadU_format,
                loadU_dir,
                update_method,
                βeff,
                firstlearn,
                log_dir,
                logfile,
                load_fp,
                measuredir,
                integratedFermionAction,
                training_data_name
            )

        end

        function Params(params_set::Params_set)
            return Params(params_set.system,params_set.actions,params_set.md,params_set.cg,params_set.wilson,params_set.staggered,params_set.measurement)
        end

    end

    

    function make_parametersdict(p::T) where T <: Any
        pnames = fieldnames(T)
        pdict = Dict()
        for i=1:length(pnames)
            pdict[String(pnames[i])] = getfield(p,pnames[i])
        end
        return pdict,pnames
    end

    function print_parameters_file(p)
        filename = p.logfile*"_parameters.jl"# "parameters_used.jl"
        fp = open(filename,"w")
        println(fp,"# - - parameters - - - - - - - - - - - ")
        pdict,pnames = make_parametersdict(p)
        for param in pdict
            if typeof(param[2]) == String
                println(fp,"$(param[1]) = \"$(param[2])\"")
            else
                println(fp,"$(param[1]) = $(param[2])")
            end
        end
        println(fp,"# - - - - - - - - - - - - - - - - - - -")
        close(fp)
        return
    end

    function print_parameters(params_set::Params_set,p)
        #(system,actions,md,cg,wilson,staggered,measurement)
        filename = p.logfile*"_parameters.jl"
        #logfile = pwd()*"/"*log_dir*"/"*system["logfile"]
        fp = open(filename,"w")
        fp2 = p.load_fp
        println(fp,"# - - parameters - - - - - - - - - - - ")
        println(fp2,"# - - parameters - - - - - - - - - - - ")
        println("# - - parameters - - - - - - - - - - - ")
        
        for (name,key) in params_set.system
            if typeof(key) == String
                println("system[\"$name\"] = \"$key\"")
                println(fp,"system[\"$name\"] = \"$key\"")
                println(fp2,"system[\"$name\"] = \"$key\"")
            else
                println("system[\"$name\"] = $key")
                println(fp,"system[\"$name\"] = $key")
                println(fp2,"system[\"$name\"] = $key")
            end
        end

        for (name,key) in params_set.actions
            if typeof(key) == String
                println("actions[\"$name\"] = \"$key\"")
                println(fp,"actions[\"$name\"] = \"$key\"")
                println(fp2,"actions[\"$name\"] = \"$key\"")
            else
                println("actions[\"$name\"] = $key")
                println(fp,"actions[\"$name\"] = $key")
                println(fp2,"actions[\"$name\"] = $key")
            end
        end

        for (name,key) in params_set.md
            if typeof(key) == String
                println("md[\"$name\"] = \"$key\"")
                println(fp,"md[\"$name\"] = \"$key\"")
                println(fp2,"md[\"$name\"] = \"$key\"")
            else
                println("md[\"$name\"] = $key")
                println(fp,"md[\"$name\"] = $key")
                println(fp2,"md[\"$name\"] = $key")
            end
        end

        for (name,key) in params_set.cg
            if typeof(key) == String
                println("cg[\"$name\"] = \"$key\"")
                println(fp,"cg[\"$name\"] = \"$key\"")
                println(fp2,"cg[\"$name\"] = \"$key\"")
            else
                println("cg[\"$name\"] = $key")
                println(fp,"cg[\"$name\"] = $key")
                println(fp2,"cg[\"$name\"] = $key")
            end
        end


        for (name,key) in params_set.wilson
            if typeof(key) == String
                println("wilson[\"$name\"] = \"$key\"")
                println(fp,"wilson[\"$name\"] = \"$key\"")  
                println(fp2,"wilson[\"$name\"] = \"$key\"")  
            else
                println("wilson[\"$name\"] = $key")
                println(fp,"wilson[\"$name\"] = $key")
                println(fp2,"wilson[\"$name\"] = $key")
            end
        end


        for (name,key) in params_set.staggered
            if typeof(key) == String
                println("staggered[\"$name\"] = \"$key\"")
                println(fp,"staggered[\"$name\"] = \"$key\"")
                println(fp2,"staggered[\"$name\"] = \"$key\"")
            else
                println("staggered[\"$name\"] = $key")
                println(fp,"staggered[\"$name\"] = $key")
                println(fp2,"staggered[\"$name\"] = $key")
            end
        end

        for (name,key) in params_set.measurement
            if typeof(key) == String
                println("measurement[\"$name\"] = \"$key\"")
                println(fp,"measurement[\"$name\"] = \"$key\"")
                println(fp2,"measurement[\"$name\"] = \"$key\"")
            else
                println("measurement[\"$name\"] = $key")
                println(fp,"measurement[\"$name\"] = $key")
                println(fp2,"measurement[\"$name\"] = $key")
            end
        end

        println("# - - - - - - - - - - - - - - - - - - -")
        println(fp,"# - - - - - - - - - - - - - - - - - - -")
        println(fp2,"# - - - - - - - - - - - - - - - - - - -")
        close(fp)

        println("""
        # Your parameters were written in $filename
        # If you want to do the simulation with same parameters, 
        # Just do 
        # julia run.jl $filename
        """)
        flush(stdout)


    end

    function print_parameters(filename,params_set::Params_set)
        #(system,actions,md,cg,wilson,staggered,measurement)
        fp = open(filename,"w")
        println(fp,"# - - parameters - - - - - - - - - - - ")
        #println("# - - parameters - - - - - - - - - - - ")
        for (name,key) in params_set.system
            if typeof(key) == String
                #println("system[\"$name\"] = \"$key\"")
                println(fp,"system[\"$name\"] = \"$key\"")
            else
                #println("system[\"$name\"] = $key")
                println(fp,"system[\"$name\"] = $key")
            end
        end

        for (name,key) in params_set.actions
            if typeof(key) == String
                #println("actions[\"$name\"] = \"$key\"")
                println(fp,"actions[\"$name\"] = \"$key\"")
            else
                #println("actions[\"$name\"] = $key")
                println(fp,"actions[\"$name\"] = $key")
            end
        end

        for (name,key) in params_set.md
            if typeof(key) == String
                #println("md[\"$name\"] = \"$key\"")
                println(fp,"md[\"$name\"] = \"$key\"")
            else
                #println("md[\"$name\"] = $key")
                println(fp,"md[\"$name\"] = $key")
            end
        end

        for (name,key) in params_set.cg
            if typeof(key) == String
                #println("cg[\"$name\"] = \"$key\"")
                println(fp,"cg[\"$name\"] = \"$key\"")
            else
                #println("cg[\"$name\"] = $key")
                println(fp,"cg[\"$name\"] = $key")
            end
        end

        for (name,key) in params_set.wilson
            if typeof(key) == String
                #println("wilson[\"$name\"] = \"$key\"")
                println(fp,"wilson[\"$name\"] = \"$key\"")  
            else
                #println("wilson[\"$name\"] = $key")
                println(fp,"wilson[\"$name\"] = $key")
            end
        end


        for (name,key) in params_set.staggered
            if typeof(key) == String
                #println("staggered[\"$name\"] = \"$key\"")
                println(fp,"staggered[\"$name\"] = \"$key\"")
            else
                #println("staggered[\"$name\"] = $key")
                println(fp,"staggered[\"$name\"] = $key")
            end
        end

        for (name,key) in params_set.measurement
            if typeof(key) == String
                #println("measurement[\"$name\"] = \"$key\"")
                println(fp,"measurement[\"$name\"] = \"$key\"")
            else
                #println("measurement[\"$name\"] = $key")
                println(fp,"measurement[\"$name\"] = $key")
            end
        end

        #println("# - - - - - - - - - - - - - - - - - - -")
        println(fp,"# - - - - - - - - - - - - - - - - - - -")
        close(fp)

        println("""
        # Your parameters were written in $filename
        # If you want to do the simulation with same parameters, 
        # Just do 
        # julia run.jl $filename
        """)
        flush(stdout)


    end


    function print_parameters(p)
        println("# - - parameters - - - - - - - - - - - ")
        
        pdict,pnames = make_parametersdict(p)
        for param in pdict
            if typeof(param[2]) == String
                println("$(param[1]) = \"$(param[2])\"")
            else
                println("$(param[1]) = $(param[2])")
            end
        end
        println("# - - - - - - - - - - - - - - - - - - -")
        
        print_parameters_file(p)
        println("""
        # Your parameters were written in parameters_used.jl
        # If you want to do the simulation with same parameters, 
        # Just do 
        # julia run.jl parameters_used.jl
        """)
        flush(stdout)
        return
    end



    #=
    system  = System_params()
    md = MDparams()
    fermimon = Wilson_parames()
    measurement = Measurement_parames()
    =#

    #=
    if length(ARGS) > 0
        include(pwd()*"/"*ARGS[1])
    end
    =#

    function defaultparameters()
        return Params_set(system,actions,md,cg,wilson,staggered,measurement)
    end

    #display(system)
    #exit()

    function parameterloading(system,actions,md,cg,wilson,staggered,measurement)
        param_set = Params_set(system,actions,md,cg,wilson,staggered,measurement)
        if param_set.system["Dirac_operator"] == nothing
            param_set.system["quench"] = true
        end
        p = Params(param_set)
        Random.seed!(p.randomseed)

        print_parameters(param_set,p)
        return p
    end

    
    function parameterloading(param_set::Params_set)
        if param_set.system["Dirac_operator"] == nothing
            param_set.system["quench"] = true
        end
        p = Params(param_set)
        Random.seed!(p.randomseed)
        
        print_parameters(param_set,p)
        return p
    end


    function parameterloading()
        param_set = defaultparameters()

        if param_set.system["Dirac_operator"] == nothing
            param_set.system["quench"] = true
        end
        p = Params(param_set)
        Random.seed!(p.randomseed)

        print_parameters(param_set)
        return p


    end



    



end

#using .System_parameters
#System_parameters.parameterloading(ARGS[1])