module Actions
    import ..Wilsonloops:Wilson_loop_set,make_staples,Wilson_loop_set,
            make_cloverloops
    import ..SUN_generator:Generator

    abstract type GaugeActionParam end

    abstract type FermiActionParam end

    mutable struct GaugeActionParam_standard <: GaugeActionParam
        β::Float64
        NTRACE::Int64
    end

    
    mutable struct GaugeActionParam_autogenerator <: GaugeActionParam
        βs::Array{Float64,1}
        β::Float64
        numactions::Int64
        NTRACE::Int64
        loops::Array{Wilson_loop_set,1}
        staples::Array{Array{Wilson_loop_set,1},1}
        #couplinglist::Array{String,1}
    end

    





    """
    ```Setup_Gauge_action(β;actiontype="standard")````
    - β: Coupling value

    Set up the information about the Gauge action.
    You can set the coupling value β.

    Now only SU(3) case is supported.
    """
    function Setup_Gauge_action(β;actiontype="standard")
        if actiontype == "standard"
            NTRACE =3
            return GaugeActionParam_standard(β,NTRACE)
        end
    end

    Base.@kwdef struct FermiActionParam_Wilson <: FermiActionParam
        hop::Float64  = 0.141139#Hopping parameter
        r::Float64  = 1#Wilson term
        #eps::Float64 = 1e-8
        eps::Float64 = 1e-19
        Dirac_operator::String = "Wilson"
        MaxCGstep::Int64 = 3000 
        quench::Bool = false
    end

    const Clover_coefficient  = 1.5612

    Base.@kwdef struct FermiActionParam_WilsonClover <: FermiActionParam
        hop::Float64  = 0.141139#Hopping parameter
        r::Float64  = 1#Wilson term
        #eps::Float64 = 1e-8
        eps::Float64 = 1e-19
        Dirac_operator::String = "WilsonClover"
        MaxCGstep::Int64 = 3000 
        Clover_coefficient::Float64 = Clover_coefficient 
        #CloverFμν::Array{ComplexF64,4}
        internal_flags::Array{Bool,1}
        inn_table::Array{Int64,3}
        _ftmp_vectors::Array{Array{ComplexF64,3},1}
        _is1::Array{Int64,1}
        _is2::Array{Int64,1}
        quench::Bool = false
        SUNgenerator::Union{Nothing,Generator}
        _cloverloops::Array{Wilson_loop_set,2}
    end

    Base.@kwdef struct FermiActionParam_Staggered <: FermiActionParam
        mass::Float64 = 0.5
        eps::Float64 = 1e-19
        Dirac_operator::String = "Staggered"
        MaxCGstep::Int64 = 3000 
        quench::Bool = false
        Nf::Int8 = 4
    end


    function GaugeActionParam_autogenerator(βs,loops,NC)#,couplinglist)
        @assert length(βs) == length(loops) "The number of loops should be the number of βs!"
        numactions = length(loops)
        β = βs[1]
        staples = Array{Array{Wilson_loop_set,1},1}(undef,numactions)
        for (i,loop) in enumerate(loops)
            staples[i] = make_staples(loop)
            #=
            println("$i-th action")
            display(loop)
            println("$i-th action's staple")
            for mu=1:4
                println("$mu -direction")
                display(staples[i][mu])
            end
            =#
        end

        return GaugeActionParam_autogenerator(βs,β,numactions,NC,loops,staples)#,couplinglist)
    end

    function show_parameters_action(fparam::FermiActionParam_Staggered)
        println("#--------------------------------------------------")
        println("#Fermion Action parameters")
        println("#type: ",typeof(fparam))
        println("mass = ",fparam.mass)
        println("eps = ",fparam.eps)
        println("Dirac_operator = ",fparam.Dirac_operator)
        println("MaxCGstep = ",fparam.MaxCGstep)
        println("quench = ",fparam.quench)
        println("Nf = ",fparam.Nf)
        println("#--------------------------------------------------")
    end

    function show_parameters_action(fparam::FermiActionParam_Wilson)
        println("#--------------------------------------------------")
        println("#Fermion Action parameters")
        println("#type: ",typeof(fparam))
        println("hop = ",fparam.hop)
        println("r = ",fparam.r)
        println("eps = ",fparam.eps)
        println("Dirac_operator = ",fparam.Dirac_operator)
        println("MaxCGstep = ",fparam.MaxCGstep)
        println("quench = ",fparam.quench)
        println("#--------------------------------------------------")
    end

    function show_parameters_action(fparam::FermiActionParam_WilsonClover)
        println("#--------------------------------------------------")
        println("#Fermion Action parameters")
        println("#type: ",typeof(fparam))
        println("hop = ",fparam.hop)
        println("r = ",fparam.r)
        println("eps = ",fparam.eps)
        println("Dirac_operator = ",fparam.Dirac_operator)
        println("MaxCGstep = ",fparam.MaxCGstep)
        println("Clover_coefficient = ",fparam.Clover_coefficient)
        println("quench = ",fparam.quench)
        println("#--------------------------------------------------")
    end

    function show_parameters_action(gparam::GaugeActionParam_standard)
        println("#--------------------------------------------------")
        println("#Gauge Action parameters")
        println("#type: ",typeof(gparam))
        println("β = ",gparam.β)
        println("NC = ",gparam.NTRACE)
        println("#--------------------------------------------------")
    end

    function show_parameters_action(gparam::GaugeActionParam_autogenerator)
        println("#--------------------------------------------------")
        println("#Gauge Action parameters")
        println("#type: ",typeof(gparam))
        println("#Num. of action terms: ",gparam.numactions)
        println("#coupling βs = ",gparam.βs)
        #=
        println("#actions that we consider: ")
        for (i,loop) in enumerate(gparam.loops)
            println("#------------------------------")
            println("#$i-th action term consists of")
            display(loop)
            println("#------------------------------")
        end
        =#
        println("NC = ",gparam.NTRACE)        
        println("#--------------------------------------------------")
    end

    

    function FermiActionParam_WilsonClover(hop,r,eps,MaxCGstep,NV,Clover_coefficient,NC;quench=false)
        #CloverFμν = zeros(ComplexF64,3,3,NV,6)
        inn_table= zeros(Int64,NV,4,2)
        internal_flags = zeros(Bool,2)
        _ftmp_vectors = Array{Array{ComplexF64,3},1}(undef,6)
        _is1 = zeros(Int64,NV)
        _is2 = zeros(Int64,NV)
        for i=1:6
            _ftmp_vectors[i] = zeros(ComplexF64,3,NV,4)
        end
        #return FermiActionParam_WilsonClover(hop,r,eps,"WilsonClover",MaxCGstep,Clover_coefficient,CloverFμν,inn_table,internal_flags,_ftmp_vectors,_is1,_is2,quench)
        if NC ≥ 4
            SUNgenerator = Generator(NC)
        else
            SUNgenerator = nothing
        end

        _cloverloops = Array{Wilson_loop_set,2}(undef,3,4)
        for μ=1:3
            for ν=μ+1:4
                _cloverloops[μ,ν] = make_cloverloops(μ,ν)
            end
        end

        return FermiActionParam_WilsonClover(hop,r,eps,"WilsonClover",MaxCGstep,Clover_coefficient,inn_table,internal_flags,_ftmp_vectors,_is1,_is2,quench,
                        SUNgenerator,_cloverloops)
    end

    function FermiActionParam_Wilson(hop,r,eps,MaxCGstep;quench=false)
        return FermiActionParam_Wilson(hop,r,eps,"Wilson",MaxCGstep,quench)
    end

    function FermiActionParam_Wilson(hop,r,eps,Dirac,MaxCGstep;quench=false)
        return FermiActionParam_Wilson(hop,r,eps,Dirac,MaxCGstep,quench)
    end

    """
    ```Setup_Fermi_action(Dirac_operator= "Wilson")```

    Set up the information about the Fermion action.

    Now only WilsonFermion case is supported.

    # For example
    ```julia
        fparam = Setup_Fermi_action()
    ```

    The default values are 

    ```julia
    hop::Float64  = 0.141139
    r::Float64  = 1
    eps::Float64 = 1e-19
    Dirac_operator::String = "Wilson"
    MaxCGstep::Int64 = 3000
    ```

    
    - hop : hopping parameter
    - r : Wilson term
    - eps : convergence criteria in the CG method
    - MaxCGstep : maximum number of the CG steps

    If you want to change the parameters for the Wilson Fermions, 
    please do as follows.

    ```julia
        fparam = FermiActionParam_Wilson(hop,r,eps,MaxCGstep)
    ```

    """
    function Setup_Fermi_action(Dirac_operator= "Wilson")
        if Dirac_operator == "Wilson"
            return FermiActionParam_Wilson()
        end
    end




end