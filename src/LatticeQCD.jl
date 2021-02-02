module LatticeQCD
    using Requires

    include("./gaugefields/SUN_generator.jl")
    include("./output/verboseprint.jl")
    include("./fermions/cgmethod.jl")

    

    include("./autostaples/wilsonloops.jl")
    
    include("./system/system_parameters.jl")
    #include("parallel.jl")
    #include("site.jl")
    include("./system/rand.jl")
    include("./actions/actions.jl")
    include("./gaugefields/gaugefields.jl")
    
    #include("gaugefields.jl")
    include("./fermions/fermionfields.jl")
    include("./liealgebra/liealgebrafields.jl")

    include("./rationalapprox/rationalapprox.jl")

    
    include("./fermions/clover.jl")
    
    include("./fermions/diracoperator.jl")

    include("./output/io.jl")
    include("./output/ildg_format.jl")

    

    include("./system/LTK_universe.jl")
    include("./gaugefields/smearing.jl")


    include("./output/print_config.jl")


    
    
    #include("cg.jl")


    include("./measurements/measurements.jl")
    include("./heatbath/heatbath.jl")
    include("./md/md.jl")
    include("./system/wizard.jl")

    include("./SLMC/SLMC.jl")



    include("./system/mainrun.jl")
    include("./output/analyze.jl")

    
    

    
    
    function __init__()
        @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin 
            include("./demo/demo.jl")
            import .Demo:demo
            export demo
            
            import .Analyze:plot_plaquette,plot_polyakov,plot_plaq_and_poly
            export plot_plaquette,plot_polyakov,plot_plaq_and_poly
        end
    end
    
    
    

    import .LTK_universe:Universe,show_parameters,make_WdagWmatrix,calc_Action,set_β!,set_βs!
    import .Actions:Setup_Gauge_action,Setup_Fermi_action,GaugeActionParam_autogenerator
    import .Measurements:calc_plaquette,measure_correlator,Measurement,calc_polyakovloop,measure_chiral_cond,calc_topological_charge,
                measurements,Measurement_set
    import  .MD:md_initialize!,MD_parameters_standard,md!,metropolis_update!,construct_MD_parameters
    import .System_parameters:Params,print_parameters,parameterloading,Params_set#,parameterloading2
    import .Print_config:write_config
    import .Smearing:gradientflow!
    import .ILDG_format:ILDG,load_gaugefield
    import .Heatbath:heatbath!
    import .Wilsonloops:make_plaq
    import .IOmodule:saveU,loadU,loadU!
    import .Wizard:run_wizard
    import .Mainrun:run_LQCD
    import .RationalApprox:calc_exactvalue,calc_Anϕ,calc_det
    #,run_LQCD!

    import .Analyze:analyze,get_plaquette,get_polyakov,get_trjs
    
    
    #import .Fermionfields:make_WdagWmatrix
    

    export Setup_Gauge_action,Setup_Fermi_action,GaugeActionParam_autogenerator
    export Universe,set_β!,set_βs!
    export calc_plaquette,calc_polyakovloop,calc_topological_charge
    export md_initialize!,MD_parameters_standard,md!,metropolis_update!,construct_MD_parameters
    export show_parameters
    export Params,print_parameters,parameterloading,Params_set#,parameterloading2
    export measure_correlator,measure_chiral_cond,Measurement,measurements,Measurement_set
    export gradientflow!
    export ILDG,load_gaugefield
    export make_WdagWmatrix
    export heatbath!
    export make_plaq
    export calc_Action
    export calc_topological_charge
    export saveU,loadU,loadU!
    export run_LQCD,run_LQCD!

    export write_config
    export run_wizard
    export analyze,get_plaquette,get_polyakov,get_trjs





end
