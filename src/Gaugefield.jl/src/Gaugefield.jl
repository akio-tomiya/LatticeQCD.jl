module Gaugefield
    #include("../../Wilsonloops/src/Wilsonloops.jl")
    using Requires

    
    include("./output/verboseprint.jl")
    include("./SUN_generator.jl")
    include("./autostaples/wilsonloops.jl")
    include("./AbstractGaugefields.jl") 
    include("./output/io.jl")
    include("./output/ildg_format.jl")
    include("./autostaples/Loops.jl")
    include("./smearing/Abstractsmearing.jl")

    function __init__()
        @require MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195" begin   
            import .AbstractGaugefields_module:identityGaugefields_4D_wing_mpi,
                    Gaugefields_4D_wing_mpi,calc_rank_and_indices,barrier,comm,
                    setvalue!
        end
    end



    import Wilsonloop:loops_staple_prime
    

    
    import .AbstractGaugefields_module:AbstractGaugefields,identitymatrix,Abstractfields,
                                        shift_U,construct_staple!,set_wing_U!,
                                        calculate_Plaquette,substitute_U!,calculate_Polyakov_loop,construct_gauges,
                                        Gaugefields_4D_wing,
                                        identityGaugefields_4D_wing,
                                        add_force!,exp_aF_U!,clear_U!,add_U!,exptU!,
                                        Traceless_antihermitian!,Traceless_antihermitian,Generator,
                                        Staggered_Gaugefields,staggered_U,
                                        Traceless_antihermitian_add!,
                                        IdentityGauges,RandomGauges,Oneinstanton,
                                        construct_Λmatrix_forSTOUT!,
                                        evaluate_gaugelinks_evenodd!,
                                        map_U!

                                        
    import .Loops_module:Loops,evaluate_loops,calc_large_wilson_loop!,evaluate_loops!
    import .Wilsonloops_module:Wilson_loop_set,make_staples,Wilson_loop_set,
                make_cloverloops,Tensor_derivative_set, make_loops,
                make_plaq_staple,make_links,make_plaq,
                make_loopforactions,make_plaqloops,make_rectloops,make_polyakovloops,
                make_plaq_staple_prime,
                calc_coordinate,make_plaq_staple_prime,calc_shift,
                Tensor_wilson_lines_set,Tensor_wilson_lines,Tensor_derivative_set,
                get_leftstartposition,get_rightstartposition,Wilson_loop,calc_loopset_μν_name,
                make_originalactions_fromloops
    import .AbstractGaugefields_module:TA_Gaugefields,initialize_TA_Gaugefields
    import .Abstractsmearing_module:Abstractsmearing,Nosmearing,Stoutsmearing,calc_smearedU,
            construct_smearing,Gradientflow,get_tempG,flow!,get_eps,back_prop,CovNeuralnet
    #import .Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3
    import .Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3,println_verbose2,println_verbose1,
    print_verbose1,print_verbose2,print_verbose3
    import .ILDG_format:ILDG,load_gaugefield,load_gaugefield!,save_binarydata
    import .IOmodule:saveU,loadU,loadU!
    import .SUN_generator:Generator
    
end