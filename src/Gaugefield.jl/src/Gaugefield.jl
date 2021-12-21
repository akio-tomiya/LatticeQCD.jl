module Gaugefield
    include("./output/verboseprint.jl")
    include("./SUN_generator.jl")
    include("./autostaples/wilsonloops.jl")
    include("./AbstractGaugefields.jl") 
    include("./autostaples/Loops.jl")
    include("./smearing/Abstractsmearing.jl")
    

    
    import .AbstractGaugefields_module:AbstractGaugefields,identitymatrix,Abstractfields,
                                        shift_U,construct_staple!,set_wing_U!,
                                        calculate_Plaquette,substitute_U!,calculate_Polyakov_loop,construct_gauges,
                                        Gaugefields_4D_wing_mpi,identityGaugefields_4D_wing_mpi,
                                        calc_rank_and_indices,barrier,comm,setvalue!,
                                        Gaugefields_4D_wing,
                                        identityGaugefields_4D_wing,
                                        add_force!,exp_aF_U!,clear_U!,add_U!,
                                        Traceless_antihermitian!,Traceless_antihermitian,Generator,
                                        Staggered_Gaugefields,staggered_U

                                        
    import .Loops_module:Loops,evaluate_loops,calc_large_wilson_loop!
    import .Wilsonloops:Wilson_loop_set,make_staples,Wilson_loop_set,
                make_cloverloops,Tensor_derivative_set, make_loops,
                make_plaq_staple,make_links,make_plaq,
                make_loopforactions,make_plaqloops,make_rectloops,make_polyakovloops,
                make_plaq_staple_prime,
                calc_coordinate,make_plaq_staple_prime,calc_shift,
                Tensor_wilson_lines_set,Tensor_wilson_lines,Tensor_derivative_set,
                get_leftstartposition,get_rightstartposition,Wilson_loop,calc_loopset_μν_name,
                make_originalactions_fromloops
    import .AbstractGaugefields_module:TA_Gaugefields,initialize_TA_Gaugefields
    import .Abstractsmearing_module:Abstractsmearing,Nosmearing,Stoutsmearing,calc_smearedU,construct_smearing,Gradientflow,get_tempG,flow!,get_eps
    #import .Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3
    import .Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3,println_verbose2,println_verbose1,
    print_verbose1,print_verbose2,print_verbose3
    
end