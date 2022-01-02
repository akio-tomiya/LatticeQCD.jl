module Fermionfield_LQCD
    import ..Gaugefield
    include("./cgmethod.jl")

    include("./rhmc/AlgRemez.jl")
    include("./rhmc/rhmc.jl")
    include("./Fermionaction.jl")
    include("./AbstractFermionfields.jl")
    include("./clover.jl")
    
    include("./WilsonFermion/clover.jl")
    include("./Dirac_operators.jl")
    include("./rationalapprox/rationalapprox.jl")
    
    import .FermionAction_module:FermiActionParam_Wilson,FermiActionParam_Staggered,FermiActionParam_WilsonClover,
    FermiActionParam
    import .AbstractFermionfields_module:Fermionfields,Z4_distribution_fermion!,
            set_wing_fermion!,clear_fermion!,Z4_distribution_fermion!,
            add_fermion!,get_origin,AbstractFermionfields,StaggeredFermion,
            gauss_distribution_fermion!,shift_fermion
    import .Dirac_operators_module:Dirac_operator,DdagD_operator,Adjoint_Wilson_operator,Wilson_operator,
                                    WilsonClover_operator

    import .CGmethods:bicg,cg,shiftedcg
    import .Rhmc:get_order,get_β,get_α,get_α0,get_β_inverse,get_α_inverse,get_α0_inverse,RHMC
    import .RationalApprox:AlgRemez_coeffs,calc_exactvalue,calc_Anϕ
    import .Clover:Make_CloverFμν!,dSclover!
    #import .Clover:Make_CloverFμν
    
end