module Fermionfield_LQCD
    import ..Gaugefield
    include("./cgmethod.jl")
    include("./rhmc/AlgRemez.jl")
    include("./rhmc/rhmc.jl")
    include("./AbstractFermionfields.jl")
    include("./Fermionaction.jl")
    include("./Dirac_operators.jl")

    import .AbstractFermionfields_module:Fermionfields,Z4_distribution_fermion!,set_wing_fermion!,clear_fermion!,Z4_distribution_fermion!,add_fermion!
    import .Dirac_operators_module:Dirac_operator 
    import .FermionAction_module:FermiActionParam_Wilson,FermiActionParam_Staggered,FermiActionParam_WilsonClover,
                FermiActionParam
    import .CGmethods:bicg,cg
    
end