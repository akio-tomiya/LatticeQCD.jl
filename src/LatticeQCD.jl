module LatticeQCD
using Requires

const LatticeQCDversion = pkgversion(LatticeQCD)

include("./mpi/simpleprint.jl")
#include("./SLMC/logdet.jl")
include("./system/parameter_structs.jl")

#include("./rhmc/AlgRemez.jl")
#include("./rhmc/rhmc.jl")

#include("./gaugefields/SUN_generator.jl")
#include("./output/verboseprint.jl")
#include("./fermions/cgmethod.jl")


#include("./autostaples/wilsonloops.jl")

include("./system/transform_oldinputfile.jl")
include("./system/system_parameters.jl")
include("./system/parameters_TOML.jl")

include("./system/universe.jl")
include("./md/AbstractMD.jl")
include("./updates/AbstractUpdate.jl")

include("./measurements/measurement_parameters_set.jl")
include("./measurements/Measurement_set.jl")
#include("./measurements/AbstractMeasurement.jl")

#include("parallel.jl")
#include("site.jl")
#include("./system/rand.jl")
#include("./actions/actions.jl")
#include("./gaugefields/gaugefields.jl")

#include("gaugefields.jl")
#include("./fermions/AbstractFermion.jl")
#include("./fermions/WilsonFermion.jl")
#include("./fermions/DomainwallFermion.jl")
#include("./fermions/StaggeredFermion.jl")
#include("./fermions/fermionfields.jl")
#include("./liealgebra/liealgebrafields.jl")

#include("./rationalapprox/rationalapprox.jl")


#include("./fermions/clover.jl")

#include("./fermions/diracoperator.jl")
#include("./fermions/misc.jl")


#include("./output/io.jl")
#include("./output/ildg_format.jl")
#include("./output/bridge_format.jl")



#include("./system/LTK_universe.jl")
#include("./gaugefields/smearing.jl")


#include("./output/print_config.jl")




#include("cg.jl")


#include("./measurements/measurements.jl")
#include("./heatbath/heatbath.jl")
#include("./md/md.jl")
include("./system/wizard.jl")

#include("./SLMC/SLMC.jl")



#include("./system/mainrun.jl")

#include("./output/analyze.jl")






function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        include("./demo/demo.jl")
        import .Demo: demo
        export demo

        #import .Analyze: plot_plaquette, plot_polyakov, plot_plaq_and_poly
        #export plot_plaquette, plot_polyakov, plot_plaq_and_poly
    end

    @require MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195" begin
        include("./mpi/mpimodule.jl")
        import .MPImodules: get_myrank, get_nprocs, println_rank0, set_PEs, get_PEs
        export get_myrank, get_nprocs, println_rank0, set_PEs, get_PEs
    end
end

include("./system/lqcd.jl")


#import .LTK_universe:
#    Universe, show_parameters, make_WdagWmatrix, calc_Action, set_β!, set_βs!
#import .Actions: Setup_Gauge_action, Setup_Fermi_action, GaugeActionParam_autogenerator
#import .Measurements:
#    calc_plaquette,
#    measure_correlator,
#    Measurement,
#    calc_polyakovloop,
#    measure_chiral_cond,
#    calc_topological_charge,
#measurements,
#    Measurement_set
#import .MD:
#    md_initialize!, MD_parameters_standard, md!, metropolis_update!, construct_MD_parameters
import .System_parameters: Params
#import .Print_config: write_config
#import .Smearing: gradientflow!
#import .ILDG_format: ILDG, load_gaugefield
#import .Heatbath: heatbath!
#import .Wilsonloops: make_plaq
#import .IOmodule: saveU, loadU, loadU!
import .Wizard: run_wizard
#import .Mainrun: run_LQCD
#import .RationalApprox: calc_exactvalue, calc_Anϕ, calc_det
#,run_LQCD!


#import .Analyze:
#   analyze,
#    get_plaquette,
#    get_polyakov,
#    get_plaquette_average,
#    get_polyakov_average,
#    get_trjs
import .LQCD: run_LQCD_file, run_LQCD##


#import .Fermionfields:make_WdagWmatrix


#export Setup_Gauge_action, Setup_Fermi_action, GaugeActionParam_autogenerator
#export Universe, set_β!, set_βs!
#export calc_plaquette, calc_polyakovloop, calc_topological_charge
#export md_initialize!,
#    MD_parameters_standard, md!, metropolis_update!, construct_MD_parameters
#export show_parameters
export Params
#export measure_correlator, measure_chiral_cond, Measurement, measurements, Measurement_set
#export gradientflow!
#export ILDG, load_gaugefield
#export make_WdagWmatrix
#export heatbath!
#export make_plaq
#export calc_Action
#export calc_topological_charge
#export saveU, loadU, loadU!
export run_LQCD, run_LQCD!

#export write_config
export run_wizard
export analyze,
    get_plaquette, get_polyakov, get_plaquette_average, get_polyakov_average, get_trjs

export run_LQCD_file



end
