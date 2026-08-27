module LatticeQCD
using Requires

const LatticeQCDversion = pkgversion(LatticeQCD)

include("./communication.jl")
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
include("./system/lqcd_config.jl")
include("./system/simulation.jl")
include("./measurements/measurement_plan.jl")
include("./system/simulation_session.jl")
include("./system/simulation_input.jl")
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
import .LQCDCommunication: get_myrank, get_nprocs, set_PEs, get_PEs
import .Simpleprint: println_rank0
import .LQCDConfig_module:
    LatticeConfig,
    AbstractGaugeInitializationConfig,
    ColdStartConfig,
    HotStartConfig,
    FileStartConfig,
    InstantonConfig,
    EmbeddedInstantonConfig,
    GaugeConfig,
    GaugeActionTermConfig,
    GaugeActionConfig,
    AbstractDiracOperatorConfig,
    WilsonDiracConfig,
    WilsonCloverDiracConfig,
    StaggeredDiracConfig,
    HISQDiracConfig,
    DomainwallDiracConfig,
    MobiusDomainwallDiracConfig,
    FermionSolverConfig,
    AbstractFermionSmearingConfig,
    NoFermionSmearingConfig,
    StoutFermionSmearingConfig,
    FermionActionConfig,
    AbstractUpdateConfig,
    AbstractConfigurationSourceConfig,
    DirectorySourceConfig,
    ManifestSourceConfig,
    ConfigurationSequenceConfig,
    AbstractMDIntegratorConfig,
    QPQConfig,
    PQPConfig,
    ForceGroupConfig,
    LeapfrogConfig,
    SextonWeingartenConfig,
    MDConfig,
    RandomStreamConfig,
    GaussianMomentumConfig,
    RankZeroMetropolisConfig,
    PseudofermionRefreshConfig,
    HMCConfig,
    SLHMCConfig,
    HeatbathConfig,
    update_config,
    md_trajectory_length,
    show_config,
    LQCDConfig
import .Simulation_module:
    AbstractConfiguration,
    GaugeConfiguration,
    GaugeFermionConfiguration,
    GaugefieldsEnvironment,
    HMCUpdater,
    HMCState,
    HMCUpdateResult,
    SLHMCUpdater,
    SLHMCUpdateResult,
    PseudofermionRefreshRuntime,
    HeatbathUpdater,
    HeatbathState,
    HeatbathUpdateResult,
    ConfigurationSequenceUpdater,
    ConfigurationSequenceState,
    ConfigurationSequenceUpdateResult,
    Simulation,
    build_configuration,
    build_gauge_action,
    build_fermion_action,
    build_simulation,
    copy_configuration!,
    save_configuration,
    load_configuration!,
    metropolis_rule,
    apply_metropolis_decision!,
    has_next_configuration,
    current_configuration_path,
    update!
import .MeasurementPlan_module:
    AbstractObservableConfig,
    PlaquetteObservableConfig,
    PolyakovLoopObservableConfig,
    TopologicalChargeObservableConfig,
    WilsonLoopObservableConfig,
    EnergyDensityObservableConfig,
    AbstractMeasurementFermionConfig,
    WilsonMeasurementFermionConfig,
    StaggeredMeasurementFermionConfig,
    MeasurementSolverConfig,
    AbstractMeasurementSmearingConfig,
    NoMeasurementSmearingConfig,
    StoutMeasurementSmearingConfig,
    ChiralCondensateObservableConfig,
    PionCorrelatorObservableConfig,
    PeriodicSchedule,
    is_due,
    ScheduledMeasurementConfig,
    MeasurementPlan,
    NoGradientFlowMeasurementConfig,
    GradientFlowMeasurementConfig,
    MeasurementProgram,
    measurement_config,
    measurement_plan,
    gradient_flow_measurement_config,
    measurement_program,
    observable_name,
    ScheduledMeasurementRuntime,
    MeasurementRuntime,
    NoGradientFlowMeasurementRuntime,
    GradientFlowMeasurementRuntime,
    MeasurementProgramRuntime,
    build_measurement_runtime,
    MeasurementPoint,
    MeasurementRecord,
    measure_now!,
    measure_due!,
    NoMeasurementSink,
    FunctionMeasurementSink,
    SimulationRunner,
    SimulationStepResult,
    build_simulation_runner,
    step!,
    run!
import .SimulationSession_module:
    SimulationSchedule,
    NoConfigurationOutput,
    JLD2ConfigurationOutput,
    BridgeTextConfigurationOutput,
    ILDGConfigurationOutput,
    NoCheckpointOutput,
    JLD2CheckpointOutput,
    OutputConfig,
    SimulationSpec,
    ValidationIssue,
    validate,
    simulation_schedule,
    output_config,
    AbstractSimulationEvent,
    RunStarted,
    ThermalizationStepFinished,
    TrajectoryFinished,
    MeasurementFinished,
    ConfigurationSaved,
    CheckpointSaved,
    CheckpointLoaded,
    ConfigurationLoaded,
    SimulationRunSummary,
    RunStopped,
    RunFinished,
    RunFailed,
    NoSimulationEventSink,
    ConsoleSimulationEventSink,
    FunctionSimulationEventSink,
    RecordingSimulationEventSink,
    CompositeSimulationEventSink,
    SimulationSessionState,
    SimulationSession,
    SessionStepResult,
    build_simulation_session,
    emit_event!,
    configuration_output_path,
    checkpoint_output_path,
    save_checkpoint,
    load_checkpoint!,
    is_finished,
    is_running,
    request_stop!,
    run_summary
import .SimulationInput_module:
    legacy_simulation_values,
    simulation_spec_from_legacy_toml,
    SIMULATION_SPEC_FORMAT,
    SIMULATION_SPEC_SCHEMA_VERSION,
    simulation_spec_dictionary,
    simulation_spec_from_toml,
    parse_simulation_spec,
    load_simulation_spec,
    write_simulation_spec
#import .Print_config: write_config
#import .Smearing: gradientflow!
#import .ILDG_format: ILDG, load_gaugefield
#import .Heatbath: heatbath!
#import .Wilsonloops: make_plaq
#import .IOmodule: saveU, loadU, loadU!
import .Wizard: run_wizard, run_wizardv2, run_wizard_legacy
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
export get_myrank, get_nprocs, println_rank0, set_PEs, get_PEs
export LatticeConfig,
    AbstractGaugeInitializationConfig,
    ColdStartConfig,
    HotStartConfig,
    FileStartConfig,
    InstantonConfig,
    EmbeddedInstantonConfig,
    GaugeConfig,
    GaugeActionTermConfig,
    GaugeActionConfig,
    AbstractDiracOperatorConfig,
    WilsonDiracConfig,
    WilsonCloverDiracConfig,
    StaggeredDiracConfig,
    HISQDiracConfig,
    DomainwallDiracConfig,
    MobiusDomainwallDiracConfig,
    FermionSolverConfig,
    AbstractFermionSmearingConfig,
    NoFermionSmearingConfig,
    StoutFermionSmearingConfig,
    FermionActionConfig,
    AbstractUpdateConfig,
    AbstractConfigurationSourceConfig,
    DirectorySourceConfig,
    ManifestSourceConfig,
    ConfigurationSequenceConfig,
    AbstractMDIntegratorConfig,
    QPQConfig,
    PQPConfig,
    ForceGroupConfig,
    LeapfrogConfig,
    SextonWeingartenConfig,
    MDConfig,
    RandomStreamConfig,
    GaussianMomentumConfig,
    RankZeroMetropolisConfig,
    PseudofermionRefreshConfig,
    HMCConfig,
    SLHMCConfig,
    HeatbathConfig,
    update_config,
    md_trajectory_length,
    show_config,
    LQCDConfig
export AbstractConfiguration,
    GaugeConfiguration,
    GaugeFermionConfiguration,
    GaugefieldsEnvironment,
    HMCUpdater,
    HMCState,
    HMCUpdateResult,
    SLHMCUpdater,
    SLHMCUpdateResult,
    PseudofermionRefreshRuntime,
    HeatbathUpdater,
    HeatbathState,
    HeatbathUpdateResult,
    ConfigurationSequenceUpdater,
    ConfigurationSequenceState,
    ConfigurationSequenceUpdateResult,
    Simulation,
    build_configuration,
    build_gauge_action,
    build_fermion_action,
    build_simulation,
    copy_configuration!,
    save_configuration,
    load_configuration!,
    metropolis_rule,
    apply_metropolis_decision!,
    has_next_configuration,
    current_configuration_path,
    update!
export AbstractObservableConfig,
    PlaquetteObservableConfig,
    PolyakovLoopObservableConfig,
    TopologicalChargeObservableConfig,
    WilsonLoopObservableConfig,
    EnergyDensityObservableConfig,
    AbstractMeasurementFermionConfig,
    WilsonMeasurementFermionConfig,
    StaggeredMeasurementFermionConfig,
    MeasurementSolverConfig,
    AbstractMeasurementSmearingConfig,
    NoMeasurementSmearingConfig,
    StoutMeasurementSmearingConfig,
    ChiralCondensateObservableConfig,
    PionCorrelatorObservableConfig,
    PeriodicSchedule,
    is_due,
    ScheduledMeasurementConfig,
    MeasurementPlan,
    NoGradientFlowMeasurementConfig,
    GradientFlowMeasurementConfig,
    MeasurementProgram,
    measurement_config,
    measurement_plan,
    gradient_flow_measurement_config,
    measurement_program,
    observable_name,
    ScheduledMeasurementRuntime,
    MeasurementRuntime,
    NoGradientFlowMeasurementRuntime,
    GradientFlowMeasurementRuntime,
    MeasurementProgramRuntime,
    build_measurement_runtime,
    MeasurementPoint,
    MeasurementRecord,
    measure_now!,
    measure_due!,
    NoMeasurementSink,
    FunctionMeasurementSink,
    SimulationRunner,
    SimulationStepResult,
    build_simulation_runner,
    step!,
    run!
export SimulationSchedule,
    NoConfigurationOutput,
    JLD2ConfigurationOutput,
    BridgeTextConfigurationOutput,
    ILDGConfigurationOutput,
    NoCheckpointOutput,
    JLD2CheckpointOutput,
    OutputConfig,
    SimulationSpec,
    ValidationIssue,
    validate,
    simulation_schedule,
    output_config,
    AbstractSimulationEvent,
    RunStarted,
    ThermalizationStepFinished,
    TrajectoryFinished,
    MeasurementFinished,
    ConfigurationSaved,
    CheckpointSaved,
    CheckpointLoaded,
    ConfigurationLoaded,
    SimulationRunSummary,
    RunStopped,
    RunFinished,
    RunFailed,
    NoSimulationEventSink,
    ConsoleSimulationEventSink,
    FunctionSimulationEventSink,
    RecordingSimulationEventSink,
    CompositeSimulationEventSink,
    SimulationSessionState,
    SimulationSession,
    SessionStepResult,
    build_simulation_session,
    emit_event!,
    configuration_output_path,
    checkpoint_output_path,
    save_checkpoint,
    load_checkpoint!,
    is_finished,
    is_running,
    request_stop!,
    run_summary
export legacy_simulation_values,
    simulation_spec_from_legacy_toml,
    SIMULATION_SPEC_FORMAT,
    SIMULATION_SPEC_SCHEMA_VERSION,
    simulation_spec_dictionary,
    simulation_spec_from_toml,
    parse_simulation_spec,
    load_simulation_spec,
    write_simulation_spec
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
export run_wizard, run_wizardv2, run_wizard_legacy
export analyze,
    get_plaquette, get_polyakov, get_plaquette_average, get_polyakov_average, get_trjs

export run_LQCD_file



end
