using LatticeQCD
import Gaugefields
import LatticeDiracOperators

const REFERENCE_KIND = get(ENV, "LQCD_REFERENCE_FERMION", "HISQ")
const IS_HISQ = REFERENCE_KIND == "HISQ"
const IS_MOBIUS = REFERENCE_KIND == "MobiusDomainwall"
IS_HISQ || IS_MOBIUS || error(
    "LQCD_REFERENCE_FERMION must be HISQ or MobiusDomainwall",
)
const LATTICE = IS_HISQ ? (4, 4, 4, 4) : (2, 2, 2, 2)
const COLORS = IS_HISQ ? 3 : 2

function force_fingerprint(force, mu)
    values = ComplexF64[]
    for t in 1:LATTICE[4], z in 1:LATTICE[3]
        for y in 1:LATTICE[2], x in 1:LATTICE[1]
            for column in 1:COLORS, row in 1:COLORS
                push!(values, force[mu][row, column, x, y, z, t])
            end
        end
    end
    return (
        sum(real, values),
        sum(imag, values),
        sum(i * real(values[i]) for i in eachindex(values)),
        sum(i * imag(values[i]) for i in eachindex(values)),
        sum(abs2, values),
    )
end

lattice = LatticeConfig(LATTICE)
gauge_config = GaugeConfig(
    COLORS,
    IS_HISQ ? 3 : 1,
    "cold",
    nothing,
    0x1234,
)
gauge_action = GaugeActionConfig(
    GaugeActionTermConfig(:gauge_plaquette, "plaquette", 1.9),
)
operator = IS_HISQ ?
    HISQDiracConfig(0.5, -0.083) :
    MobiusDomainwallDiracConfig(0.1, -1.0, 2, 2.0, 1.0)
solver = FermionSolverConfig(1e-10, 2_000, 0, (1, 1, 1, -1))
fermion_config = FermionActionConfig(
    :fermion_1,
    operator,
    IS_HISQ ? 4 : 2,
    solver,
    NoFermionSmearingConfig(),
)
md = MDConfig(
    0.001,
    1,
    LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge, :fermion_1)),
)
update = HMCConfig(
    md,
    GaussianMomentumConfig(1.0, RandomStreamConfig(0x5678, :momentum)),
    RankZeroMetropolisConfig(RandomStreamConfig(0x9abc, :metropolis)),
    (PseudofermionRefreshConfig(
        :fermion_1,
        RandomStreamConfig(0xdef0, :pseudofermion),
        1,
    ),),
)
input = LQCDConfig(
    lattice,
    gauge_config,
    gauge_action,
    (fermion_config,),
    update,
)
environment = GaugefieldsEnvironment(
    process_grid=(1, 1, 1, 1),
    communicator=Gaugefields.SerialCommunicator(),
    element_type=ComplexF64,
    verbose=0,
)
configuration = build_configuration(input, environment)
gauge = configuration.gauge
field = configuration.fermions.fermion_1

typed_action = build_fermion_action(fermion_config, gauge, field)
direct_parameters = Dict{String,Any}(
    "Dirac_operator" => REFERENCE_KIND,
    "eps_CG" => 1e-10,
    "MaxCGstep" => 2_000,
    "verbose_level" => 0,
    "boundarycondition" => [1, 1, 1, -1],
)
if IS_HISQ
    direct_parameters["mass"] = 0.5
    direct_parameters["naik_epsilon"] = -0.083
else
    direct_parameters["mass"] = 0.1
    direct_parameters["M"] = -1.0
    direct_parameters["L5"] = 2
    direct_parameters["b"] = 2.0
    direct_parameters["c"] = 1.0
end
direct_dirac = LatticeDiracOperators.Dirac_operator(
    gauge,
    field,
    direct_parameters,
)
direct_action = LatticeDiracOperators.FermiAction(
    direct_dirac,
    Dict{String,Any}("Nf" => (IS_HISQ ? 4 : 2)),
)

provider = LatticeDiracOperators.PseudofermionMDAction(typed_action, field)
noise = similar(field)
LatticeDiracOperators.refresh_pseudofermion!(
    provider,
    gauge,
    noise;
    seed=0x314159,
    sweep=7,
    subgroup=1,
)

for (name, action) in (("LatticeQCD", typed_action), ("direct-LDO", direct_action))
    value = LatticeDiracOperators.evaluate_FermiAction(action, gauge, field)
    force = LatticeDiracOperators.calc_UdSfdU(action, gauge, field)
    println("ACTION implementation=", name, " value=", value)
    for mu in 1:4
        values = force_fingerprint(force, mu)
        println(
            "FORCE implementation=", name,
            " mu=", mu,
            " sum_re=", values[1],
            " sum_im=", values[2],
            " weighted_re=", values[3],
            " weighted_im=", values[4],
            " norm2=", values[5],
        )
    end
end
