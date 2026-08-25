using LinearAlgebra
using LatticeQCD
import Gaugefields

BLAS.set_num_threads(1)

const MODE = get(ENV, "LQCD_BENCH_MODE", "typed_lm")
const SAMPLES = parse(Int, get(ENV, "LQCD_BENCH_SAMPLES", "9"))
const WARMUP = parse(Int, get(ENV, "LQCD_BENCH_WARMUP", "3"))
const SAMPLE_SECONDS = parse(
    Float64,
    get(ENV, "LQCD_BENCH_SAMPLE_SECONDS", "0.35"),
)

struct BenchmarkCase
    name::String
    lattice::NTuple{4,Int}
    colors::Int
    beta::Float64
    kind::Symbol
    even_odd::Bool
    overrelaxation_steps::Int
    md_steps::Int
    md_step_size::Float64
end

const CASES = (
    BenchmarkCase(
        "heatbath_eo_su2_8x8x8x8",
        (8, 8, 8, 8),
        2,
        1.9,
        :heatbath,
        true,
        0,
        0,
        0.0,
    ),
    BenchmarkCase(
        "heatbath_eo_su3_8x8x8x8",
        (8, 8, 8, 8),
        3,
        5.7,
        :heatbath,
        true,
        0,
        0,
        0.0,
    ),
    BenchmarkCase(
        "heatbath_general_su3_8x8x8x8",
        (8, 8, 8, 8),
        3,
        5.7,
        :heatbath,
        false,
        0,
        0,
        0.0,
    ),
    BenchmarkCase(
        "heatbath_eo_su3_or3_8x8x8x8",
        (8, 8, 8, 8),
        3,
        5.7,
        :heatbath,
        true,
        3,
        0,
        0.0,
    ),
    BenchmarkCase(
        "gauge_hmc_su2_6x6x6x6_md10",
        (6, 6, 6, 6),
        2,
        1.9,
        :hmc,
        false,
        0,
        10,
        0.02,
    ),
    BenchmarkCase(
        "gauge_hmc_su3_6x6x6x6_md10",
        (6, 6, 6, 6),
        3,
        5.7,
        :hmc,
        false,
        0,
        10,
        0.02,
    ),
)

function legacy_fixture(case::BenchmarkCase)
    gauge = Gaugefields.Initialize_Gaugefields(
        case.colors,
        1,
        case.lattice...;
        condition="cold",
        verbose_level=0,
    )
    action = Gaugefields.GaugeAction(gauge)
    loops = Gaugefields.make_loops_fromname("plaquette"; Dim=4)
    append!(loops, loops')
    push!(action, case.beta / 2, loops)

    updater = if case.kind == :heatbath
        LatticeQCD.AbstractUpdate_module.Updatemethod(
            gauge,
            action,
            "Heatbath",
            true;
            isevenodd=case.even_odd,
            β=case.beta,
            ITERATION_MAX=100_000,
            numOR=case.overrelaxation_steps,
            useOR=case.overrelaxation_steps > 0,
        )
    else
        LatticeQCD.AbstractUpdate_module.Updatemethod(
            gauge,
            action,
            "HMC",
            true,
            case.md_step_size,
            case.md_steps;
            QPQ=true,
        )
    end
    run_update!() = LatticeQCD.AbstractUpdate_module.update!(updater, gauge)
    return run_update!, gauge
end

function typed_fixture(case::BenchmarkCase, backend)
    lattice = LatticeConfig(case.lattice)
    gauge = GaugeConfig(case.colors, 1, "cold", nothing, 0x1234)
    action = GaugeActionConfig(
        GaugeActionTermConfig(:gauge_plaquette, "plaquette", case.beta),
    )
    update = if case.kind == :heatbath
        HeatbathConfig(
            case.even_odd,
            100_000,
            case.overrelaxation_steps,
            RandomStreamConfig(0x5678, :heatbath),
        )
    else
        integrator = LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge))
        md = MDConfig(case.md_step_size, case.md_steps, integrator)
        momentum = GaussianMomentumConfig(
            1.0,
            RandomStreamConfig(0x5678, :momentum),
        )
        acceptance = RankZeroMetropolisConfig(
            RandomStreamConfig(0x9abc, :metropolis),
        )
        HMCConfig(md, momentum, acceptance)
    end
    input = LQCDConfig(lattice, gauge, action, update)
    environment = GaugefieldsEnvironment(backend=backend, verbose=0)
    simulation = build_simulation(input, environment)
    run_update!() = LatticeQCD.update!(simulation)
    return run_update!, simulation.configuration.gauge
end

function make_fixture(case::BenchmarkCase)
    if MODE == "legacy"
        return legacy_fixture(case)
    elseif MODE == "typed_lm"
        return typed_fixture(case, Gaugefields.LatticeMatricesBackend())
    elseif MODE == "typed_legacy"
        return typed_fixture(case, Gaugefields.LegacyBackend())
    end
    error("LQCD_BENCH_MODE must be legacy, typed_lm, or typed_legacy")
end

function median_value(values)
    ordered = sort(values)
    middle = length(ordered) ÷ 2
    return isodd(length(ordered)) ?
           ordered[middle + 1] :
           (ordered[middle] + ordered[middle + 1]) / 2
end

function percentile(values, fraction)
    ordered = sort(values)
    index = clamp(ceil(Int, fraction * length(ordered)), 1, length(ordered))
    return ordered[index]
end

function calibrated_repetitions(run_update!)
    elapsed = @elapsed run_update!()
    return clamp(ceil(Int, SAMPLE_SECONDS / max(elapsed, eps())), 1, 100)
end

function normalized_plaquette(gauge, case::BenchmarkCase)
    raw = Gaugefields.calculate_Plaquette(
        gauge,
        similar(gauge[1]),
        similar(gauge[1]),
    )
    return real(raw) / (6 * prod(case.lattice) * case.colors)
end

function run_case(case::BenchmarkCase)
    run_update!, gauge = make_fixture(case)
    for _ in 1:WARMUP
        run_update!()
    end
    repetitions = calibrated_repetitions(run_update!)
    times = Float64[]
    for _ in 1:SAMPLES
        GC.gc()
        start = time_ns()
        for _ in 1:repetitions
            run_update!()
        end
        push!(times, (time_ns() - start) / 1.0e9 / repetitions)
    end
    allocation = @allocated run_update!()
    println(
        "RESULT",
        " mode=", MODE,
        " case=", case.name,
        " samples=", SAMPLES,
        " repetitions=", repetitions,
        " median_s=", median_value(times),
        " min_s=", minimum(times),
        " p90_s=", percentile(times, 0.9),
        " allocation_bytes=", allocation,
        " plaquette=", normalized_plaquette(gauge, case),
    )
    flush(stdout)
    return nothing
end

println(
    "ENV",
    " mode=", MODE,
    " julia=", VERSION,
    " latticeqcd=", Base.pkgversion(LatticeQCD),
    " gaugefields=", Base.pkgversion(Gaugefields),
    " threads=", Threads.nthreads(),
    " blas_threads=", BLAS.get_num_threads(),
)
for case in CASES
    if MODE == "typed_legacy" && case.kind == :hmc
        println(
            "SKIP",
            " mode=", MODE,
            " case=", case.name,
            " reason=legacy_backend_cannot_honor_explicit_momentum_seed",
        )
    else
        run_case(case)
    end
end
