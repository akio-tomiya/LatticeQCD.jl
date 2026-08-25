using LatticeQCD
using TOML

const PS = LatticeQCD.Parameter_structs
const WIZARD = LatticeQCD.Wizard
const OUTPUT_DIRECTORY = joinpath(@__DIR__, "wizard_fermion")

const FERMION_VARIANTS = (
    (id="wilson", operator="Wilson", flavors=2),
    (id="wilson_clover", operator="WilsonClover", flavors=2),
    (id="staggered_nf1", operator="Staggered", flavors=1),
    (id="staggered_nf2", operator="Staggered", flavors=2),
    (id="staggered_nf3", operator="Staggered", flavors=3),
    (id="staggered_nf4", operator="Staggered", flavors=4),
    (id="staggered_nf8", operator="Staggered", flavors=8),
    (id="hisq_nf1", operator="HISQ", flavors=1),
    (id="hisq_nf2", operator="HISQ", flavors=2),
    (id="hisq_nf3", operator="HISQ", flavors=3),
    (id="hisq_nf4", operator="HISQ", flavors=4),
    (id="hisq_nf8", operator="HISQ", flavors=8),
    (id="domainwall", operator="Domainwall", flavors=2),
    (id="mobius_domainwall", operator="MobiusDomainwall", flavors=2),
)

function fermion_parameters(variant)
    if variant.operator == "Wilson"
        return PS.Wilson_parameters()
    elseif variant.operator == "WilsonClover"
        return PS.Wilson_parameters(
            Dirac_operator="WilsonClover",
            hasclover=true,
            Clover_coefficient=1.5612,
        )
    elseif variant.operator == "Staggered"
        return PS.Staggered_parameters(Nf=variant.flavors)
    elseif variant.operator == "HISQ"
        return PS.HISQ_parameters(
            Nf=variant.flavors,
            naik_epsilon=-0.083,
        )
    elseif variant.operator == "Domainwall"
        return PS.Domainwall_parameters()
    elseif variant.operator == "MobiusDomainwall"
        return PS.MobiusDomainwall_parameters(b=2.0, c=1.0)
    end
    error("unsupported Wizard fermion $(variant.operator)")
end

function wizard_sections(
    colors::Int,
    variant,
    stout::Bool,
    sexton::Bool,
    update_method::String,
)
    beta = colors == 3 ? 5.7 : 2.7
    physical = PS.Print_Physical_parameters(
        L=variant.operator == "HISQ" ? [4, 4, 4, 4] : [2, 2, 2, 2],
        NC=colors,
        β=beta,
        Nthermalization=0,
        Nsteps=2,
        initial="hot",
        update_method=update_method,
        Nwing=variant.operator == "HISQ" ? 3 : 1,
    )
    fermions = PS.Print_Fermions_parameters(
        quench=false,
        Dirac_operator=variant.operator,
        smearing_for_fermion=stout ? "stout" : "nothing",
        stout_numlayers=stout ? 1 : nothing,
        stout_ρ=stout ? [0.1] : nothing,
        stout_loops=stout ? ["plaquette"] : nothing,
    )
    control = PS.Print_System_control_parameters(
        randomseed=111,
        verboselevel=0,
        log_dir="./logs",
        logfile="wizard-fermion.log",
        measurement_basedir="./measurements",
        measurement_dir="wizard-fermion",
        saveU_format="nothing",
    )
    hmc = PS.Print_HMCrelated_parameters(
        Δτ=0.001,
        MDsteps=1,
        SextonWeingargten=sexton,
        N_SextonWeingargten=2,
        eps=1e-10,
        MaxCGstep=2_000,
        QPQ=true,
    )
    measurement = PS.Measurement_parameterset(
        measurement_methods=PS.Measurement_parameters[PS.Plaq_parameters()],
    )
    gradient = PS.Print_Gradientflow_parameters(hasgradientflow=false)
    gradient_measurement = PS.Measurement_parameterset()
    return (;
        physical,
        fermions,
        fermion_parameters=fermion_parameters(variant),
        control,
        hmc,
        measurement,
        gradient,
        gradient_measurement,
        slhmc_beta=beta - 0.5,
    )
end

function printable_dictionary(sections)
    return WIZARD.wizard_parameter_dictionary(
        sections.physical,
        sections.fermions,
        sections.fermion_parameters,
        sections.control,
        sections.hmc,
        sections.measurement,
        sections.gradient,
        sections.gradient_measurement,
        slhmc_beta=sections.slhmc_beta,
    )
end

function validate_parameter_file(
    path::String,
    colors,
    variant,
    stout,
    sexton,
    update_method,
)
    parameters = TOML.parsefile(path)
    physical = parameters["Physical setting"]
    fermions = parameters["Physical setting(fermions)"]
    hmc = parameters["HMC related"]
    measurements = parameters["Measurement set"]["measurement_methods"]
    gradient = parameters["gradientflow_measurements"]

    get(physical, "NC", 3) == colors || error("$path has the wrong NC")
    physical["initial"] == "hot" || error("$path is not a hot start")
    physical["update_method"] == update_method || error(
        "$path has the wrong update method",
    )
    if update_method == "SLHMC"
        sections_beta = (colors == 3 ? 5.7 : 2.7) - 0.5
        parameters["SLHMC related"]["βeff"] == sections_beta || error(
            "$path has the wrong SLHMC effective beta",
        )
    end
    fermions["quench"] == false || error("$path is quenched")
    fermions["Dirac_operator"] == variant.operator || error(
        "$path has the wrong Dirac operator",
    )
    if variant.operator in ("Staggered", "HISQ")
        get(fermions, "Nf", 4) == variant.flavors || error(
            "$path has the wrong Nf",
        )
    end
    if variant.operator == "HISQ"
        get(physical, "Nwing", 0) == 3 || error(
            "$path has the wrong HISQ halo width",
        )
        fermions["naik_epsilon"] == -0.083 || error(
            "$path has the wrong Naik correction",
        )
        stout && error("$path applies unsupported outer stout smearing to HISQ")
    end
    if variant.operator == "WilsonClover"
        fermions["Clover_coefficient"] == 1.5612 || error(
            "$path has the wrong cSW",
        )
    end
    if variant.operator == "MobiusDomainwall"
        fermions["b"] == 2.0 || error("$path has the wrong Möbius b")
        fermions["c"] == 1.0 || error("$path has the wrong Möbius c")
    end
    (get(fermions, "smearing_for_fermion", "nothing") == "stout") == stout || error(
        "$path has the wrong smearing setting",
    )
    get(hmc, "SextonWeingargten", false) == sexton || error(
        "$path has the wrong Sexton-Weingarten setting",
    )
    Set(keys(measurements)) == Set(["Plaquette"]) || error(
        "$path does not contain exactly one Plaquette measurement",
    )
    gradient["hasgradientflow"] == false || error("$path enables flow")
    isempty(gradient["measurements_for_flow"]) || error(
        "$path contains a gradient-flow measurement",
    )
    return nothing
end

function write_parameter_file(
    filename,
    sections,
    colors,
    variant,
    stout,
    sexton,
    update_method,
)
    path = joinpath(OUTPUT_DIRECTORY, filename)
    open(path, "w") do io
        TOML.print(io, printable_dictionary(sections))
    end
    validate_parameter_file(
        path,
        colors,
        variant,
        stout,
        sexton,
        update_method,
    )
    return path
end

function main()
    mkpath(OUTPUT_DIRECTORY)
    filenames = String[]
    for variant in FERMION_VARIANTS
        colors_values = variant.operator == "HISQ" ? (3,) : (2, 3)
        stout_values = variant.operator == "HISQ" ? (false,) : (false, true)
        for colors in colors_values
            for stout in stout_values
                for sexton in (false, true)
                    for update_method in ("HMC", "SLHMC")
                        parts = String["su$colors"]
                        update_method == "SLHMC" && push!(parts, "slhmc")
                        append!(parts, [
                            variant.id,
                            stout ? "stout" : "thin",
                            sexton ? "sw" : "leapfrog",
                        ])
                        filename = join(parts, "_") * ".toml"
                        sections = wizard_sections(
                            colors,
                            variant,
                            stout,
                            sexton,
                            update_method,
                        )
                        write_parameter_file(
                            filename,
                            sections,
                            colors,
                            variant,
                            stout,
                            sexton,
                            update_method,
                        )
                        push!(filenames, filename)
                    end
                end
            end
        end
    end

    length(filenames) == 164 || error(
        "expected 164 Wizard fermion cases, generated $(length(filenames))",
    )
    existing = sort(filter(
        filename -> endswith(filename, ".toml"),
        readdir(OUTPUT_DIRECTORY),
    ))
    sort(filenames) == existing || error(
        "the output directory contains stale or missing TOML files",
    )
    println("generated $(length(filenames)) fermion Wizard parameter files")
    println("output: $OUTPUT_DIRECTORY")
    return nothing
end

main()
