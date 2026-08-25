using LatticeQCD
using TOML

const PS = LatticeQCD.Parameter_structs
const WIZARD = LatticeQCD.Wizard
const OUTPUT_DIRECTORY = joinpath(@__DIR__, "wizard_gauge_only")

const LOAD_FORMATS = (
    (id="jld2", value="JLD", extension="jld2"),
    (id="ildg", value="ILDG", extension="ildg"),
    (id="bridge", value="BridgeText", extension="txt"),
)

const SAVE_FORMATS = (
    (id="none", value="nothing"),
    (id="jld2", value="JLD"),
    (id="ildg", value="ILDG"),
    (id="bridge", value="BridgeText"),
)

const UPDATE_VARIANTS = (
    (id="heatbath_or", method="Heatbath", use_or=true, sexton=false),
    (id="heatbath_no_or", method="Heatbath", use_or=false, sexton=false),
    (id="hmc_no_sw", method="HMC", use_or=false, sexton=false),
    (id="hmc_sw", method="HMC", use_or=false, sexton=true),
)

function initialization_variants(colors::Int)
    variants = Any[
        (id="cold", kind=:cold, format=nothing),
        (id="hot", kind=:hot, format=nothing),
    ]
    for format in LOAD_FORMATS
        push!(variants, (
            id="file_$(format.id)",
            kind=:file,
            format,
        ))
    end
    if colors == 2
        push!(variants, (id="instanton", kind=:instanton, format=nothing))
    end
    push!(variants, (
        id="embedded_instanton",
        kind=:embedded_instanton,
        format=nothing,
    ))
    return variants
end

function wizard_sections(colors::Int)
    physical = PS.Print_Physical_parameters()
    physical.L = [4, 4, 4, 4]
    physical.NC = colors
    physical.β = colors == 3 ? 5.7 : 2.7

    fermions = PS.Print_Fermions_parameters()
    fermions.quench = true
    fermions.Dirac_operator = "nothing"
    fermion_parameters = PS.Quench_parameters()

    control = PS.Print_System_control_parameters()
    control.randomseed = 111
    control.verboselevel = 2

    hmc = PS.Print_HMCrelated_parameters()
    measurement = PS.Measurement_parameterset()
    measurement.measurement_methods = PS.Measurement_parameters[
        PS.Plaq_parameters(),
    ]
    gradient = PS.Print_Gradientflow_parameters()
    gradient.hasgradientflow = false
    gradient_measurement = PS.Measurement_parameterset()

    return (;
        physical,
        fermions,
        fermion_parameters,
        control,
        hmc,
        measurement,
        gradient,
        gradient_measurement,
    )
end

function apply_initialization!(sections, initialization)
    physical = sections.physical
    control = sections.control
    if initialization.kind === :cold
        physical.initial = "cold"
    elseif initialization.kind === :hot
        physical.initial = "hot"
    elseif initialization.kind === :instanton
        physical.initial = "one instanton"
    elseif initialization.kind === :embedded_instanton
        physical.initial = "embedded instanton"
    elseif initialization.kind === :file
        format = initialization.format
        control.loadU_format = format.value
        physical.initial = "./confs/conf_00000001.$(format.extension)"
        physical.initialtrj = 1
    else
        error("unsupported initialization $(initialization.kind)")
    end
    return sections
end

function apply_update!(sections, update)
    physical = sections.physical
    hmc = sections.hmc
    physical.Nthermalization = 0
    physical.Nsteps = 101
    physical.update_method = update.method
    physical.useOR = update.use_or
    physical.numOR = update.use_or ? 3 : 0

    hmc.MDsteps = 20
    hmc.Δτ = 1 / hmc.MDsteps
    hmc.SextonWeingargten = update.sexton
    hmc.N_SextonWeingargten = 2
    return sections
end

function apply_output!(sections, save_format)
    header = WIZARD.make_headername(
        sections.physical,
        sections.fermions,
        sections.fermion_parameters,
    )
    control = sections.control
    control.measurement_basedir = "./measurements"
    control.measurement_dir = header
    control.log_dir = "./logs"
    control.logfile = "$header.txt"
    control.saveU_format = save_format.value
    if save_format.value != "nothing"
        control.saveU_every = 10
        control.saveU_dir = "./confs_$header"
    end
    return sections
end

function printable_dictionary(sections)
    parameters = Dict{String,Any}(
        "Physical setting" => PS.struct2dict(sections.physical),
        "Physical setting(fermions)" => merge(
            PS.struct2dict(sections.fermions),
            PS.struct2dict(sections.fermion_parameters),
        ),
        "System Control" => PS.struct2dict(sections.control),
        "HMC related" => PS.struct2dict(sections.hmc),
        "Measurement set" => PS.struct2dict(sections.measurement),
        "gradientflow_measurements" => merge(
            PS.struct2dict(sections.gradient),
            PS.struct2dict(sections.gradient_measurement),
        ),
    )

    if sections.physical.update_method == "Heatbath"
        delete!(parameters["HMC related"], "Δτ")
        delete!(parameters["HMC related"], "MDsteps")
    end
    PS.remove_default_values!(parameters)
    gradient = parameters["gradientflow_measurements"]
    gradient["measurements_for_flow"] = deepcopy(
        gradient["measurement_methods"],
    )
    delete!(gradient, "measurement_methods")
    return parameters
end

function validate_parameter_file(path::String)
    parameters = TOML.parsefile(path)
    fermions = parameters["Physical setting(fermions)"]
    fermions["Dirac_operator"] == "nothing" || error(
        "$path is not a gauge-only parameter file",
    )
    measurements = parameters["Measurement set"]["measurement_methods"]
    Set(keys(measurements)) == Set(["Plaquette"]) || error(
        "$path does not contain exactly one Plaquette measurement",
    )
    measurements["Plaquette"]["methodname"] == "Plaquette" || error(
        "$path contains an invalid Plaquette measurement",
    )
    gradient = parameters["gradientflow_measurements"]
    !gradient["hasgradientflow"] || error("$path enables gradient flow")
    isempty(gradient["measurements_for_flow"]) || error(
        "$path contains gradient-flow measurements",
    )
    return nothing
end

function write_parameter_file(filename::String, parameters)
    path = joinpath(OUTPUT_DIRECTORY, filename)
    open(path, "w") do io
        TOML.print(io, parameters)
    end
    validate_parameter_file(path)
    return path
end

function generate_update_cases!(filenames::Vector{String}, colors::Int)
    for initialization in initialization_variants(colors)
        for update in UPDATE_VARIANTS
            for save_format in SAVE_FORMATS
                sections = wizard_sections(colors)
                apply_initialization!(sections, initialization)
                apply_update!(sections, update)
                apply_output!(sections, save_format)
                filename = join((
                    "su$colors",
                    initialization.id,
                    update.id,
                    "save_$(save_format.id)",
                ), "_") * ".toml"
                write_parameter_file(
                    filename,
                    printable_dictionary(sections),
                )
                push!(filenames, filename)
            end
        end
    end
    return filenames
end

function generate_fileloading_cases!(filenames::Vector{String}, colors::Int)
    for load_format in LOAD_FORMATS
        for list_mode in ("all", "list")
            for save_format in SAVE_FORMATS
                sections = wizard_sections(colors)
                sections.physical.update_method = "Fileloading"
                sections.control.loadU_format = load_format.value
                sections.control.loadU_dir = "./confs"
                sections.control.loadU_fromfile = list_mode == "list"
                if list_mode == "list"
                    sections.control.loadU_filename = "filelist.txt"
                end
                apply_output!(sections, save_format)
                filename = join((
                    "su$colors",
                    "fileloading_$(load_format.id)",
                    list_mode,
                    "save_$(save_format.id)",
                ), "_") * ".toml"
                write_parameter_file(
                    filename,
                    printable_dictionary(sections),
                )
                push!(filenames, filename)
            end
        end
    end
    return filenames
end

function main()
    mkpath(OUTPUT_DIRECTORY)
    filenames = String[]
    for colors in (3, 2)
        generate_update_cases!(filenames, colors)
        generate_fileloading_cases!(filenames, colors)
    end

    length(filenames) == 256 || error(
        "expected 256 Wizard cases, generated $(length(filenames))",
    )
    length(unique(filenames)) == length(filenames) || error(
        "generated duplicate filenames",
    )
    existing = sort(filter(
        name -> endswith(name, ".toml"),
        readdir(OUTPUT_DIRECTORY),
    ))
    sort(filenames) == existing || error(
        "the output directory contains stale or missing TOML files",
    )

    println("generated $(length(filenames)) gauge-only Wizard parameter files")
    println("output: $OUTPUT_DIRECTORY")
    return nothing
end

main()
