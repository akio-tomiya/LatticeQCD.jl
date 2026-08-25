using LatticeQCD
using TOML

const PS = LatticeQCD.Parameter_structs
const WIZARD = LatticeQCD.Wizard
const OUTPUT_DIRECTORY = joinpath(@__DIR__, "wizard_fermion_io")

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

function wizard_fermion_io_sections(colors, load_format, save_format)
    physical = PS.Print_Physical_parameters(
        L=[2, 2, 2, 2],
        NC=colors,
        β=colors == 3 ? 5.7 : 2.7,
        Nthermalization=0,
        Nsteps=1,
        initial="./confs/wizard_fermion_source.$(load_format.extension)",
        initialtrj=1,
        update_method="HMC",
        Nwing=1,
    )
    fermions = PS.Print_Fermions_parameters(
        quench=false,
        Dirac_operator="Wilson",
        smearing_for_fermion="nothing",
    )
    control = PS.Print_System_control_parameters(
        randomseed=111,
        verboselevel=0,
        log_dir="./logs",
        logfile="wizard-fermion-io.log",
        measurement_basedir="./measurements",
        measurement_dir="wizard-fermion-io",
        loadU_format=load_format.value,
        saveU_format=save_format.value,
        saveU_dir=save_format.value == "nothing" ? "" : "./confs_output",
        saveU_every=1,
    )
    hmc = PS.Print_HMCrelated_parameters(
        Δτ=0.001,
        MDsteps=1,
        SextonWeingargten=false,
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
        fermion_parameters=PS.Wilson_parameters(),
        control,
        hmc,
        measurement,
        gradient,
        gradient_measurement,
    )
end

function wizard_fermion_io_dictionary(sections)
    return WIZARD.wizard_parameter_dictionary(
        sections.physical,
        sections.fermions,
        sections.fermion_parameters,
        sections.control,
        sections.hmc,
        sections.measurement,
        sections.gradient,
        sections.gradient_measurement,
    )
end

function validate_wizard_fermion_io(path, colors, load_format, save_format)
    parameters = TOML.parsefile(path)
    physical = parameters["Physical setting"]
    fermions = parameters["Physical setting(fermions)"]
    control = parameters["System Control"]
    @assert get(physical, "NC", 3) == colors
    @assert physical["update_method"] == "HMC"
    @assert fermions["Dirac_operator"] == "Wilson"
    @assert fermions["quench"] == false
    @assert control["loadU_format"] == load_format.value
    @assert get(control, "saveU_format", "nothing") == save_format.value
    @assert Set(keys(
        parameters["Measurement set"]["measurement_methods"],
    )) == Set(["Plaquette"])
    @assert !parameters["gradientflow_measurements"]["hasgradientflow"]
    return nothing
end

function main()
    mkpath(OUTPUT_DIRECTORY)
    filenames = String[]
    for colors in (2, 3)
        for load_format in LOAD_FORMATS
            for save_format in SAVE_FORMATS
                filename = join((
                    "su$colors",
                    "load_$(load_format.id)",
                    "save_$(save_format.id)",
                ), "_") * ".toml"
                sections = wizard_fermion_io_sections(
                    colors,
                    load_format,
                    save_format,
                )
                path = joinpath(OUTPUT_DIRECTORY, filename)
                open(path, "w") do io
                    TOML.print(io, wizard_fermion_io_dictionary(sections))
                end
                validate_wizard_fermion_io(
                    path,
                    colors,
                    load_format,
                    save_format,
                )
                push!(filenames, filename)
            end
        end
    end

    length(filenames) == 24 || error(
        "expected 24 Wizard fermion I/O cases, generated $(length(filenames))",
    )
    existing = sort(filter(
        filename -> endswith(filename, ".toml"),
        readdir(OUTPUT_DIRECTORY),
    ))
    sort(filenames) == existing || error(
        "the output directory contains stale or missing TOML files",
    )
    println("generated $(length(filenames)) fermion I/O Wizard parameter files")
    println("output: $OUTPUT_DIRECTORY")
    return nothing
end

main()
