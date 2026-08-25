import ..Parameter_structs:
    Plaq_parameters,
    Poly_parameters,
    ChiralCondensate_parameters,
    TopologicalCharge_parameters,
    Pion_parameters,
    Wilson_loop_parameters,
    Energy_density_parameters

@enum WizardV2Page begin
    WizardV2ModePage
    WizardV2LatticePage
    WizardV2SourcePage
    WizardV2FermionPage
    WizardV2UpdatePage
    WizardV2MeasurementPage
    WizardV2GradientFlowPage
    WizardV2OutputPage
    WizardV2ReviewPage
end

@enum WizardV2Action begin
    WizardV2Next
    WizardV2BackAction
    WizardV2QuitAction
    WizardV2Save
    WizardV2Jump
end

struct WizardV2PageResult
    action::WizardV2Action
    target::Union{Nothing,WizardV2Page}
end

WizardV2PageResult(action::WizardV2Action) = WizardV2PageResult(action, nothing)

abstract type AbstractWizardV2PromptResult end

struct WizardV2Answer{T} <: AbstractWizardV2PromptResult
    value::T
end

struct WizardV2Back <: AbstractWizardV2PromptResult end
struct WizardV2Quit <: AbstractWizardV2PromptResult end

abstract type AbstractWizardV2UI end
struct TerminalWizardV2UI <: AbstractWizardV2UI end

"""A small non-interactive UI used to exercise Wizard v2 in tests."""
mutable struct ScriptedWizardV2UI <: AbstractWizardV2UI
    answers::Vector{Any}
    prompts::Vector{String}
end

ScriptedWizardV2UI(answers) = ScriptedWizardV2UI(Any[answers...], String[])

function wizard_v2_scripted_answer!(ui::ScriptedWizardV2UI, prompt)
    push!(ui.prompts, prompt)
    isempty(ui.answers) && error("no scripted Wizard v2 answer remains for: $prompt")
    answer = popfirst!(ui.answers)
    answer === :back && return WizardV2Back()
    answer === :quit && return WizardV2Quit()
    return WizardV2Answer(answer)
end

function wizard_v2_choice(
    ::TerminalWizardV2UI,
    prompt,
    options::Vector{String};
    default::Int=1,
    allow_back::Bool=true,
)
    displayed = copy(options)
    back_index = 0
    if allow_back
        push!(displayed, "← Back")
        back_index = length(displayed)
    end
    push!(displayed, "Quit wizard")
    quit_index = length(displayed)

    println(prompt)
    selected = request(
        RadioMenu(displayed; charset=:unicode, scroll_wrap=true);
        cursor=clamp(default, 1, length(options)),
    )
    selected == -1 && return WizardV2Quit()
    allow_back && selected == back_index && return WizardV2Back()
    selected == quit_index && return WizardV2Quit()
    return WizardV2Answer(selected)
end

function wizard_v2_choice(
    ui::ScriptedWizardV2UI,
    prompt,
    options::Vector{String};
    default::Int=1,
    allow_back::Bool=true,
)
    result = wizard_v2_scripted_answer!(ui, prompt)
    result isa WizardV2Answer || return result
    selected = result.value
    selected isa Integer || error("scripted choice for $prompt must be an integer")
    1 <= selected <= length(options) || error(
        "scripted choice $selected is outside 1:$(length(options)) for $prompt",
    )
    return WizardV2Answer(Int(selected))
end

function wizard_v2_multiselect(
    ::TerminalWizardV2UI,
    prompt,
    options::Vector{String};
    selected=Set{Int}(),
)
    displayed = copy(options)
    push!(displayed, "← Back")
    back_index = length(displayed)
    push!(displayed, "Quit wizard")
    quit_index = length(displayed)

    println(prompt)
    choices = request(MultiSelectMenu(
        displayed;
        charset=:unicode,
        selected=Set{Int}(selected),
    ))
    back_index in choices && return WizardV2Back()
    quit_index in choices && return WizardV2Quit()
    return WizardV2Answer(Set(filter(index -> index <= length(options), choices)))
end

function wizard_v2_multiselect(
    ui::ScriptedWizardV2UI,
    prompt,
    options::Vector{String};
    selected=Set{Int}(),
)
    result = wizard_v2_scripted_answer!(ui, prompt)
    result isa WizardV2Answer || return result
    choices = Set(Int(index) for index in result.value)
    all(index -> 1 <= index <= length(options), choices) || error(
        "a scripted selection is outside 1:$(length(options)) for $prompt",
    )
    return WizardV2Answer(choices)
end

function wizard_v2_value(
    ::TerminalWizardV2UI,
    ::Type{T},
    prompt;
    default,
    valid=(_ -> true),
    validation_message="Invalid value.",
) where {T}
    while true
        raw = String(Base.prompt(
            "$prompt  [:back / :quit]",
            default=string(default),
        ))
        command = lowercase(strip(raw))
        command == ":back" && return WizardV2Back()
        command == ":quit" && return WizardV2Quit()

        value = if T === String
            raw
        else
            tryparse(T, strip(raw))
        end
        if value === nothing || !valid(value)
            println(validation_message)
            continue
        end
        return WizardV2Answer(value)
    end
end

function wizard_v2_value(
    ui::ScriptedWizardV2UI,
    ::Type{T},
    prompt;
    default,
    valid=(_ -> true),
    validation_message="Invalid value.",
) where {T}
    result = wizard_v2_scripted_answer!(ui, prompt)
    result isa WizardV2Answer || return result
    raw = result.value
    value = if T === String
        String(raw)
    elseif raw isa T
        raw
    elseif raw isa AbstractString
        tryparse(T, raw)
    else
        try
            convert(T, raw)
        catch
            nothing
        end
    end
    value === nothing && error("$validation_message Prompt: $prompt")
    valid(value) || error("$validation_message Prompt: $prompt")
    return WizardV2Answer(value)
end

function wizard_v2_page_result(result::AbstractWizardV2PromptResult)
    result isa WizardV2Back && return WizardV2PageResult(WizardV2BackAction)
    result isa WizardV2Quit && return WizardV2PageResult(WizardV2QuitAction)
    error("an answer cannot be converted to a page navigation result")
end

macro wizard_v2_answer(variable, expression)
    return quote
        local wizard_v2_prompt_result = $(esc(expression))
        if !(wizard_v2_prompt_result isa WizardV2Answer)
            return wizard_v2_page_result(wizard_v2_prompt_result)
        end
        $(esc(variable)) = wizard_v2_prompt_result.value
    end
end

@enum WizardV2FermionKind begin
    WizardV2Quenched
    WizardV2Wilson
    WizardV2WilsonClover
    WizardV2Staggered
    WizardV2HISQ
    WizardV2Domainwall
    WizardV2MobiusDomainwall
end

"""Concrete storage for every editable Wizard fermion form."""
mutable struct WizardV2FermionChoices{Q,W,C,S,H,D,M}
    selected::WizardV2FermionKind
    quench::Q
    wilson::W
    wilson_clover::C
    staggered::S
    hisq::H
    domainwall::D
    mobius_domainwall::M
end

WizardV2FermionChoices() = WizardV2FermionChoices(
    WizardV2Quenched,
    Quench_parameters(),
    Wilson_parameters(),
    Wilson_parameters(
        Dirac_operator="WilsonClover",
        hasclover=true,
        Clover_coefficient=1.5612,
    ),
    Staggered_parameters(),
    HISQ_parameters(),
    Domainwall_parameters(),
    MobiusDomainwall_parameters(),
)

function active_fermion_parameters(choices::WizardV2FermionChoices)
    choices.selected == WizardV2Quenched && return choices.quench
    choices.selected == WizardV2Wilson && return choices.wilson
    choices.selected == WizardV2WilsonClover && return choices.wilson_clover
    choices.selected == WizardV2Staggered && return choices.staggered
    choices.selected == WizardV2HISQ && return choices.hisq
    choices.selected == WizardV2Domainwall && return choices.domainwall
    return choices.mobius_domainwall
end

function select_fermion_parameters!(choices::WizardV2FermionChoices, value)
    if value isa Quench_parameters
        choices.quench = value
        choices.selected = WizardV2Quenched
    elseif value isa Wilson_parameters
        if value.Dirac_operator == "WilsonClover" || value.hasclover
            value.Dirac_operator = "WilsonClover"
            value.hasclover = true
            choices.wilson_clover = value
            choices.selected = WizardV2WilsonClover
        else
            value.Dirac_operator = "Wilson"
            value.hasclover = false
            choices.wilson = value
            choices.selected = WizardV2Wilson
        end
    elseif value isa Staggered_parameters
        choices.staggered = value
        choices.selected = WizardV2Staggered
    elseif value isa HISQ_parameters
        choices.hisq = value
        choices.selected = WizardV2HISQ
    elseif value isa Domainwall_parameters
        choices.domainwall = value
        choices.selected = WizardV2Domainwall
    elseif value isa MobiusDomainwall_parameters
        choices.mobius_domainwall = value
        choices.selected = WizardV2MobiusDomainwall
    else
        throw(ArgumentError(
            "unsupported Wizard v2 fermion parameter type $(typeof(value))",
        ))
    end
    return value
end

mutable struct WizardV2Draft{F<:WizardV2FermionChoices}
    mode::Wizardmode
    filename::String
    isfileloading::Bool
    physicalparams::Print_Physical_parameters
    fermionparams::Print_Fermions_parameters
    fermion_choices::F
    controlparams::Print_System_control_parameters
    hmcparams::Print_HMCrelated_parameters
    slhmc_beta::Float64
    cg::ConjugateGradient
    measurement::Measurement_parameterset
    gradient_params::Print_Gradientflow_parameters
    measurement_gradientflow::Measurement_parameterset
    last_header::String
    completed_pages::Set{WizardV2Page}
end

function Base.getproperty(draft::WizardV2Draft, name::Symbol)
    name === :fermion_parameters && return active_fermion_parameters(
        getfield(draft, :fermion_choices),
    )
    return getfield(draft, name)
end

function Base.setproperty!(draft::WizardV2Draft, name::Symbol, value)
    name === :fermion_parameters && return select_fermion_parameters!(
        getfield(draft, :fermion_choices),
        value,
    )
    return setfield!(draft, name, value)
end

function WizardV2Draft()
    return WizardV2Draft(
        simple,
        "my_parameters.toml",
        false,
        Print_Physical_parameters(),
        Print_Fermions_parameters(),
        WizardV2FermionChoices(),
        Print_System_control_parameters(),
        Print_HMCrelated_parameters(),
        5.7,
        ConjugateGradient(),
        Measurement_parameterset(),
        Print_Gradientflow_parameters(),
        Measurement_parameterset(),
        "",
        Set{WizardV2Page}(),
    )
end

function copy_wizard_v2_draft!(destination::WizardV2Draft, source::WizardV2Draft)
    for name in fieldnames(WizardV2Draft)
        setfield!(destination, name, deepcopy(getfield(source, name)))
    end
    return destination
end

function reset_wizard_v2_draft!(draft::WizardV2Draft, mode::Wizardmode, filename)
    replacement = WizardV2Draft()
    replacement.mode = mode
    replacement.filename = filename
    copy_wizard_v2_draft!(draft, replacement)
    return draft
end

wizard_v2_choice_index(value, values; default=1) =
    something(findfirst(==(value), values), default)

function edit_wizard_v2_mode!(ui, draft)
    working = deepcopy(draft)
    @wizard_v2_answer mode_index wizard_v2_choice(
        ui,
        "Choose Wizard mode",
        ["simple", "expert"];
        default=Int(working.mode),
        allow_back=false,
    )
    selected_mode = Wizardmode(mode_index)

    @wizard_v2_answer filename wizard_v2_value(
        ui,
        String,
        "Name of the parameter file";
        default=working.filename,
        valid=value -> !isempty(strip(value)),
        validation_message="The parameter filename must not be empty.",
    )
    if selected_mode != working.mode
        reset_wizard_v2_draft!(working, selected_mode, filename)
    else
        working.mode = selected_mode
        working.filename = filename
    end
    copy_wizard_v2_draft!(draft, working)
    return WizardV2PageResult(WizardV2Next)
end

function edit_wizard_v2_lattice!(ui, draft)
    working = deepcopy(draft)
    physical = working.physicalparams
    control = working.controlparams
    isexpert = working.mode == expert

    if isexpert
        sizes = copy(physical.L)
        labels = ("Nx", "Ny", "Nz", "Nt")
        for index in eachindex(sizes)
            @wizard_v2_answer extent wizard_v2_value(
                ui,
                Int64,
                labels[index];
                default=sizes[index],
                valid=value -> value > 0,
                validation_message="Lattice extents must be positive integers.",
            )
            sizes[index] = extent
        end
        physical.L = sizes

        previous_colors = physical.NC
        @wizard_v2_answer color_index wizard_v2_choice(
            ui,
            "Choose a gauge group",
            ["SU(3)", "SU(2)"];
            default=previous_colors == 2 ? 2 : 1,
        )
        physical.NC = color_index == 1 ? 3 : 2
        if physical.NC != previous_colors
            physical.β = physical.NC == 3 ? 5.7 : 2.7
        end

        @wizard_v2_answer randomseed wizard_v2_value(
            ui,
            Int64,
            "Random seed";
            default=control.randomseed,
        )
        control.randomseed = randomseed

        @wizard_v2_answer verboselevel wizard_v2_value(
            ui,
            Int64,
            "Verbose level";
            default=WizardV2LatticePage in working.completed_pages ?
                    control.verboselevel : 2,
            valid=value -> 1 <= value <= 3,
            validation_message="Verbose level must be between 1 and 3.",
        )
        control.verboselevel = verboselevel
    else
        if physical.NC != 3
            physical.NC = 3
            physical.β = 5.7
        end
        @wizard_v2_answer spatial_size wizard_v2_value(
            ui,
            Int64,
            "Spatial lattice size";
            default=physical.L[1],
            valid=value -> value > 0,
            validation_message="The spatial lattice size must be positive.",
        )
        @wizard_v2_answer temporal_size wizard_v2_value(
            ui,
            Int64,
            "Temporal lattice size";
            default=physical.L[4],
            valid=value -> value > 0,
            validation_message="The temporal lattice size must be positive.",
        )
        physical.L = [spatial_size, spatial_size, spatial_size, temporal_size]
    end

    @wizard_v2_answer beta wizard_v2_value(
        ui,
        Float64,
        "Beta";
        default=physical.β,
        valid=value -> value > 0,
        validation_message="Beta must be positive.",
    )
    physical.β = beta

    copy_wizard_v2_draft!(draft, working)
    return WizardV2PageResult(WizardV2Next)
end

const WIZARD_V2_LOAD_FORMATS = ["JLD", "ILDG", "BridgeText"]

function wizard_v2_load_format_index(format)
    return wizard_v2_choice_index(format, WIZARD_V2_LOAD_FORMATS; default=1)
end

function edit_wizard_v2_source!(ui, draft)
    working = deepcopy(draft)
    physical = working.physicalparams
    control = working.controlparams
    was_fileloading = working.isfileloading

    @wizard_v2_answer loading_index wizard_v2_choice(
        ui,
        "Only measure configurations from a directory?",
        ["No", "Yes"];
        default=working.isfileloading ? 2 : 1,
    )
    working.isfileloading = loading_index == 2
    if working.isfileloading != was_fileloading
        delete!(working.completed_pages, WizardV2FermionPage)
        delete!(working.completed_pages, WizardV2UpdatePage)
    end

    if working.isfileloading
        physical.update_method = "Fileloading"
        physical.initial = "cold"
        physical.initialtrj = 1
        physical.Nthermalization = 0
        physical.Nsteps = 100
        physical.useOR = false
        physical.numOR = 0
        working.fermionparams = Print_Fermions_parameters()
        working.fermionparams.Dirac_operator = "nothing"
        working.fermionparams.quench = true
        working.fermion_parameters = Quench_parameters()
        working.cg = ConjugateGradient()
        working.hmcparams = Print_HMCrelated_parameters()

        @wizard_v2_answer format_index wizard_v2_choice(
            ui,
            "Configuration format for loading",
            WIZARD_V2_LOAD_FORMATS;
            default=wizard_v2_load_format_index(control.loadU_format),
        )
        control.loadU_format = WIZARD_V2_LOAD_FORMATS[format_index]

        @wizard_v2_answer load_directory wizard_v2_value(
            ui,
            String,
            "Loading directory";
            default=isempty(control.loadU_dir) ? "./confs" : control.loadU_dir,
            valid=value -> !isempty(strip(value)),
            validation_message="The loading directory must not be empty.",
        )
        control.loadU_dir = load_directory

        @wizard_v2_answer list_index wizard_v2_choice(
            ui,
            "Which configurations do you use?",
            [
                "All configurations in the directory",
                "Configurations written in a list",
            ];
            default=control.loadU_fromfile ? 2 : 1,
        )
        control.loadU_fromfile = list_index == 2
        if control.loadU_fromfile
            @wizard_v2_answer list_filename wizard_v2_value(
                ui,
                String,
                "Name of the configuration list";
                default=isempty(control.loadU_filename) ?
                        "filelist.txt" : control.loadU_filename,
                valid=value -> !isempty(strip(value)),
                validation_message="The list filename must not be empty.",
            )
            control.loadU_filename = list_filename
        else
            control.loadU_filename = ""
        end
    else
        if was_fileloading
            physical.update_method = "HMC"
            physical.Nthermalization = 0
            physical.Nsteps = 100
            physical.useOR = false
            physical.numOR = 0
            working.fermionparams = Print_Fermions_parameters()
            working.fermion_parameters = Quench_parameters()
            working.cg = ConjugateGradient()
            working.hmcparams = Print_HMCrelated_parameters()
        end
        control.loadU_dir = ""
        control.loadU_fromfile = false
        control.loadU_filename = ""

        initialization_options = physical.NC == 2 ?
            [
                "cold start",
                "hot start",
                "start from a file",
                "one instanton",
                "SU(2) embedded instanton",
            ] :
            [
                "cold start",
                "hot start",
                "start from a file",
                "SU(2) embedded instanton",
            ]
        current_initial = if physical.initial == "cold"
            1
        elseif physical.initial == "hot"
            2
        elseif physical.initial == "one instanton"
            physical.NC == 2 ? 4 : 3
        elseif physical.initial == "embedded instanton"
            length(initialization_options)
        else
            3
        end
        @wizard_v2_answer initial_index wizard_v2_choice(
            ui,
            "Choose initial configurations",
            initialization_options;
            default=current_initial,
        )

        if initial_index == 1
            physical.initial = "cold"
            physical.initialtrj = 1
            control.loadU_format = nothing
        elseif initial_index == 2
            physical.initial = "hot"
            physical.initialtrj = 1
            control.loadU_format = nothing
        elseif initial_index == 3
            @wizard_v2_answer format_index wizard_v2_choice(
                ui,
                "Configuration format for loading",
                WIZARD_V2_LOAD_FORMATS;
                default=wizard_v2_load_format_index(control.loadU_format),
            )
            control.loadU_format = WIZARD_V2_LOAD_FORMATS[format_index]
            extension = get_filename_extension(Fileformat(format_index))
            default_path = startswith(physical.initial, "cold") ||
                           startswith(physical.initial, "hot") ||
                           physical.initial in
                           ("one instanton", "embedded instanton") ?
                           "./confs/conf_00000001.$extension" : physical.initial
            @wizard_v2_answer initial_path wizard_v2_value(
                ui,
                String,
                "Initial configuration filename";
                default=default_path,
                valid=value -> !isempty(strip(value)),
                validation_message="The initial configuration filename must not be empty.",
            )
            physical.initial = initial_path
            @wizard_v2_answer initialtrj wizard_v2_value(
                ui,
                Int64,
                "Start trajectory number";
                default=physical.initialtrj,
                valid=value -> value >= 0,
                validation_message="The trajectory number must be nonnegative.",
            )
            physical.initialtrj = initialtrj
        elseif physical.NC == 2 && initial_index == 4
            physical.initial = "one instanton"
            physical.initialtrj = 1
            control.loadU_format = nothing
        else
            physical.initial = "embedded instanton"
            physical.initialtrj = 1
            control.loadU_format = nothing
        end
    end

    copy_wizard_v2_draft!(draft, working)
    return WizardV2PageResult(WizardV2Next)
end

function configure_wizard_v2_cg!(ui, cg::ConjugateGradient)
    @wizard_v2_answer eps wizard_v2_value(
        ui,
        Float64,
        "Relative error in CG loops";
        default=cg.eps,
        valid=value -> value > 0,
        validation_message="The CG relative error must be positive.",
    )
    @wizard_v2_answer max_steps wizard_v2_value(
        ui,
        Int64,
        "Maximum iteration steps in CG loops";
        default=cg.MaxCGstep,
        valid=value -> value > 0,
        validation_message="The maximum CG step count must be positive.",
    )
    cg.eps = eps
    cg.MaxCGstep = max_steps
    return nothing
end

function configure_wizard_v2_wilson!(ui, parameters, cg)
    @wizard_v2_answer hop wizard_v2_value(
        ui,
        Float64,
        "Hopping parameter kappa";
        default=parameters.hop,
        valid=value -> value > 0,
        validation_message="Kappa must be positive.",
    )
    parameters.hop = hop
    return configure_wizard_v2_cg!(ui, cg)
end

function configure_wizard_v2_wilson_clover!(ui, parameters, cg)
    @wizard_v2_answer hop wizard_v2_value(
        ui,
        Float64,
        "Hopping parameter kappa";
        default=parameters.hop,
        valid=value -> value > 0,
        validation_message="Kappa must be positive.",
    )
    @wizard_v2_answer clover_coefficient wizard_v2_value(
        ui,
        Float64,
        "Clover coefficient cSW";
        default=parameters.Clover_coefficient,
        valid=isfinite,
        validation_message="cSW must be finite.",
    )
    parameters.Dirac_operator = "WilsonClover"
    parameters.hasclover = true
    parameters.hop = hop
    parameters.Clover_coefficient = clover_coefficient
    return configure_wizard_v2_cg!(ui, cg)
end

const WIZARD_V2_NF_VALUES = [2, 3, 4, 8, 1]

function configure_wizard_v2_staggered!(ui, parameters, cg)
    @wizard_v2_answer mass wizard_v2_value(
        ui,
        Float64,
        "Staggered fermion mass";
        default=parameters.mass,
        valid=value -> value > 0,
        validation_message="The staggered fermion mass must be positive.",
    )
    parameters.mass = mass
    @wizard_v2_answer nf_index wizard_v2_choice(
        ui,
        "Number of flavors (tastes)",
        [
            "2 (RHMC)",
            "3 (RHMC)",
            "4 (HMC)",
            "8 (HMC)",
            "1 (RHMC)",
        ];
        default=wizard_v2_choice_index(parameters.Nf, WIZARD_V2_NF_VALUES; default=1),
    )
    parameters.Nf = WIZARD_V2_NF_VALUES[nf_index]
    return configure_wizard_v2_cg!(ui, cg)
end

function configure_wizard_v2_hisq!(ui, parameters, cg)
    @wizard_v2_answer mass wizard_v2_value(
        ui,
        Float64,
        "HISQ fermion mass";
        default=parameters.mass,
        valid=value -> value > 0,
        validation_message="The HISQ mass must be positive.",
    )
    @wizard_v2_answer naik_epsilon wizard_v2_value(
        ui,
        Float64,
        "Naik correction epsilon_N";
        default=parameters.naik_epsilon,
        valid=isfinite,
        validation_message="The Naik correction must be finite.",
    )
    parameters.mass = mass
    parameters.naik_epsilon = naik_epsilon
    @wizard_v2_answer nf_index wizard_v2_choice(
        ui,
        "Number of HISQ flavors (tastes)",
        [
            "2 (RHMC)",
            "3 (RHMC)",
            "4 (HMC)",
            "8 (HMC)",
            "1 (RHMC)",
        ];
        default=wizard_v2_choice_index(
            parameters.Nf,
            WIZARD_V2_NF_VALUES;
            default=3,
        ),
    )
    parameters.Nf = WIZARD_V2_NF_VALUES[nf_index]
    return configure_wizard_v2_cg!(ui, cg)
end

function configure_wizard_v2_domainwall!(ui, parameters, cg; use_legacy_defaults=false)
    @wizard_v2_answer size5 wizard_v2_value(
        ui,
        Int64,
        "Size of the extra dimension L5";
        default=parameters.N5,
        valid=value -> value > 0,
        validation_message="L5 must be positive.",
    )
    @wizard_v2_answer domainwall_m5 wizard_v2_value(
        ui,
        Float64,
        "Domain-wall M";
        default=parameters.M,
        valid=value -> value < 0,
        validation_message="Domain-wall M must be negative.",
    )
    @wizard_v2_answer mass wizard_v2_value(
        ui,
        Float64,
        "Domain-wall physical mass";
        default=use_legacy_defaults ? 0.25 : parameters.m,
    )
    parameters.N5 = size5
    parameters.M = domainwall_m5
    parameters.m = mass
    return configure_wizard_v2_cg!(ui, cg)
end

function configure_wizard_v2_mobius_domainwall!(
    ui,
    parameters,
    cg;
    use_legacy_defaults=false,
)
    @wizard_v2_answer size5 wizard_v2_value(
        ui,
        Int64,
        "Size of the extra dimension L5";
        default=parameters.N5,
        valid=value -> value > 0,
        validation_message="L5 must be positive.",
    )
    @wizard_v2_answer domainwall_m5 wizard_v2_value(
        ui,
        Float64,
        "Möbius domain-wall M";
        default=parameters.M,
        valid=value -> value < 0,
        validation_message="Möbius domain-wall M must be negative.",
    )
    @wizard_v2_answer mass wizard_v2_value(
        ui,
        Float64,
        "Möbius domain-wall physical mass";
        default=use_legacy_defaults ? 0.25 : parameters.m,
        valid=isfinite,
        validation_message="The physical mass must be finite.",
    )
    @wizard_v2_answer b wizard_v2_value(
        ui,
        Float64,
        "Möbius coefficient b";
        default=parameters.b,
        valid=isfinite,
        validation_message="The Möbius coefficient b must be finite.",
    )
    @wizard_v2_answer c wizard_v2_value(
        ui,
        Float64,
        "Möbius coefficient c";
        default=parameters.c,
        valid=isfinite,
        validation_message="The Möbius coefficient c must be finite.",
    )
    parameters.N5 = size5
    parameters.M = domainwall_m5
    parameters.m = mass
    parameters.b = b
    parameters.c = c
    return configure_wizard_v2_cg!(ui, cg)
end

function configure_wizard_v2_stout!(ui, stout::Stout_parameters)
    previous = Dict{String,Float64}()
    if stout.stout_loops !== nothing
        for (loop, coefficient) in zip(stout.stout_loops, stout.ρ)
            previous[loop] = coefficient
        end
    end
    selected = Set(
        index for (index, loop) in enumerate(kindsof_loops) if haskey(previous, loop)
    )
    @wizard_v2_answer choices wizard_v2_multiselect(
        ui,
        "Loops used in stout smearing",
        collect(kindsof_loops);
        selected=selected,
    )
    loops = String[]
    coefficients = Float64[]
    for index in sort!(collect(choices))
        loop = kindsof_loops[index]
        @wizard_v2_answer coefficient wizard_v2_value(
            ui,
            Float64,
            "Coefficient rho for $loop";
            default=get(previous, loop, 0.1),
        )
        push!(loops, loop)
        push!(coefficients, coefficient)
    end
    stout.stout_loops = loops
    stout.ρ = coefficients
    return nothing
end

function configure_wizard_v2_primary_smearing!(ui, working)
    fermions = working.fermionparams
    current_stout = fermions.smearing_for_fermion == "stout"
    @wizard_v2_answer smearing_index wizard_v2_choice(
        ui,
        "Smearing for the fermion action",
        ["No smearing", "stout smearing"];
        default=current_stout ? 2 : 1,
    )
    if smearing_index == 1
        fermions.smearing_for_fermion = "nothing"
        fermions.stout_numlayers = nothing
        fermions.stout_ρ = nothing
        fermions.stout_loops = nothing
    else
        stout = Stout_parameters(
            numlayers=something(fermions.stout_numlayers, 1),
            ρ=isnothing(fermions.stout_ρ) ?
              Float64[] : copy(fermions.stout_ρ),
            stout_loops=isnothing(fermions.stout_loops) ?
                        nothing : copy(fermions.stout_loops),
        )
        result = configure_wizard_v2_stout!(ui, stout)
        result === nothing || return result
        fermions.smearing_for_fermion = "stout"
        fermions.stout_numlayers = stout.numlayers
        fermions.stout_ρ = stout.ρ
        fermions.stout_loops = stout.stout_loops
    end
    return nothing
end

function edit_wizard_v2_fermion!(ui, draft)
    working = deepcopy(draft)
    if working.mode == simple
        if !(working.fermion_parameters isa Wilson_parameters) ||
           working.fermionparams.quench ||
           working.fermionparams.Dirac_operator != "Wilson"
            working.fermion_parameters = Wilson_parameters()
            working.cg = ConjugateGradient()
            working.fermionparams = Print_Fermions_parameters()
        end
        working.fermionparams.quench = false
        working.fermionparams.Dirac_operator = "Wilson"
        @wizard_v2_answer hop wizard_v2_value(
            ui,
            Float64,
            "Hopping parameter kappa";
            default=working.fermion_parameters.hop,
            valid=value -> value > 0,
            validation_message="Kappa must be positive.",
        )
        working.fermion_parameters.hop = hop
    else
        fermion_kinds = working.physicalparams.NC == 3 ?
            (
                WizardV2Quenched,
                WizardV2Wilson,
                WizardV2WilsonClover,
                WizardV2Staggered,
                WizardV2HISQ,
                WizardV2Domainwall,
                WizardV2MobiusDomainwall,
            ) :
            (
                WizardV2Quenched,
                WizardV2Wilson,
                WizardV2WilsonClover,
                WizardV2Staggered,
                WizardV2Domainwall,
                WizardV2MobiusDomainwall,
            )
        labels = working.physicalparams.NC == 3 ?
            [
                "Nothing (quenched approximation)",
                "Wilson fermion (2-flavor)",
                "Wilson--clover fermion (2-flavor)",
                "Staggered fermion",
                "HISQ fermion (SU(3))",
                "Domain-wall fermion (experimental)",
                "Möbius domain-wall fermion (experimental)",
            ] :
            [
                "Nothing (quenched approximation)",
                "Wilson fermion (2-flavor)",
                "Wilson--clover fermion (2-flavor)",
                "Staggered fermion",
                "Domain-wall fermion (experimental)",
                "Möbius domain-wall fermion (experimental)",
            ]
        current_kind = if working.fermionparams.quench
            WizardV2Quenched
        elseif working.fermionparams.Dirac_operator == "Wilson"
            WizardV2Wilson
        elseif working.fermionparams.Dirac_operator == "WilsonClover"
            WizardV2WilsonClover
        elseif working.fermionparams.Dirac_operator == "Staggered"
            WizardV2Staggered
        elseif working.fermionparams.Dirac_operator == "HISQ"
            WizardV2HISQ
        elseif working.fermionparams.Dirac_operator == "Domainwall"
            WizardV2Domainwall
        else
            WizardV2MobiusDomainwall
        end
        current_type = something(findfirst(==(current_kind), fermion_kinds), 1)
        @wizard_v2_answer fermion_index wizard_v2_choice(
            ui,
            "Choose a dynamical fermion",
            labels;
            default=current_type,
        )
        selected_kind = fermion_kinds[fermion_index]

        if selected_kind != current_kind
            working.fermionparams = Print_Fermions_parameters()
            working.cg = ConjugateGradient()
        end
        if selected_kind == WizardV2Quenched
            working.fermionparams.quench = true
            working.fermionparams.Dirac_operator = "nothing"
            working.fermion_parameters = Quench_parameters()
        elseif selected_kind == WizardV2Wilson
            working.fermionparams.quench = false
            working.fermionparams.Dirac_operator = "Wilson"
            if !(working.fermion_parameters isa Wilson_parameters) ||
               working.fermion_parameters.Dirac_operator != "Wilson" ||
               working.fermion_parameters.hasclover
                working.fermion_parameters = Wilson_parameters()
            end
            result = configure_wizard_v2_wilson!(
                ui,
                working.fermion_parameters,
                working.cg,
            )
            result === nothing || return result
        elseif selected_kind == WizardV2WilsonClover
            working.fermionparams.quench = false
            working.fermionparams.Dirac_operator = "WilsonClover"
            if !(working.fermion_parameters isa Wilson_parameters) ||
               working.fermion_parameters.Dirac_operator != "WilsonClover"
                working.fermion_parameters = Wilson_parameters(
                    Dirac_operator="WilsonClover",
                    hasclover=true,
                    Clover_coefficient=1.5612,
                )
            end
            result = configure_wizard_v2_wilson_clover!(
                ui,
                working.fermion_parameters,
                working.cg,
            )
            result === nothing || return result
        elseif selected_kind == WizardV2Staggered
            working.fermionparams.quench = false
            working.fermionparams.Dirac_operator = "Staggered"
            working.fermion_parameters isa Staggered_parameters ||
                (working.fermion_parameters = Staggered_parameters())
            result = configure_wizard_v2_staggered!(
                ui,
                working.fermion_parameters,
                working.cg,
            )
            result === nothing || return result
        elseif selected_kind == WizardV2HISQ
            working.fermionparams.quench = false
            working.fermionparams.Dirac_operator = "HISQ"
            working.fermion_parameters isa HISQ_parameters ||
                (working.fermion_parameters = HISQ_parameters())
            result = configure_wizard_v2_hisq!(
                ui,
                working.fermion_parameters,
                working.cg,
            )
            result === nothing || return result
        elseif selected_kind == WizardV2Domainwall
            working.fermionparams.quench = false
            working.fermionparams.Dirac_operator = "Domainwall"
            working.fermion_parameters isa Domainwall_parameters ||
                (working.fermion_parameters = Domainwall_parameters())
            result = configure_wizard_v2_domainwall!(
                ui,
                working.fermion_parameters,
                working.cg,
                use_legacy_defaults=(
                    selected_kind != current_kind ||
                    !(WizardV2FermionPage in working.completed_pages)
                ),
            )
            result === nothing || return result
        else
            working.fermionparams.quench = false
            working.fermionparams.Dirac_operator = "MobiusDomainwall"
            working.fermion_parameters isa MobiusDomainwall_parameters ||
                (working.fermion_parameters = MobiusDomainwall_parameters())
            result = configure_wizard_v2_mobius_domainwall!(
                ui,
                working.fermion_parameters,
                working.cg,
                use_legacy_defaults=(
                    selected_kind != current_kind ||
                    !(WizardV2FermionPage in working.completed_pages)
                ),
            )
            result === nothing || return result
        end
    end

    if working.fermionparams.Dirac_operator == "HISQ"
        # HISQ already contains its two-level Fat7/U(3)/Lepage/Naik link
        # construction. The Wizard does not add a second, outer stout layer.
        working.fermionparams.smearing_for_fermion = "nothing"
        working.fermionparams.stout_numlayers = nothing
        working.fermionparams.stout_ρ = nothing
        working.fermionparams.stout_loops = nothing
    elseif !working.fermionparams.quench
        result = configure_wizard_v2_primary_smearing!(ui, working)
        result === nothing || return result
    else
        working.fermionparams.smearing_for_fermion = "nothing"
    end

    working.hmcparams.eps = working.cg.eps
    working.hmcparams.MaxCGstep = working.cg.MaxCGstep
    copy_wizard_v2_draft!(draft, working)
    return WizardV2PageResult(WizardV2Next)
end

function edit_wizard_v2_update!(ui, draft)
    working = deepcopy(draft)
    physical = working.physicalparams
    hmc = working.hmcparams
    isexpert = working.mode == expert
    isrevisit = WizardV2UpdatePage in working.completed_pages

    if isexpert
        if working.fermionparams.quench
            @wizard_v2_answer update_index wizard_v2_choice(
                ui,
                "Choose an update method",
                ["Heatbath", "Hybrid Monte Carlo"];
                default=isrevisit ?
                        (physical.update_method == "Heatbath" ? 1 : 2) : 1,
            )
            physical.update_method = update_index == 1 ? "Heatbath" : "HMC"
        else
            @wizard_v2_answer update_index wizard_v2_choice(
                ui,
                "Choose an update method",
                ["Hybrid Monte Carlo", "Self-learning Hybrid Monte Carlo"];
                default=physical.update_method == "SLHMC" ? 2 : 1,
            )
            physical.update_method = update_index == 1 ? "HMC" : "SLHMC"
            if physical.update_method == "SLHMC"
                @wizard_v2_answer effective_beta wizard_v2_value(
                    ui,
                    Float64,
                    "Effective beta used by the SLHMC MD action";
                    default=isrevisit ? working.slhmc_beta : physical.β,
                    valid=isfinite,
                    validation_message="The SLHMC effective beta must be finite.",
                )
                working.slhmc_beta = effective_beta
            end
        end
    else
        physical.update_method = "HMC"
    end

    if physical.update_method == "Heatbath"
        hmc = Print_HMCrelated_parameters()
        working.hmcparams = hmc
        @wizard_v2_answer overrelax_index wizard_v2_choice(
            ui,
            "Use the overrelaxation method?",
            ["true", "false"];
            default=isrevisit ? (physical.useOR ? 1 : 2) : 1,
        )
        physical.useOR = overrelax_index == 1
        if physical.useOR
            @wizard_v2_answer num_or wizard_v2_value(
                ui,
                Int64,
                "Number of overrelaxation updates";
                default=physical.numOR == 0 ? 3 : physical.numOR,
                valid=value -> value > 0,
                validation_message="The number of overrelaxation updates must be positive.",
            )
            physical.numOR = num_or
        else
            physical.numOR = 0
        end
    else
        physical.useOR = false
        physical.numOR = 0
    end

    if isexpert
        @wizard_v2_answer thermalization wizard_v2_value(
            ui,
            Int64,
            "Number of thermalization steps";
            default=physical.Nthermalization,
            valid=value -> value >= 0,
            validation_message="Thermalization steps must be nonnegative.",
        )
        physical.Nthermalization = thermalization
    else
        physical.Nthermalization = 0
    end

    @wizard_v2_answer nsteps wizard_v2_value(
        ui,
        Int64,
        "Number of total trajectories after thermalization";
        default=isrevisit ? physical.Nsteps : 100 + physical.initialtrj,
        valid=value -> value > 0,
        validation_message="The number of trajectories must be positive.",
    )
    physical.Nsteps = nsteps

    if isexpert && physical.update_method != "Heatbath"
        @wizard_v2_answer mdsteps wizard_v2_value(
            ui,
            Int64,
            "MD steps";
            default=hmc.MDsteps,
            valid=value -> value > 0,
            validation_message="MD steps must be positive.",
        )
        default_delta = hmc.MDsteps == mdsteps ? hmc.Δτ : 1 / mdsteps
        @wizard_v2_answer delta_tau wizard_v2_value(
            ui,
            Float64,
            "Delta tau";
            default=default_delta,
            valid=value -> value > 0,
            validation_message="Delta tau must be positive.",
        )
        hmc.MDsteps = mdsteps
        hmc.Δτ = delta_tau

        @wizard_v2_answer sexton_index wizard_v2_choice(
            ui,
            "Use the Sexton-Weingarten multi-time-scale method?",
            ["false", "true"];
            default=hmc.SextonWeingargten ? 2 : 1,
        )
        hmc.SextonWeingargten = sexton_index == 2
        if hmc.SextonWeingargten
            @wizard_v2_answer sexton_steps wizard_v2_value(
                ui,
                Int64,
                "Number of Sexton-Weingarten steps";
                default=hmc.N_SextonWeingargten,
                valid=value -> value > 0,
                validation_message="Sexton-Weingarten steps must be positive.",
            )
            hmc.N_SextonWeingargten = sexton_steps
        else
            hmc.N_SextonWeingargten = 2
        end
    elseif !isexpert
        defaults = MD()
        hmc.MDsteps = defaults.MDsteps
        hmc.Δτ = defaults.Δτ
        hmc.SextonWeingargten = defaults.SextonWeingargten
        hmc.N_SextonWeingargten = defaults.N_SextonWeingargten
    end
    hmc.eps = working.cg.eps
    hmc.MaxCGstep = working.cg.MaxCGstep

    copy_wizard_v2_draft!(draft, working)
    return WizardV2PageResult(WizardV2Next)
end

const WIZARD_V2_MEASUREMENT_NAMES = [
    "Plaquette",
    "Polyakov_loop",
    "Topological_charge",
    "Chiral_condensate",
    "Pion_correlator",
    "Wilson_loop",
    "Energy_density",
]

function wizard_v2_measurement_option(method)
    return something(
        findfirst(==(method.methodname), WIZARD_V2_MEASUREMENT_NAMES),
        0,
    )
end

function wizard_v2_new_measurement(option::Int)
    option == Int(Plaquette) && return Plaq_parameters()
    option == Int(Polyakov_loop) && return Poly_parameters()
    option == Int(Topological_charge) && return TopologicalCharge_parameters()
    option == Int(Chiral_condensate) && return ChiralCondensate_parameters()
    option == Int(Pion_correlator) && return Pion_parameters()
    option == Int(Wilson_loop) && return Wilson_loop_parameters()
    option == Int(Energy_density) && return Energy_density_parameters()
    error("unsupported measurement option $option")
end

function wizard_v2_existing_measurements(measurement_set)
    return Dict(
        wizard_v2_measurement_option(method) => method for
        method in measurement_set.measurement_methods if
        wizard_v2_measurement_option(method) != 0
    )
end

function configure_wizard_v2_measurement_smearing!(ui, method)
    current_stout = method.smearing_for_fermion == "stout"
    @wizard_v2_answer smearing_index wizard_v2_choice(
        ui,
        "Smearing for this measurement",
        ["No smearing", "stout smearing"];
        default=current_stout ? 2 : 1,
    )
    if smearing_index == 1
        method.smearing_for_fermion = "nothing"
        defaults = typeof(method)()
        for field in (:stout_numlayers, :stout_ρ, :stout_loops)
            hasfield(typeof(method), field) || continue
            setfield!(method, field, deepcopy(getfield(defaults, field)))
        end
    else
        stout = Stout_parameters()
        if current_stout
            stout.numlayers = method.stout_numlayers
            stout.ρ = copy(method.stout_ρ)
            stout.stout_loops = copy(method.stout_loops)
        end
        result = configure_wizard_v2_stout!(ui, stout)
        result === nothing || return result
        method.smearing_for_fermion = "stout"
        method.stout_numlayers = stout.numlayers
        method.stout_ρ = stout.ρ
        method.stout_loops = stout.stout_loops
    end
    return nothing
end

function configure_wizard_v2_measurement!(ui, method, option, lattice; isnew=false)
    default_interval = isnew ? 1 : method.measure_every
    @wizard_v2_answer interval wizard_v2_value(
        ui,
        Int64,
        "Measure $(WIZARD_V2_MEASUREMENT_NAMES[option]) every";
        default=default_interval,
        valid=value -> value > 0,
        validation_message="Measurement intervals must be positive.",
    )
    method.measure_every = interval

    if option == Int(Wilson_loop)
        r_default = isnew ? min(lattice[1], lattice[2], lattice[3]) ÷ 2 : method.Rmax
        t_default = isnew ? lattice[4] ÷ 2 : method.Tmax
        @wizard_v2_answer rmax wizard_v2_value(
            ui,
            Int64,
            "Maximum R for the RxT Wilson loop";
            default=r_default,
            valid=value -> value >= 0,
            validation_message="Rmax must be nonnegative.",
        )
        @wizard_v2_answer tmax wizard_v2_value(
            ui,
            Int64,
            "Maximum T for the RxT Wilson loop";
            default=t_default,
            valid=value -> value >= 0,
            validation_message="Tmax must be nonnegative.",
        )
        method.Rmax = rmax
        method.Tmax = tmax
    elseif option == Int(Chiral_condensate)
        @wizard_v2_answer mass wizard_v2_value(
            ui,
            Float64,
            "Mass for the chiral-condensate measurement";
            default=method.mass,
        )
        @wizard_v2_answer eps wizard_v2_value(
            ui,
            Float64,
            "Relative error in measurement CG loops";
            default=method.eps,
            valid=value -> value > 0,
            validation_message="The CG relative error must be positive.",
        )
        @wizard_v2_answer max_steps wizard_v2_value(
            ui,
            Int64,
            "Maximum iteration steps in measurement CG loops";
            default=method.MaxCGstep,
            valid=value -> value > 0,
            validation_message="The maximum CG step count must be positive.",
        )
        method.mass = mass
        method.Nf = 4
        method.eps = eps
        method.MaxCGstep = max_steps
        result = configure_wizard_v2_measurement_smearing!(ui, method)
        result === nothing || return result
    elseif option == Int(Pion_correlator)
        @wizard_v2_answer fermion_index wizard_v2_choice(
            ui,
            "Fermion type for the pion-correlator measurement",
            ["Standard Wilson fermion", "Staggered fermion"];
            default=method.fermiontype == "Staggered" ? 2 : 1,
        )
        cg = ConjugateGradient(eps=method.eps, MaxCGstep=method.MaxCGstep)
        if fermion_index == 1
            method.fermiontype = "Wilson"
            parameters = Wilson_parameters()
            result = configure_wizard_v2_wilson!(ui, parameters, cg)
            result === nothing || return result
        else
            method.fermiontype = "Staggered"
            parameters = Staggered_parameters()
            result = configure_wizard_v2_staggered!(ui, parameters, cg)
            result === nothing || return result
        end
        # Keep compatibility with the original Wizard, which records the
        # selected fermion type and solver values but leaves the nested
        # `fermion_parameters` field at its constructor default.
        method.eps = cg.eps
        method.MaxCGstep = cg.MaxCGstep
        result = configure_wizard_v2_measurement_smearing!(ui, method)
        result === nothing || return result
    end
    return nothing
end

function configure_wizard_v2_measurement_set!(ui, measurement_set, lattice, choices)
    existing = wizard_v2_existing_measurements(measurement_set)
    configured = Measurement_parameters[]
    for option in sort!(collect(choices))
        isnew = !haskey(existing, option)
        method = isnew ? wizard_v2_new_measurement(option) : deepcopy(existing[option])
        result = configure_wizard_v2_measurement!(
            ui,
            method,
            option,
            lattice;
            isnew=isnew,
        )
        result === nothing || return result
        push!(configured, method)
    end
    measurement_set.measurement_methods = configured
    return nothing
end

function edit_wizard_v2_measurements!(ui, draft)
    working = deepcopy(draft)
    existing = wizard_v2_existing_measurements(working.measurement)
    choices = if working.mode == simple
        Set([Int(Plaquette), Int(Polyakov_loop), Int(Pion_correlator)])
    else
        @wizard_v2_answer selected wizard_v2_multiselect(
            ui,
            "Select measurement methods",
            WIZARD_V2_MEASUREMENT_NAMES;
            selected=Set(keys(existing)),
        )
        selected
    end
    result = configure_wizard_v2_measurement_set!(
        ui,
        working.measurement,
        working.physicalparams.L,
        choices,
    )
    result === nothing || return result
    copy_wizard_v2_draft!(draft, working)
    return WizardV2PageResult(WizardV2Next)
end

function edit_wizard_v2_gradient_flow!(ui, draft)
    working = deepcopy(draft)
    gradient = working.gradient_params
    if working.mode == simple
        gradient.hasgradientflow = false
        working.measurement_gradientflow = Measurement_parameterset()
    else
        @wizard_v2_answer flow_index wizard_v2_choice(
            ui,
            "Perform measurements with gradient flow?",
            ["No", "Yes"];
            default=gradient.hasgradientflow ? 2 : 1,
        )
        gradient.hasgradientflow = flow_index == 2
        if gradient.hasgradientflow
            @wizard_v2_answer eps_flow wizard_v2_value(
                ui,
                Float64,
                "Gradient-flow time step";
                default=gradient.eps_flow,
                valid=value -> value > 0,
                validation_message="The gradient-flow time step must be positive.",
            )
            @wizard_v2_answer numflow wizard_v2_value(
                ui,
                Int64,
                "Number of gradient-flow updates";
                default=gradient.numflow,
                valid=value -> value >= 0,
                validation_message="The number of flow updates must be nonnegative.",
            )
            gradient.eps_flow = eps_flow
            gradient.numflow = numflow

            existing = wizard_v2_existing_measurements(
                working.measurement_gradientflow,
            )
            @wizard_v2_answer choices wizard_v2_multiselect(
                ui,
                "Select measurements performed during gradient flow",
                WIZARD_V2_MEASUREMENT_NAMES;
                selected=Set(keys(existing)),
            )
            result = configure_wizard_v2_measurement_set!(
                ui,
                working.measurement_gradientflow,
                working.physicalparams.L,
                choices,
            )
            result === nothing || return result
        else
            working.measurement_gradientflow = Measurement_parameterset()
        end
    end
    copy_wizard_v2_draft!(draft, working)
    return WizardV2PageResult(WizardV2Next)
end

function refresh_wizard_v2_automatic_names!(draft)
    header = make_headername(
        draft.physicalparams,
        draft.fermionparams,
        draft.fermion_parameters,
    )
    previous = draft.last_header
    control = draft.controlparams

    if isempty(control.measurement_dir) || control.measurement_dir == previous
        control.measurement_dir = header
    end
    if isempty(control.logfile) || control.logfile == "$previous.txt"
        control.logfile = "$header.txt"
    end
    if isempty(control.saveU_dir) || control.saveU_dir == "./confs_$previous"
        control.saveU_dir = "./confs_$header"
    end
    draft.last_header = header
    return header
end

const WIZARD_V2_SAVE_FORMATS = ["nothing", "JLD", "ILDG", "BridgeText"]

function edit_wizard_v2_output!(ui, draft)
    working = deepcopy(draft)
    control = working.controlparams
    header = refresh_wizard_v2_automatic_names!(working)

    if !isempty(working.measurement.measurement_methods)
        @wizard_v2_answer measurement_basedir wizard_v2_value(
            ui,
            String,
            "Base directory for measurements";
            default=isempty(control.measurement_basedir) ?
                    "./measurements" : control.measurement_basedir,
            valid=value -> !isempty(strip(value)),
            validation_message="The measurement base directory must not be empty.",
        )
        @wizard_v2_answer measurement_dir wizard_v2_value(
            ui,
            String,
            "Measurement directory inside $measurement_basedir";
            default=control.measurement_dir,
            valid=value -> !isempty(strip(value)),
            validation_message="The measurement directory must not be empty.",
        )
        control.measurement_basedir = measurement_basedir
        control.measurement_dir = measurement_dir
    else
        control.measurement_basedir = "./measurements"
        control.measurement_dir = header
    end

    @wizard_v2_answer log_dir wizard_v2_value(
        ui,
        String,
        "Log directory";
        default=isempty(control.log_dir) ? "./logs" : control.log_dir,
        valid=value -> !isempty(strip(value)),
        validation_message="The log directory must not be empty.",
    )
    @wizard_v2_answer logfile wizard_v2_value(
        ui,
        String,
        "Log filename";
        default=control.logfile,
        valid=value -> !isempty(strip(value)),
        validation_message="The log filename must not be empty.",
    )
    control.log_dir = log_dir
    control.logfile = logfile

    if working.mode == expert && !working.isfileloading
        @wizard_v2_answer save_index wizard_v2_choice(
            ui,
            "Configuration format for saving",
            ["no save", "JLD", "ILDG", "Text format (BridgeText)"];
            default=wizard_v2_choice_index(
                something(control.saveU_format, "nothing"),
                WIZARD_V2_SAVE_FORMATS;
                default=1,
            ),
        )
        control.saveU_format = WIZARD_V2_SAVE_FORMATS[save_index]
        if control.saveU_format != "nothing"
            @wizard_v2_answer save_every wizard_v2_value(
                ui,
                Int64,
                "Save a configuration every";
                default=WizardV2OutputPage in working.completed_pages ?
                        control.saveU_every : 10,
                valid=value -> value > 0,
                validation_message="The configuration save interval must be positive.",
            )
            @wizard_v2_answer save_dir wizard_v2_value(
                ui,
                String,
                "Configuration saving directory";
                default=control.saveU_dir,
                valid=value -> !isempty(strip(value)),
                validation_message="The configuration saving directory must not be empty.",
            )
            control.saveU_every = save_every
            control.saveU_dir = save_dir
        else
            control.saveU_every = 1
            control.saveU_dir = ""
        end
    else
        control.saveU_format = nothing
        control.saveU_every = 1
        control.saveU_dir = ""
    end

    copy_wizard_v2_draft!(draft, working)
    return WizardV2PageResult(WizardV2Next)
end

function wizard_v2_route(draft)
    route = WizardV2Page[
        WizardV2ModePage,
        WizardV2LatticePage,
        WizardV2SourcePage,
    ]
    if !draft.isfileloading
        push!(route, WizardV2FermionPage)
        push!(route, WizardV2UpdatePage)
    end
    append!(route, [
        WizardV2MeasurementPage,
        WizardV2GradientFlowPage,
        WizardV2OutputPage,
        WizardV2ReviewPage,
    ])
    return route
end

const WIZARD_V2_PAGE_LABELS = Dict(
    WizardV2ModePage => "Mode and filename",
    WizardV2LatticePage => "Lattice and gauge settings",
    WizardV2SourcePage => "Configuration source",
    WizardV2FermionPage => "Fermion settings",
    WizardV2UpdatePage => "Update settings",
    WizardV2MeasurementPage => "Measurements",
    WizardV2GradientFlowPage => "Gradient flow",
    WizardV2OutputPage => "Output settings",
    WizardV2ReviewPage => "Review",
)

function print_wizard_v2_page_header(::TerminalWizardV2UI, page, draft)
    route = wizard_v2_route(draft)
    index = findfirst(==(page), route)
    position = index === nothing ? "?" : string(index)
    println(
        "\n--- Wizard v2 [$position/$(length(route))]: ",
        WIZARD_V2_PAGE_LABELS[page],
        " ---",
    )
end

print_wizard_v2_page_header(::ScriptedWizardV2UI, page, draft) = nothing

function wizard_v2_parameter_dictionary(draft::WizardV2Draft)
    return wizard_parameter_dictionary(
        draft.physicalparams,
        draft.fermionparams,
        draft.fermion_parameters,
        draft.controlparams,
        draft.hmcparams,
        draft.measurement,
        draft.gradient_params,
        draft.measurement_gradientflow,
        slhmc_beta=draft.slhmc_beta,
    )
end

function print_wizard_v2_review(::TerminalWizardV2UI, draft)
    println("\nReview generated TOML")
    println("---------------------")
    io = IOBuffer()
    TOML.print(io, wizard_v2_parameter_dictionary(draft))
    print(String(take!(io)))
    println("---------------------")
end

print_wizard_v2_review(::ScriptedWizardV2UI, draft) = nothing

function edit_wizard_v2_review!(ui, draft)
    print_wizard_v2_review(ui, draft)
    editable_pages = wizard_v2_route(draft)[1:end-1]
    options = ["Save TOML"]
    append!(options, ["Edit: $(WIZARD_V2_PAGE_LABELS[page])" for page in editable_pages])
    @wizard_v2_answer review_index wizard_v2_choice(
        ui,
        "Confirm parameters or edit a section",
        options;
        default=1,
    )
    review_index == 1 && return WizardV2PageResult(WizardV2Save)
    return WizardV2PageResult(WizardV2Jump, editable_pages[review_index-1])
end

function edit_wizard_v2_page!(ui, page, draft)
    page == WizardV2ModePage && return edit_wizard_v2_mode!(ui, draft)
    page == WizardV2LatticePage && return edit_wizard_v2_lattice!(ui, draft)
    page == WizardV2SourcePage && return edit_wizard_v2_source!(ui, draft)
    page == WizardV2FermionPage && return edit_wizard_v2_fermion!(ui, draft)
    page == WizardV2UpdatePage && return edit_wizard_v2_update!(ui, draft)
    page == WizardV2MeasurementPage && return edit_wizard_v2_measurements!(ui, draft)
    page == WizardV2GradientFlowPage && return edit_wizard_v2_gradient_flow!(ui, draft)
    page == WizardV2OutputPage && return edit_wizard_v2_output!(ui, draft)
    page == WizardV2ReviewPage && return edit_wizard_v2_review!(ui, draft)
    error("unsupported Wizard v2 page $page")
end

function wizard_v2_next_page(draft, current)
    route = wizard_v2_route(draft)
    index = findfirst(==(current), route)
    index === nothing && error("$current is not active in the current Wizard v2 route")
    index == length(route) && return current
    return route[index+1]
end

function wizard_v2_write_and_build_spec(draft)
    parameters = write_wizard_parameter_file(
        draft.filename,
        draft.physicalparams,
        draft.fermionparams,
        draft.fermion_parameters,
        draft.controlparams,
        draft.hmcparams,
        draft.measurement,
        draft.gradient_params,
        draft.measurement_gradientflow,
        slhmc_beta=draft.slhmc_beta,
    )
    return simulation_spec_from_legacy_toml(parameters)
end

function print_wizard_v2_completion(::TerminalWizardV2UI, draft)
    println("""
    --------------------------------------------------------------------------------
    run_wizard is done.

    The returned value is a typed SimulationSpec. To run it, use

        session = build_simulation(spec, GaugefieldsEnvironment())
        run!(session)

    The legacy file runner also remains available:

        run_LQCD("$(draft.filename)")

    The output parameter file is $(draft.filename).
    --------------------------------------------------------------------------------
    """)
end

print_wizard_v2_completion(::ScriptedWizardV2UI, draft) = nothing

function run_wizard(ui::AbstractWizardV2UI)
    ui isa TerminalWizardV2UI && print_wizard_logo(stdout)
    draft = WizardV2Draft()
    history = WizardV2Page[]
    page = WizardV2ModePage

    while true
        print_wizard_v2_page_header(ui, page, draft)
        result = edit_wizard_v2_page!(ui, page, draft)
        if result.action == WizardV2Next
            push!(draft.completed_pages, page)
            push!(history, page)
            page = wizard_v2_next_page(draft, page)
        elseif result.action == WizardV2BackAction
            isempty(history) || (page = pop!(history))
        elseif result.action == WizardV2Jump
            route = wizard_v2_route(draft)
            target_index = findfirst(==(result.target), route)
            target_index === nothing && error(
                "$(result.target) is not editable in the current Wizard route",
            )
            history = copy(route[1:target_index-1])
            page = result.target
        elseif result.action == WizardV2Save
            spec = wizard_v2_write_and_build_spec(draft)
            print_wizard_v2_completion(ui, draft)
            return spec
        elseif result.action == WizardV2QuitAction
            return nothing
        end
    end
end

"""
    run_wizard()

Create a LatticeQCD TOML parameter file using the navigable terminal Wizard.
Every prompt supports returning to the previous page, and the review page can
jump back to any active section. The returned value is a typed
`SimulationSpec`; no `Params`, logfile, or output directory is constructed.

Use `run_wizard_legacy()` for the original Wizard and its historical `Params`
return value.
"""
run_wizard() = run_wizard(TerminalWizardV2UI())

"""Backward-compatible alias for the typed [`run_wizard`](@ref)."""
const run_wizardv2 = run_wizard
