module LegacyInputRegressions
using LatticeQCD, Gaugefields, QCDMeasurements, Test, TOML
const root = mktempdir()
function document(name)
    dir=joinpath(root,name)
    return Dict{String,Any}(
        "Physical setting"=>Dict{String,Any}("L"=>[2,2,2,2],"NC"=>2,"Nwing"=>1,"β"=>1.9,
            "initial"=>"cold","update_method"=>"Heatbath","Nsteps"=>3,"initialtrj"=>1,"Nthermalization"=>3),
        "Physical setting(fermions)"=>Dict("Dirac_operator"=>"nothing"),
        "HMC related"=>Dict{String,Any}(),
        "System Control"=>Dict{String,Any}("logfile"=>"run.log","log_dir"=>dir,
            "measurement_basedir"=>dir,"measurement_dir"=>"measurements","verboselevel"=>1,
            "hasgradientflow"=>false,"randomseed"=>1234),
        "Measurement set"=>Dict{String,Any}("measurement_methods"=>Dict("Plaquette"=>Dict("methodname"=>"Plaquette","measure_every"=>1))),
        "gradientflow_measurements"=>Dict{String,Any}("numflow"=>2,"Nflow"=>1,"eps_flow"=>0.01,
            "measurements_for_flow"=>Dict("Polyakov_loop"=>Dict("methodname"=>"Polyakov_loop","measure_every"=>1))))
end
function save_doc(d,name)
    file=joinpath(root,name*".toml"); open(io->TOML.print(io,d),file,"w"); file
end
parse_old(d)=LatticeQCD.Parameters_TOML.construct_Params_from_TOML(d)

@testset "Legacy input rejects ambiguous keys before I/O" begin
    for reader in (parse_old, simulation_spec_from_legacy_toml)
        for kind in ("duplicate", "unknown", "collision", "section", "table", "flow_conflict")
            d = document(kind)
            if kind == "duplicate"
                d["Measurement set"]["measurement_basedir"] = "other"
            elseif kind == "unknown"
                d["Physical setting"]["Nstepz"] = 99
            elseif kind == "collision"
                d["System Control"]["Nsteps"] = 7
            elseif kind == "section"
                d["Physical seting"] = Dict("NC" => 2)
            elseif kind == "flow_conflict"
                d["gradientflow_measurements"]["hasgradientflow"] = true
            else
                d["Physical setting"] = 3
            end
            @test_throws ArgumentError reader(d)
            @test !ispath(joinpath(root, kind))
        end
        d = document("equal_flow")
        d["gradientflow_measurements"]["hasgradientflow"] = false
        # Validation alone is free of I/O, even when reporting legacy duplicates.
        @test_logs (:warn, r"duplicate legacy parameter.*hasgradientflow") LatticeQCD.Parameter_structs.validate_legacy_parameters(d)
        @test !ispath(joinpath(root, "equal_flow"))
    end
end

data_lines(path) = filter(l -> !isempty(l) && !startswith(l, "#"), readlines(path))
trajectories(path) = [parse(Int, first(split(l))) for l in data_lines(path)]

@testset "Legacy runner respects measurement controls" begin
    for (name, thermal, flow, expected) in (
        ("thermal", 3, false, [3]),
        ("flow", 3, true, [3]),
        ("unthermalized", 0, false, [0, 1, 2, 3]),
        ("allthermal", 4, true, Int[]),
    )
        d = document(name)
        d["Physical setting"]["Nthermalization"] = thermal
        d["System Control"]["hasgradientflow"] = flow
        file = save_doc(d, name)
        trace = joinpath(root, name * "-stdout.log")
        plaq = open(trace, "w") do io
            redirect_stdout(io) do
                run_LQCD_file(file)
            end
        end
        @test isfinite(plaq)
        @test trajectories(joinpath(root, name, "measurements", "Plaquette.txt")) == expected
        flowfile = joinpath(root, name, "measurements", "Polyakov_loop.txt")
        if flow
            # Data must be readable immediately on return (the flow stream closes).
            lines = data_lines(flowfile)
            @test length(lines) == (thermal == 3 ? 2 : 0)
            @test all(l -> first(split(l)) == "3", lines)
        else
            @test !isfile(flowfile)
            @test !occursin("# poly", read(trace, String))
        end
    end
    # Nflow integration steps must agree with the same number of one-step flows.
    final_flow = Float64[]
    for (name, nflow, samples) in (("two_steps", 2, 2), ("one_step", 1, 4))
        d = document(name)
        d["System Control"]["hasgradientflow"] = true
        d["gradientflow_measurements"]["Nflow"] = nflow
        d["gradientflow_measurements"]["numflow"] = samples
        run_LQCD_file(save_doc(d, name))
        line = last(data_lines(joinpath(root, name, "measurements", "Polyakov_loop.txt")))
        columns = split(line)
        @test parse(Float64, columns[3]) ≈ 0.04
        push!(final_flow, parse(Float64, columns[4]))
    end
    @test final_flow[1] ≈ final_flow[2] atol=1e-13
    file = save_doc(document("modern"), "modern")
    summary = run!(build_simulation(load_simulation_spec(file)); verbose=false)
    @test trajectories(joinpath(root, "modern", "measurements", "plaquette.txt")) == [3]
end

@testset "Compatibility measurement collection" begin
    U = Initialize_Gaugefields(2, 1, 2, 2, 2, 2; condition="cold", verbose_level=0)
    m = LatticeQCD.Measurements_set(U, root, Dict[])
    @test isempty(m.baremesurement_indices)
    @test isempty(m.flowmeasurement_indices)
    methods = Dict[
        Dict("methodname" => "Wilson_loop", "Tmax" => 1, "Rmax" => 1),
        Dict("methodname" => "Polyakov_loop", "bare measure" => false, "flow measure" => true),
    ]
    m = LatticeQCD.Measurements_set(U, root, methods)
    @test m.baremesurement_indices == [1]
    @test m.flowmeasurement_indices == [2]
    @test m.measurements[1] isa Wilson_loop_measurement
    @test m.measurements[1].filename == joinpath(root, "Wilson_loop.txt")
    @test m.measurements[2].filename == joinpath(root, "Polyakov_loop.txt")
    # Exercise the repaired Wilson loop, including its initialized index list.
    @test isnothing(LatticeQCD.measure(m, 1, U))
    @test m.measurements[1].outputvalues ≈ ones(1, 1)
end
end
