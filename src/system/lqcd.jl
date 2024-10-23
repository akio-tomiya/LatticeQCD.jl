module LQCD
import ..Universe_module: Univ
import ..Transform_oldinputfile: transform_to_toml
import ..Parameters_TOML: construct_Params_from_TOML
import ..AbstractUpdate_module: Updatemethod, update!
import Gaugefields:
    Gradientflow,
    println_verbose_level1,
    get_myrank,
    flow!,
    save_binarydata,
    save_textdata,
    saveU
using QCDMeasurements
#import ..AbstractMeasurement_module:Measurement_methods,
#calc_measurement_values,measure,Plaquette_measurement,get_temporary_gaugefields
import QCDMeasurements: measure, Plaquette_measurement#, get_temporary_gaugefields
import ..LatticeQCD: Measurement_methods, calc_measurement_values, LatticeQCDversion


using Gaugefields
using LatticeDiracOperators
using Wilsonloop
using InteractiveUtils
using Dates
using Random
import ..Simpleprint: println_rank0


function run_LQCD(filenamein::String; MPIparallel=false)
    plaq = run_LQCD_file(filenamein, MPIparallel=MPIparallel)
    return plaq
end

function printTOMLfile(univ, filename)
    data = readlines(filename)
    for d in data
        println_verbose_level1(univ, d)
    end
end

function run_LQCD_file(filenamein::String; MPIparallel=false)
    if MPIparallel == true
        println_rank0("test")
    end

    filename_head = splitext(filenamein)[1]
    ext = splitext(filenamein)[end]
    filename = filename_head * ".toml"
    if ext == ".jl"
        transform_to_toml(filenamein)
        println_rank0("input file $filenamein is transformed to $filename")
    elseif ext == ".toml"
    else
        @error "$filenamein is not supported. use a TOML format."
    end

    parameters = construct_Params_from_TOML(filename)
    Random.seed!(parameters.randomseed)


    univ = Univ(parameters)

    println_verbose_level1(univ, "# ", pwd())

    println_verbose_level1(univ, "# ", Dates.now())
    io = IOBuffer()

    InteractiveUtils.versioninfo(io)
    versioninfo = String(take!(io))
    println_verbose_level1(univ, versioninfo)

    versionstring = """
    LatticeQCD $(LatticeQCDversion)
    LatticeDiracOperators $(pkgversion(LatticeDiracOperators))
    Gaugefields $(pkgversion(Gaugefields))
    Wilsonloop $(pkgversion(Wilsonloop))
    """
    println_verbose_level1(univ, versionstring)


    println_verbose_level1(univ, "# Input TOML file------------------")
    printTOMLfile(univ, filename)
    println_verbose_level1(univ, "#----------------------------------")


    updatemethod = Updatemethod(parameters, univ)



    eps_flow = parameters.eps_flow
    numflow = parameters.numflow
    Nflow = parameters.Nflow
    dτ = Nflow * eps_flow
    gradientflow = Gradientflow(univ.U, Nflow=1, eps=eps_flow)

    measurements =
        Measurement_methods(univ.U, parameters.measuredir, parameters.measurement_methods)
    plaq_m = Plaquette_measurement(univ.U, printvalues=false)
    i_plaq = 0
    for i = 1:measurements.num_measurements
        if isa(measurements.measurements[i], Plaquette_measurement)
            i_plaq = i
            plaq_m = measurements.measurements[i]
        end
    end
    #if i_plaq == 0
    #    plaq_m = Plaquette_measurement(univ.U, printvalues=false)
    #end

    measurements_for_flow =
        Measurement_methods(univ.U, parameters.measuredir, parameters.measurements_for_flow)


    calc_measurement_values(measurements, 0, univ.U)
    savedata = Savedata(
        parameters.saveU_format,
        parameters.saveU_dir,
        parameters.saveU_every,
        parameters.update_method,
        univ.U,
    )

    value, runtime_all = @timed begin

        numaccepts = 0
        for itrj = parameters.initialtrj:parameters.Nsteps
            println_verbose_level1(univ, "# itrj = $itrj")
            accepted, runtime = @timed update!(updatemethod, univ.U)

            println_verbose_level1(univ, "Update: Elapsed time $runtime [s]")
            if accepted
                numaccepts += 1
            end
            save_gaugefield(savedata, univ.U, itrj)

            measurestrings = calc_measurement_values(measurements, itrj, univ.U)
            #println(measurestrings)
            #println("$(univ.verbose_print.fp)")
            for st in measurestrings
                println(univ.verbose_print.fp, st)
            end



            Usmr = deepcopy(univ.U)
            for istep = 1:numflow
                τ = istep * dτ
                flow!(Usmr, gradientflow)
                additional_string = "$itrj $istep $τ "
                for i = 1:measurements_for_flow.num_measurements
                    interval = measurements_for_flow.intervals[i]
                    if istep % interval == 0
                        measure(
                            measurements_for_flow.measurements[i],
                            Usmr,
                            additional_string=additional_string,
                        )
                    end
                end
            end

            println_verbose_level1(
                univ,
                "Acceptance $numaccepts/$itrj : $(round(numaccepts*100/itrj)) %",
            )
            flush(univ.verbose_print.fp)

        end
        #println(stdout,univ.verbose_print.fp)
        #println("close file")
        close(univ.verbose_print.fp)
        for meas in measurements.measurements
            close(meas.verbose_print.fp)
        end
        #println(stdout,univ.verbose_print.fp)
    end

    println_verbose_level1(univ, "Total Elapsed time $(runtime_all) [s]")


    temps = QCDMeasurements.get_temporary_gaugefields(plaq_m)
    plaq = real(calculate_Plaquette(univ.U, temps[1], temps[2]) * plaq_m.factor)
    return plaq

end


mutable struct Savedata
    issaved::Bool
    saveU_format::Union{Nothing,String}
    saveU_dir::String
    saveU_every::Int64
    itrjsavecount::Int64


    function Savedata(saveU_format, saveU_dir, saveU_every, update_method, U)
        itrjsavecount = 0

        if saveU_format != nothing && update_method != "Fileloading"
            itrj = 0
            itrjstring = lpad(itrj, 8, "0")
            println_verbose_level1(
                U[1],
                "save gaugefields U every $(saveU_every) trajectory",
            )
            #println("save gaugefields U every $(parameters.saveU_every) trajectory")   
            issaved = true
        else
            issaved = false
        end

        return new(issaved, saveU_format, saveU_dir, saveU_every, itrjsavecount)
    end
end

function save_gaugefield(savedata::Savedata, U, itrj)
    if savedata.issaved == false
        return
    end

    if itrj % savedata.saveU_every == 0
        savedata.itrjsavecount += 1
        itrjstring = lpad(itrj, 8, "0")
        if savedata.saveU_format == "JLD"
            filename = savedata.saveU_dir * "/conf_$(itrjstring).jld2"
            saveU(filename, U)
        elseif savedata.saveU_format == "ILDG"
            filename = savedata.saveU_dir * "/conf_$(itrjstring).ildg"
            save_binarydata(U, filename)
        elseif savedata.saveU_format == "BridgeText"
            filename = savedata.saveU_dir * "/conf_$(itrjstring).txt"
            save_textdata(U, filename)
        else
            error("$(savedata.saveU_format) is not supported")
        end
    end
end




end
