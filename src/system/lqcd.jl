module LQCD
import ..Universe_module: Univ
import ..Transform_oldinputfile: transform_to_toml
import ..Parameters_TOML: construct_Params_from_TOML
import ..AbstractUpdate_module: Updatemethod, update!
import Gaugefields: Gradientflow, println_verbose_level1, get_myrank,flow!,
save_binarydata,save_textdata,saveU
import ..AbstractMeasurement_module:Measurement_methods,
calc_measurement_values,measure


function run_LQCD_file(filenamein::String)
    filename_head = splitext(filenamein)[1]
    ext = splitext(filenamein)[end]
    filename = filename_head * ".toml"
    if ext == ".jl"
        transform_to_toml(filenamein)
        println("input file $filenamein is transformed to $filename")
    elseif ext == ".toml"
    else
        @error "$filenamein is not supported. use a TOML format."
    end

    parameters = construct_Params_from_TOML(filename)

    univ = Univ(parameters)

    updatemethod = Updatemethod(parameters,univ)

    

    eps_flow = parameters.eps_flow 
    numflow = parameters.numflow
    Nflow = parameters.Nflow
    dτ = Nflow * eps_flow
    gradientflow = Gradientflow(univ.U, Nflow = 1, eps = eps_flow)  

    measurements = Measurement_methods(univ.U, parameters.measuredir, parameters.measurement_methods)
    measurements_for_flow = Measurement_methods(univ.U, parameters.measuredir, parameters.measurements_for_flow)


    calc_measurement_values(measurements,0, univ.U)
    savedata = Savedata(parameters.saveU_format,parameters.saveU_dir,
    parameters.saveU_every,parameters.update_method,univ.U)



    for itrj = parameters.initialtrj:parameters.Nsteps
        println_verbose_level1(univ.U[1], "# itrj = $itrj")
        @time update!(updatemethod, univ.U)
        save_gaugefield(savedata,univ.U,itrj)

        calc_measurement_values(measurements,itrj, univ.U)


        Usmr = deepcopy(univ.U)
        for istep = 1:numflow
            τ = istep * dτ
            flow!(Usmr, gradientflow)
            additional_string = "$istep $τ "
            for i =1:measurements_for_flow.num_measurements
                interval = measurements_for_flow.intervals[i]
                if istep % interval == 0
                    measure(measurements_for_flow.measurements[i], itrj, Usmr, additional_string = additional_string)
                end
            end
        end
    end

end


mutable struct Savedata
    issaved::Bool
    saveU_format::Union{Nothing,String}
    saveU_dir::String
    saveU_every::Int64
    itrjsavecount::Int64


    function Savedata(saveU_format,saveU_dir,saveU_every,update_method,U)
        itrjsavecount = 0

        if saveU_format != nothing && update_method != "Fileloading"
            itrj = 0
            itrjstring = lpad(itrj, 8, "0")
            println_verbose_level1(U[1], "save gaugefields U every $(saveU_every) trajectory")
            #println("save gaugefields U every $(parameters.saveU_every) trajectory")   
            issaved = true
        else
            issaved = false
        end

        return new(issaved,saveU_format,saveU_dir,saveU_every,itrjsavecount)
    end
end

function save_gaugefield(savedata::Savedata,U,itrj)
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