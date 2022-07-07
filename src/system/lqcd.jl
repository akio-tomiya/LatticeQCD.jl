module LQCD
import ..Universe_module: Univ
import ..Transform_oldinputfile: transform_to_toml
import ..Parameters_TOML: construct_Params_from_TOML
import ..AbstractUpdate_module: Updatemethod, update!
import Gaugefields: Gradientflow, println_verbose_level1, get_myrank,flow!
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

    for itrj = parameters.initialtrj:parameters.Nsteps
        println_verbose_level1(univ.U[1], "# itrj = $itrj")
        @time update!(updatemethod, univ.U)
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





end