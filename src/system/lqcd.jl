module LQCD
import ..Universe_module: Univ
import ..Transform_oldinputfile: transform_to_toml
import ..Parameters_TOML: construct_Params_from_TOML
import ..AbstractUpdate_module: Updatemethod, update!
import Gaugefields: Gradientflow, println_verbose_level1, get_myrank
import ..AbstractMeasurement_module:Measurement_methods,
calc_measurement_values


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

    calc_measurement_values(measurements,0, univ.U)

    for itrj = parameters.initialtrj:parameters.Nsteps
        println_verbose_level1(univ.U[1], "# itrj = $itrj")
        @time update!(updatemethod, univ.U)
        calc_measurement_values(measurements,itrj, univ.U)
    end

end





end