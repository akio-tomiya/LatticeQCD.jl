
import ..Parameter_structs:construct_Measurement_parameters_from_dict,Measurement_parameters


struct Measurement_methods
    measurement_parameters_set::Vector{Measurement_parameters}
    measurements::Vector{AbstractMeasurement}
    num_measurements::Int64
    intervals::Vector{Int64}
end

function calc_measurement_values(m::Measurement_methods,itrj,U;additional_string ="")
    for i =1:m.num_measurements
        interval = m.intervals[i]
        if itrj % interval == 0
            measure(m.measurements[i], itrj, U, additional_string = additional_string)
        end
    end
end

function prepare_measurement(U,measurement_parameters::T,filename) where T
    if T == Plaq_parameters
        measurement = Plaquette_measurement(U,measurement_parameters,filename)
    elseif T == Poly_parameters
        measurement = Polyakov_measurement(U,measurement_parameters,filename)
    elseif T == TopologicalCharge_parameters
        measurement = Topological_charge_measurement(U,measurement_parameters,filename)
    elseif T == ChiralCondensate_parameters
        measurement = Chiral_condensate_measurement(U,measurement_parameters,filename)
    elseif T == Pion_parameters
        measurement = Pion_correlator_measurement(U,measurement_parameters,filename)
    elseif T == Energy_density_parameters
        measurement = Energy_density_measurement(U,measurement_parameters,filename)
    else
        error(T, " is not supported in measurements")
    end
    return measurement
end

function Measurement_methods(U, measurement_dir, measurement_methods::T) where T <: Vector{Dict}
    #println( measurement_methods)
    nummeasurements = length(measurement_methods)
    measurements = Vector{AbstractMeasurement}(undef, nummeasurements)
    measurement_parameters_set = Vector{Measurement_parameters}(undef,nummeasurements)
    intervals = zeros(Int64,nummeasurements)

    for (i,method) in enumerate(measurement_methods)
        measurement_parameters = construct_Measurement_parameters_from_dict(method)
        #println(measurement_parameters)
        intervals[i] = measurement_parameters.measure_every
        
        filename = measurement_dir*measurement_parameters.methodname*".txt"
        measurements[i] = prepare_measurement(U,measurement_parameters,filename)
        measurement_parameters_set[i] = deepcopy(measurement_parameters)
        
    end

    for i=1:nummeasurements
        println(measurement_parameters_set[i] )
    end

    return Measurement_methods(measurement_parameters_set,measurements,nummeasurements,intervals)
end
