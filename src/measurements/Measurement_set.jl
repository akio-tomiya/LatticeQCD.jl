struct Measurements_set
    nummeasurements::Int64
    measurements::Array{AbstractMeasurement,1}
    measurement_methods::Array{Dict,1}
    methodnames::Array{String,1}
    baremesurement_indices::Array{Int64,1}
    flowmeasurement_indices::Array{Int64,1}

    function Measurements_set(U,measurement_dir,measurement_methods)
        nummeasurements = length(measurement_methods)
        measurements = Array{AbstractMeasurement,1}(undef,nummeasurements)
        methodnames = Array{String,1}(undef,nummeasurements)
        baremesurement_indices = Int64[]
        flowmeasurement_indices = Int64[]

        for i=1:nummeasurements
            method = measurement_methods[i]
            methodnames[i] = method["methodname"]

            if method["methodname"] == "Plaquette"
                measurements[i] = Plaquette_measurement(U,filename=measurement_dir*"/Plaquette.txt")

                if haskey(method,"bare measure")
                    baremeasure = method["bare measure"]
                else
                    baremeasure = true
                end
                if haskey(method,"flow measure")
                    flowmeasure = method["flow measure"]
                    
                else
                    flowmeasure = true
                end

            elseif method["methodname"] == "Polyakov_loop"
                measurements[i] = Polyakov_measurement(U,filename=measurement_dir*"/Polyakov_loop.txt.txt")
                if haskey(method,"bare measure")
                    baremeasure = method["bare measure"]
                else
                    baremeasure = true
                end
                if haskey(method,"flow measure")
                    flowmeasure = method["flow measure"]
                    
                else
                    flowmeasure = true
                end

            elseif method["methodname"] == "Topological_charge"
                if haskey(method,"TC_methods")
                    TC_methods = method["TC_methods"]
                else
                    TC_methods = ["plaquette","clover"]
                end
                measurements[i] = Topological_charge_measurement(U, filename=measurement_dir*"/Topological_charge.txt",
                                             TC_methods = TC_methods)
                if haskey(method,"bare measure")
                    baremeasure = method["bare measure"]
                else
                    baremeasure = true
                end
                if haskey(method,"flow measure")
                    flowmeasure = method["flow measure"]
                    
                else
                    flowmeasure = true
                end
            elseif method["methodname"] == "Energy_density"
                measurements[i] = Energy_density_measurement(U, filename=measurement_dir*"/Energy_density.txt")
                if haskey(method,"bare measure")
                    baremeasure = method["bare measure"]
                else
                    baremeasure = true
                end
                if haskey(method,"flow measure")
                    flowmeasure = method["flow measure"]
                    
                else
                    flowmeasure = true
                end
                #=
            elseif method["methodname"] == "Chiral_condensate" 
                measurements[i] = Measure_chiral_condensate(measurement_dir*"/Chiral_condensate.txt",univ.U,method)
            elseif method["methodname"] == "Pion_correlator" 
                measurements[i] = Measure_Pion_correlator(measurement_dir*"/Pion_correlator.txt",univ.U,method)
                =#
            else
                error("$(method["methodname"]) is not supported")
            end

            if baremeasure
                push!(baremesurement_indices,i)
            end

            if flowmeasure
                push!(flowmeasurement_indices,i)
            end

        end


        return new( nummeasurements,#::Int64
                    measurements,#::Array{AbstractMeasurement,1}
                    measurement_methods,#::Array{Dict,1}
                    methodnames,#::Array{String,1}
                    baremesurement_indices,#::Array{Int64,1}
                    flowmeasurement_indices,#::Array{Int64,1}

        )
    end


end


function measure(m::Measurements_set,itrj,U)
    for i in m.baremesurement_indices
        measure(m.measurements[i],itrj,U)
    end
end

function measure_withflow()
end