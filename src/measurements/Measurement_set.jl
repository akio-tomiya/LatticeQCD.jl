function set_parameter_default(method,key,defaultvalue)
    if haskey(method,key)
        return method[key]
    else
        return defaultvalue
    end
end

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

                baremeasure = set_parameter_default(method,"bare measure",true)
                flowmeasure = set_parameter_default(method,"flow measure",true)

            elseif method["methodname"] == "Polyakov_loop"
                measurements[i] = Polyakov_measurement(U,filename=measurement_dir*"/Polyakov_loop.txt.txt")
                baremeasure = set_parameter_default(method,"bare measure",true)
                flowmeasure = set_parameter_default(method,"flow measure",true)
            elseif method["methodname"] == "Topological_charge"
                TC_methods = set_parameter_default(method,"TC_methods",["plaquette","clover"])
                measurements[i] = Topological_charge_measurement(U, filename=measurement_dir*"/Topological_charge.txt",
                                             TC_methods = TC_methods)
                baremeasure = set_parameter_default(method,"bare measure",true)
                flowmeasure = set_parameter_default(method,"flow measure",true)

            elseif method["methodname"] == "Energy_density"
                measurements[i] = Energy_density_measurement(U, filename=measurement_dir*"/Energy_density.txt")
                baremeasure = set_parameter_default(method,"bare measure",true)
                flowmeasure = set_parameter_default(method,"flow measure",true)
                
            elseif method["methodname"] == "Chiral_condensate" 
                baremeasure = set_parameter_default(method,"bare measure",true)
                flowmeasure = set_parameter_default(method,"flow measure",false)
                fermiontype = set_parameter_default(method,"fermiontype","Staggered")
                mass = set_parameter_default(method,"mass",0.1)
                Nf = set_parameter_default(method,"Nf",2)
                κ = set_parameter_default(method,"hop",0.141139)
                r = set_parameter_default(method,"r",1)
                M = set_parameter_default(method,"Domainwall_M",-1)
                L5 = set_parameter_default(method,"Domainwall_L5",2)
                mass = set_parameter_default(method,"Domainwall_m",mass)
                BoundaryCondition = set_parameter_default(method,"BoundaryCondition",nothing)
                eps_CG = set_parameter_default(method,"eps",1e-14)
                MaxCGstep = set_parameter_default(method,"MaxCGstep",3000)
                Nr = set_parameter_default(method,"Nr",10)



                measurements[i] = Chiral_condensate_measurement(U,
                                        filename = filename=measurement_dir*"/Chiral_condensate.txt.txt",
                                        fermiontype = fermiontype,
                                        mass = mass,
                                        Nf = Nf,
                                        κ = κ,
                                        r = r,
                                        L5 = L5,
                                        M = M,
                                        eps_CG = eps_CG,
                                        MaxCGstep = MaxCGstep,
                                        BoundaryCondition = BoundaryCondition,
                                        Nr = Nr
                ) 

                #measurements[i] = Measure_chiral_condensate(measurement_dir*"/Chiral_condensate.txt",univ.U,method)
                
            #=
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
    additional_string = "0 0.0 "
    for i in m.baremesurement_indices
        measure(m.measurements[i],itrj,U,additional_string= additional_string)
    end
end

function measure_withflow(m::Measurements_set,itrj,U,smearing::Abstractsmearing,numstep,dτ)
    Usmr = deepcopy(U)
    for istep = 1:numstep
        τ = istep*dτ
        flow!(Usmr,smearing)
        additional_string = "$istep $τ "
        for i in m.flowmeasurement_indices
            measure(m.measurements[i],itrj,Usmr,additional_string= additional_string)
        end
    end
end