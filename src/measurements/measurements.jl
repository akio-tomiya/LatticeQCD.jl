module Measurements
    using LinearAlgebra
    import ..LTK_universe:Universe
    import ..Gaugefields:GaugeFields,set_wing!,substitute!,
            make_staple!,calc_Plaq!,SU3GaugeFields,
            SU2GaugeFields,SU3GaugeFields_1d,SU2GaugeFields_1d,
            GaugeFields_1d,calc_Polyakov,calc_Plaq,calc_Plaq_notrace_1d,SUn,SU2,SU3,TA,add!,
            SUNGaugeFields,SUNGaugeFields_1d,
            Loops,evaluate_loops!,evaluate_loops
    import ..Fermionfields:clear!,FermionFields,WilsonFermion
    import ..Fermionfields:Z4_distribution_fermi!,gauss_distribution_fermi!,set_wing_fermi!
    import ..CGmethods:bicg
    #import ..CGfermion:cg0!
    import ..System_parameters:Params
    import ..Actions:FermiActionParam_Wilson,FermiActionParam_Staggered,FermiActionParam_WilsonClover
    import ..Diracoperators:Dirac_operator
    import ..Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3,println_verbose2,println_verbose1,
            print_verbose1,print_verbose2,print_verbose3
    import ..Smearing:gradientflow!
    import ..Wilsonloops:Wilson_loop,Wilson_loop_set

    #=
    abstract type MeasureMethod end
    
    Base.@kwdef mutable struct Pion_correlator <: MeasureMethod
        measure_every::Int64 = 10
        Dirac_operator::String = "Wilson"
    end
    =#



    struct Measurement{Fermi,GaugeP,FermiP,Gauge_temp}
        gparam::GaugeP
        fparam::Union{Nothing,FermiP}

        _temporal_gauge::Array{Gauge_temp,1}
        _temporal_fermi1::Union{Nothing,Array{Fermi,1}}
        _temporal_fermi2::Union{Nothing,Array{Fermi,1}}

        function Measurement(p::Params,univ::Universe)
            if p.Dirac_operator == p.Dirac_operator_measurement
                return Measurement(univ,Dirac_operator=univ.Dirac_operator)
            else
                if p.Dirac_operator_measurement == nothing
                    fparam = nothing
                else
                    if p.Dirac_operator_measurement == "Wilson"
                        fparam = FermiActionParam_Wilson(p.hop_measurement,p.r_measurement,p.eps_measurement,
                                                            p.Dirac_operator_measurement,p.MaxCGstep_measurement,p.quench)
                    elseif p.Dirac_operator_measurement == "Staggered"
                        fparam = FermiActionParam_Staggered(p.mass_measurement,p.eps_measurement,
                                                            p.Dirac_operator_measurement,p.MaxCGstep_measurement,p.quench,p.Nf)
                    end
                end
                return Measurement(univ,fparam;Dirac_operator=p.Dirac_operator_measurement)
            end
        end


        function Measurement(univ::Universe,gparam,fparam;Dirac_operator=nothing)
            NX = univ.NX
            NY = univ.NY
            NZ = univ.NZ
            NT = univ.NT
            NC = univ.NC

            if NC == 3
                U = Array{SU3GaugeFields,1}(undef,4)
                _temporal_gauge = Array{SU3GaugeFields_1d,1}(undef,4)
            elseif NC == 2
                U = Array{SU2GaugeFields,1}(undef,4)
                _temporal_gauge = Array{SU2GaugeFields_1d,1}(undef,4)
            elseif NC ≥ 4
                U = Array{SUNGaugeFields,1}(undef,4)
                _temporal_gauge = Array{SUNGaugeFields_1d,1}(undef,4)
            
            end

            for i=1:length(_temporal_gauge)
                _temporal_gauge[i] = GaugeFields_1d(NC,NX,NY,NZ,NT) #similar(U[1])
            end

            if Dirac_operator == nothing
                _temporal_fermi1 = nothing
                _temporal_fermi2 = nothing
                φ = nothing
            else
                φ = FermionFields(NC,NX,NY,NZ,NT,fparam,univ.BoundaryCondition)

                _temporal_fermi1 = Array{FermionFields,1}(undef,4)
                _temporal_fermi2 = Array{FermionFields,1}(undef,4)
                for i=1:length(_temporal_fermi1)
                    _temporal_fermi1[i] = similar(φ)
                end
                for i=1:length(_temporal_fermi2)
                    _temporal_fermi2[i] = similar(φ)
                end
            end

            Fermi = typeof(φ)
            GaugeP = typeof(gparam)
            FermiP = typeof(fparam)
            Gauge_temp = eltype(_temporal_gauge)

            return new{Fermi,GaugeP,FermiP,Gauge_temp}(gparam,fparam,_temporal_gauge,_temporal_fermi1,_temporal_fermi2)

        end

        function Measurement(univ::Universe;Dirac_operator=nothing)
            return Measurement(univ::Universe,univ.gparam,univ.fparam,Dirac_operator=Dirac_operator)
        end

        function Measurement(univ::Universe,fparam;Dirac_operator=nothing)
            return Measurement(univ::Universe,univ.gparam,fparam,Dirac_operator=Dirac_operator)
        end
    end

    defaultmeasures = Array{Dict,1}(undef,2)
    for i=1:length(defaultmeasures)
        defaultmeasures[i] = Dict()
    end
    defaultmeasures[1]["methodname"] = "Plaquette"
    defaultmeasures[1]["measure_every"] = 1
    defaultmeasures[1]["fermiontype"] = nothing
    defaultmeasures[2]["methodname"] = "Polyakov_loop"
    defaultmeasures[2]["measure_every"] = 1
    defaultmeasures[2]["fermiontype"] = nothing
    #=
    defaultmeasures[3]["methodname"] = "Chiral_cond" 
    defaultmeasures[3]["measure_every"] = 10
    defaultmeasures[3]["fermiontype"] = "Staggered"
    
    defaultmeasures[4]["methodname"] = "Pion_correlator" 
    defaultmeasures[4]["measure_every"] = 20
    defaultmeasures[4]["fermiontype"] = "Wilson"

    defaultmeasures[5]["methodname"] = "Topologicalcharge"
    defaultmeasures[5]["measure_every"] = 5
    defaultmeasures[5]["fermiontype"] = nothing
    defaultmeasures[5]["numflow"]  = 10
    =#

    struct Measurement_set
        nummeasurement::Int64
        fermions::Array{Measurement,1}
        measurement_methods::Array{Dict,1}
        measurementfps::Array{IOStream,1}

        function Measurement_set(univ::Universe,measurement_dir;measurement_methods=defaultmeasures)
            nummeasurement = length(measurement_methods)
            fermions = Array{Measurement,1}(undef,nummeasurement)
            measurementfps = Array{IOStream,1}(undef,nummeasurement)
            for i=1:nummeasurement
 
                method = measurement_methods[i]
                fermiontype = method["fermiontype"] 

                if method["methodname"] == "Plaquette"
                    measurementfps[i] = open(measurement_dir*"/Plaquette.txt","w")
                elseif method["methodname"] == "Polyakov_loop"
                    measurementfps[i] = open(measurement_dir*"/Polyakov_loop.txt","w")
                elseif method["methodname"] == "Topological_charge"
                    measurementfps[i] = open(measurement_dir*"/Topological_charge.txt","w")
                elseif method["methodname"] == "Chiral_condensate" 
                    measurementfps[i] = open(measurement_dir*"/Chiral_condensate.txt","w")
                elseif method["methodname"] == "Pion_correlator" 
                    measurementfps[i] = open(measurement_dir*"/Pion_correlator.txt","w")
                else
                    error("$(method["methodname"]) is not supported in measurement functions")
                end

                if haskey(method,"eps")
                    eps = method["eps"]
                else
                    eps = 1e-16
                end

                if haskey(method,"MaxCGstep")
                    MaxCGstep = method["MaxCGstep"]
                else
                    MaxCGstep = 5000
                end
                

                if fermiontype  == "Wilson"

                    if haskey(method,"hop")
                        hop = method["hop"]
                    else
                        hop = 0.141139
                        println("Warning. hop = $hop, Default value is used in measurement $(method["methodname"])")
                    end

                    if haskey(method,"r")
                        r = method["r"]
                    else
                        r = 1
                    end

                    

                    quench = false
                    fparam = FermiActionParam_Wilson(hop,r,eps,fermiontype,MaxCGstep,quench)
                elseif fermiontype == "WilsonClover"
                    if haskey(method,"hop")
                        hop = method["hop"]
                    else
                        hop = 0.141139
                        println("Warning. hop = $hop, Default value is used in measurement $(method["methodname"])")
                    end

                    if haskey(method,"r")
                        r = method["r"]
                    else
                        r = 1
                    end

                    if haskey(method,"Clover_coefficient")
                        Clover_coefficient = method["Clover_coefficient"]
                    else
                        Clover_coefficient = 1.5612
                    end


                    NV = univ.NV
                    inn_table= zeros(Int64,NV,4,2)
                    internal_flags = zeros(Bool,2)
                    _ftmp_vectors = Array{Array{ComplexF64,3},1}(undef,6)
                    for i=1:6
                        _ftmp_vectors[i] = zeros(ComplexF64,univ.NC,NV,4)
                    end

                    _is1 = zeros(Int64,NV)
                    _is2 = zeros(Int64,NV)

                    quench = false
                    fparam = FermiActionParam_WilsonClover(hop,r,eps,fermiontype,MaxCGstep,Clover_coefficient,
                                    internal_flags,inn_table,_ftmp_vectors,_is1,_is2,
                                    quench)
                elseif fermiontype == "Staggered"
                    if haskey(method,"mass")
                        mass = method["mass"]
                    else
                        mass = 0.5
                        println("Warning. mass = $mass, Default value is used in measurement $(method["methodname"])")
                    end

                    if haskey(method,"Nf")
                        Nf = method["Nf"]
                    else                        
                        error("Nf should be set if you want to use the staggered fermion in measurements")
                        #println("Warning. mass = $hop, Default value is used")
                    end

                    quench = false
                    
                    #println("Measurement_set::mass_measurement = $(p.mass_measurement)")
                    fparam = FermiActionParam_Staggered(mass,eps,fermiontype,MaxCGstep,quench,Nf)
                                
                elseif fermiontype == nothing
                    fparam = nothing
                else
                    error("$fermiontype is not supported. use Wilson or Staggered")
                end

                fermions[i] = Measurement(univ,fparam;Dirac_operator=fermiontype )
            end
            return new(nummeasurement,fermions,measurement_methods,measurementfps)
        end
    end



    function measurements(itrj,U,univ,measset::Measurement_set;verbose = Verbose_2())
        plaq = 0.0
        poly = 0.0 + 0.0im
        for i = 1:measset.nummeasurement
            method = measset.measurement_methods[i]
            measfp = measset.measurementfps[i]
            #println(method)
            #println(method["measure_every"])
            if itrj % method["measure_every"] == 0
                println_verbose1(verbose,"-----------------")
                if method["methodname"] == "Plaquette"
                    plaq = calc_plaquette(U)
                    println_verbose1(verbose,"$itrj $plaq # plaq")
                    println(measfp,"$itrj $plaq # plaq")
                elseif method["methodname"] == "Polyakov_loop"
                    poly = calc_Polyakov(U)
                    println_verbose1(verbose,"$itrj $(real(poly)) $(imag(poly)) # poly")
                    println(measfp,"$itrj $(real(poly)) $(imag(poly)) # poly")
                elseif method["methodname"] == "Topological_charge"
                    Nflowsteps = method["Nflowsteps"]#50
                    eps_flow =  method["eps_flow"] #0.01

                    println_verbose2(verbose,"# epsilon for the Wilson flow is $eps_flow")
                    Usmr = deepcopy(U)
                    W1 = deepcopy(univ.U)
                    W2 = deepcopy(univ.U)

                    temp_UμνTA = Array{GaugeFields_1d,2}(undef,4,4)

                    τ = 0.0
                    plaq = calc_plaquette(Usmr)
                    #Q = calc_topological_charge(Usmr)
                    # sign of topological charge defined to be positive for one-instanton.
                    Qplaq = calc_topological_charge_plaq(Usmr,temp_UμνTA)
                    Qclover= calc_topological_charge_clover(Usmr)
                    Qimproved= calc_topological_charge_improved(Usmr,temp_UμνTA,Qclover)
                    #println_verbose1(verbose,"$itrj $τ $plaq $(real(Qplaq)) $(real(Qclover)) #flow itrj flowtime plaq Qplaq Qclov")
                    #println(measfp,"$itrj $τ $plaq $(real(Qplaq)) $(real(Qclover)) #flow itrj flowtime plaq Qplaq Qclov")
                    println_verbose1(verbose,"$itrj $τ $plaq $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #flow itrj flowtime plaq Qplaq Qclov Qimproved")
                    println(measfp,"$itrj $τ $plaq $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #flow itrj flowtime plaq Qplaq Qclov Qimproved")
                    flush(stdout)
                    smearing_type = "gradient_flow"
                    if smearing_type == "gradient_flow"
                        for iflow = 1:method["numflow"]#5000 # eps=0.01: t_flow = 50
                            gradientflow!(Usmr,univ,W1,W2,Nflowsteps,eps_flow)
                            plaq = calc_plaquette(Usmr)
                            Qplaq = calc_topological_charge_plaq(Usmr,temp_UμνTA)
                            Qclover= calc_topological_charge_clover(Usmr,temp_UμνTA)
                            Qimproved= calc_topological_charge_improved(Usmr,temp_UμνTA,Qclover)
                            #@time Q = calc_topological_charge(Usmr)
                            τ = iflow*eps_flow*Nflowsteps
                            println_verbose1(verbose,"$itrj $τ $plaq $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #flow itrj flowtime plaq Qplaq Qclov Qimproved")
                            println(measfp,"$itrj $τ $plaq $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #flow itrj flowtime plaq Qplaq Qclov Qimproved")
                            if iflow%10 == 0
                                flush(stdout)
                            end
                        end
                    elseif smearing_type == "APE" # TODO this should be commonized to gradient flow
                        α = 6/13 # magic number from hep-lat/9907019 for reproduction. This should be an input parameter
                        for iflow = 1:method["numflow"]
                            calc_fatlink_APE!(Usmr,Usmr,α,α,normalize_method="special unitary")
                            plaq = calc_plaquette(Usmr)
                            Qplaq = calc_topological_charge_plaq(Usmr,temp_UμνTA)
                            Qclover= calc_topological_charge_clover(Usmr,temp_UμνTA)
                            Qimproved= calc_topological_charge_improved(Usmr,temp_UμνTA,Qclover)
                            τ = iflow
                            println_verbose1(verbose,"$itrj $τ $plaq $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #ape itrj flowtime plaq Qplaq Qclov Qimproved")
                            println(measfp,"$itrj $τ $plaq $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #ape itrj flowtime plaq Qplaq Qclov Qimproved")
                            if iflow%10 == 0
                                flush(stdout)
                            end
                        end
                        ρ = α/13 
                        Usmr = deepcopy(U)
                        for iflow = 1:method["numflow"]
                            calc_stout!(Usmr,Usmr,ρ)
                            plaq = calc_plaquette(Usmr)
                            Qplaq = calc_topological_charge_plaq(Usmr,temp_UμνTA)
                            Qclover= calc_topological_charge_clover(Usmr,temp_UμνTA)
                            Qimproved= calc_topological_charge_improved(Usmr,temp_UμνTA,Qclover)
                            τ = iflow
                            println_verbose1(verbose,"$itrj $τ $plaq $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #stout itrj flowtime plaq Qplaq Qclov Qimproved")
                            println(measfp,"$itrj $τ $plaq $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #stout itrj flowtime plaq Qplaq Qclov Qimproved")
                            if iflow%10 == 0
                                flush(stdout)
                            end
                        end
                    else
                        error("Invalid smearing_type = $gradient_flow")
                    end
                elseif method["methodname"] == "Chiral_condensate" 
                    #fermiontype = method["fermiontype"]
                    measure_chiral_cond(univ,measset.fermions[i],itrj,measfp,verbose)
                elseif method["methodname"] == "Pion_correlator" 
                    #fermiontype = method["fermiontype"]
                    #calc_pion_correlator(univ,measset.fermions[i])
                    measure_correlator(univ,measset.fermions[i],itrj,measfp,verbose)
                else
                    error("$(method["methodname"]) is not supported")
                end
                println_verbose1(verbose,"-----------------")
                flush(measfp)
            end
        end
        return plaq,poly
    end




    function calc_polyakovloop(univ::Universe)
        poly = calc_Polyakov(univ.U)

        return poly

    end

    function calc_factor_plaq(U)
        factor = 2/(U[1].NV*4*3*U[1].NC)
    end

    function calc_plaquette(U::Array{T,1}) where T <: GaugeFields
        plaq = 0
        factor = calc_factor_plaq(U)
        plaq = calc_Plaq(U)*factor
        return real(plaq)
    end

    function calc_plaquette(univ::Universe,U::Array{T,1}) where T <: GaugeFields
        plaq = 0
        temps = univ._temporal_gauge
        factor = calc_factor_plaq(U)
        plaq = calc_Plaq!(U,temps)*factor
        return real(plaq)
    end

    function calc_topological_charge(U::Array{GaugeFields{S},1}) where S <: SUn
        return calc_topological_charge_plaq(U)
    end

    function calc_topological_charge_clover(U::Array{GaugeFields{S},1}) where S <: SUn
        UμνTA = Array{GaugeFields_1d,2}(undef,4,4)
        Q = calc_topological_charge_clover(U,UμνTA)
        return Q
    end

    function calc_topological_charge_clover(U::Array{GaugeFields{S},1},temp_UμνTA) where S <: SUn
        UμνTA = temp_UμνTA
        numofloops = calc_UμνTA!(UμνTA,"clover",U)
        Q = calc_Q(UμνTA,numofloops,U)
        return Q
    end

    function calc_topological_charge_improved(U::Array{GaugeFields{S},1},temp_UμνTA) where S <: SUn
        UμνTA = temp_UμνTA
        numofloops = calc_UμνTA!(UμνTA,"clover",U)
        Qclover = calc_Q(UμνTA,numofloops,U)
        Q = calc_topological_charge_improved(U,UμνTA,Qclover)
        return Q
    end

    function calc_topological_charge_improved(U::Array{GaugeFields{S},1},temp_UμνTA,Qclover) where S <: SUn
        UμνTA = temp_UμνTA
        #numofloops = calc_UμνTA!(UμνTA,"clover",U)
        #Qclover = calc_Q(UμνTA,numofloops,U)

        numofloops = calc_UμνTA!(UμνTA,"rect",U)
        Qrect = 2*calc_Q(UμνTA,numofloops,U)
        c1 = -1/12
        c0 = 5/3
        Q = c0*Qclover + c1*Qrect
        return Q
    end


    function calc_topological_charge_plaq(U::Array{GaugeFields{S},1},temp_UμνTA) where S <: SUn
        UμνTA = temp_UμνTA
        numofloops = calc_UμνTA!(UμνTA,"plaq",U)
        Q = calc_Q(UμνTA,numofloops,U)
        return Q
    end

    function calc_topological_charge_plaq(U::Array{GaugeFields{S},1}) where S <: SUn
        UμνTA = Array{GaugeFields_1d,2}(undef,4,4)
        return calc_topological_charge_plaq(U,UμνTA )
    end

    function calc_loopset_μν(name)
        loops_μν= Array{Wilson_loop_set,2}(undef,4,4)
        if name == "plaq"
            numofloops = 1
            for μ=1:4
                for ν=1:4
                    if ν == μ
                        continue
                    end
                    loops = Wilson_loop_set()
                    loop = Wilson_loop([(μ,1),(ν,1),(μ,-1),(ν,-1)])
                    push!(loops,loop)
                    loops_μν[μ,ν] = loops
                end
            end
        elseif name == "clover"
            numofloops = 4
            for μ=1:4
                for ν=1:4
                    if ν == μ
                        continue
                    end
                    loops = Wilson_loop_set()

                    loop_righttop = Wilson_loop([(μ,1),(ν,1),(μ,-1),(ν,-1)])
                    loop_lefttop = Wilson_loop([(ν,1),(μ,-1),(ν,-1),(μ,1)])
                    loop_rightbottom = Wilson_loop([(ν,-1),(μ,1),(ν,1),(μ,-1)])
                    loop_leftbottom= Wilson_loop([(μ,-1),(ν,-1),(μ,1),(ν,1)])
                    push!(loops,loop_righttop)
                    push!(loops,loop_lefttop)
                    push!(loops,loop_rightbottom)
                    push!(loops,loop_leftbottom)

                    loops_μν[μ,ν] = loops
                end
            end
        elseif name == "rect"
            numofloops = 8
            for μ=1:4
                for ν=1:4
                    if ν == μ
                        continue
                    end
                    loops = Wilson_loop_set()

                    loop_righttop = Wilson_loop([(μ,2),(ν,1),(μ,-2),(ν,-1)])
                    loop_lefttop = Wilson_loop([(ν,1),(μ,-2),(ν,-1),(μ,2)])
                    loop_rightbottom = Wilson_loop([(ν,-1),(μ,2),(ν,1),(μ,-2)])
                    loop_leftbottom= Wilson_loop([(μ,-2),(ν,-1),(μ,2),(ν,1)])
                    push!(loops,loop_righttop)
                    push!(loops,loop_lefttop)
                    push!(loops,loop_rightbottom)
                    push!(loops,loop_leftbottom)

                    loop_righttop = Wilson_loop([(μ,1),(ν,2),(μ,-1),(ν,-2)])
                    loop_lefttop = Wilson_loop([(ν,2),(μ,-1),(ν,-2),(μ,1)])
                    loop_rightbottom = Wilson_loop([(ν,-2),(μ,1),(ν,2),(μ,-1)])
                    loop_leftbottom= Wilson_loop([(μ,-1),(ν,-2),(μ,1),(ν,2)])
                    push!(loops,loop_righttop)
                    push!(loops,loop_lefttop)
                    push!(loops,loop_rightbottom)
                    push!(loops,loop_leftbottom)

                    loops_μν[μ,ν] = loops
                end
            end
        else
            error("$name is not supported")
        end
        return loops_μν,numofloops
    end



    function calc_UμνTA!(temp_UμνTA,loops_μν,U)
        UμνTA = temp_UμνTA
        for μ=1:4
            for ν=1:4
                if ν == μ
                    continue
                end
                loopset = Loops(U,loops_μν[μ,ν])
                UμνTA[μ,ν] = evaluate_loops(loopset,U)
                UμνTA[μ,ν] = TA(UμνTA[μ,ν])
            end
        end
        return 
    end

    function calc_UμνTA!(temp_UμνTA,name::String,U)
        loops_μν,numofloops = calc_loopset_μν(name)
        calc_UμνTA!(temp_UμνTA,loops_μν,U)
        return numofloops
    end

    #=
    implementation of topological charge is based on
    https://arxiv.org/abs/1509.04259
    =#
    function calc_Q(UμνTA,numofloops,U::Array{GaugeFields{S},1}) where S <: SUn
        if S == SU3
            NC = 3
        elseif S == SU2
            NC = 2
        else
            NC = U[1].NC
            #error("NC != 2,3 is not supported")
        end

        Q = 0.0

        ε(μ,ν,ρ,σ) = epsilon_tensor(μ,ν,ρ,σ)  
        NV=UμνTA[1,2].NV
        for n=1:NV
            for μ=1:4
                for ν=1:4
                    if ν == μ
                        continue
                    end
                    for ρ =1:4
                        for σ=1:4
                            if ρ == σ
                                continue
                            end
                            
                            s = 0im
                            for i=1:NC
                                for j=1:NC
                                    s += UμνTA[μ,ν][i,j,n]*UμνTA[ρ,σ][j,i,n]
                                end
                            end

                            Q += ε(μ,ν,ρ,σ)*s/numofloops^2 #*tr(tmp1*tmp2) 
                        end
                    end
                end
            end
        end
        return -Q/(32*(π^2))
    end



    

    function TAm(M)
        NC = size(M)[1]
        AM = (M - M')/2
        t = tr(AM)/NC
        for i=1:NC
            AM[i,i] -= t
        end
        return AM
    end
    
    function calc_plaquette(univ::Universe)
        return calc_plaquette(univ,univ.U)
    end

    function spincolor(ic,is,univ::Universe)
        return ic-1 + (is-1)*univ.NC + 1
    end

    function calc_chiral_cond(univ::Universe,meas,measfp,Nr = 10, verbose = Verbose_2())
        #error("calc_chiral_cond is not implemented!")
        #
        # pbp = (1/Nr) Σ_i p_i
        # p_i = r_i^\dag xi_i
        # xi_i = D^{-1} r_i   # D xi = r : r is a random veccor
        #
        pbp = 0.0
        # setup a massive Dirac operator
        M = Dirac_operator(univ.U,meas._temporal_fermi2[1],meas.fparam)
        for ir=1:Nr
            r = similar(meas._temporal_fermi2[1]) 
            p = similar(r) 
            clear!(p)
            #gauss_distribution_fermi!(r,univ.ranf)
            Z4_distribution_fermi!(r)
            #set_wing_fermi!(r) 
            bicg(p,M,r,eps=meas.fparam.eps,maxsteps = meas.fparam.MaxCGstep,verbose = verbose) # solve Mp=b, we get p=M^{-1}b
            tmp = r*p # hemitian inner product
            println_verbose2(verbose,"$(real(tmp)/univ.NV) # chiralcond $ir")
            println(measfp,"$(real(tmp)/univ.NV) # chiralcond $ir")
            pbp+=tmp 
        end
        return real(pbp/Nr)/univ.NV
    end
    function measure_chiral_cond(univ::Universe,meas::Measurement,itrj,measfp,verbose = Verbose_2())
        Nr = 10
        pbp = calc_chiral_cond(univ,meas,measfp,Nr,verbose)
        println_verbose1(verbose,"$itrj $pbp # pbp")
        flush(stdout)
    end

    function calc_quark_propagators_point_source_each(M,meas,i,NC,verbose)
        # calculate D^{-1} for a given source at the origin.
        # Nc*4 elements has to be gathered.
        b = similar(meas._temporal_fermi2[1]) # source is allocated
        p = similar(b) # sink is allocated (propagator to sink position)
        #k = meas._temporal_fermi2[2]
        clear!(b)
        is = ((i-1) % 4)+1 # spin index 
        ic = ((i-is) ÷ 4)+ 1 # color index
        println_verbose1(verbose,"$ic $is")
        #println("$ic $is")
        b[ic,1,1,1,1,is]=1 # source at the origin

        @time bicg(p,M,b) # solve Mp=b, we get p=M^{-1}b

        #@time cg0!(k,b,1, univ.U, meas._temporal_gauge, meas._temporal_fermi1, meas.fparam) # k[x] = M^{-1}b[0]
        println_verbose1(verbose,"Hadron spectrum: Inversion $(i)/$(NC*4) is done")
        #println("Hadron spectrum: Inversion $(i)/$(NC*4) is done")
        flush(stdout)
        return p
    end

    function calc_quark_propagators_point_source(D,meas,NC,verbose)
        # D^{-1} for each spin x color element
        propagators = map(i -> calc_quark_propagators_point_source_each(D,meas,i,NC,verbose),1:NC*4)
        return propagators
    end

    function calc_pion_correlator(univ::Universe,meas::Measurement,verbose)
        println_verbose2(verbose,"Hadron spectrum started")
        #println("Hadron spectrum started")
        #U = univ._temporal_gauge
        #if univ.NC != 3
        #    error("not implemented yet for calc_pion_correlator")
        #end
        #b = meas._temporal_fermi2[1] # source allocate
        #k = meas._temporal_fermi2[2] # sink allocate

        # setup a massive Dirac operator
        M = Dirac_operator(univ.U,meas._temporal_fermi2[1],meas.fparam)

        # Allocate the Wilson matrix = S
        S = zeros(ComplexF64, (univ.NX,univ.NY,univ.NZ,univ.NT, 4*univ.NC, 4*univ.NC) )

        # calculate quark propagators from a point source at he origin
        propagators = calc_quark_propagators_point_source(M,meas,univ.NC,verbose)

        #ctr = 0 # a counter
        for ic=1:univ.NC
            for is=1:4
                icum = (ic-1)*4+is
                
                #=
                clear!(b)
                b[ic,1,1,1,1,is]=1 # source at the origin
                @time cg0!(k,b,1, univ.U, meas._temporal_gauge, meas._temporal_fermi1, meas.fparam) # k[x] = M^{-1}b[0]
                =#
                propagator = propagators[icum]
                α0=spincolor(ic,is,univ) # source(color-spinor) index
                #
                
                # reconstruction
                for t=1:univ.NT
                    for z=1:univ.NZ
                        for y=1:univ.NY
                            for x=1:univ.NX
                                for ic2=1:univ.NC 
                                    for is2=1:4 # 4 is the number of spinor index in 4d.
                                        β=spincolor(ic2,is2,univ)
                                        S[x,y,z,t,α0,β]+= propagator[ic,x,y,z,t,is]
                                    end
                                end
                            end
                        end
                    end
                end
                # end for the substitution
                
                #ctr+=1
            end 
        end
        # contruction end.

        println_verbose2(verbose,"Hadron spectrum: Reconstruction")
        #println("Hadron spectrum: Reconstruction")
        Cpi = zeros( univ.NT )
        # Construct Pion propagator 
        for t=1:univ.NT
            tmp = 0.0+0.0im
            for z=1:univ.NZ
                for y=1:univ.NY
                    for x=1:univ.NX
                        for ic=1:univ.NC 
                            for is=1:4 # 4 is the number of spinor index in 4d.
                                α=spincolor(ic,is,univ)
                                for ic2=1:univ.NC 
                                    for is2=1:4 # 4 is the number of spinor index in 4d.
                                        β=spincolor(ic2,is2,univ)
                                        tmp += S[x,y,z,t,α,β] * S[x,y,z,t,α,β]'#inner product.
                                        # complex conjugate = g5 S g5.
                                    end
                                end
                                # complex conjugate = g5 S g5.
                            end
                        end
                    end
                end
            end
            Cpi[t] = real(tmp)
        end
        #println(typeof(verbose),"\t",verbose)
        println_verbose2(verbose,"Hadron spectrum end")
        #println("Hadron spectrum end")
        return Cpi
    end
    function measure_correlator(univ::Universe,meas::Measurement,itrj,measfp,verbose)
        C = calc_pion_correlator(univ,meas,verbose)
        print_verbose1(verbose,"$itrj ")
        print(measfp,"$itrj ")
        for it=1:length(C)
            cc = C[it]
            print_verbose1(verbose,"$cc ")
            print(measfp,"$cc ")
        end
        println_verbose1(verbose,"#pioncorrelator")
        println(measfp,"#pioncorrelator")
    end
    #topological charge
    function epsilon_tensor(mu::Int,nu::Int,rho::Int,sigma::Int) 
        sign=1 # (3) 1710.09474 extended epsilon tensor
        if mu < 0
            sign*=-1
            mu=-mu
        end
        if nu < 0
            sign*=-1
            nu=-nu
        end
        if rho < 0
            sign*=-1
            rh=-rho
        end
        if sigma < 0
            sign*=-1
            sigma=-sigma
        end
        epsilon = zeros(Int,4,4,4,4)
        epsilon[ 1, 2, 3, 4 ] = 1
        epsilon[ 1, 2, 4, 3 ] = -1
        epsilon[ 1, 3, 2, 4 ] = -1
        epsilon[ 1, 3, 4, 2 ] = 1
        epsilon[ 1, 4, 2, 3 ] = 1
        epsilon[ 1, 4, 3, 2 ] = -1
        epsilon[ 2, 1, 3, 4 ] = -1
        epsilon[ 2, 1, 4, 3 ] = 1
        epsilon[ 2, 3, 1, 4 ] = 1
        epsilon[ 2, 3, 4, 1 ] = -1
        epsilon[ 2, 4, 1, 3 ] = -1
        epsilon[ 2, 4, 3, 1 ] = 1
        epsilon[ 3, 1, 2, 4 ] = 1
        epsilon[ 3, 1, 4, 2 ] = -1
        epsilon[ 3, 2, 1, 4 ] = -1
        epsilon[ 3, 2, 4, 1 ] = 1
        epsilon[ 3, 4, 1, 2 ] = 1
        epsilon[ 3, 4, 2, 1 ] = -1
        epsilon[ 4, 1, 2, 3 ] = -1
        epsilon[ 4, 1, 3, 2 ] = 1
        epsilon[ 4, 2, 1, 3 ] = 1
        epsilon[ 4, 2, 3, 1 ] = -1
        epsilon[ 4, 3, 1, 2 ] = -1
        epsilon[ 4, 3, 2, 1 ] = 1
        return epsilon[mu,nu,rho,sigma]*sign
    end

end