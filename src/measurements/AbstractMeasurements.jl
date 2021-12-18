module AbstractMeasurement_module
    import ..Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3,println_verbose2,println_verbose1,
            print_verbose1,print_verbose2,print_verbose3
    import ..Gaugefield:AbstractGaugefields
    import ..Actions:GaugeActionParam_autogenerator,GaugeActionParam
    import ..Actions:FermiActionParam_Wilson,FermiActionParam_Staggered,FermiActionParam_WilsonClover,
                FermiActionParam

    abstract type AbstractMeasurement end

    function measure(measurement::M,itrj,U::Array{T,1};verbose = Verbose_2()) where {M <: AbstractMeasurement,T <: AbstractGaugefields}
        error("measure is not implemented in $(typeof(measurement))")
    end

    function set_params(dict,string,default)
        if haskey(dict,string)
            return dict[string] 
        else
            return default
        end
    end


    function set_Wilson(params)
        if haskey(params,"hop")
            hop = params["hop"]
        else
            hop = 0.141139
            println("Warning. hop = $hop, Default value is used in measurement $(params["paramsname"])")
        end

        if haskey(params,"r")
            r = params["r"]
        else
            r = 1
        end

        quench = false

        smearing_for_fermion = set_params(params,"smearing_for_fermion","nothing")
        stout_numlayers = set_params(params,"stout_numlayers",nothing)
        stout_ρ = set_params(params,"stout_ρ",nothing)
        stout_loops = set_params(params,"stout_loops",nothing)


        if smearing_for_fermion == "nothing"
            fparam = FermiActionParam_Wilson(hop,r,eps,fermiontype,MaxCGstep,quench)
        else

            L = (univ.NX,univ.NY,univ.NZ,univ.NT)
            fparam = FermiActionParam_Wilson(hop,r,eps,fermiontype,MaxCGstep,quench,
                                                smearingparameters = "stout",
                                                loops_list = stout_loops,
                                                coefficients  = stout_ρ,
                                                numlayers = stout_numlayers,
                                                L = L)
        end
        return fparam
    end

    function set_WilsonClover(params)
        if haskey(params,"hop")
            hop = params["hop"]
        else
            hop = 0.141139
            println("Warning. hop = $hop, Default value is used in measurement $(params["paramsname"])")
        end

        if haskey(params,"r")
            r = params["r"]
        else
            r = 1
        end

        if haskey(params,"Clover_coefficient")
            Clover_coefficient = params["Clover_coefficient"]
        else
            Clover_coefficient = 1.5612
        end


        NV = U[1].NV
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

        smearing_for_fermion = set_params(params,"smearing_for_fermion","nothing")
        stout_numlayers = set_params(params,"stout_numlayers",nothing)
        stout_ρ = set_params(params,"stout_ρ",nothing)
        stout_loops = set_params(params,"stout_loops",nothing)

        if smearing_for_fermion == "nothing"
            FermiActionParam_WilsonClover(hop,r,eps,fermiontype,MaxCGstep,Clover_coefficient,
                        internal_flags,inn_table,_ftmp_vectors,_is1,_is2,
                        quench)
        else
            error("stout for WilsonClover is not supported yet!")
            L = (univ.NX,univ.NY,univ.NZ,univ.NT)
            FermiActionParam_WilsonClover(hop,r,eps,fermiontype,MaxCGstep,Clover_coefficient,
                        internal_flags,inn_table,_ftmp_vectors,_is1,_is2,
                        quench,
                        smearingparameters = "stout",
                        loops_list = stout_loops,
                        coefficients  = stout_ρ,
                        numlayers = stout_numlayers,
                        L = L)
        end
        return fparam
    end

    function set_Staggered(params)

        if haskey(params,"mass")
            mass = params["mass"]
        else
            mass = 0.5
            println("Warning. mass = $mass, Default value is used in measurement $(params["paramsname"])")
        end

        if haskey(params,"Nf")
            Nf = params["Nf"]
        else                        
            error("Nf should be set if you want to use the staggered fermion in measurements")
            #println("Warning. mass = $hop, Default value is used")
        end

        quench = false

        smearing_for_fermion = set_params(params,"smearing_for_fermion","nothing")
        stout_numlayers = set_params(params,"stout_numlayers",nothing)
        stout_ρ = set_params(params,"stout_ρ",nothing)
        stout_loops = set_params(params,"stout_loops",nothing)

        if smearing_for_fermion == "nothing"
            fparam = FermiActionParam_Staggered(mass,eps,fermiontype,MaxCGstep,quench,Nf)
        else
            L = (univ.NX,univ.NY,univ.NZ,univ.NT)
            fparam = FermiActionParam_Staggered(mass,eps,fermiontype,MaxCGstep,quench,Nf,
                        smearingparameters = "stout",
                        loops_list = stout_loops,
                        coefficients  = stout_ρ,
                        numlayers = stout_numlayers,
                        L = L)
        end

        
        #println("Measurement_set::mass_measurement = $(p.mass_measurement)")
        #fparam = FermiActionParam_Staggered(mass,eps,fermiontype,MaxCGstep,quench,Nf)
        return fparam
    end

end