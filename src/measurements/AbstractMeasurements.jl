module AbstractMeasurement_module
    import ..Gaugefield:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3,println_verbose2,println_verbose1,
                print_verbose1,print_verbose2,print_verbose3
    import ..Gaugefield:AbstractGaugefields
    import ..Actions:GaugeActionParam_autogenerator,GaugeActionParam
    #import ..Actions:FermiActionParam_Wilson,FermiActionParam_Staggered,FermiActionParam_WilsonClover,
    #            FermiActionParam
    import ..Fermionfield_LQCD:Fermionfields,FermiActionParam_Wilson,FermiActionParam_Staggered,FermiActionParam_WilsonClover,
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

    function set_Fermiontype(U,params)
        if haskey(params,"fermiontype")
            fermiontype = params["fermiontype"] 
        else
            error("fermiontype should be defined")
        end

        if fermiontype  == "Wilson"
            fparam,fermion = set_Wilson(U,params)
        elseif fermiontype == "WilsonClover"
            fparam,fermion = set_WilsonClover(U,params)
        elseif fermiontype == "Staggered"
            fparam,fermion = set_Staggered(U,params)
        elseif fermiontype == nothing
            fparam = nothing
            fermion = nothing
        else
            error("$fermiontype is not supported. use Wilson or Staggered")
        end

        return fparam,fermion

    end


    function set_Wilson(U,params)
        fermiontype = "Wilson"
        if haskey(params,"MaxCGstep")
            MaxCGstep = params["MaxCGstep"]
        else
            MaxCGstep = 5000
        end
        if haskey(params,"eps")
            eps = params["eps"]
        else
            eps = 1e-16
        end



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

            L = (U[1].NX,U[1].NY,U[1].NZ,U[1].NT)
            fparam = FermiActionParam_Wilson(hop,r,eps,fermiontype,MaxCGstep,quench,
                                                smearingparameters = "stout",
                                                loops_list = stout_loops,
                                                coefficients  = stout_ρ,
                                                numlayers = stout_numlayers,
                                                L = L)
        end



        return fparam,fermion
    end

    function set_WilsonClover(U,params)
        fermiontype = "WilsonClover"
        if haskey(params,"MaxCGstep")
            MaxCGstep = params["MaxCGstep"]
        else
            MaxCGstep = 5000
        end
        if haskey(params,"eps")
            eps = params["eps"]
        else
            eps = 1e-16
        end



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
            _ftmp_vectors[i] = zeros(ComplexF64,U[1].NC,NV,4)
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
            L = (U[1].NX,U[1].NY,U[1].NZ,U[1].NT)
            FermiActionParam_WilsonClover(hop,r,eps,fermiontype,MaxCGstep,Clover_coefficient,
                        internal_flags,inn_table,_ftmp_vectors,_is1,_is2,
                        quench,
                        smearingparameters = "stout",
                        loops_list = stout_loops,
                        coefficients  = stout_ρ,
                        numlayers = stout_numlayers,
                        L = L)
        end
        return fparam,fermion
    end

    function set_Staggered(U,params)
        fermiontype = "Staggered"
        if haskey(params,"MaxCGstep")
            MaxCGstep = params["MaxCGstep"]
        else
            error("MaxCGstep is not set in measurement $(params["paramsname"])")
        end
        if haskey(params,"eps")
            eps = params["eps"]
        else
            error("eps is not set in measurement $(params["paramsname"])")
        end



        if haskey(params,"mass")
            mass = params["mass"]
        else
            error("mass is not set in measurement $(params["paramsname"])")
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
            L = (U[1].NX,U[1].NY,U[1].NZ,U[1].NT)
            fparam = FermiActionParam_Staggered(mass,eps,fermiontype,MaxCGstep,quench,Nf,
                        smearingparameters = "stout",
                        loops_list = stout_loops,
                        coefficients  = stout_ρ,
                        numlayers = stout_numlayers,
                        L = L)
        end


        NC,_,NN... = size(U[1]) 
        ϕ = Fermionfields(params,NC,U[1].NDW,fermiontype,NN...)

        _temporal_fermions = Array{typeof(ϕ),1}(undef,4)
        for i=1:length(_temporal_fermions)
            _temporal_fermions[i] = similar(ϕ)
        end

        
        #println("Measurement_set::mass_measurement = $(p.mass_measurement)")
        #fparam = FermiActionParam_Staggered(mass,eps,fermiontype,MaxCGstep,quench,Nf)
        return fparam,_temporal_fermions
    end

end