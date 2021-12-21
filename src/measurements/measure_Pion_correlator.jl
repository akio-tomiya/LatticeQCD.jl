module Measure_Pion_correlator_module
    using LinearAlgebra
    #using InteractiveUtils
    import ..AbstractMeasurement_module:AbstractMeasurement,measure,set_params,
            set_Fermiontype
    import ..Gaugefield:AbstractGaugefields,
                                        Traceless_antihermitian!,Traceless_antihermitian,
                                        initialize_TA_Gaugefields
    import ..Gaugefield:calculate_Plaquette,calculate_Polyakov_loop,substitute_U!

    import ..Gaugefield:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3,println_verbose2,println_verbose1,
            print_verbose1,print_verbose2,print_verbose3
    import ..Gaugefield:Wilson_loop,Wilson_loop_set,calc_loopset_μν_name            
    import ..Gaugefield:Loops,evaluate_loops,TA_Gaugefields,calc_smearedU
    import ..Fermionfield_LQCD:Dirac_operator
    import ..Smearing:gradientflow!,calc_stout!,calc_fatlink_APE!,calc_stout,calc_fatlink_APE,calc_multihit!
    #import ..Actions:GaugeActionParam_autogenerator,GaugeActionParam
    import ..Fermionfield_LQCD:FermiActionParam_Wilson,FermiActionParam_Staggered,FermiActionParam_WilsonClover,
                FermiActionParam,clear_fermion!,Z4_distribution_fermion!,bicg
    #import ..Actions:FermiActionParam_Wilson,FermiActionParam_Staggered,FermiActionParam_WilsonClover,
    #            FermiActionParam




    mutable struct Measure_Pion_correlator{T,FP,Ftype} <: AbstractMeasurement
        filename::String
        fp::IOStream
        tempU::Array{T,1}
        printvalues::Bool
        fparam::FP
        _temporal_fermions::Array{Ftype,1}        

        function Measure_Pion_correlator(filename,
                    U::Array{T,1},params;printvalues = true) where T
            fp = open(filename,"w")
            
            tempU = Array{T,1}(undef,3)
            for i=1:3
                tempU[i] = similar(U[1])
            end

            if haskey(params,"BoundaryCondition")
            else
                params["BoundaryCondition"] = [1,1,1,-1]
            end


            fparam,_temporal_fermions = set_Fermiontype(U,params)
            FP = typeof(fparam)
            Ftype = eltype(_temporal_fermions)

            m = new{T,FP,Ftype}(filename,fp,tempU,printvalues,fparam,_temporal_fermions)

            finalizer(m) do m
                close(m.fp)
            end
            return m
        end

    end


    function measure(m::M,itrj,U::Array{<: AbstractGaugefields{NC,Dim},1};verbose = Verbose_2()) where {M <: Measure_Pion_correlator,NC,Dim}
        error("not implemented")
        return 
    end



end