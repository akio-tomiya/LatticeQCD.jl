module Measure_chiral_condensate_module
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




    mutable struct Measure_chiral_condensate{T,FP,Ftype} <: AbstractMeasurement
        filename::String
        fp::IOStream
        tempU::Array{T,1}
        printvalues::Bool
        fparam::FP
        _temporal_fermions::Array{Ftype,1}        

        function Measure_chiral_condensate(filename,
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


    function measure(m::M,itrj,U::Array{<: AbstractGaugefields{NC,Dim},1};verbose = Verbose_2()) where {M <: Measure_chiral_condensate,NC,Dim}
        Nr = 10
        pbp = calc_chiral_cond(m,itrj,U,Nr,verbose)
        if m.printvalues
            println_verbose1(verbose,"$itrj $pbp # pbp Nr=$Nr")
            println(m.fp,"$itrj $pbp # pbp Nr=$Nr")
            flush(stdout)
        end
        #error("not implemented")
        return 
    end

    function calc_chiral_cond(m::Me,itrj,Uin::Array{<: AbstractGaugefields{NC,Dim},1}, Nr = 10, verbose = Verbose_2()) where {Me <: Measure_chiral_condensate,NC,Dim}
        #(univ::Universe,meas,measfp,itrj, Nr = 10, verbose = Verbose_2())
        # pbp = (1/Nr) Σ_i p_i
        # p_i = r_i^\dag xi_i
        # xi_i = D^{-1} r_i   # D xi = r : r is a random veccor
        #
        Nfbase = ifelse( m.fparam.Dirac_operator == "Staggered",4,1)
        Nf = m.fparam.Nf
        factor = Nf/Nfbase
        #
        println_verbose2(verbose,"Chiral condensate for Nf = $(Nf), Dirac_operator = $(m.fparam.Dirac_operator), factor=$factor is multiplied.")
        #
        pbp = 0.0
        # setup a massive Dirac operator
        println(m.fparam.smearing)
        U,_... = calc_smearedU(Uin,m.fparam.smearing)
        #U,_... = calc_smearingU(U,m.fparam.smearing)
        M = Dirac_operator(U,m._temporal_fermions[1],m.fparam)
        #M = Dirac_operator(univ.U,meas._temporal_fermi2[1],meas.fparam)
        r = similar(m._temporal_fermions[1]) 
        p = similar(r) 
        for ir=1:Nr
            clear_fermion!(p)
            #gauss_distribution_fermi!(r,univ.ranf)
            Z4_distribution_fermion!(r)
            #set_wing_fermi!(r) 
            @time bicg(p,M,r,eps=m.fparam.eps,maxsteps = m.fparam.MaxCGstep,verbose = verbose) # solve Mp=b, we get p=M^{-1}b
            tmp = dot(r,p) # hermitian inner product
            if m.printvalues
                println_verbose2(verbose,"# $itrj $ir $(real(tmp)/U[1].NV) # itrj irand chiralcond")
                println(m.fp,"# $itrj $ir $(real(tmp)/U[1].NV) # itrj irand chiralcond")
            end
            pbp+=tmp
        end
        return real(pbp/Nr)/U[1].NV * factor

    end

end