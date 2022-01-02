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
    #import ..Smearing:gradientflow!,calc_stout!,calc_fatlink_APE!,calc_stout,calc_fatlink_APE,calc_multihit!
    #import ..Actions:GaugeActionParam_autogenerator,GaugeActionParam
    import ..Fermionfield_LQCD:FermiActionParam_Wilson,FermiActionParam_Staggered,FermiActionParam_WilsonClover,
                FermiActionParam,clear_fermion!,Z4_distribution_fermion!,bicg,get_origin
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
        measure_correlator(U,m,itrj,verbose)
        #error("not implemented")
        return 
    end

    function measure_correlator(Uin::Array{<: AbstractGaugefields{NC,Dim},1},m::M,itrj,verbose) where {M <: Measure_Pion_correlator,NC,Dim}
        C = calc_pion_correlator(Uin,m,verbose)
        print_verbose1(verbose,"$itrj ")
        print(m.fp,"$itrj ")
        for it=1:length(C)
            cc = C[it]
            print_verbose1(verbose,"$cc ")
            print(m.fp,"$cc ")
        end
        println_verbose1(verbose,"#pioncorrelator")
        println(m.fp,"#pioncorrelator")
    end

    function calc_pion_correlator(Uin::Array{<: AbstractGaugefields{NC,Dim},1},m::Me,verbose) where {Me <: Measure_Pion_correlator,NC,Dim}
        println_verbose2(verbose,"Hadron spectrum started")
        #println("Hadron spectrum started")
        #U = univ._temporal_gauge
        #if univ.NC != 3
        #    error("not implemented yet for calc_pion_correlator")
        #end
        #b = meas._temporal_fermi2[1] # source allocate
        #k = meas._temporal_fermi2[2] # sink allocate

        # setup a massive Dirac operator
        U,_... = calc_smearedU(Uin,m.fparam.smearing)
        M = Dirac_operator(U,m._temporal_fermions[1],m.fparam)
        #M = Dirac_operator(univ.U,meas._temporal_fermi2[1],meas.fparam)

        # Allocate the Wilson matrix = S
        Nspinor = ifelse( m.fparam.Dirac_operator == "Staggered" ,1,4)
        _,_,NN... = size(Uin[1])
        S = zeros(ComplexF64, (NN..., Nspinor*NC, Nspinor*NC) )

        # calculate quark propagators from a point source at he origin
        propagators = calc_quark_propagators_point_source(M,m,NC,verbose)



        #ctr = 0 # a counter
        for ic=1:NC
            for is=1:Nspinor
                icum = (ic-1)*Nspinor+is
                #=
                clear!(b)
                b[ic,1,1,1,1,is]=1 # source at the origin
                @time cg0!(k,b,1, univ.U, meas._temporal_gauge, meas._temporal_fermi1, meas.fparam) # k[x] = M^{-1}b[0]
                =#
                propagator = propagators[icum]
                α0=spincolor(ic,is,NC) # source(color-spinor) index
                # reconstruction
                if Dim==4
                    for t=1:NN[4]
                        for z=1:NN[3]
                            for y=1:NN[2]
                                for x=1:NN[1]
                                    for ic2=1:NC 
                                        for is2=1:Nspinor # Nspinor is the number of spinor index in 4d.
                                            β=spincolor(ic2,is2,NC)
                                            S[x,y,z,t,α0,β]+= propagator[ic,x,y,z,t,is]
                                        end
                                    end
                                end
                            end
                        end
                    end
                else
                    error("Dim = $Dim is not supported")
                end
                # end for the substitution
                
                
                #ctr+=1
            end 
        end
        # contruction end.

        println_verbose2(verbose,"Hadron spectrum: Reconstruction")
        #println("Hadron spectrum: Reconstruction")
        Cpi = zeros( NN[end] )
        #Cpi = zeros( univ.NT )
        # Construct Pion propagator 
        if Dim==4
            for t=1:NN[4]
                tmp = 0.0+0.0im
                for z=1:NN[3]
                    for y=1:NN[2]
                        for x=1:NN[1]
                            for ic=1:NC 
                                for is=1:Nspinor # Nspinor is the number of spinor index in 4d.
                                    α=spincolor(ic,is,NC)
                                    for ic2=1:NC 
                                        for is2=1:Nspinor # Nspinor is the number of spinor index in 4d.
                                            β=spincolor(ic2,is2,NC)
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
                # staggered Pion correlator relies on https://itp.uni-frankfurt.de/~philipsen/theses/breitenfelder_ba.pdf (3.33)
                # we adopt ignoreing the staggering factor. See detail above reference.
                ksfact = 1.0 #ifelse( meas.fparam.Dirac_operator == "Staggered" , (-1)^(t-1) * 64, 1)
                Cpi[t] = real(tmp)*ksfact
            end
        end
        
        #println(typeof(verbose),"\t",verbose)
        println_verbose2(verbose,"Hadron spectrum end")
        #println("Hadron spectrum end")
        return Cpi
    end

    function spincolor(ic,is,NC)
        return ic-1 + (is-1)*NC + 1
    end

    function calc_quark_propagators_point_source(D,meas,NC,verbose)
        # D^{-1} for each spin x color element
        Nspinor = ifelse( meas.fparam.Dirac_operator == "Staggered" ,1,4)
        propagators = map(i -> calc_quark_propagators_point_source_each(D,meas,i,NC,verbose),1:NC*Nspinor)
        return propagators
    end



    function calc_quark_propagators_point_source_each(M,meas,i,NC,verbose)
        # calculate D^{-1} for a given source at the origin.
        # Nc*Ns (Ns: dim of spinor, Wilson=4, ks=1) elements has to be gathered.
        # staggered Pion correlator relies on https://itp.uni-frankfurt.de/~philipsen/theses/breitenfelder_ba.pdf (3.33)
        b = similar(meas._temporal_fermions[1]) # source is allocated
        p = similar(b) # sink is allocated (propagator to sink position)
        #k = meas._temporal_fermi2[2]
        clear_fermion!(b)
        Nspinor = ifelse( meas.fparam.Dirac_operator == "Staggered" ,1,4)
        is = ((i-1) % Nspinor)+1 # spin index   
        ic = ((i-is) ÷ Nspinor)+ 1 # color index
        println_verbose1(verbose,"$ic $is")
        #println("$ic $is")
        iorigin = get_origin(b)
        b[ic,iorigin...,is]=1 # source at the origin

        @time bicg(p,M,b,eps=meas.fparam.eps,maxsteps = meas.fparam.MaxCGstep,verbose = verbose) # solve Mp=b, we get p=M^{-1}b

        #@time cg0!(k,b,1, univ.U, meas._temporal_gauge, meas._temporal_fermi1, meas.fparam) # k[x] = M^{-1}b[0]
        println_verbose1(verbose,"Hadron spectrum: Inversion $(i)/$(NC*Nspinor) is done")
        #println("Hadron spectrum: Inversion $(i)/$(NC*4) is done")
        flush(stdout)
        return p
    end

end