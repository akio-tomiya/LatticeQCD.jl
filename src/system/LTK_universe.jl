module LTK_universe
    using LinearAlgebra
    using SparseArrays
    #export Universe
    

    import ..Gaugefields:GaugeFields,GaugeFields_1d,
                        IdentityGauges,RandomGauges,
                        set_wing!,
                        substitute!,
                        calc_GaugeAction,
                        SU3GaugeFields,SU3GaugeFields_1d,
                        SU2GaugeFields,SU2GaugeFields_1d,
                        SUNGaugeFields,SUNGaugeFields_1d,
                        Oneinstanton,
                        evaluate_wilson_loops!,
                        U1GaugeFields,U1GaugeFields_1d,
                        apply_smearing,calc_smearingU
    import ..Gaugefields
                        
    import ..Fermionfields:FermionFields,WilsonFermion,StaggeredFermion,substitute_fermion!,gauss_distribution_fermi!,set_wing_fermi!
    import ..Fermionfields
    import ..Actions:GaugeActionParam,FermiActionParam,
                Setup_Gauge_action,Setup_Fermi_action,
                GaugeActionParam_standard,FermiActionParam_Wilson,
                show_parameters_action,
                FermiActionParam_WilsonClover,
                FermiActionParam_Staggered,
                GaugeActionParam_autogenerator,SmearingParam_multi,SmearingParam_single,Nosmearing,
                FermiActionParam_Domainwall
    import ..LieAlgebrafields:LieAlgebraFields,clear!,add_gaugeforce!
    import ..LieAlgebrafields
    import ..Rand:Random_LCGs
    import ..System_parameters:Params
    import ..Diracoperators:DdagD_operator,Domainwall_operator,make_densematrix,
                D5DW_Domainwall_operator
    import ..Diracoperators
    import ..Wilsonloops:make_loopforactions,Wilson_loop_set,make_originalactions_fromloops,
                make_cloverloops
    import ..Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1
    import ..IOmodule:loadU
    import ..SUN_generator:Generator

    import ..RationalApprox:calc_exactvalue,calc_Anϕ
    import ..ILDG_format:ILDG,load_gaugefield,load_gaugefield!,save_binarydata
    import ..Othermethods:tdlogdet
    import ..CGmethods:bicg,cg,shiftedcg
    import ..Rhmc:get_order,get_β,get_α,get_α0,get_β_inverse,get_α_inverse,get_α0_inverse
    #import ..Smearing:apply_smearing



    """
    Your universe is described in this type.
    """
    struct Universe{Gauge,Lie,Fermi,GaugeP,FermiP,Gauge_temp}
        NX::Int64
        NY::Int64
        NZ::Int64
        NT::Int64
        NV::Int64
        NC::Int64
        Nwing::Int64
        Dirac_operator::Union{Nothing,String}
        U::Array{Gauge,1}
        Uold::Array{Gauge,1}
        p::Array{Lie,1}
        φ::Union{Nothing,Fermi}
        ξ::Union{Nothing,Fermi}
        η::Union{Nothing,Fermi}
        gparam::GaugeP
        fparam::Union{Nothing,FermiP}
        BoundaryCondition::Array{Int8,1}
        initial::String
        quench::Bool
        NDFALG::Int64

        
        
        

        _temporal_gauge::Array{Gauge_temp,1}
        _temporal_fermi::Union{Nothing,Array{Fermi,1}}
        _temporal_algebra::Array{Lie,1}
        ranf::Random_LCGs
        verboselevel::Int8
        kind_of_verboselevel::Verbose_level

    end





    function set_β!(univ::Universe,β)
        if typeof(univ.gparam) == GaugeActionParam_standard 
            univ.gparam.β = β
            #univ.gparam = GaugeActionParam_standard(β,univ.gparam.NTRACE)
        elseif typeof(univ.gparam) == GaugeActionParam_autogenerator
            #univ.gparam =GaugeActionParam_autogenerator(univ.gparam.βs * (univ.gparam.β/β),β,univ.gparam.numactions,univ.gparam.NTRACE,univ.gparam.loops,univ.gparam.staples)
            univ.gparam.β = β
            univ.gparam.βs .= 0#univ.gparam.βs * (univ.gparam.β/β)
            univ.gparam.βs[1] = β
        else
            error("$(typeof(univ.gparam)) is not supported!")
        end
    end

    function set_β!(univ::Universe,βs::Array{T,1}) where T <: Number
        if typeof(univ.gparam) == GaugeActionParam_standard 
            univ.gparam.β = βs[1]
            #univ.gparam = GaugeActionParam_standard(β,univ.gparam.NTRACE)
        elseif typeof(univ.gparam) == GaugeActionParam_autogenerator
            @assert length(univ.gparam.βs) == length(βs) "univ.gparam.βs = $(univ.gparam.βs) and βs = $βs"
            #univ.gparam =GaugeActionParam_autogenerator(univ.gparam.βs * (univ.gparam.β/β),β,univ.gparam.numactions,univ.gparam.NTRACE,univ.gparam.loops,univ.gparam.staples)
            univ.gparam.β = βs[1]
            univ.gparam.βs .= βs
            #univ.gparam.βs = univ.gparam.βs * (univ.gparam.β/β)
        else
            error("$(typeof(univ.gparam)) is not supported!")
        end
    end

    function get_β(univ::Universe)
        if typeof(univ.gparam) == GaugeActionParam_standard 
            return univ.gparam.β
        elseif typeof(univ.gparam) == GaugeActionParam_autogenerator
            return univ.gparam.β
        else
            error("$(typeof(univ.gparam)) is not supported!")
        end
    end

    function set_βs!(univ::Universe,βs)
        if typeof(univ.gparam) == GaugeActionParam_standard 
            univ.gparam.β = βs[1]
        elseif typeof(univ.gparam) == GaugeActionParam_autogenerator
            univ.gparam.β = βs[1]
            univ.gparam.βs = βs[:]
        else
            error("$(typeof(univ.gparam)) is not supported!")
        end
    end


    include("default.jl")

    function Universe()
        file = "default.jl"
        return Universe(file)
    end

    function Universe(p::Params)
        L = p.L
        β = p.β
        NTRACE = p.NC
        if p.use_autogeneratedstaples 
            #loops = make_loopforactions(p.couplinglist,L)

            #println(loops)
            loops = make_originalactions_fromloops(p.coupling_loops)
            
            gparam = GaugeActionParam_autogenerator(p.couplingcoeff,loops,p.NC,p.couplinglist)
        else
            gparam =  GaugeActionParam_standard(β,NTRACE)
        end

        if p.Dirac_operator == nothing
            fparam = nothing
        else
            if p.Dirac_operator == "Wilson"
                if p.smearing_for_fermion == "nothing"
                    fparam = FermiActionParam_Wilson(p.hop,p.r,p.eps,p.Dirac_operator,p.MaxCGstep,p.quench)
                else
                    fparam = FermiActionParam_Wilson(p.hop,p.r,p.eps,p.Dirac_operator,p.MaxCGstep,p.quench,
                                                        smearingparameters = "stout",
                                                        loops_list = p.stout_loops,
                                                        coefficients  = p.stout_ρ,
                                                        numlayers = p.stout_numlayers,
                                                        L = p.L)
                end
            elseif p.Dirac_operator == "WilsonClover"
                #if p.NC == 2
                #    error("You use NC = 2. But WilsonClover Fermion is not supported in SU2 gauge theory yet. Use NC = 3.")
                #end
                NV = prod(p.L)
                #CloverFμν = zeros(ComplexF64,p.NC,p.NC,NV,6)
                inn_table= zeros(Int64,NV,4,2)
                internal_flags = zeros(Bool,2)
                _ftmp_vectors = Array{Array{ComplexF64,3},1}(undef,6)
                for i=1:6
                    _ftmp_vectors[i] = zeros(ComplexF64,p.NC,NV,4)
                end

                _is1 = zeros(Int64,NV)
                _is2 = zeros(Int64,NV)

                if p.NC ≥ 4
                    SUNgenerator = Generator(p.NC)
                else
                    SUNgenerator = nothing
                end

                _cloverloops = Array{Wilson_loop_set,2}(undef,3,4)
                for μ=1:3
                    for ν=μ+1:4
                        _cloverloops[μ,ν] = make_cloverloops(μ,ν)
                    end
                end


                #fparam = FermiActionParam_WilsonClover(p.hop,p.r,p.eps,p.Dirac_operator,p.MaxCGstep,p.Clover_coefficient,CloverFμν,
                #                internal_flags,inn_table,_ftmp_vectors,_is1,_is2,
                #                p.quench)
                if p.smearing_for_fermion == "nothing"
                    fparam = FermiActionParam_WilsonClover(p.hop,p.r,p.eps,p.Dirac_operator,p.MaxCGstep,p.Clover_coefficient,
                                    internal_flags,inn_table,_ftmp_vectors,_is1,_is2,
                                    p.quench,SUNgenerator,_cloverloops)
                else
                    error("stout for Wilson with clover is not supported yet!")
                    fparam = FermiActionParam_WilsonClover(p.hop,p.r,p.eps,p.Dirac_operator,p.MaxCGstep,p.Clover_coefficient,
                                    internal_flags,inn_table,_ftmp_vectors,_is1,_is2,
                                    p.quench,SUNgenerator,_cloverloops,
                                    smearingparameters = "stout",
                                    loops_list = p.stout_loops,
                                    coefficients  = p.stout_ρ,
                                    numlayers = p.stout_numlayers,
                                    L = p.L)
                end
            elseif p.Dirac_operator == "Staggered"
                if p.smearing_for_fermion == "nothing"
                    fparam = FermiActionParam_Staggered(p.mass,p.eps,p.Dirac_operator,p.MaxCGstep,p.quench,p.Nf)
                else
                    fparam = FermiActionParam_Staggered(p.mass,p.eps,p.Dirac_operator,p.MaxCGstep,p.quench,p.Nf,
                                                            smearingparameters = "stout",
                                                            loops_list = p.stout_loops,
                                                            coefficients  = p.stout_ρ,
                                                            numlayers = p.stout_numlayers,
                                                            L = p.L
                                                        )
                end
            elseif p.Dirac_operator == "Domainwall"
                if p.smearing_for_fermion == "nothing"
                    fparam = FermiActionParam_Domainwall(p.Domainwall_N5,p.Domainwall_r,p.Domainwall_M,p.Domainwall_m,p.Domainwall_ωs,
                    p.Domainwall_b,p.Domainwall_c,
                    p.eps,p.Dirac_operator,
                    p.MaxCGstep,p.quench)
                else
                    fparam = FermiActionParam_Domainwall(p.Domainwall_N5,p.Domainwall_r,p.Domainwall_M,p.Domainwall_m,p.Domainwall_ωs,
                                                            p.Domainwall_b,p.Domainwall_c,
                                                            p.eps,p.Dirac_operator,
                                                            p.MaxCGstep,p.quench,
                                                            smearingparameters = "stout",
                                                            loops_list = p.stout_loops,
                                                            coefficients  = p.stout_ρ,
                                                            numlayers = p.stout_numlayers,
                                                            L = p.L
                                                        )
                end
                #error("Not implemented! here is in LTK_universe.jl")
            else
                error(p.Dirac_operator," is not supported!")
            end
        end

        univ = Universe(L,gparam,p.Nwing,fparam,p.BoundaryCondition,p.initial,p.NC,p.verboselevel,p.load_fp,p.loadU_format)

    end




    """
    ```Universe(file)```
    - file: file name of the input file.

    Make your universe. The input file is loaded.

    Undefined parameters in your input file are defined with the default values.

    The default values are as follows.

    ```julia
        L = (4,4,4,4)
        β = 6
        NTRACE = 3
        #gparam = Setup_Gauge_action(β)
        gparam =  GaugeActionParam_standard(β,NTRACE)

        BoundaryCondition=[1,1,1,-1]
        Nwing = 1
        initial="cold"
        NC =3


        hop= 0.141139 #Hopping parameter
        r= 1 #Wilson term
        eps= 1e-19
        Dirac_operator= "Wilson"
        MaxCGstep= 3000

        fparam = FermiActionParam_Wilson(hop,r,eps,Dirac_operator,MaxCGstep)

    ```
    """
    function Universe(file)
        if file != "default.jl"
            include(pwd()*"/"*file)
        end
        univ = Universe(L,gparam,Nwing,fparam,BoundaryCondition,initial,NC)

    end

    function show_parameters(univ::Universe)
        println("""
            L = $((univ.NX,univ.NY,univ.NZ,univ.NT))
            """)

        show_parameters_action(univ.gparam)

        println("""
            BoundaryCondition= $(univ.BoundaryCondition)
            Nwing = $(univ.Nwing)
            initial= $(univ.initial)
            NC =$(univ.NC)
            """)

        if univ.fparam != nothing
            show_parameters_action(univ.fparam)
        end
    end

    """
    ```Universe(L::Tuple,gparam::GaugeActionParam,initial="cold",fparam=nothing)```

    - L: system size (NX,NY,NZ,NT)
    - gparam: parameters for gauge actions
    - [initial]: initial Gauge configuration
    - [fparam]: parameters for fermion actions
    """
    function Universe(L::Tuple,gparam::GaugeActionParam,initial="cold",fparam=nothing)
        Nwing = 1
        BoundaryCondition=[1,1,1,-1]
        NC =3
        Universe(L,gparam,Nwing,fparam,BoundaryCondition,initial,NC)
    end


    function Universe(L,gparam,Nwing,fparam,BoundaryCondition,initial,NC,verboselevel,load_fp,loadU_format)
        #(L::Tuple,gparam::GaugeActionParam;
        #   Nwing = 1,fparam=nothing,
        #  BoundaryCondition=[1,1,1,-1],initial="cold",NC =3)

        if verboselevel == 1
            kind_of_verboselevel = Verbose_1(load_fp)
        elseif verboselevel == 2
            kind_of_verboselevel = Verbose_2(load_fp)
        elseif verboselevel == 3
            kind_of_verboselevel = Verbose_3(load_fp)
        end


        ranf = Random_LCGs(1)

        NX = L[1]
        NY = L[2]
        NZ = L[3]
        NT = L[4]
        NV = NX*NY*NZ*NT
        NDFALG  = ifelse(NC==1,1,NC^2-1)

        num_tempfield_g = 5
        num_tempfield_f = 4

        if fparam != nothing && fparam.Dirac_operator == "WilsonClover"
            numbasis = NC^2-1
            num_tempfield_g += 4+4+2numbasis
        elseif fparam != nothing # && fparam.Dirac_operator == "Staggered"
            num_tempfield_f += 4
        end

        if fparam != nothing && fparam.Dirac_operator == "Staggered"
            if fparam.Nf != 4 && fparam.Nf != 8
                N_action = get_order(fparam.rhmc_action)
                N_MD = get_order(fparam.rhmc_MD)
                num_tempfield_f += maximum((N_action,N_MD))+1
            end
        end

        if fparam != nothing && fparam.smearing != nothing
            num_tempfield_g += 4
        end

        if typeof(gparam) == GaugeActionParam_autogenerator
            num_tempfield_g += 1
        end
        

        if NC == 3
            U = Array{SU3GaugeFields,1}(undef,4)
            _temporal_gauge = Array{SU3GaugeFields_1d,1}(undef,num_tempfield_g)
        elseif NC == 2
            U = Array{SU2GaugeFields,1}(undef,4)
            _temporal_gauge = Array{SU2GaugeFields_1d,1}(undef,num_tempfield_g)
        elseif NC ≥ 4
            U = Array{SUNGaugeFields{NC},1}(undef,4)
            _temporal_gauge = Array{SUNGaugeFields_1d{NC},1}(undef,num_tempfield_g)
        elseif NC == 1
            U = Array{U1GaugeFields,1}(undef,4)
            _temporal_gauge = Array{U1GaugeFields_1d,1}(undef,num_tempfield_g)
        end
        
#        Uold = Array{GaugeFields,1}(undef,4)

        p = Array{LieAlgebraFields,1}(undef,4)
        for μ=1:4
            p[μ] = LieAlgebraFields(NC,NX,NY,NZ,NT)
        end

        if initial == "cold"
            println(".....  Cold start")
            for μ=1:4
                U[μ] = IdentityGauges(NC,NX,NY,NZ,NT,Nwing)
            end
        elseif initial == "hot"
            println(".....  Hot start")
            for μ=1:4
                U[μ] = RandomGauges(NC,NX,NY,NZ,NT,Nwing)
            end
        elseif initial == "one instanton"
            @assert NC == 2 "From one instanton start, NC should be 2!"
            U = Oneinstanton(NC,NX,NY,NZ,NT,Nwing)

        else #if initial == "file"
            println(".....  File start")
            println("File name is $initial")
            if loadU_format == "JLD"
                U = loadU(initial,NX,NY,NZ,NT,NC)
            elseif loadU_format == "ILDG"
                ildg = ILDG(initial)
                i = 1
                for μ=1:4
                    U[μ] = IdentityGauges(NC,NX,NY,NZ,NT,Nwing)
                end
                load_gaugefield!(U,i,ildg,L,NC)
            else
                error("loadU_format should be JLD or ILDG")
            end

            #error("not supported yet.")
        end

        set_wing!(U)
        Uold = similar(U)
        substitute!(Uold,U)

        
        #=
        for μ=1:4
            set_wing!(U[μ])
            Uold[μ] = similar(U[μ])
            substitute!(Uold[μ],U[μ])
        end
        =#

        
        #_temporal_gauge = Array{GaugeFields,1}(undef,4)
        for i=1:length(_temporal_gauge)
            _temporal_gauge[i] = GaugeFields_1d(NC,NX,NY,NZ,NT) #similar(U[1])
            #_temporal_gauge[i] = similar(U[1])
        end

        _temporal_algebra = Array{LieAlgebraFields,1}(undef,1)
        for i=1:length(_temporal_algebra)
            _temporal_algebra[i] = similar(p[1])
        end


        if fparam == nothing
            Dirac_operator = nothing
            φ = nothing
            η = nothing
            ξ = nothing
            _temporal_fermi = nothing
            quench = true
        else
            Dirac_operator = fparam.Dirac_operator
            φ = FermionFields(NC,NX,NY,NZ,NT,fparam,BoundaryCondition)
            η = similar(φ)
            ξ = similar(φ)
            _temporal_fermi = Array{FermionFields,1}(undef,num_tempfield_f )
            for i=1:length(_temporal_fermi)
                _temporal_fermi[i] = similar(φ)
            end
            quench = fparam.quench


            if Dirac_operator =="Domainwall" 
                gauss_distribution_fermi!(ξ)

                #A = D5DW_Domainwall_operator(U,ξ,fparam)
                A = Domainwall_operator(U,ξ,fparam)
                fn = "_LS8_3333_hot.txt"
                #
                
                A_wilson = make_densematrix(A.D5DW.wilsonoperator)
                e,v = eigen(A_wilson)
                fp = open("eigenvalues_wilson"*fn,"w")
                for ene in e
                    println(fp,real(ene),"\t",imag(ene))
                end
                close(fp)

                A_dense = make_densematrix(A.D5DW)
                AdA = A_dense'*A_dense

                e,v = eigen(AdA)
                fp = open("eigenvalues_AdA"*fn,"w")
                for ene in e
                    println(fp,real(ene),"\t",imag(ene))
                end
                close(fp)
                #error("dd")

                 
                
                e,v = eigen(A_dense)
                fp = open("eigenvalues_domain"*fn,"w")
                for ene in e
                    println(fp,real(ene),"\t",imag(ene))
                end
                close(fp)
                #error("dd")

                 
                A_dense_PV = make_densematrix(A.D5DW_PV)
                M = A_dense*inv(A_dense_PV)
                #e,v = eigen(A_dense)
                e,v = eigen(M)
                fp = open("eigenvalues_M"*fn,"w")
                for ene in e
                    println(fp,real(ene),"\t",imag(ene))
                end
                close(fp)
                #println(e)
                error("stop!")
                bicg(η,A,ξ)
                error("stop!")
            end

        end

        

        Gauge =eltype(U)
        Fermi = typeof(φ)
        Lie = eltype(p)
        GaugeP = typeof(gparam)
        FermiP = typeof(fparam)
        Gauge_temp = eltype(_temporal_gauge)



        return Universe{Gauge,Lie,Fermi,GaugeP,FermiP,Gauge_temp}(
            NX,
            NY,
            NZ,
            NT,
            NV,
            NC,
            Nwing,
            Dirac_operator,
            U,
            Uold,
            p,
            φ,
            ξ,
            η,
            gparam,
            fparam,
            BoundaryCondition,
            initial,
            quench,
            NDFALG,
            _temporal_gauge,
            _temporal_fermi,
            _temporal_algebra,
            ranf,
            verboselevel,
            kind_of_verboselevel
        )


    end

    function gauss_distribution(univ::Universe,nv) 
        variance = 1
        nvh = div(nv,2)
        granf = zeros(Float64,nv)
        for i=1:nvh
            rho = sqrt(-2*log(univ.ranf())*variance)
            theta = 2pi*univ.ranf()
            granf[i] = rho*cos(theta)
            granf[i+nvh] = rho*sin(theta)
        end
        if 2*nvh == nv
            return granf
        end

        granf[nv] = sqrt(-2*log(univ.ranf())*variance) * cos(2pi*univ.ranf())
        return granf
    end






    function expF_U!(U::Array{T,1},F::Array{N,1},Δτ,univ::Universe) where {T<: GaugeFields, N <: LieAlgebraFields} 
        LieAlgebrafields.expF_U!(U,F,Δτ,univ._temporal_gauge,univ._temporal_algebra[1])
    end

    function calc_gaugeforce!(F::Array{N,1},univ::Universe) where N <: LieAlgebraFields
        clear!(F)
        calc_gaugeforce!(F,univ.U,univ)
        #add_gaugeforce!(F,univ.U,univ._temporal_gauge,univ._temporal_algebra[1]) 
        return
    end

    function calc_gaugeforce!(F::Array{N,1},U::Array{T,1},univ::Universe) where {N<: LieAlgebraFields, T<: GaugeFields} 
        clear!(F)
        calc_gaugeforce!(F,U,univ,univ.gparam)
        #add_gaugeforce!(F,U,univ._temporal_gauge,univ._temporal_algebra[1],univ.gparam) 
        return
    end

    function calc_gaugeforce!(F::Array{N,1},U::Array{T,1},univ::Universe,gparam::GP) where {N<: LieAlgebraFields, 
                                                                                            T<: GaugeFields,
                                                                                            GP <: GaugeActionParam_autogenerator} 
        clear!(F)
        add_gaugeforce!(F,U,univ._temporal_gauge,univ._temporal_algebra[1],gparam) 
        return
    end

    function calc_gaugeforce!(F::Array{N,1},U::Array{T,1},univ::Universe,gparam::GP) where {N<: LieAlgebraFields, 
                                                                                            T<: GaugeFields,
                                                                                            GP} 
        clear!(F)
        add_gaugeforce!(F,U,univ._temporal_gauge,univ._temporal_algebra[1]) 
        return
    end

    function Gaugefields.calc_GaugeAction(univ::Universe)
        Sg,plaq = calc_GaugeAction(univ.U,univ.gparam,univ._temporal_gauge)
        return real(Sg),real(plaq)
    end

    function calc_Action(univ::Universe)
        Sg,plaq = calc_GaugeAction(univ.U,univ.gparam,univ._temporal_gauge)
        
        SP = univ.p*univ.p/2
        #println("ek = $SP")
        S = Sg + SP
        #println("eym = $plaq")
        #println("e = $S")

        if univ.quench == false
        #if univ.fparam != nothing
            Sf = univ.η*univ.η
            S += Sf
        end
        return real(S),real(plaq)
    end


    function calc_IntegratedFermionAction(univ::Universe;debug = false,M=16,m=100,nc=20,nonc=false)
        println("Making W^+W matrix...")
        NX = univ.U[1].NX
        NY = univ.U[1].NY
        NZ = univ.U[1].NZ
        NT = univ.U[1].NT
        NC = univ.U[1].NC
        T = eltype(univ._temporal_fermi)
        if T == StaggeredFermion
            NG = 1
        else
            NG = 4
        end
        Nsize = NX*NY*NZ*NT*NC*NG

        WdagW = zeros(ComplexF64,Nsize,Nsize)
        @time  make_WdagWmatrix!(WdagW,univ)
        #@time WdagW = make_WdagWmatrix(univ)
        Sfnew = calc_IntegratedFermionAction(univ,WdagW,debug = debug,M=M,m=m,nc=nc,nonc=nonc)
        return Sfnew
    end


    function calc_IntegratedFermionAction(univ::Universe,WdagW;debug = false,M=16,m=100,nc=20,nonc=false)
        if debug

            #WdagW = DdagD_operator(univ.U,univ.η,univ.fparam)
            #println("Calculating approximated logdet")
            #@time Sfnew = -tdlogdet(WdagW,M,m,nc=nc,nonc=nonc)
            #exit()

            WdagWsp = sparse(WdagW)
            N,_ = size(WdagW)
            tempvecs = Array{Array{ComplexF64,1}}(undef,3)
            for i=1:3
                tempvecs[i] = zeros(ComplexF64,N)
            end
            println("Calculating approximated logdet")

            #@time Sfnew = -tdlogdet(WdagWsp,M2,m2,tempvecs,nc=nc2,nonc=nonc)

            #M2 = 100
            #m2 = 100
            #nc2 = 100
            @time Sfnew = -tdlogdet(WdagWsp,M,m,tempvecs,nc=nc,nonc=nonc)

            #exit()
            #=


            M = 8*2
            #M = 4#8*2
            m = 100
            NX = univ.U[1].NX
            #m = 20
            N,_ = size(WdagW)
            tempvecs = Array{Array{ComplexF64,1}}(undef,3)
            for i=1:3
                tempvecs[i] = zeros(ComplexF64,N)
            end

            nc = 20
            nonc = false

            @time ldet = logdet(WdagW)
            if nonc
                fp = open("deltavalue_m$(m)_nonc_NX$(NX)_6.dat","w")
            else
                fp = open("deltavalue_m$(m)_nc$(nc)_NX$(NX)_6.dat","w")
            end

            for M in [4,8,16,32,64,128]
                if nonc
                    filename="delta_$(NX)$(NX)$(NX)$(NX)_$(M)_m$(m)_nonc_6.dat"
                else
                    filename="delta_$(NX)$(NX)$(NX)$(NX)_$(M)_m$(m)_nc$(nc)_6.dat"
                end
                a = tdlogdet(WdagWsp,M,m,tempvecs,filename=filename,nc=nc,nonc=nonc)
                println("$M $a $ldet")
                println(fp,"$M $a $ldet")
            end
            close(fp)

            exit()
            @time Sfnew = -tdlogdet(WdagWsp,M,m,tempvecs)
            =#
        else
            println("Calculating logdet")
            @time Sfnew = -real(logdet(WdagW))
        end

       

        if univ.Dirac_operator == "Staggered" 
            #if univ.fparam.Nf == 4
            #    Sfnew /= 2
            #end
            Sfnew /= (8/univ.fparam.Nf)
        end
        return Sfnew
    end

    struct Wilsonloops_actions{Gauge_temp}
        loops::Array{Wilson_loop_set,1}
        loopaction::Gauge_temp
        temps::Array{Gauge_temp,1}
        numloops::Int64
        β::Float64
        NC::Float64
        couplinglist::Array{String,1}

        function Wilsonloops_actions(univ)
            couplinglist = ["plaq","rect","polyx","polyy","polyz","polyt"]
            L = (univ.U[1].NX,univ.U[1].NY,univ.U[1].NZ,univ.U[1].NT)
            loops = make_loopforactions(couplinglist,L)
            loopaction = similar(univ._temporal_gauge[1])
            sutype = typeof(loopaction)
            temps = Array{sutype,1}(undef,3)
            for i=1:3
                temps[i] = similar(loopaction)
            end
            numloops = length(couplinglist)
            β = univ.gparam.β
            NC = univ.NC

            return new{sutype}(loops,loopaction,temps,numloops,β,NC,couplinglist)
        end
    end


    function calc_looptrvalues(w::Wilsonloops_actions,univ)
        numloops = w.numloops
        trs = zeros(ComplexF64,numloops)
        for i=1:numloops
            evaluate_wilson_loops!(w.loopaction,w.loops[i],univ.U,w.temps[1:3])
            trs[i] = tr(w.loopaction)
        end
        return trs
    end

    function calc_looptrvalues_site(w::Wilsonloops_actions,univ)
        return 
        numloops = w.numloops
        NX = univ.U[1].NX
        NY = univ.U[1].NY
        NZ = univ.U[1].NZ
        NT = univ.U[1].NT
        NC = univ.NC
        V = zeros(ComplexF64,NC,NC)

        trs = zeros(ComplexF64,numloops,NX,NY,NZ,NT)
        for it = 1:NT
            for iz = 1:NZ
                for iy = 1:NY
                    for ix=1:NX
                        for i=1:numloops
                            evaluate_wilson_loops!(V,w.loops[i],univ.U,ix,iy,iz,it)
                            trs[i,ix,iy,iz,it] = tr(V)
                        end
                    end
                end
            end
        end
        
        return trs
    end

    function get_looptrvalues(w::Wilsonloops_actions)
        return w.trs
    end

    function calc_IntegratedSf(w::Wilsonloops_actions,univ)
        Sf = calc_IntegratedFermionAction(univ)
        return Sf
    end

    function calc_gaugeSg(w::Wilsonloops_actions,univ)
        trs = calc_looptrvalues(w,univ)
        Sg = (-w.β/w.NC)*trs[1]
        return real(Sg),trs
    end

    function calc_trainingdata(w::Wilsonloops_actions,univ)
        Sf = calc_IntegratedSf(w,univ)
        Sg,trs = calc_gaugeSg(w,univ)
        return trs,Sg,Sf
    end

    function make_WdagWmatrix(univ::Universe)
        return make_WdagWmatrix(univ.U,univ._temporal_fermi,univ.fparam)
    end

    function make_Wmatrix(univ::Universe)
        return make_Wmatrix(univ.U,univ._temporal_fermi,univ.fparam)
    end


    function make_WdagWmatrix!(WdagW,univ::Universe)
        make_WdagWmatrix!(WdagW,univ.U,univ._temporal_fermi,univ.fparam)
        return
    end

    function make_Wmatrix!(WdagW,univ::Universe)
        make_Wmatrix!(WdagW,univ.U,univ._temporal_fermi,univ.fparam)
        return
    end

    function make_WdagWmatrix(univ::Universe,U::Array{T,1}) where T <: GaugeFields
        return make_WdagWmatrix(U,univ._temporal_fermi,univ.fparam)
    end

    function make_Wmatrix(univ::Universe,U::Array{T,1}) where T <: GaugeFields
        return make_Wmatrix(U,univ._temporal_fermi,univ.fparam)
    end

    function make_WdagWmatrix(U::Array{G,1},temps::Array{T,1},fparam) where {G <: GaugeFields,T <:FermionFields}
        x0 = temps[7]
        xi = temps[8]
        return make_WdagWmatrix(U,x0,xi,fparam) 
    end

    function make_Wmatrix(U::Array{G,1},temps::Array{T,1},fparam) where {G <: GaugeFields,T <:FermionFields}
        x0 = temps[7]
        xi = temps[8]
        return make_Wmatrix(U,x0,xi,fparam) 
    end


    function make_WdagWmatrix!(WdagW,U::Array{G,1},temps::Array{T,1},fparam) where {G <: GaugeFields,T <:FermionFields}
        x0 = temps[7]
        xi = temps[8]
        make_WdagWmatrix!(WdagW,U,x0,xi,fparam) 
        return
    end

    function make_Wmatrix!(WdagW,U::Array{G,1},temps::Array{T,1},fparam) where {G <: GaugeFields,T <:FermionFields}
        x0 = temps[7]
        xi = temps[8]
        make_Wmatrix!(WdagW,U,x0,xi,fparam) 
        return
    end

    function make_WdagWmatrix(U::Array{G,1},x0::T,xi::T,fparam) where {G <: GaugeFields,T <:FermionFields}
        NX = x0.NX
        NY = x0.NY
        NZ = x0.NZ
        NT = x0.NT
        NC = x0.NC
        if T == StaggeredFermion
            NG = 1
        else
            NG = 4
        end
        Nsize = NX*NY*NZ*NT*NC*NG
        
        WdagW = zeros(ComplexF64,Nsize,Nsize)
        make_WdagWmatrix!(WdagW,U,x0,xi,fparam)
        return WdagW
    end

    function make_Wmatrix(U::Array{G,1},x0::T,xi::T,fparam) where {G <: GaugeFields,T <:FermionFields}
        NX = x0.NX
        NY = x0.NY
        NZ = x0.NZ
        NT = x0.NT
        NC = x0.NC
        if T == StaggeredFermion
            NG = 1
        else
            NG = 4
        end
        Nsize = NX*NY*NZ*NT*NC*NG
        
        W = zeros(ComplexF64,Nsize,Nsize)
        make_Wmatrix!(W,U,x0,xi,fparam)
        return W
    end
    
    function make_WdagWmatrix!(WdagW,U::Array{G,1},x0::T,xi::T,fparam) where {G <: GaugeFields,T <:FermionFields}
        
        Fermionfields.clear!(x0)
        NX = x0.NX
        NY = x0.NY
        NZ = x0.NZ
        NT = x0.NT
        NC = x0.NC
        if T == StaggeredFermion
            NG = 1
        else
            NG = 4
        end
        Nsize = NX*NY*NZ*NT*NC*NG
        #WdagW = zeros(ComplexF64,Nsize,Nsize)
        WdagWoperator = DdagD_operator(U,x0,fparam)
        
        
        j = 0
        for α=1:NG
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            for ic=1:NC
                                Fermionfields.clear!(x0)
                                j += 1
                                x0[ic,ix,iy,iz,it,α] = 1
                                #set_wing_fermi!(x0)
                                
                                #WdagWx!(xi,U,x0,temps,fparam)
                                #println(xi*xi)

                                #Wx!(xi,U,x0,temps,fparam)
                                #println(xi*xi)

                                #Wx!(xi,U,x0,temps,fparam,(ix,iy,iz,it,α))
                                #println(xi*xi)
                                
                                mul!(xi,WdagWoperator,x0,(ix,iy,iz,it,α))
                                #WdagWx!(xi,U,x0,temps,fparam,(ix,iy,iz,it,α))
                                #println("fast ",xi*xi)
                                #exit()
                                #display(xi)
                                #exit()
                                substitute_fermion!(WdagW,j,xi)
                        
                            end
                        end
                    end
                end
            end            
        end


        #return WdagW

    end

        
    function make_Wmatrix!(WdagW,U::Array{G,1},x0::T,xi::T,fparam) where {G <: GaugeFields,T <:FermionFields}
        
        Fermionfields.clear!(x0)
        NX = x0.NX
        NY = x0.NY
        NZ = x0.NZ
        NT = x0.NT
        NC = x0.NC
        if T == StaggeredFermion
            NG = 1
        else
            NG = 4
        end
        Nsize = NX*NY*NZ*NT*NC*NG
        #WdagW = zeros(ComplexF64,Nsize,Nsize)
        WdagWoperator = Diracoperators.Dirac_operator(U,x0,fparam)
        
        
        j = 0
        for α=1:NG
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            for ic=1:NC
                                Fermionfields.clear!(x0)
                                j += 1
                                x0[ic,ix,iy,iz,it,α] = 1
                                #set_wing_fermi!(x0)
                                
                                #WdagWx!(xi,U,x0,temps,fparam)
                                #println(xi*xi)

                                #Wx!(xi,U,x0,temps,fparam)
                                #println(xi*xi)

                                #Wx!(xi,U,x0,temps,fparam,(ix,iy,iz,it,α))
                                #println(xi*xi)
                                
                                #mul!(xi,WdagWoperator,x0,(ix,iy,iz,it,α))
                                mul!(xi,WdagWoperator,x0)
                                #WdagWx!(xi,U,x0,temps,fparam,(ix,iy,iz,it,α))
                                #println("fast ",xi*xi)
                                #exit()
                                #display(xi)
                                #exit()
                                substitute_fermion!(WdagW,j,xi)
                        
                            end
                        end
                    end
                end
            end            
        end


        #return WdagW

    end


    
    function construct_fermion_gauss_distribution!(univ::Universe{Gauge,Lie,Fermi,GaugeP,FermiP,Gauge_temp})  where {Gauge,Lie,Fermi,GaugeP,FermiP,Gauge_temp}
        gauss_distribution_fermi!(univ.η,univ.ranf)
    end

    function construct_fermion_gauss_distribution!(univ::Universe{Gauge,Lie,Fermi,GaugeP,FermiP,Gauge_temp})  where {Gauge,Lie,Fermi <: StaggeredFermion,GaugeP,FermiP,Gauge_temp}
        gauss_distribution_fermi!(univ.η,univ.ranf)
        if univ.fparam.Nf == 4
            U,_... = calc_smearingU(univ.U,univ.fparam.smearing)

            evensite = false
            W = Diracoperators.Dirac_operator(U,univ.η,univ.fparam)
            mul!(univ.φ,W',univ.η)

            #Wdagx!(univ.φ,univ.U,univ.η,univ._temporal_fermi,univ.fparam)
            Fermionfields.clear!(univ.φ,evensite)

            bicg(univ.η,W',univ.φ,eps = univ.fparam.eps,maxsteps= univ.fparam.MaxCGstep)
        end
    end

    function construct_fermionfield_φ!(univ::Universe{Gauge,Lie,Fermi,GaugeP,FermiP,Gauge_temp})  where {Gauge,Lie,Fermi,GaugeP,FermiP,Gauge_temp}
        U,_... = calc_smearingU(univ.U,univ.fparam.smearing)


        W = Diracoperators.Dirac_operator(U,univ.η,univ.fparam)
        mul!(univ.φ,W',univ.η)
        set_wing_fermi!(univ.φ)
    end

    function construct_fermionfield_φ!(univ::Universe{Gauge,Lie,Fermi,GaugeP,FermiP,Gauge_temp})  where {Gauge,Lie,Fermi <: StaggeredFermion,GaugeP,FermiP,Gauge_temp}
        U,_... = calc_smearingU(univ.U,univ.fparam.smearing)
        
        if univ.fparam.Nf == 4 || univ.fparam.Nf == 8
        
            W = Diracoperators.Dirac_operator(U,univ.η,univ.fparam)
            mul!(univ.φ,W',univ.η)
        else
            WdagW = Diracoperators.DdagD_operator(U,univ.η,univ.fparam)
            N = get_order(univ.fparam.rhmc_action)

            x = univ.φ
            vec_x = univ._temporal_fermi[end-N+1:end]
            for j=1:N
                Fermionfields.clear!(vec_x[j])
            end

            # eta =  (MdagM)^{-alpha/2} phi -> phi = (MdagM)^{alpha/2} eta
            #(MdagM)^{alpha/2} eta ~ α0 eta + sum_k αk (MdagM + βk)^{-1} eta
            vec_β = get_β(univ.fparam.rhmc_action)
            vec_α = get_α(univ.fparam.rhmc_action)
            α0 = get_α0(univ.fparam.rhmc_action)
            shiftedcg(vec_x,vec_β,x,WdagW,univ.η,eps = univ.fparam.eps,maxsteps= univ.fparam.MaxCGstep)
            Fermionfields.clear!(univ.φ)
            Fermionfields.add!(univ.φ,α0,univ.η)
            for j=1:N
                αk = vec_α[j]
                Fermionfields.add!(univ.φ,αk,vec_x[j])
            end


        end

        set_wing_fermi!(univ.φ)
    end

    
    function construct_fermionfield_η!(univ::Universe{Gauge,Lie,Fermi,GaugeP,FermiP,Gauge_temp})  where {Gauge,Lie,Fermi,GaugeP,FermiP,Gauge_temp}
        U,_... = calc_smearingU(univ.U,univ.fparam.smearing)


        W = Diracoperators.Dirac_operator(U,univ.η,univ.fparam)
        bicg(univ.η,W',univ.φ,eps = univ.fparam.eps,maxsteps= univ.fparam.MaxCGstep,verbose = univ.kind_of_verboselevel)
    end

    function construct_fermionfield_η!(univ::Universe{Gauge,Lie,Fermi,GaugeP,FermiP,Gauge_temp})  where {Gauge,Lie,Fermi <: StaggeredFermion,GaugeP,FermiP,Gauge_temp} # eta = Wdag^{-1} phi
        U,_... = calc_smearingU(univ.U,univ.fparam.smearing)


        if univ.fparam.Nf == 4 || univ.fparam.Nf == 8
            W = Diracoperators.Dirac_operator(U,univ.η,univ.fparam)
            bicg(univ.η,W',univ.φ,eps = univ.fparam.eps,maxsteps= univ.fparam.MaxCGstep,verbose = univ.kind_of_verboselevel)
        else
            WdagW = Diracoperators.DdagD_operator(U,univ.η,univ.fparam)
            N = get_order(univ.fparam.rhmc_action)

            x = univ.η
            vec_x = univ._temporal_fermi[end-N+1:end]
            for j=1:N
                Fermionfields.clear!(vec_x[j])
            end

            # eta =  (MdagM)^{-alpha/2} phi -> phi = (MdagM)^{alpha/2} eta
            #(MdagM)^{alpha/2} eta ~ α0 eta + sum_k αk (MdagM + βk)^{-1} eta
            vec_β = get_β_inverse(univ.fparam.rhmc_action)
            vec_α = get_α_inverse(univ.fparam.rhmc_action)
            α0 = get_α0_inverse(univ.fparam.rhmc_action)
            shiftedcg(vec_x,vec_β,x,WdagW,univ.φ,eps = univ.fparam.eps,maxsteps= univ.fparam.MaxCGstep)
            Fermionfields.clear!(univ.η)
            Fermionfields.add!(univ.η,α0,univ.φ)
            for j=1:N
                αk = vec_α[j]
                Fermionfields.add!(univ.η,αk,vec_x[j])
            end
        end
    end





end