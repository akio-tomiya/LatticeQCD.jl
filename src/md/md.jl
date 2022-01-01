module MD
    using LinearAlgebra
    using Base.Threads
    using Random,Distributions

    import ..Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3


    import ..Actions:GaugeActionParam_standard,
        GaugeActionParam_autogenerator
    import ..Gaugefield:Stoutsmearing,Nosmearing
                        #FermiActionParam_WilsonClover,FermiActionParam_Wilson,
                        
    import ..LTK_universe:Universe,calc_Action,gauss_distribution,make_WdagWmatrix,set_β!,
                calc_IntegratedFermionAction,construct_fermion_gauss_distribution!,
                construct_fermionfield_φ!,construct_fermionfield_η!
    import ..Gaugefields:GaugeFields,substitute!,SU3GaugeFields,set_wing!,
                            projlink!,make_staple_double!,
                            GaugeFields_1d,SU3GaugeFields_1d,SU2GaugeFields_1d,calc_Plaq_notrace_1d,
                            calc_GaugeAction,apply_smearing,calc_smearingU
    import ..Gaugefields
    import ..LieAlgebrafields:gauss_distribution_lie!,LieAlgebraFields,SU3AlgebraFields,SU2AlgebraFields,
                                Gauge2Lie!,add!,add_gaugeforce!,expA!,stoutfource
    import ..Fermionfields:gauss_distribution_fermi!,set_wing_fermi!,Wdagx!,vvmat!,
            FermionFields,fermion_shift!,WilsonFermion, fermion_shiftB!,
            StaggeredFermion,Wx!,clear!,WdagWx!,substitute_fermion!
    #using ..Fermionfields
    import ..Heatbath:heatbath!
    import ..Gaugefield:AbstractGaugefields,exptU!,substitute_U!,set_wing_U!,
                        add_force!,TA_Gaugefields,staggered_U,clear_U!,
                        Traceless_antihermitian_add!,add_U!,back_prop,CovNeuralnet
    import ..Fermionfield_LQCD:AbstractFermionfields,
                                FermiActionParam_WilsonClover,FermiActionParam_Wilson,
                                clear_fermion!,set_wing_fermion!,shift_fermion



    #import ..CGfermion:cg0!,cg0_WdagW!,shiftedcg0_WdagW!
    import ..Clover:Make_CloverFμν!,dSclover!
    import ..System_parameters:Params

    import ..Fermionfield_LQCD:Wilson_operator,Adjoint_Wilson_operator,WilsonClover_operator,
                Dirac_operator,DdagD_operator
    import ..Fermionfield_LQCD:bicg,cg,shiftedcg
    import ..RationalApprox:calc_exactvalue,calc_Anϕ
    import ..Fermionfield_LQCD:get_order,get_β,get_α,get_α0,get_β_inverse,get_α_inverse,get_α0_inverse


    abstract type MD_parameters end

    struct MD_parameters_SextonWeingargten  <: MD_parameters
        Δτ ::Float64
        MDsteps::Int64 
        N_SextonWeingargten::Int64
        QPQ::Bool

        function MD_parameters_SextonWeingargten(Δτ,MDsteps,N_SextonWeingargten;QPQ=true)
            return new(Δτ,MDsteps,N_SextonWeingargten,QPQ)
        end
    end

    struct MD_parameters_standard  <: MD_parameters
        Δτ ::Float64
        MDsteps::Int64 
        QPQ::Bool

        function MD_parameters_standard(Δτ,MDsteps;QPQ=true)
            return new(Δτ,MDsteps,QPQ)
        end
    end

    struct MD_parameters_IntegratedHMC  <: MD_parameters
        Δτ ::Float64
        MDsteps::Int64 
    end

    struct MD_parameters_IntegratedHB  <: MD_parameters
    end

    mutable struct MD_parameters_SLMC  <: MD_parameters
        βeff::Union{Float64,Array{Float64,1}}
        debug::Bool
        function MD_parameters_SLMC(βeff)
            debug = false
            #debug = true
            return new(βeff,debug)
        end
    end

    mutable struct MD_parameters_SLHMC  <: MD_parameters
        Δτ ::Float64
        MDsteps::Int64 
        βeff::Union{Float64,Array{Float64,1}}
    end

    function calc_factor(univ::Universe,mdparams::MD_parameters)
        return -univ.gparam.β*mdparams.Δτ/ (2*univ.NC)
    end

    

    function construct_MD_parameters(p::Params)
        #fac = - p.Δτ * p.β / (2*p.NC)
        if p.SextonWeingargten 
            if p.quench == true
                error("The quench update does not need the SextonWeingargten method. Put SextonWeingargten = false")
            end
            #return MD_parameters_SextonWeingargten(p.Δτ,p.MDsteps,fac,p.β,p.N_SextonWeingargten)
            return MD_parameters_SextonWeingargten(p.Δτ,p.MDsteps,p.N_SextonWeingargten)
        elseif p.SextonWeingargten ==false && p.update_method == "IntegratedHMC"
            if p.quench == false
                error("quench = false. The IntegratedHMC needs the quench update. Put update_method != IntegratedHMC or quench = true")
            end
            #return MD_parameters_IntegratedHMC(p.Δτ,p.MDsteps,fac,p.β)
            return MD_parameters_IntegratedHMC(p.Δτ,p.MDsteps)
        elseif p.SextonWeingargten ==false && p.update_method == "SLHMC"
            if p.quench == false
                error("quench = false. The SLHMC needs the quench update. Put update_method != SLHMC or quench = true")
            end
            #return MD_parameters_IntegratedHMC(p.Δτ,p.MDsteps,fac,p.β)
            return MD_parameters_SLHMC(p.Δτ,p.MDsteps,p.βeff)  
        elseif  p.SextonWeingargten ==false && p.update_method == "IntegratedHB"
            if p.quench == false
                error("quench = false. The IntegratedHB needs the quench update. Put update_method != IntegratedHB or quench = true")
            end
            #return MD_parameters_IntegratedHMC(p.Δτ,p.MDsteps,fac,p.β)
            return MD_parameters_IntegratedHB()
        elseif p.SextonWeingargten ==false && p.update_method == "SLMC"
            if p.quench == false
                error("quench = false. The SLMC needs the quench update. Put update_method != SLMC or quench = true")
            end
            #return MD_parameters_IntegratedHMC(p.Δτ,p.MDsteps,fac,p.β)
            return MD_parameters_SLMC(p.βeff)  
        else
            #return MD_parameters_standard(p.Δτ,p.MDsteps,fac,p.β)
            return MD_parameters_standard(p.Δτ,p.MDsteps)
        end
    end

    function MD_parameters_standard(gparam::GaugeActionParam_standard)
        #βMD = gparam.β
        MDsteps = 10
        Δτ = 0.1
        
        #fac = - Δτ * βMD / (2*gparam.NTRACE)
        return MD_parameters_standard(gparam,Δτ,MDsteps)#,βMD)
    end

    function MD_parameters_standard(gparam::GaugeActionParam_standard,Δτ,MDsteps,βMD)
        #fac = - Δτ * βMD / (2*gparam.NTRACE)
        return MD_parameters_standard(Δτ,MDsteps)#,fac,βMD)
    end



    function md!(univ::Universe,mdparams::MD_parameters_standard)
        #Sold = md_initialize!(univ::Universe)
        QPQ = mdparams.QPQ #true
        if QPQ 
            for istep=1:mdparams.MDsteps
                #println("istep = $istep")
                updateU!(univ,mdparams,0.5)

                updateP!(univ,mdparams,1.0)


                if univ.quench == false
                #if univ.fparam != nothing
                    updateP_fermi!(univ,mdparams,1.0)
                end
                updateU!(univ,mdparams,0.5)

                #Sold,Sg = calc_Action(univ)
                #println("$istep Sold,plaq = ",Sold,"\t",Sg)

                #exit()
            end
        else
            for istep=1:mdparams.MDsteps
                #println("istep = $istep")
                updateP!(univ,mdparams,0.5)


                if univ.quench == false
                #if univ.fparam != nothing
                    updateP_fermi!(univ,mdparams,0.5)
                end


                updateU!(univ,mdparams,1.0)

                updateP!(univ,mdparams,0.5)


                if univ.quench == false
                #if univ.fparam != nothing
                    updateP_fermi!(univ,mdparams,0.5)
                end
                


            end
        end

        if univ.quench == false
            construct_fermionfield_η!(univ)
        end


        Snew,plaq = calc_Action(univ)
        return Snew
    end

    function md!(univ::Universe,mdparams::MD_parameters_SextonWeingargten)
        #Sold = md_initialize!(univ::Universe)
        QPQ = mdparams.QPQ
        if QPQ != true
            for istep=1:mdparams.MDsteps
                
                if univ.quench == false
                    updateP_fermi!(univ,mdparams,0.5)
                end
                for istep_SW=1:mdparams.N_SextonWeingargten
                    #println("istep_SW = $(istep_SW)")
                    updateP!(univ,mdparams,0.5/mdparams.N_SextonWeingargten)
                    updateU!(univ,mdparams,1/mdparams.N_SextonWeingargten)
                    updateP!(univ,mdparams,0.5/mdparams.N_SextonWeingargten)
                end
                if univ.quench == false
                    updateP_fermi!(univ,mdparams,0.5)
                end
            end
        elseif QPQ
            @assert mdparams.N_SextonWeingargten % 2 == 0 "mdparams.N_SextonWeingargten % 2 should be 0!"

            for istep=1:mdparams.MDsteps

                for istep_SW=1:div(mdparams.N_SextonWeingargten,2)
                    updateU!(univ,mdparams,0.5/mdparams.N_SextonWeingargten)
                    updateP!(univ,mdparams,1/mdparams.N_SextonWeingargten)
                    updateU!(univ,mdparams,0.5/mdparams.N_SextonWeingargten)
                end

                if univ.quench == false
                    updateP_fermi!(univ,mdparams,1.0)
                end

                for istep_SW=1:div(mdparams.N_SextonWeingargten,2)
                    updateU!(univ,mdparams,0.5/mdparams.N_SextonWeingargten)
                    updateP!(univ,mdparams,1/mdparams.N_SextonWeingargten)
                    updateU!(univ,mdparams,0.5/mdparams.N_SextonWeingargten)
                end
            
            end
        end

        if univ.quench == false
            construct_fermionfield_η!(univ)
        end

        Snew,plaq = calc_Action(univ)

        return Snew
    end

    function md!(univ::Universe,Sfold,Sgold,mdparams::MD_parameters_IntegratedHMC)
        @assert univ.quench == true "quench should be true!"
        if Sfold == nothing       
            
            WdagW = make_WdagWmatrix(univ)
            #e,_ = eigen(WdagW)
            #println(e)
            Sfold = -real(logdet(WdagW))
            #Uold = deepcopy(univ.U)
            if univ.Dirac_operator == "Staggered" 
                if univ.fparam.Nf == 4
                    Sfold /= 2
                end
            end
            
        end

        for istep=1:mdparams.MDsteps
            updateU!(univ,mdparams,0.5)

            updateP!(univ,mdparams,1.0)

            updateU!(univ,mdparams,0.5)
        end


        Sgnew,plaq = calc_Action(univ)
        println("Making W^+W matrix...")
        @time WdagW = make_WdagWmatrix(univ)
        println("Calculating logdet")
        @time Sfnew = -real(logdet(WdagW))
        if univ.Dirac_operator == "Staggered" 
            if univ.fparam.Nf == 4
                Sfnew /= 2
            end
        end


            
        #

        #println("Snew,plaq = ",Snew,"\t",plaq)
        println("Sgnew,Sfnew,Sgold,Sfold: ", Sgnew,"\t",Sfnew,"\t",Sgold,"\t",Sfold)
        return Sgnew,Sfnew,Sgold,Sfold
    end

    function md!(univ::Universe,Sfold,Sgold,mdparams::MD_parameters_IntegratedHB)
        @assert univ.quench == true "quench should be true!"
        @assert univ.NC == 2 "Only SU(2) is supported now."
        if Sfold == nothing       
            @time Sfold = calc_IntegratedFermionAction(univ)    
             
            #=
            WdagW = make_WdagWmatrix(univ)
            #e,_ = eigen(WdagW)
            #println(e)
            Sfold = -real(logdet(WdagW))
            #Uold = deepcopy(univ.U)
            if univ.Dirac_operator == "Staggered" 
                if univ.fparam.Nf == 4
                    Sfold /= 2
                end
            end
            =#
        end

        Sgold,plaq = calc_GaugeAction(univ)

        Sgeffold = Sgold

        heatbath!(univ)


        Sgnew,plaq = calc_GaugeAction(univ)

        Sgeffnew = Sgnew
        #=
        println("Making W^+W matrix...")
        @time WdagW = make_WdagWmatrix(univ)
        println("Calculating logdet")
        @time Sfnew = -real(logdet(WdagW))
        if univ.Dirac_operator == "Staggered" 
            if univ.fparam.Nf == 4
                Sfnew /= 2
            end
        end
        =#
        @time Sfnew = calc_IntegratedFermionAction(univ)


        
        #

        #println("Snew,plaq = ",Snew,"\t",plaq)
        println("Sgnew,Sfnew,Sgold,Sfold: ", Sgnew,"\t",Sfnew,"\t",Sgold,"\t",Sfold)
        return Sgnew,Sfnew,Sgeffnew,Sgold,Sfold,Sgeffold
    end

    function md!(univ::Universe,Sfold,Sgold,mdparams::MD_parameters_SLMC)
        @assert univ.quench == true "quench should be true!"



        
        #@assert univ.NC == 2 "Only SU(2) is supported now."
        if Sfold == nothing       
            @time Sfold = calc_IntegratedFermionAction(univ,debug=mdparams.debug)
            #=

            WdagW = make_WdagWmatrix(univ)
            #e,_ = eigen(WdagW)
            #println(e)
            Sfold = -real(logdet(WdagW))
            #Uold = deepcopy(univ.U)
            if univ.Dirac_operator == "Staggered" 
                if univ.fparam.Nf == 4
                    Sfold /= 2
                end
            end
            =#
        end

        β = univ.gparam.β
        Sgold,plaq = calc_GaugeAction(univ)

        set_β!(univ,mdparams.βeff)
        Sgeffold,plaq = calc_GaugeAction(univ)

        heatbath!(univ)
        #heatbath!(univ)

        Sgeffnew,plaq = calc_GaugeAction(univ)
        set_β!(univ,β)

        Sgnew,plaq = calc_GaugeAction(univ)



        @time Sfnew = calc_IntegratedFermionAction(univ,debug=mdparams.debug)
        #println("Sfnew comparison: $Sfnew,$Sfnew2")

        

        
        
        #

        #println("Snew,plaq = ",Snew,"\t",plaq)
        println("Sgnew,Sfnew,Sgeffnew,Sgold,Sfold,Sgeffold: ", Sgnew,"\t",Sfnew,"\t",Sgeffnew,"\t",Sgold,"\t",Sfold,"\t",Sgeffold)
        return Sgnew,Sfnew,Sgeffnew,Sgold,Sfold,Sgeffold
    end

    function md!(univ::Universe,Sfold,Sgold,mdparams::MD_parameters_SLHMC)
        @assert univ.quench == true "quench should be true!"
        if Sfold == nothing  
            

            

            WdagW = make_WdagWmatrix(univ)
            #e,_ = eigen(WdagW)
            
            Sfold = -real(logdet(WdagW))
            #Uold = deepcopy(univ.U)
            if univ.Dirac_operator == "Staggered" 
                if univ.fparam.Nf == 4
                    Sfold /= 2
                end
            end
        end

        β = univ.gparam.β
        set_β!(univ,mdparams.βeff)

        for istep=1:mdparams.MDsteps
            updateU!(univ,mdparams,0.5)

            updateP!(univ,mdparams,1.0)

            updateU!(univ,mdparams,0.5)
        end

        set_β!(univ,β)


        Sgnew,plaq = calc_Action(univ)
        println("Making W^+W matrix...")
        @time WdagW = make_WdagWmatrix(univ)
        println("Calculating logdet")

        
        @time Sfnew = -real(logdet(WdagW))
        if univ.Dirac_operator == "Staggered" 
            if univ.fparam.Nf == 4
                Sfnew /= 2
            end
        end


        #
        #println("Snew,plaq = ",Snew,"\t",plaq)
        println("Sgnew,Sfnew,Sgold,Sfold: ", Sgnew,"\t",Sfnew,"\t",Sgold,"\t",Sfold)
        return Sgnew,Sfnew,Sgold,Sfold,plaq
    end

    function md_initialize!(univ::Universe)
        substitute!(univ.Uold,univ.U)
        NV = univ.NV
        NDFALG = univ.NDFALG

        for μ=1:4
            pwork = gauss_distribution(univ,NV*NDFALG)
            substitute!(univ.p[μ],pwork)
            #display(univ.p[μ])
        end
        #gauss_distribution_lie!(univ.p)

        if univ.quench == false 
            #if univ.fparam != nothing 
            construct_fermion_gauss_distribution!(univ) #generate η
            #println(univ.η*univ.η)
            construct_fermionfield_φ!(univ)  #generate φ

        end

        Sold,Sg = calc_Action(univ)

        #println("-----------------------------------------------")
        #println("Sold,plaq in init = ",Sold,"\t",Sg)
        #println("-----------------------------------------------")
        #exit()
        return Sold
    end

    function metropolis_update!(univ::Universe,Sold,Snew)
        accept = exp(Sold - Snew) >= univ.ranf()
        println("Snew,Sold $Snew $Sold")
        println("Sold,Snew,Diff,accept: $Sold $Snew $(Sold-Snew) $accept")
        if accept
        else
            substitute!(univ.U,univ.Uold)
        end
        return accept
    end

    function updateP_fermi!(univ::Universe,mdparams::MD_parameters,τ)

        updateP_fermi!(univ.η,univ.φ,univ.ξ,univ.fparam,
            univ.p,mdparams,τ,univ.U,
            univ._temporal_gauge,univ._temporal_algebra,
            univ._temporal_fermi,kind_of_verboselevel =univ.kind_of_verboselevel
        )

        return
    end






    function updateP_fermi_fromX!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,U::Array{T,1},
        temps::Array{T_1d,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2()
        ) where {F <: WilsonFermion, T<: GaugeFields,N<: LieAlgebraFields,T_1d <: GaugeFields_1d} 
        temp0_f = temps_fermi[1] #F_field
        temp1_f = temps_fermi[2] #F_field
        temp2_g = temps[1] #G_field1
        temp3_g = temps[2] #G_field1
        c = temp_a[1]
        NV = temp2_g.NV
        


        W = Dirac_operator(U,φ,fparam)
        mul!(Y,W,X)
        set_wing_fermi!(Y)
        
        

        for μ=1:4
            #!  Construct U(x,mu)*P1

            # U_{k,μ} X_{k+μ}
            fermion_shift!(temp0_f,U,μ,X)
            # (r-γ_μ) U_{k,μ} X_{k+μ}
            mul!(temp0_f,view(X.rminusγ,:,:,μ),temp0_f)
            # κ (r-γ_μ) U_{k,μ} X_{k+μ}
            mul!(temp1_f,X.hopp[μ],temp0_f)

            # κ ((r-γ_μ) U_{k,μ} X_{k+μ}) ⊗ Y_k
            vvmat!(temp2_g,temp1_f,Y,1)


            #.....   Projection onto Lie Algebra   .....
            projlink!(temp3_g,temp2_g)

            Gauge2Lie!(c,temp3_g)


            #...  p(new) = p(old) + fac * c  .....
            add!(p[μ],τ*mdparams.Δτ,c)



            #!  Construct P2*U_adj(x,mu)
            # Y_{k+μ}^dag U_{k,μ}^dag
            fermion_shiftB!(temp0_f,U,-μ,Y)
            # Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)
            mul!(temp0_f,temp0_f,view(X.rplusγ,:,:,μ))

            # κ Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)
            mul!(temp1_f,X.hopm[μ],temp0_f)

            # X_k ⊗ κ Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)
            vvmat!(temp2_g,X,temp1_f,2)


            #.....   Projection onto Lie Algebra   .....
            projlink!(temp3_g,temp2_g)
            Gauge2Lie!(c,temp3_g)

            add!(p[μ],-τ*mdparams.Δτ,c)

            if typeof(fparam) == FermiActionParam_WilsonClover
                dSclover!(c,μ,X,Y,U,fparam,temps)
                add!(p[μ],-τ*mdparams.Δτ,c)
            end


        end

        

    end

    function updateP_fermi_fromX_smearing!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,U::Array{T,1},Uout_multi,dSdU,Uin,
        temps::Array{T_1d,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2()
        ) where {F <: WilsonFermion, T<: GaugeFields,N<: LieAlgebraFields,T_1d <: GaugeFields_1d} 
        temp0_f = temps_fermi[1] #F_field
        temp1_f = temps_fermi[2] #F_field
        temp2_g = temps[1] #G_field1
        temp3_g = temps[2] #G_field1
        c = temp_a[1]
        NV = temp2_g.NV
        


        W = Dirac_operator(U,φ,fparam)
        mul!(Y,W,X)
        set_wing_fermi!(Y)
        
        

        for μ=1:4
            #!  Construct U(x,mu)*P1

            # U_{k,μ} X_{k+μ}
            fermion_shift!(temp0_f,U,μ,X)
            # (r-γ_μ) U_{k,μ} X_{k+μ}
            mul!(temp0_f,view(X.rminusγ,:,:,μ),temp0_f)
            # κ (r-γ_μ) U_{k,μ} X_{k+μ}
            mul!(temp1_f,X.hopp[μ],temp0_f)

            # κ ((r-γ_μ) U_{k,μ} X_{k+μ}) ⊗ Y_k
            vvmat!(temp2_g,temp1_f,Y,1)
            mul!(dSdU[μ],U[μ]',temp2_g) #additional term


            #!  Construct P2*U_adj(x,mu)
            # Y_{k+μ}^dag U_{k,μ}^dag
            fermion_shiftB!(temp0_f,U,-μ,Y)
            # Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)
            mul!(temp0_f,temp0_f,view(X.rplusγ,:,:,μ))

            # κ Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)
            mul!(temp1_f,X.hopm[μ],temp0_f)

            # X_k ⊗ κ Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)
            vvmat!(temp2_g,X,temp1_f,2)

            mul!(temp3_g,U[μ]',temp2_g)
            #Gaugefields.add!(dSdU[μ],temp3_g)
            Gaugefields.muladd!(dSdU[μ],-1,temp3_g)


            if typeof(fparam) == FermiActionParam_WilsonClover
                dSclover!(c,μ,X,Y,U,fparam,temps)
                Gaugefields.muladd!(dSdU[μ],-1,c)
                #add!(p[μ],-τ*mdparams.Δτ,c)
            end



        end

        if typeof(fparam.smearing) <: SmearingParam_single
            dSdUnew,_ = stoutfource(dSdU,Uin,fparam.smearing) 
        elseif typeof(fparam.smearing) <: SmearingParam_multi
            dSdUnew,_ = stoutfource(dSdU,Uout_multi,Uin,fparam.smearing) 
        elseif typeof(fparam.smearing) <: Nosmearing
            dSdUnew = dSdU
        else
            error("$(typeof(fparam.smearing)) is not supported")
        end


        for μ=1:4

            mul!(temp2_g,Uin[μ],dSdUnew[μ])

            #.....   Projection onto Lie Algebra   .....
            projlink!(temp3_g,temp2_g)
            Gauge2Lie!(c,temp3_g)


            #...  p(new) = p(old) + fac * c  .....
            add!(p[μ],τ*mdparams.Δτ,c)

        end

        

    end

    function  updateP_fermi_fromX!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,U::Array{T,1},
        temps::Array{T_1d,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2(),coeff=1
        ) where {F <: StaggeredFermion, T<: GaugeFields,N<: LieAlgebraFields,T_1d <: GaugeFields_1d} 

        temp0_f = temps_fermi[1] #F_field
        temp1_f = temps_fermi[2] #F_field
        temp2_g = temps[1] #G_field1
        temp3_g = temps[2] #G_field1
        c = temp_a[1]
        NV = temp2_g.NV

        W = Dirac_operator(U,φ,fparam)
        mul!(Y,W,X)

        for μ=1:4
            #!  Construct U(x,mu)*P1

            # U_{k,μ} X_{k+μ}
            fermion_shift!(temp0_f,U,μ,X)


            # κ ((r-γ_μ) U_{k,μ} X_{k+μ}) ⊗ Y_k
            vvmat!(temp2_g,temp0_f,Y,1)


            #.....   Projection onto Lie Algebra   .....
            projlink!(temp3_g,temp2_g)

            Gauge2Lie!(c,temp3_g)


            #...  p(new) = p(old) + fac * c  .....
            add!(p[μ],-0.5*τ*mdparams.Δτ*coeff,c)



            #!  Construct P2*U_adj(x,mu)
            # Y_{k+μ}^dag U_{k,μ}^dag
            fermion_shiftB!(temp0_f,U,-μ,Y)
            # Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)

            # X_k ⊗ κ Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)
            vvmat!(temp2_g,X,temp0_f,2)


            #.....   Projection onto Lie Algebra   .....
            projlink!(temp3_g,temp2_g)
            Gauge2Lie!(c,temp3_g)

            add!(p[μ],-0.5*τ*mdparams.Δτ*coeff,c)


        end
        return

    end

    function  updateP_fermi_fromX_smearing!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,U::Array{T,1},Uout_multi,dSdU,Uin,
        temps::Array{T_1d,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2(),coeff=1
        ) where {F <: StaggeredFermion, T<: GaugeFields,N<: LieAlgebraFields,T_1d <: GaugeFields_1d} 

        temp0_f = temps_fermi[1] #F_field
        temp1_f = temps_fermi[2] #F_field
        temp2_g = temps[1] #G_field1
        temp3_g = temps[2] #G_field1
        c = temp_a[1]
        NV = temp2_g.NV

        W = Dirac_operator(U,φ,fparam)
        mul!(Y,W,X)

        for μ=1:4
            #!  Construct U(x,mu)*P1

            # U_{k,μ} X_{k+μ}
            fermion_shift!(temp0_f,U,μ,X)


            # κ ((r-γ_μ) U_{k,μ} X_{k+μ}) ⊗ Y_k
            vvmat!(temp2_g,temp0_f,Y,1)
            mul!(dSdU[μ],U[μ]',temp2_g) #additional term


            #!  Construct P2*U_adj(x,mu)
            # Y_{k+μ}^dag U_{k,μ}^dag
            fermion_shiftB!(temp0_f,U,-μ,Y)

            # X_k ⊗ κ Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)
            #vvmat!(UdSdU2[μ],X,temp0_f,2)
            #if fparam.smearing != nothing
            vvmat!(temp2_g,X,temp0_f,2)

            mul!(temp3_g,U[μ]',temp2_g)
            Gaugefields.add!(dSdU[μ],temp3_g)
        end

        
        if typeof(fparam.smearing) <: SmearingParam_single
            dSdUnew,_ = stoutfource(dSdU,Uin,fparam.smearing) 
        elseif typeof(fparam.smearing) <: SmearingParam_multi
            dSdUnew,_ = stoutfource(dSdU,Uout_multi,Uin,fparam.smearing) 
        else
            error("$(typeof(fparam.smearing)) is not supported")
        end
        

        for μ=1:4
            mul!(temp2_g,Uin[μ],dSdUnew[μ])
            projlink!(temp3_g,temp2_g)


            Gauge2Lie!(c,temp3_g)


            add!(p[μ],-0.5*τ*mdparams.Δτ*coeff,c)
        end
        return

    end

    function updateP_fermi!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,Uin::Array{T,1},
        temps::Array{T2,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2()
        ) where {F <: AbstractFermionfields, T<: AbstractGaugefields,N<: TA_Gaugefields,T2 <: AbstractGaugefields} 

        if Y.Dirac_operator == "Staggered"
            updateP_fermi_Staggered!(Y,φ,X,fparam,
                p,mdparams,τ,Uin,
                temps,temp_a,temps_fermi;kind_of_verboselevel = kind_of_verboselevel)
        elseif findfirst("Wilson",Y.Dirac_operator) != nothing #Dirac_operator == "Wilson"
            updateP_fermi_Wilson!(Y,φ,X,fparam,
                p,mdparams,τ,Uin,
                temps,temp_a,temps_fermi;kind_of_verboselevel = kind_of_verboselevel)
        else
            error("Dirac_operator = $(Y.Dirac_operator) is not supported")
        end
    end

    function updateP_fermi_Staggered!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,Uin::Array{T,1},
        temps::Array{T2,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2()
        ) where {F <: AbstractFermionfields, T<: AbstractGaugefields,N<: TA_Gaugefields,T2 <: AbstractGaugefields} 
        temp0_f = temps_fermi[1] #F_field
        temp1_f = temps_fermi[2] #F_field
        temp2_g = temps[1] #G_field1
        temp3_g = temps[2] #G_field1
        c = temp_a[1]
        NV = temp2_g.NV

        U,Uout_multi,dSdU = calc_smearingU(Uin,fparam.smearing,calcdSdU = true,temps = temps)
        #println("Uin ",Uin[1][1,1,1,1,1,1])
        #println("U ",U[1][1,1,1,1,1,1])

        WdagW = DdagD_operator(U,φ,fparam)

        if fparam.Nf == 4 || fparam.Nf == 8
            #X = (D^dag D)^(-1) ϕ 
            #
            cg(X,WdagW,φ,eps = fparam.eps,maxsteps= fparam.MaxCGstep,verbose = kind_of_verboselevel)
            set_wing_fermion!(X)
            

            if fparam.smearing != nothing && typeof(fparam.smearing) != Nosmearing
                updateP_fermi_fromX_smearing!(Y,φ,X,fparam,
                p,mdparams,τ,U,Uout_multi,dSdU,Uin,
                temps,temp_a,temps_fermi,kind_of_verboselevel = kind_of_verboselevel)
            else
                updateP_fermi_fromX!(Y,φ,X,fparam,
                p,mdparams,τ,U,
                temps,temp_a,temps_fermi,kind_of_verboselevel = kind_of_verboselevel)
            end
        else
            N_MD = get_order(fparam.rhmc_MD)
            #numtemps_fermi = length(temps_fermi)
            x = temps_fermi[end-N_MD]
            vec_x = temps_fermi[end-N_MD+1:end]
            for j=1:N_MD
                clear_fermion!(vec_x[j])
            end
            vec_β = get_β_inverse(fparam.rhmc_MD)
            vec_α = get_α_inverse(fparam.rhmc_MD)
            α0 = get_α0_inverse(fparam.rhmc_MD)
            #println("φ ",φ[1,1,1,1,1,1])
            shiftedcg(vec_x,vec_β,x,WdagW,φ,eps = fparam.eps,maxsteps= fparam.MaxCGstep)
            for j=1:N_MD
                set_wing_fermion!(vec_x[j])
                #println("X ",vec_x[j][1,1,1,1,1,1])
                if fparam.smearing != nothing && typeof(fparam.smearing) != Nosmearing
                    updateP_fermi_fromX_smearing!(Y,φ,vec_x[j],fparam,
                        p,mdparams,τ,U,Uout_multi,dSdU,Uin,
                        temps,temp_a,temps_fermi,kind_of_verboselevel = kind_of_verboselevel,coeff=vec_α[j])
                else
                    updateP_fermi_fromX!(Y,φ,vec_x[j],fparam,
                        p,mdparams,τ,U,
                        temps,temp_a,temps_fermi,kind_of_verboselevel = kind_of_verboselevel,coeff=vec_α[j])
                end
            end
        end



    end



    function updateP_fermi_Wilson!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,Uin::Array{T,1},
        temps::Array{T2,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2()
        ) where {F <: AbstractFermionfields, T<: AbstractGaugefields,N<: TA_Gaugefields,T2 <: AbstractGaugefields} 
        temp0_f = temps_fermi[1] #F_field
        temp1_f = temps_fermi[2] #F_field
        temp2_g = temps[1] #G_field1
        temp3_g = temps[2] #G_field1
        c = temp_a[1]
        NV = temp2_g.NV

        U,Uout_multi,dSdU = calc_smearingU(Uin,fparam.smearing,calcdSdU = true,temps = temps)
        WdagW = DdagD_operator(U,φ,fparam)
        cg(X,WdagW,φ,eps = fparam.eps,maxsteps= fparam.MaxCGstep,verbose = kind_of_verboselevel)
        set_wing_fermion!(X)

        if fparam.smearing != nothing && typeof(fparam.smearing) != Nosmearing
            updateP_fermi_fromX_smearing!(Y,φ,X,fparam,
            p,mdparams,τ,U,Uout_multi,dSdU,Uin,
            temps,temp_a,temps_fermi,kind_of_verboselevel = kind_of_verboselevel)
        else
            updateP_fermi_fromX!(Y,φ,X,fparam,
            p,mdparams,τ,U,
            temps,temp_a,temps_fermi,kind_of_verboselevel = kind_of_verboselevel
            )
        end

    end

    function  updateP_fermi_fromX!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,U::Array{T,1},
        temps::Array{T2,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2(),coeff=1
        ) where {F <: AbstractFermionfields, T<: AbstractGaugefields,N<: TA_Gaugefields,T2 <: AbstractGaugefields} 

        if Y.Dirac_operator == "Staggered"
            updateP_fermi_fromX_Staggered!(Y,φ,X,fparam,
                p,mdparams,τ,U,
                temps,temp_a,temps_fermi;kind_of_verboselevel = kind_of_verboselevel,coeff=coeff
                )
        elseif findfirst("Wilson",Y.Dirac_operator) != nothing #Dirac_operator == "Wilson"
            updateP_fermi_fromX_Wilson!(Y,φ,X,fparam,
                p,mdparams,τ,U,
                temps,temp_a,temps_fermi;kind_of_verboselevel = kind_of_verboselevel,coeff=coeff
                )
        else
            error("Dirac_operator = $(Y.Dirac_operator) is not supported")
        end
    end

    function  updateP_fermi_fromX_Staggered!(Y::F,φ::F,X::AbstractFermionfields{NC,Dim},fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,U::Array{T,1},
        temps::Array{T2,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2(),coeff=1
        ) where {NC,Dim,F <: AbstractFermionfields, T<: AbstractGaugefields,N<: TA_Gaugefields,T2 <: AbstractGaugefields} 
        temp0_f = temps_fermi[1] #F_field
        temp1_f = temps_fermi[2] #F_field
        temp2_g = temps[1] #G_field1
        temp3_g = temps[2] #G_field1
        c = temp_a[1]
        NV = temp2_g.NV

        W = Dirac_operator(U,φ,fparam)
        mul!(Y,W,X)

        for μ=1:Dim
            #!  Construct U(x,mu)*P1

            # U_{k,μ} X_{k+μ}
            Xplus = shift_fermion(X,μ)
            Us = staggered_U(U[μ],μ)
            mul!(temp0_f,Us,Xplus)

            #U_{k,μ} X_{k+μ}) ⊗ Y_k
            mul!(temp2_g,temp0_f,Y') 
            #...  p(new) = p(old) +factor * temp2_g
            Traceless_antihermitian_add!(p[μ],-0.5*τ*mdparams.Δτ*coeff,temp2_g)



            #!  Construct P2*U_adj(x,mu)
            # Y_{k+μ}^dag U_{k,μ}^dag
            Yplus = shift_fermion(Y,μ)
            mul!(temp0_f,Yplus',Us')

            #X_k ⊗ Y_{k+μ}^dag U_{k,μ}^dag
            mul!(temp2_g,X,temp0_f) 

            mul!(temp3_g,U[μ]',temp2_g)
            #...  p(new) = p(old) + fac * c  .....
            Traceless_antihermitian_add!(p[μ],-0.5*τ*mdparams.Δτ*coeff,temp2_g)

        end
    end

    function  updateP_fermi_fromX_Wilson!(Y::F,φ::F,X::AbstractFermionfields{NC,Dim},fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,U::Array{T,1},
        temps::Array{T2,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2(),coeff=1
        ) where {NC,Dim,F <: AbstractFermionfields, T<: AbstractGaugefields,N<: TA_Gaugefields,T2 <: AbstractGaugefields} 
        temp0_f = temps_fermi[1] #F_field
        temp1_f = temps_fermi[2] #F_field
        temp2_g = temps[1] #G_field1
        temp3_g = temps[2] #G_field1
        c = temp_a[1]
        NV = temp2_g.NV

        W = Dirac_operator(U,φ,fparam)
        mul!(Y,W,X)
        set_wing_fermion!(Y)

        for μ=1:Dim
            #!  Construct U(x,mu)*P1

            # U_{k,μ} X_{k+μ}
            Xplus = shift_fermion(X,μ)
            mul!(temp0_f,U[μ],Xplus)
            
            
            # (r-γ_μ) U_{k,μ} X_{k+μ}
            mul!(temp1_f,view(X.rminusγ,:,:,μ),temp0_f)
            
            # κ (r-γ_μ) U_{k,μ} X_{k+μ}
            mul!(temp0_f,X.hopp[μ],temp1_f)

            # κ ((r-γ_μ) U_{k,μ} X_{k+μ}) ⊗ Y_k
            mul!(temp2_g,temp0_f,Y') 
            #vvmat!(temp2_g,temp1_f,Y,1)


            Traceless_antihermitian_add!(p[μ],τ*mdparams.Δτ,temp2_g)

            #!  Construct P2*U_adj(x,mu)
            # Y_{k+μ}^dag U_{k,μ}^dag
            Yplus = shift_fermion(Y,μ)
            mul!(temp0_f,Yplus',U[μ]')

            # Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)
            mul!(temp1_f,temp0_f,view(X.rplusγ,:,:,μ))
 
            # κ Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)
            mul!(temp0_f,X.hopm[μ],temp1_f)

            # X_k ⊗ κ Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)
            mul!(temp2_g,X,temp0_f) 

            Traceless_antihermitian_add!(p[μ],-τ*mdparams.Δτ,temp2_g)

            if typeof(fparam) == FermiActionParam_WilsonClover
                error("not implemented yet.")
                dSclover!(c,μ,X,Y,U,fparam,temps)
                add!(p[μ],-τ*mdparams.Δτ,c)
            end

        end

    end

    function  updateP_fermi_fromX_smearing!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,U::Array{T,1},Uout_multi,dSdU,Uin,
        temps::Array{T2,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2(),coeff=1
        ) where {F <: AbstractFermionfields, T<: AbstractGaugefields,N<: TA_Gaugefields,T2 <: AbstractGaugefields}
        
        if Y.Dirac_operator == "Staggered"
            updateP_fermi_fromX_smearing_Staggered!(Y,φ,X,fparam,
                p,mdparams,τ,U,Uout_multi,dSdU,Uin,
                temps,temp_a,temps_fermi;kind_of_verboselevel = kind_of_verboselevel,coeff=coeff)
        elseif findfirst("Wilson",Y.Dirac_operator) != nothing #Dirac_operator == "Wilson"
            updateP_fermi_fromX_smearing_Wilson!(Y,φ,X,fparam,
                p,mdparams,τ,U,Uout_multi,dSdU,Uin,
                temps,temp_a,temps_fermi;kind_of_verboselevel = kind_of_verboselevel,coeff=coeff)
        else
            error("Dirac_operator = $(Y.Dirac_operator) is not supported")
        end

    end

    function  updateP_fermi_fromX_smearing_Staggered!(Y::F,φ::F,X::AbstractFermionfields{NC,Dim},fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,U::Array{T,1},Uout_multi,dSdU,Uin,
        temps::Array{T2,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2(),coeff=1
        ) where {NC,Dim,F <: AbstractFermionfields, T<: AbstractGaugefields,N<: TA_Gaugefields,T2 <: AbstractGaugefields}
        temp0_f = temps_fermi[1] #F_field
        temp1_f = temps_fermi[2] #F_field
        temp2_g = temps[1] #G_field1
        temp3_g = temps[2] #G_field1
        c = temp_a[1]
        NV = temp2_g.NV

        W = Dirac_operator(U,φ,fparam)
        mul!(Y,W,X)

        for μ=1:Dim
            #!  Construct U(x,mu)*P1

            # U_{k,μ} X_{k+μ}
            Xplus = shift_fermion(X,μ)
            Us = staggered_U(U[μ],μ)
            mul!(temp0_f,Us,Xplus)
            println(temp0_f[1,1,1,1,1,1])

            #U_{k,μ} X_{k+μ}) ⊗ Y_k
            mul!(temp2_g,temp0_f,Y') 
            #println(temp2_g[1,1,1,1,1,1])
            mul!(dSdU[μ],U[μ]',temp2_g) #additional term

            #!  Construct P2*U_adj(x,mu)
            # Y_{k+μ}^dag U_{k,μ}^dag
            Yplus = shift_fermion(Y,μ)
            mul!(temp0_f,Yplus',Us')

            #X_k ⊗ Y_{k+μ}^dag U_{k,μ}^dag
            mul!(temp2_g,X,temp0_f) 

            mul!(temp3_g,U[μ]',temp2_g)
            add_U!(dSdU[μ],temp3_g)
        end

        #println(dSdU[1][1,1,1,1,1,1])

        if typeof(fparam.smearing) <: CovNeuralnet
            dSdUnew = back_prop(dSdU,fparam.smearing,Uout_multi,Uin)
        else
            error("$(typeof(fparam.smearing)) is not supported")
        end

        println("dSdUnew[1][1,1,1,1,1,1] = ", dSdUnew[1][1,1,1,1,1,1])
        println("dSdUnew[2][1,1,1,1,1,1] = ", dSdUnew[2][1,1,1,1,1,1])
        println("dSdUnew[4][1,1,1,1,1,1] = ", dSdUnew[4][:,:,1,1,1,1])
        
        #=

        if typeof(fparam.smearing) <: SmearingParam_single
            dSdUnew,_ = stoutfource(dSdU,Uin,fparam.smearing) 
        elseif typeof(fparam.smearing) <: SmearingParam_multi
            dSdUnew,_ = stoutfource(dSdU,Uout_multi,Uin,fparam.smearing) 
        else
            error("$(typeof(fparam.smearing)) is not supported")
        end
        =#

        println("p[1][1,1,1,1,1] = ",p[1][1,1,1,1,1])
        println("p[2][:,2,2,2,2] = ",p[2][:,2,2,2,2])
        println("-0.5*τ*mdparams.Δτ*coeff ", -0.5*τ*mdparams.Δτ*coeff)
        println("-0.5 $τ $(mdparams.Δτ) $coeff ")



        for μ=1:Dim
            #println("Uin[μ] = ",Uin[μ][:,:,1,1,1,1])
            mul!(temp2_g,Uin[μ],dSdUnew[μ])
            #println("temp2_g = ",temp2_g[:,:,1,1,1,1])

            Traceless_antihermitian_add!(p[μ],-0.5*τ*mdparams.Δτ*coeff,temp2_g)
        end

        println("p[1][1,1,1,1,1] = ",p[1][1,1,1,1,1])
        println("p[2][:,2,2,2,2] = ",p[2][:,2,2,2,2])

        error("stga")
        return


    end

    function  updateP_fermi_fromX_smearing_Wilson!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,U::Array{T,1},Uout_multi,dSdU,Uin,
        temps::Array{T2,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2(),coeff=1
        ) where {F <: AbstractFermionfields, T<: AbstractGaugefields,N<: TA_Gaugefields,T2 <: AbstractGaugefields}
        error("updateP_fermi_fromX_smearing_Wilson! is not implemented")
    end

    function updateP_fermi!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,Uin::Array{T,1},
        temps::Array{T_1d,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2()
        ) where {F <: WilsonFermion, T<: GaugeFields,N<: LieAlgebraFields,T_1d <: GaugeFields_1d} 
        temp0_f = temps_fermi[1] #F_field
        temp1_f = temps_fermi[2] #F_field
        temp2_g = temps[1] #G_field1
        temp3_g = temps[2] #G_field1
        c = temp_a[1]
        NV = temp2_g.NV

        U,Uout_multi,dSdU = calc_smearingU(Uin,fparam.smearing,calcdSdU = true,temps = temps)


        WdagW = DdagD_operator(U,φ,fparam)
        cg(X,WdagW,φ,eps = fparam.eps,maxsteps= fparam.MaxCGstep,verbose = kind_of_verboselevel)
        set_wing_fermi!(X)

        if fparam.smearing != nothing && typeof(fparam.smearing) != Nosmearing
            updateP_fermi_fromX_smearing!(Y,φ,X,fparam,
            p,mdparams,τ,U,Uout_multi,dSdU,Uin,
            temps,temp_a,temps_fermi,kind_of_verboselevel = kind_of_verboselevel)
        else
            

            updateP_fermi_fromX!(Y,φ,X,fparam,
            p,mdparams,τ,U,
            temps,temp_a,temps_fermi,kind_of_verboselevel = kind_of_verboselevel
            )
        end


       
        return

    end


    function updateP_fermi!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,Uin::Array{T,1},
        temps::Array{T_1d,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2()
        ) where {F <: StaggeredFermion, T<: GaugeFields,N<: LieAlgebraFields,T_1d <: GaugeFields_1d} 
        temp0_f = temps_fermi[1] #F_field
        temp1_f = temps_fermi[2] #F_field
        temp2_g = temps[1] #G_field1
        temp3_g = temps[2] #G_field1
        c = temp_a[1]
        NV = temp2_g.NV

        U,Uout_multi,dSdU = calc_smearingU(Uin,fparam.smearing,calcdSdU = true,temps = temps)

        WdagW = DdagD_operator(U,φ,fparam)

        if fparam.Nf == 4 || fparam.Nf == 8
            #X = (D^dag D)^(-1) ϕ 
            #
            cg(X,WdagW,φ,eps = fparam.eps,maxsteps= fparam.MaxCGstep,verbose = kind_of_verboselevel)
            set_wing_fermi!(X)

            if fparam.smearing != nothing && typeof(fparam.smearing) != Nosmearing
                updateP_fermi_fromX_smearing!(Y,φ,X,fparam,
                p,mdparams,τ,U,Uout_multi,dSdU,Uin,
                temps,temp_a,temps_fermi,kind_of_verboselevel = kind_of_verboselevel)
            else
                updateP_fermi_fromX!(Y,φ,X,fparam,
                p,mdparams,τ,U,
                temps,temp_a,temps_fermi,kind_of_verboselevel = kind_of_verboselevel)
            end
        else
            N_MD = get_order(fparam.rhmc_MD)
            #numtemps_fermi = length(temps_fermi)
            x = temps_fermi[end-N_MD]
            vec_x = temps_fermi[end-N_MD+1:end]
            for j=1:N_MD
                Fermionfields.clear!(vec_x[j])
            end
            vec_β = get_β_inverse(fparam.rhmc_MD)
            vec_α = get_α_inverse(fparam.rhmc_MD)
            α0 = get_α0_inverse(fparam.rhmc_MD)
            shiftedcg(vec_x,vec_β,x,WdagW,φ,eps = fparam.eps,maxsteps= fparam.MaxCGstep)
            for j=1:N_MD
                set_wing_fermi!(vec_x[j])
                if fparam.smearing != nothing && typeof(fparam.smearing) != Nosmearing
                    updateP_fermi_fromX_smearing!(Y,φ,vec_x[j],fparam,
                        p,mdparams,τ,U,Uout_multi,dSdU,Uin,
                        temps,temp_a,temps_fermi,kind_of_verboselevel = kind_of_verboselevel,coeff=vec_α[j])
                else
                    updateP_fermi_fromX!(Y,φ,vec_x[j],fparam,
                        p,mdparams,τ,U,
                        temps,temp_a,temps_fermi,kind_of_verboselevel = kind_of_verboselevel,coeff=vec_α[j])
                end
            end
        end

        #exit()
        return

    end


    function updateP!(univ::Universe,mdparams::MD_parameters,τ)
        factor = calc_factor(univ,mdparams)
        updateP!(univ.p,factor,τ,univ.U,univ._temporal_gauge,univ._temporal_algebra[1],univ.gparam)
        return
    end

    function updateP!(p::Array{N,1},factor,τ,U::Array{T,1},temps::Array{T_1d,1},temp_a::N,gparam::GaugeActionParam_autogenerator) where {T<: GaugeFields,N<: LieAlgebraFields,T_1d <: GaugeFields_1d} 
        add_gaugeforce!(p,U,temps,temp_a,gparam,fac = τ*factor)
        return
    end

    function updateP!(p::Array{N,1},factor,τ,U::Array{T,1},temps::Array{T_1d,1},temp_a::N,gparam::GaugeActionParam_standard) where {T<: GaugeFields,N<: LieAlgebraFields,T_1d <: GaugeFields_1d} 
        add_gaugeforce!(p,U,temps,temp_a,fac = τ*factor)
        
        return
    end

    function updateU!(univ::Universe,mdparams::MD_parameters,τ)
        updateU!(univ.U,mdparams,τ,univ.p,univ._temporal_gauge,univ._temporal_algebra[1])
        return
    end

    function updateU!(U,mdparams::MD_parameters,τ,p::Array{N,1},temps,temp_a::N) where {N<: LieAlgebraFields}  
        temp1 = temps[1]
        temp2 = temps[2]
        temp3 = temps[3]
        temp4 = temps[4]

        c = temp_a

        for μ=1:4
            

            substitute!(temp1,U[μ])

            mul!(c,τ*mdparams.Δτ,p[μ])

            expA!(temp2,c,temp3,temp4)


            #mul!(temp3,temp1,U[μ])
            mul!(temp3,temp2,temp1)

            substitute!(U[μ],temp3)

            set_wing!(U[μ])

        end
        
    end

    #function updateU!(U,mdparams::MD_parameters,τ,p::Array{N,1},temps::Array{T_1d,1},temp_a::N) where {N<: LieAlgebraFields,T_1d <: GaugeFields_1d}  
    function updateU!(U,mdparams::MD_parameters,τ,p::Array{N,1},temps::Array{T_1d,1},temp_a::N) where {N,T_1d}  
        temp1 = temps[1]
        temp2 = temps[2]
        temp3 = temps[3]
        temp4 = temps[4]

        c = temp_a

        for μ=1:4

            substitute!(temp1,U[μ])

            mul!(c,τ*mdparams.Δτ,p[μ])

            expA!(temp2,c,temp3,temp4)

            
            #mul!(temp3,temp1,U[μ])
            mul!(temp3,temp2,temp1)

            substitute!(U[μ],temp3)

            set_wing!(U[μ])


        end
        
    end

    function updateP!(p::Array{N,1},factor,τ,U::Array{T,1},temps::Array{T_1d,1},temp_a::N,gparam::GaugeActionParam_autogenerator) where {T<: AbstractGaugefields,N,T_1d} 
        add_force!(p,U,temps,plaqonly = false,factor =τ*factor,staplefactors = gparam.βs./gparam.β)
        #add_gaugeforce!(p,U,temps,temp_a,gparam,fac = τ*factor)
        return
    end

    function updateP!(p::Array{N,1},factor,τ,U::Array{T,1},temps::Array{T_1d,1},temp_a::N,gparam::GaugeActionParam_standard) where {T<: AbstractGaugefields,N,T_1d} 
        add_force!(p,U,temps,plaqonly = true,factor =τ*factor)
        #add_gaugeforce!(p,U,temps,temp_a,fac = τ*factor)
        #println(p[1][1,2,1,1,1])
        return
    end

    function updateU!(U::Array{<: AbstractGaugefields{NC,Dim}},mdparams::MD_parameters,τ,p::Array{N,1},temps::Array{T_1d,1},temp_a::N) where {N,T_1d,NC,Dim}  
        temp1 = temps[1]
        temp2 = temps[2]
        expU = temps[3]
        W = temps[4]

        for μ=1:Dim
            exptU!(expU,τ*mdparams.Δτ,p[μ],[temp1,temp2])
            mul!(W,expU,U[μ])
            substitute_U!(U[μ],W)
            set_wing_U!(U[μ])
        end

        #println("U ",U[1][1,1,2,1,1,1])
        
    end




end