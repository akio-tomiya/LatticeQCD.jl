module MD
    using LinearAlgebra
    using Base.Threads

    import ..Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3


    import ..Actions:GaugeActionParam_standard,
                        FermiActionParam_WilsonClover,FermiActionParam_Wilson,
                        GaugeActionParam_autogenerator
    import ..LTK_universe:Universe,calc_Action,gauss_distribution,make_WdagWmatrix,set_β!,
                calc_IntegratedFermionAction
    import ..Gaugefields:GaugeFields,substitute!,SU3GaugeFields,set_wing!,
                            projlink!,make_staple_double!,
                            GaugeFields_1d,SU3GaugeFields_1d,SU2GaugeFields_1d,calc_Plaq_notrace_1d,
                            calc_GaugeAction
    import ..LieAlgebrafields:gauss_distribution_lie!,LieAlgebraFields,SU3AlgebraFields,SU2AlgebraFields,
                                Gauge2Lie!,add!,add_gaugeforce!,expA!
    import ..Fermionfields:gauss_distribution_fermi!,set_wing_fermi!,Wdagx!,vvmat!,
            FermionFields,fermion_shift!,WilsonFermion, fermion_shiftB!,
            StaggeredFermion,Wx!,clear!,WdagWx!,substitute_fermion!
    using ..Fermionfields
    import ..Heatbath:heatbath!


    #import ..CGfermion:cg0!,cg0_WdagW!,shiftedcg0_WdagW!
    import ..Clover:Make_CloverFμν!,dSclover!
    import ..System_parameters:Params

    import ..Diracoperators:Wilson_operator,Adjoint_Wilson_operator,WilsonClover_operator,
                Dirac_operator,DdagD_operator
    import ..CGmethods:bicg,cg,shiftedcg
    import ..RationalApprox:calc_exactvalue,calc_Anϕ

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
        #if univ.fparam != nothing
            W = Dirac_operator(univ.U,univ.η,univ.fparam)
            bicg(univ.η,W',univ.φ,eps = univ.fparam.eps,maxsteps= univ.fparam.MaxCGstep,verbose = univ.kind_of_verboselevel)
            #cg0!(univ.η,univ.φ,2,univ.U,univ._temporal_gauge,univ._temporal_fermi,univ.fparam)
            #cg0!(univ.η,univ.φ,2,univ.U,univ._temporal_fermi,univ.fparam)
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
            W = Dirac_operator(univ.U,univ.η,univ.fparam)
            bicg(univ.η,W',univ.φ,eps = univ.fparam.eps,maxsteps= univ.fparam.MaxCGstep,verbose = univ.kind_of_verboselevel)
            #cg0!(univ.η,univ.φ,2,univ.U,univ._temporal_gauge,univ._temporal_fermi,univ.fparam)
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

        function debug_gauss()
            WdagW = make_WdagWmatrix(univ)
            
            #WdagWold = make_WdagWmatrix(univ,univ.Uold)
            #ldet = logdet(WdagW)
            #ldetold = logdet(WdagWold)
            #ldetsa= ldetold -ldet
            #println(ldetsa)
            #exit()

            N,_ = size(WdagW)
            n = 8

            e,v = eigen(WdagW)
            println(e)
            exit()
            en = e.^(-1/n)
            esa = 1 .- en
            #println(esa)
            #exit()


            x = deepcopy(univ.η)
            ϕ = zeros(ComplexF64,N)
            gauss_distribution_fermi!(x,univ.ranf)
            substitute_fermion!(ϕ,x)
            ϕ2 = calc_exactvalue(-n,WdagW,ϕ) 
            println(ϕ2'*ϕ)

            WdagW2 = DdagD_operator(univ.U,univ.η,univ.fparam)
            
            ϕ4 = calc_Anϕ(-n,WdagW2,x)   
            
            println(ϕ4*x)
            ldet = logdet(WdagW)/n
            detex = det(WdagW)
            #println(detex,"\t",detex^(1/n))
            println(ldet)
            

            itemax = 1000000
            detAprrox = 0
            gaussdata = zeros(ComplexF64,n,itemax)
            #gaussdata2 = zeros(ComplexF64,itemax)
            

            
            
            fp = open("detn.dat","w")
            for ite=1:itemax
            
                rexp = 0
                k = 1
                #for k=1:n
                    gauss_distribution_fermi!(x,univ.ranf)
                    substitute_fermion!(ϕ,x)
                    vx = v'*ϕ
                    vx2 = vx .* esa
                    #ed = Diagonal(en)
                    #wda = v*ed*v'
                    #display(sum(v*ed*v'- WdagW^(1/n)))
                    #exit()


                    
                    #ϕ2 = calc_exactvalue(-n,WdagW,ϕ) 
                    #println(ϕ'*ϕ2)
                    #println("\t")
                    #println(vx'*ed*vx)
                    #println(ϕ'*WdagW^(1/n)*ϕ)
                    #exit()
                    #ϕ3 = calc_exactvalue(n,WdagWold,ϕ2) 
                    
                    #ϕWϕ = ϕ'*ϕ - ϕ'*ϕ3
                    #ϕWϕ = ϕ'*ϕ - (univ.fparam.mass^2)^(1/n)*ϕ'*ϕ2
                    

                    #Wϕ = calc_Anϕ(-n,WdagW2,x)   
                    #ϕWϕ =  x*x - (univ.fparam.mass^2)^(1/n)*(x*Wϕ)
                    #ϕWϕ2 = ϕ'*ϕ - ϕ'*ϕ2
                    ϕWϕ =vx'*vx2
                    #println(ϕWϕ,"\t",ϕWϕ2)
                    #exit()

                    #println("xWϕ = ",x*Wϕ,"\t",x*x)

                    #println(ϕWϕ)
                    gaussdata[k,ite] = ϕWϕ 
                    #gaussdata2[ite]= x*x
                    ϕWϕmax = maximum(real.(gaussdata[k,1:ite]))
                    #ϕϕmax = maximum(real.(gaussdata2[1:ite]))

                    csum = sum(exp.(gaussdata[k,1:ite] .- ϕWϕmax))
                    csumlog = ϕWϕmax + log(csum) -log(ite) #+ N*log(2) #  + N*(2/n)*log(univ.fparam.mass) -N*log(pi)
                    #csumlog = ϕWϕmax + log(csum) #-log(ite) + N*(2/n)*log(univ.fparam.mass) -log(pi)
                    rexp += csumlog
                #end
                #csumlog = ϕWϕmax + log(csum) -log(ite) + N*(2/n)*log(univ.fparam.mass)
                #detAprrox += rexp
                
                if ite % 100 == 0
                    #println("$ite ",real(rexp),"\t",ldet*n,"\t",real(ldet*n-rexp))
                    #println(fp,ite,"\t",real(rexp),"\t",real(ldet*n))
                    println("$ite ",real(rexp),"\t",ldet,"\t",real(ldet-rexp))
                    println(fp,ite,"\t",real(rexp),"\t",real(ldet))
                end
                #println("$ite ",real(rexp),"\t",ldetsa,"\t",real(ldetsa-rexp))
            end
            close(fp)


            #e,_ = eigen(WdagW)
            #println(e)
            exit()
        end
        #debug_gauss()

        
        #@assert univ.NC == 2 "Only SU(2) is supported now."
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

        β = univ.gparam.β
        Sgold,plaq = calc_GaugeAction(univ)

        set_β!(univ,mdparams.βeff)
        Sgeffold,plaq = calc_GaugeAction(univ)

        heatbath!(univ)
        #heatbath!(univ)

        Sgeffnew,plaq = calc_GaugeAction(univ)
        set_β!(univ,β)

        Sgnew,plaq = calc_GaugeAction(univ)

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

            gauss_distribution_fermi!(univ.η,univ.ranf)
            #set_wing_fermi!(univ.η) 
            
            if univ.Dirac_operator == "Staggered" 
                if univ.fparam.Nf == 4

                    function shifttest()
                        println("Test for ShiftedCD")
                        println("(WdagW+ β_i) x_i = b")
                        WdagW = DdagD_operator(univ.U,univ.η,univ.fparam)
                        x = deepcopy(univ.η)
                        N = 10
                        vec_β = rand(N)
                        vec_β[1] = 0
                        vec_β[2] = 0.1
                        vec_β[3] = 0.5
                        vec_β[4] = -0.2
                        vec_x = Array{typeof(x),1}(undef,N)
                        for j=1:N
                            vec_x[j] = similar(x)
                        end
                        vec_β[1] = 0
                        shiftedcg(vec_x,vec_β,x,WdagW,univ.η,eps = univ.fparam.eps,maxsteps= univ.fparam.MaxCGstep)

                        #shiftedcg0_WdagW!(vec_x,x,univ.η,vec_β,univ.U,univ._temporal_gauge,univ._temporal_fermi,univ.fparam)  #(WdagW + vec_beta)*x = b
                        for j=1:N
                            println("β_$j = ",vec_β[j])
                            mul!(univ.φ,WdagW,vec_x[j])
                            #WdagWx!(univ.φ,univ.U,vec_x[j],univ._temporal_fermi,univ.fparam)

                            Fermionfields.add!(univ.φ,vec_β[j],vec_x[j])
                            #println("norm; ", univ.φ*univ.φ)
                            println("residual: ", univ.φ*univ.φ-univ.η*univ.η)
                            
                        end
                    end
                    #shifttest()

                    evensite = false
                    W = Dirac_operator(univ.U,univ.η,univ.fparam)
                    mul!(univ.φ,W',univ.η)

                    #Wdagx!(univ.φ,univ.U,univ.η,univ._temporal_fermi,univ.fparam)
                    clear!(univ.φ,evensite)

                    bicg(univ.η,W',univ.φ,eps = univ.fparam.eps,maxsteps= univ.fparam.MaxCGstep)

                    #cg0!(univ.η,univ.φ,2,univ.U,univ._temporal_gauge,univ._temporal_fermi,univ.fparam)

                    

                end
            end

            #gauss_distribution_fermi!(univ.η,univ.ranf)

            
            W = Dirac_operator(univ.U,univ.η,univ.fparam)

            #if typeof(univ.fparam) == FermiActionParam_WilsonClover
            #    Make_CloverFμν!(univ.fparam,univ.U,univ._temporal_gauge)
            #end
            mul!(univ.φ,W',univ.η)

            #Wdagx!(univ.φ,univ.U,univ.η,univ._temporal_fermi,univ.fparam) #φ = Wdag*η
            set_wing_fermi!(univ.φ)
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


    function updateP_fermi!(Y::F,φ::F,X::F,fparam,
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
        

        #=
        X = (D^dag D)^(-1) ϕ = D^(-1) (D^dag)^(-1) ϕ 
        Y = D X = (D^dag)^(-1) ϕ
        
        Y = D X = (D^dag)^(-1) ϕ
        Solve D^dag Y = ϕ
        =#
        bicg(Y,W',φ,eps = fparam.eps,maxsteps= fparam.MaxCGstep,verbose = kind_of_verboselevel)
        set_wing_fermi!(Y)

        
        bicg(X,W,Y,eps = fparam.eps,maxsteps= fparam.MaxCGstep,verbose = kind_of_verboselevel)
        set_wing_fermi!(X)


        


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

    function updateP_fermi!(Y::F,φ::F,X::F,fparam,
        p::Array{N,1},mdparams::MD_parameters,τ,U::Array{T,1},
        temps::Array{T_1d,1},temp_a::Array{N,1},temps_fermi;kind_of_verboselevel = Verbose_2()
        ) where {F <: StaggeredFermion, T<: GaugeFields,N<: LieAlgebraFields,T_1d <: GaugeFields_1d} 
        temp0_f = temps_fermi[1] #F_field
        temp1_f = temps_fermi[2] #F_field
        temp2_g = temps[1] #G_field1
        temp3_g = temps[2] #G_field1
        c = temp_a[1]
        NV = temp2_g.NV

        
        #X = (D^dag D)^(-1) ϕ 
        #

        WdagW = DdagD_operator(U,φ,fparam)
        cg(X,WdagW,φ,eps = fparam.eps,maxsteps= fparam.MaxCGstep,verbose = kind_of_verboselevel)
        set_wing_fermi!(X)

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
            add!(p[μ],-0.5*τ*mdparams.Δτ,c)



            #!  Construct P2*U_adj(x,mu)
            # Y_{k+μ}^dag U_{k,μ}^dag
            fermion_shiftB!(temp0_f,U,-μ,Y)
            # Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)

            # X_k ⊗ κ Y_{k+μ}^dag U_{k,μ}^dag*(r+γ_μ)
            vvmat!(temp2_g,X,temp0_f,2)


            #.....   Projection onto Lie Algebra   .....
            projlink!(temp3_g,temp2_g)
            Gauge2Lie!(c,temp3_g)

            add!(p[μ],-0.5*τ*mdparams.Δτ,c)


        end
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

    function updateU!(U,mdparams::MD_parameters,τ,p::Array{N,1},temps::Array{T_1d,1},temp_a::N) where {N<: LieAlgebraFields,T_1d <: GaugeFields_1d}  
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




end