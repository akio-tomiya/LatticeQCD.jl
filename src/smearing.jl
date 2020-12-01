module Smearing
    import ..LTK_universe:Universe,calc_gaugeforce!,expF_U!
    import ..LieAlgebrafields:LieAlgebraFields,clear!,add!
    import ..Gaugefields:GaugeFields

    function add!(a::Array{N,1},α,b::Array{N,1}) where N <: LieAlgebraFields
        for μ=1:4
            add!(a[μ],α,b[μ]) # a_\mu+= \beta b_\mu
        end
    end

    function gradientflow!(U::Array{N,1},univ::Universe,tempW1,tempW2,Nflow::Int = 1,eps::Float64 = 0.01) where N <: GaugeFields
        #Here we use definition in 1006.4518 except for the hermiticity.
        NC = univ.NC
        NX = univ.NX
        NY = univ.NY
        NZ = univ.NZ
        NT = univ.NT

        g2 = (2*NC)/univ.gparam.β # β = 2Nc/g^2

        #W0 = deepcopy(U)
        W1 = tempW1
        W2 = tempW2

        F0 = Array{LieAlgebraFields,1}(undef,4)
        for μ=1:4
            F0[μ] = LieAlgebraFields(NC,NX,NY,NZ,NT)
        end
        F1 = Array{LieAlgebraFields,1}(undef,4)
        for μ=1:4
            F1[μ] = LieAlgebraFields(NC,NX,NY,NZ,NT)
        end
        F2 = Array{LieAlgebraFields,1}(undef,4)
        for μ=1:4
            F2[μ] = LieAlgebraFields(NC,NX,NY,NZ,NT)
        end
        Ftmp = Array{LieAlgebraFields,1}(undef,4)
        for μ=1:4
            Ftmp[μ] = LieAlgebraFields(NC,NX,NY,NZ,NT)
        end

        for istep=1:Nflow #RK4 integrator
            calc_gaugeforce!(F0,U,univ) #F
            expF_U!(W1,F0,-eps*(1/4)*g2,univ)  #exp(eps*F)*U 
            #
            calc_gaugeforce!(F1,W1,univ) #F
            clear!(Ftmp)
            add!(Ftmp,-(8/9*eps)*g2,F1)
            add!(Ftmp,(17/36*eps)*g2,F0)
            expF_U!(W2,Ftmp,1,univ)
            #
            calc_gaugeforce!(F2,W2,univ) #F
            clear!(Ftmp)
            add!(Ftmp,-(3/4*eps)*g2,F2)
            add!(Ftmp,(8/9*eps)*g2,F1)
            add!(Ftmp,-(17/36*eps)*g2,F0)
            expF_U!(U,Ftmp,1,univ)
            #
        end
        #return W1 #test
        #return W0
    end

    function gradientflow!(U::Array{N,1},univ::Universe,Nflow::Int = 1,eps::Float64 = 0.01) where N <: GaugeFields
        W1 = deepcopy(U)
        W2 = deepcopy(U)
        gradientflow!(U,univ,W1,W2,Nflow,eps)
        return 
    end
    #=
    function gradientflow_euler!(univ::Universe) #backup
        eps = 0.01
        Nflow = 1
        NC = univ.NC
        NX = univ.NX
        NY = univ.NY
        NZ = univ.NZ
        NT = univ.NT

        F = Array{LieAlgebraFields,1}(undef,4)
        for μ=1:4
            F[μ] = LieAlgebraFields(NC,NX,NY,NZ,NT)
        end

        U = similar(univ.U)

        for istep=1:Nflow #Euler integrator
            calc_gaugeforce!(F,univ) #F
            expF_U!(U,F,eps,univ)  #exp(eps*F)*U
        end
    end
    =#
end