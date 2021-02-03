module Smearing
    import ..LTK_universe:Universe,calc_gaugeforce!,expF_U!
    import ..LieAlgebrafields:LieAlgebraFields,clear!,add!
    import ..Gaugefields:GaugeFields,normalize!,SU,evaluate_wilson_loops!,TA!
    import ..Wilsonloops:Wilson_loop_set,calc_coordinate,make_plaq_staple_prime,calc_shift,make_plaq,make_plaq_staple

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

    function calc_stout(U::Array{GaugeFields{SU{NC}},1},ρ) where NC
        Uout = similar(U)
        calc_stout!(Uout,U,ρ)
        return Uout
    end

    function calc_stout!(Uout::Array{GaugeFields{SU{NC}},1},U::Array{GaugeFields{SU{NC}},1},ρ) where NC
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        Q = zeros(ComplexF64,NC,NC)
        Vtemp = zeros(ComplexF64,NC,NC)

        for μ=1:4
            loops = make_plaq_staple(μ)
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            Vtemp .= 0
                            evaluate_wilson_loops!(Vtemp,loops,U,ix,iy,iz,it)
                            TA!(Q,(ρ/im)*Vtemp[:,:]*U[μ][:,:,ix,iy,iz,it]')
                            
                            Uout[μ][:,:,ix,iy,iz,it] = exp(im*Q)*U[μ][:,:,ix,iy,iz,it]
                        end
                    end
                end
            end
        end
        return
    end

    function calc_fatlink_APE!(Uout::Array{GaugeFields{SU{NC}},1},U::Array{GaugeFields{SU{NC}},1},α,β;normalize_method= "special unitary") where NC
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        Vtemp = zeros(ComplexF64,NC,NC)
        if normalize_method == "unitary"
            nmethod(A) = A*inv(sqrt(A'*A))
        elseif normalize_method == "special unitary"
            nmethod(A) = normalize!(A)
        end

        for μ=1:4
            loops = make_plaq_staple(μ)
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            Vtemp .= 0
                            evaluate_wilson_loops!(Vtemp,loops,U,ix,iy,iz,it)
                            Vtemp[:,:] = (1-α)*U[μ][:,:,ix,iy,iz,it] .+ (β/6)*Vtemp[:,:]
                            nmethod(Vtemp)
                            Uout[μ][:,:,ix,iy,iz,it] = Vtemp[:,:]
                        end
                    end
                end
            end
        end
        return 
    end

    function calc_fatlink_APE(U::Array{GaugeFields{SU{NC}},1},α,β;normalize_method= "special unitary") where NC
        Uout = similar(U)
        calc_fatlink_APE!(Uout,U,α,β,normalize_method=normalize_method)
        return Uout
    end
end