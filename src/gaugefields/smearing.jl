module Smearing
    using LinearAlgebra, SpecialFunctions
    import ..LTK_universe:Universe,calc_gaugeforce!,expF_U!
    import ..LieAlgebrafields:LieAlgebraFields,clear!,add!
    import ..Gaugefields:GaugeFields,normalize!,SU,evaluate_wilson_loops!,TA!,set_wing!
    import ..Wilsonloops:Wilson_loop_set,calc_coordinate,make_plaq_staple_prime,calc_shift,make_plaq,make_plaq_staple
    import ..Actions:SmearingParam_single,SmearingParam_multi

    function add!(a::Array{N,1},α,b::Array{N,1}) where N <: LieAlgebraFields
        for μ=1:4
            add!(a[μ],α,b[μ]) # a_\mu+= \beta b_\mu
        end
    end

    function gradientflow!(U::Array{N,1},univ::Universe,tempW1,tempW2,Nflow::Int = 1,eps::Float64 = 0.01) where N <: GaugeFields
        # Here we use definition in 1006.4518 except for the hermiticity.
        # ref https://www2.physik.uni-bielefeld.de/fileadmin/user_upload/theory_e6/Master_Theses/Masterarbeit_LukasMazur.pdf
        NC = univ.NC
        NX = univ.NX
        NY = univ.NY
        NZ = univ.NZ
        NT = univ.NT

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

        #=
         W0 = Vt
         Z0 = Z(W0)
         W1 = exp(1/4 Z0)*W0
         Z1 = Z(W1)
         W2 = exp(8/9 * Z1 - 17/36 * Z0)*W1
         Z2 = Z(W2)
         W3 = exp(3/4* Z2 - 8/9 * Z1 +17/36 *Z0)*W2
         Vt = W3


        =#

        for istep=1:Nflow #RK4 integrator
            calc_gaugeforce!(F0,U,univ) #F
            #println(F0[1][1,1,1,1,1])
            #expF_U!(W1,F0,-eps*(1/4),univ)  #exp(eps*F)*U 
            expF_U!(U,F0,-eps*(1/4),univ)  #exp(eps*F)*U 
            #println("W1 ", U[1][1,1,1,1,1,1])
            calc_gaugeforce!(F1,U,univ) #F
            #calc_gaugeforce!(F1,W1,univ) #F
            #println("F1 ",F1[1][1,1,1,1,1,1])
            clear!(Ftmp)
            add!(Ftmp,-(8/9*eps),F1)
            #println("Ftmp ",Ftmp[1][1,1,1,1,1,1])
            add!(Ftmp,(17/36*eps),F0)
            #println("Ftmp1 ",Ftmp[1][1,1,1,1,1,1])
            expF_U!(U,Ftmp,1,univ)

            #expF_U!(W2,Ftmp,1,univ)
            #println("W2 ",W2[1][1,1,1,1,1,1])
            #
            calc_gaugeforce!(F2,U,univ) #F
            #calc_gaugeforce!(F2,W2,univ) #F
            clear!(Ftmp)
            add!(Ftmp,-(3/4*eps),F2)
            add!(Ftmp,(8/9*eps),F1)
            add!(Ftmp,-(17/36*eps),F0)
            expF_U!(U,Ftmp,1,univ)

            #println(U[1][1,1,1,1,1,1])
            #error("U")
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

    function calc_stout!(Uout::Array{GaugeFields{SU{NC}},1},ρ) where NC
        Uin = deepcopy(Uout)
        calc_stout!(Uout,Uin,ρ)
    end

    function calc_stout!(Uout::Array{GaugeFields{SU{NC}},1},U::Array{GaugeFields{SU{NC}},1},ρ) where NC
        @assert Uout != U "input U and output U should not be same!"

        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        Q = zeros(ComplexF64,NC,NC)
        C = zeros(ComplexF64,NC,NC)

        for μ=1:4
            loops = make_plaq_staple_prime(μ)
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            C .= 0
                            evaluate_wilson_loops!(C,loops,U,ix,iy,iz,it)
                            #display(C)
                            #println("\t")

                            Ω = ρ*C[:,:]*U[μ][:,:,ix,iy,iz,it]'
                            Q = (im/2)*(Ω' .- Ω) .- (im/(2NC))*tr(Ω' .- Ω)
                            Uout[μ][:,:,ix,iy,iz,it] = exp(im*Q)*U[μ][:,:,ix,iy,iz,it]
                            set_wing!(Uout[μ],ix,iy,iz,it)

                            
                        end
                    end
                end
            end
        end
        #display(Uout[1][:,:,1,1,1,1])
        return
    end

    function calc_multihit!(Uout::Array{GaugeFields{SU{NC}},1},U::Array{GaugeFields{SU{NC}},1},β) where NC
        @assert Uout != U "input U and output U should not be same!"
        @assert U[1].NC == 2 "NC=2 is required for now but $(U[1].NC)"

        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        Q = zeros(ComplexF64,NC,NC)
        A = zeros(ComplexF64,NC,NC)

        for μ=1:4
            loops = make_plaq_staple_prime(μ)
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            A .= 0
                            evaluate_wilson_loops!(A,loops,U,ix,iy,iz,it)
                            K = sqrt(abs(det(A)))

                            #display(C)
                            #println("\t")

                            #Ω = ρ*C[:,:]*U[μ][:,:,ix,iy,iz,it]'
                            #Q = (im/2)*(Ω' .- Ω) .- (im/(2NC))*tr(Ω' .- Ω)
                            Uout[μ][:,:,ix,iy,iz,it] = K*inv(A)*besseli(2, β*K)/besseli(1, β*K)
                            #0  #exp(im*Q)*U[μ][:,:,ix,iy,iz,it]
                            set_wing!(Uout[μ],ix,iy,iz,it)

                            
                        end
                    end
                end
            end
        end
        #display(Uout[1][:,:,1,1,1,1])
        return
    end

    function calc_fatlink_APE!(Uout::Array{GaugeFields{SU{NC}},1},α,β;normalize_method= "special unitary",temporal_dir_smear=true) where NC
        Uin = deepcopy(Uout)
        calc_fatlink_APE!(Uout,Uin,α,β,normalize_method=normalize_method,temporal_dir_smear=temporal_dir_smear)
    end

    function calc_fatlink_APE!(Uout::Array{GaugeFields{SU{NC}},1},U::Array{GaugeFields{SU{NC}},1},α,β;normalize_method= "special unitary",temporal_dir_smear=true) where NC
        @assert Uout != U "input U and output U should not be same!"
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        Vtemp = zeros(ComplexF64,NC,NC)
        nmethod(A) = if normalize_method == "unitary"
                        A*(A'*A)^(-1/2)#inv(sqrt(A'*A))
                    elseif normalize_method == "special unitary"
                        normalize!(A)
                        A
                    else
                        error("normalize_method should be unitary or special unitary. Now, $normalize_method")
                    end
                    

        #println(normalize_method)
        d=3
        if temporal_dir_smear
            d=4
        end
        for μ=1:d
            loops = make_plaq_staple_prime(μ)
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            Vtemp .= 0
                            evaluate_wilson_loops!(Vtemp,loops,U,ix,iy,iz,it)
                            #Vtemp = Vtemp'
                            Vtemp[:,:] = (1-α)*U[μ][:,:,ix,iy,iz,it] .+ (β/6)*Vtemp[:,:]
                            #nmethod!(Vtemp)
                            Uout[μ][:,:,ix,iy,iz,it] = nmethod(Vtemp)#Vtemp[:,:]
                            set_wing!(Uout[μ],ix,iy,iz,it)
                        end
                    end
                end
            end
        end
        return 
    end

    function calc_fatlink_APE(U::Array{GaugeFields{SU{NC}},1},α,β;normalize_method= "special unitary",temporal_dir_smear=true) where NC
        Uout = similar(U)
        calc_fatlink_APE!(Uout,U,α,β,normalize_method=normalize_method,temporal_dir_smear=temporal_dir_smear)
        return Uout
    end


    function calc_multihit!(Uout::Array{GaugeFields{SU{NC}},1},U::Array{GaugeFields{SU{NC}},1},β) where NC
        @assert Uout != U "input U and output U should not be same!"
        @assert U[1].NC == 2 "NC=2 is required for now but $(U[1].NC)"

        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        Q = zeros(ComplexF64,NC,NC)
        A = zeros(ComplexF64,NC,NC)

        for μ=1:4
            loops = make_plaq_staple_prime(μ)
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            A .= 0
                            evaluate_wilson_loops!(A,loops,U,ix,iy,iz,it)
                            K = sqrt(abs(det(A)))

                            #display(C)
                            #println("\t")

                            #Ω = ρ*C[:,:]*U[μ][:,:,ix,iy,iz,it]'
                            #Q = (im/2)*(Ω' .- Ω) .- (im/(2NC))*tr(Ω' .- Ω)
                            Uout[μ][:,:,ix,iy,iz,it] = K*inv(A)*besseli(2, β*K)/besseli(1, β*K)
                            #0  #exp(im*Q)*U[μ][:,:,ix,iy,iz,it]
                            set_wing!(Uout[μ],ix,iy,iz,it)

                            
                        end
                    end
                end
            end
        end
        #display(Uout[1][:,:,1,1,1,1])
        return
    end
end