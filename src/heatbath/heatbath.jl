module Heatbath
    using LinearAlgebra
    import ..LTK_universe:Universe
    import ..Gaugefields:GaugeFields,SU3GaugeFields,
            SU2GaugeFields,SU3GaugeFields_1d,SU2GaugeFields_1d,
            GaugeFields_1d,elementwise_tr!,set_wing!,make_staple_double!,substitute!,clear!,
            evaluate_wilson_loops!,normalize!,normalize3!,normalize2!,
            SUNGaugeFields,SUNGaugeFields_1d,normalizeN!
            #,
            #set_wing_x!,set_wing_y!,set_wing_z!,set_wing_t!
    import ..Wilsonloops:Wilson_loop,Wilson_loop_set,make_plaq_staple,make_links,make_staples,make_plaq
    import ..Actions:GaugeActionParam_standard,GaugeActionParam_autogenerator,GaugeActionParam

    function heatbath!(univ::Universe)
        heatbath!(univ.U,univ.ranf,univ.gparam,univ._temporal_gauge)
        #heatbath!(univ.U,univ.ranf,univ.gparam.β,univ._temporal_gauge)
    end

    function heatbath!(u::Array{T,1},ranf,gparam::GaugeActionParam,temps::Array{T_1d,1}) where {T <: SU2GaugeFields,T_1d <: SU2GaugeFields_1d}
        beta = gparam.β

    
        staple= temps[1]
        temp1= temps[2]
        temp2 = temps[3]
        temp3 = temps[4]



        NV = staple.NV
        ITERATION_MAX = 10^5
        
        Wnew = zeros(ComplexF64,2,2)
        NX = u[1].NX
        NY = u[1].NY
        NZ = u[1].NZ
        NT = u[1].NT

        NC = 2
        V = zeros(ComplexF64,NC,NC)
        Vtemp = zeros(ComplexF64,NC,NC)

        a = zeros(Float64,4)

        for mu=1:4
            #make_staple_double!(staple,u,mu,temp1,temp2,temp3)
            if typeof(gparam) == GaugeActionParam_standard
                loops = make_plaq_staple(mu)
            end


            
            
            #for i=1:NV
            i = 0
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix = 1:NX
                            i += 1 

                            if typeof(gparam) == GaugeActionParam_standard
                                evaluate_wilson_loops!(V,loops,u,ix,iy,iz,it)
                            elseif typeof(gparam) == GaugeActionParam_autogenerator
                                V .= 0
                                Vtemp .= 0
                                for iloop = 1:gparam.numactions
                                    evaluate_wilson_loops!(Vtemp,gparam.staples[iloop][mu],u,ix,iy,iz,it)
                                    @. V += (gparam.βs[iloop]/beta)*Vtemp
                                end
                            end
                            
                            #u[mu][:,:,ix,iy,iz,it] = SU2update(V,beta,NC,ITERATION_MAX)
                            u[mu][:,:,ix,iy,iz,it] = SU2update_KP(V,beta,NC,ITERATION_MAX)

                            set_wing!(u[mu],ix,iy,iz,it)

                        end
                    end
                end
            end
            #normalize!(u[mu])
        end

    end

    function SU2update_KP(V,beta,NC,ITERATION_MAX = 10^5)
        eps = 0.000000000001

        #R = real(sqrt(det(V)))
        #V0 = inv(V/R)
        #w = zeros(Float64,4)

        ρ0 = real(V[1,1]+V[2,2])/2
        ρ1 = -imag(V[1,2]+V[2,1])/2
        #ρ1 = imag(V[1,2]+V[2,1])/2
        ρ2 = real(V[2,1]-V[1,2])/2
        ρ3 = imag(V[2,2]-V[1,1])/2
        ρ = sqrt(ρ0^2+ρ1^2+ρ2^2+ρ3^2)
        #println("R = ",R," ρ ",ρ)
        #println("detV = , ", det(V)," ",ρ0^2+ρ1^2+ρ2^2+ρ3^2)
        V0 = inv(V/ρ)

        #
        #Nc = 2 # Since Ishikawa's book uses 1/g^2 notation.
        #k = (beta/NC)*ρ
        k = 2*(beta/NC)*ρ

        
        #k = (beta/2)*ρ

        R = rand() + eps
        Rp = rand() + eps
        X = -log(R)/k
        Xp = -log(Rp)/k
        Rpp = rand()
        C = cos(2pi*Rpp)^2
        A = X*C
        delta = Xp + A
        Rppp = rand()

        a = zeros(Float64,4)
        while(Rppp^2 > 1-0.5*delta)
            R = rand()
            Rp = rand()
            X = -log(R)/k
            Xp = -log(Rp)/k
            Rpp = rand()
            C = cos(2pi*Rpp)^2
            A = X*C
            delta = Xp + A
            Rppp = rand()
        end
        a[1] = 1-delta


        rr = sqrt(1.0-a[1]^2)
        ϕ = rand()*pi*2.0 # ϕ = [0,2pi]
        cosθ = (rand()-0.5)*2.0 # -1<cosθ<1
        sinθ = sqrt(1-cosθ^2)

        a[2]=rr*cos(ϕ)*sinθ
        a[3]=rr*sin(ϕ)*sinθ
        a[4]=rr*cosθ
        Unew = [a[1]+im*a[4] a[3]+im*a[2]
                -a[3]+im*a[2] a[1]-im*a[4]]*V0

        #=
        w[1]=(1/ρ)*( a[1]*ρ0+a[2]*ρ1+a[3]*ρ2+a[4]*ρ3)
        w[2]=(1/ρ)*(-a[1]*ρ1+a[2]*ρ0+a[3]*ρ3-a[4]*ρ2)
        w[3]=(1/ρ)*(-a[1]*ρ2-a[2]*ρ3+a[3]*ρ0+a[4]*ρ1)
        w[4]=(1/ρ)*(-a[1]*ρ3+a[2]*ρ2-a[3]*ρ1+a[4]*ρ0)

        Unew = [w[1]+im*w[4] w[3]+im*w[2]
                -w[3]+im*w[2] w[1]-im*w[4]]
        =#

        #normalize2!(Unew)
        #display(Unew)

        #α = Unew[1,1]
        #β = Unew[2,1]
        

        α = Unew[1,1]*0.5 + conj(Unew[2,2])*0.5
        β = Unew[2,1]*0.5 - conj(Unew[1,2])*0.5

        detU = abs(α)^2 + abs(β)^2
        Unew[1,1] = α/detU
        Unew[2,1]  = β/detU
        Unew[1,2] = -conj(β)/detU
        Unew[2,2] = conj(α)/detU     
        
        

        return Unew
    end


    function SU2update(V,beta,NC,ITERATION_MAX = 10^5)
        R = real(sqrt(det(V)))
        V0 = inv(V/R)

        ρ0 = real(V[1,1]+V[2,2])
        ρ1 = -imag(V[1,2]+V[2,1])
        ρ2 = real(V[2,1]-V[1,2])
        ρ3 = imag(V[2,2]-V[1,1])
        ρ = sqrt(ρ0^2+ρ1^2+ρ2^2+ρ3^2)

        #
        #Nc = 2 # Since Ishikawa's book uses 1/g^2 notation.
        k = (beta/NC)*ρ

        #k = 2beta*R


        #A = 2*sinh(k)
        emk = exp(-k)
        epk = exp(k)
        ur = 999.0
        i_count=0
        a = zeros(Float64,4)
        while(ur^2 > 1.0-a[1]^2) # rejection sampling
            s = rand()
            a[1] = log(s*epk + (1-s)*emk)/k # F.17
            #a[1] =log(A*y+B)/(k)
            ur = rand()
            i_count+=1
            if i_count> ITERATION_MAX
                error("The rejection sampling is failed after $ITERATION_MAX trials.")
            end
        end

        rr = sqrt(1.0-a[1]^2)
        ϕ = rand()*pi*2.0 # ϕ = [0,2pi]
        cosθ = (rand()-0.5)*2.0 # -1<cosθ<1
        sinθ = sqrt(1-cosθ^2)

        a[2]=rr*cos(ϕ)*sinθ
        a[3]=rr*sin(ϕ)*sinθ
        a[4]=rr*cosθ
        Unew = [a[1]+im*a[4] a[3]+im*a[2]
                -a[3]+im*a[2] a[1]-im*a[4]]*V0
        #normalize2!(Unew)
        #display(Unew)

        #α = Unew[1,1]
        #β = Unew[2,1]
        α = Unew[1,1]*0.5 + conj(Unew[2,2])*0.5
        β = Unew[2,1]*0.5 - conj(Unew[1,2])*0.5

        detU = abs(α)^2 + abs(β)^2
        Unew[1,1] = α/detU
        Unew[2,1]  = β/detU
        Unew[1,2] = -conj(β)/detU
        Unew[2,2] = conj(α)/detU                            

        return Unew
    end

    function heatbath!(u::Array{T,1},ranf,gparam,temps::Array{T_1d,1}) where {T <: SU3GaugeFields,T_1d <: SU3GaugeFields_1d}
        #println("Warning!!!!!!!!")
        #error("Heatbath update for SU(3) is not implemented")
        beta = gparam.β

    
        staple= temps[1]
        temp1= temps[2]
        temp2 = temps[3]
        temp3 = temps[4]



        NV = staple.NV
        ITERATION_MAX = 10000
        Wnew = zeros(ComplexF64,2,2)
        NX = u[1].NX
        NY = u[1].NY
        NZ = u[1].NZ
        NT = u[1].NT

        NC = 3
        V = zeros(ComplexF64,NC,NC)
        Vtemp = zeros(ComplexF64,NC,NC)

        a = zeros(Float64,4)

        #=
        for mu=1:4
            normalize!(u[mu])
            set_wing!(u[mu])
        end
        =#

        for mu=1:4
            #make_staple_double!(staple,u,mu,temp1,temp2,temp3)
            if typeof(gparam) == GaugeActionParam_standard
                loops = make_plaq_staple(mu)
            end
            


            
            
            #for i=1:NV
            i = 0
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix = 1:NX
                            i += 1 

                            if typeof(gparam) == GaugeActionParam_standard
                                V .= 0
                                evaluate_wilson_loops!(V,loops,u,ix,iy,iz,it)
                            elseif typeof(gparam) == GaugeActionParam_autogenerator
                                V .= 0
                                Vtemp .= 0
                                for iloop = 1:gparam.numactions
                                    evaluate_wilson_loops!(Vtemp,gparam.staples[iloop][mu],u,ix,iy,iz,it)
                                    @. V += (gparam.βs[iloop]/beta)*Vtemp
                                end
                            end

                            #println("#Heatbath for one SU(3) link started")
                            for l=1:3

                                UV = u[mu][:,:,ix,iy,iz,it]*V

                                if l==1
                                    n,m = 1,2
                                elseif l==2
                                    n,m = 2,3
                                else
                                    n,m = 1,3

                                end

                                S = make_submatrix(UV,n,m)
                                #gramschmidt_special!(S)
                                project_onto_SU2!(S)

                                K = SU2update_KP(S,beta,NC,ITERATION_MAX)


                                A = make_largematrix(K,n,m,NC)

                                AU = A*u[mu][:,:,ix,iy,iz,it]

                                u[mu][:,:,ix,iy,iz,it] = AU
                            end

                            AU = u[mu][:,:,ix,iy,iz,it]
                            normalize3!(AU)
                            u[mu][:,:,ix,iy,iz,it] = AU

                            set_wing!(u[mu],ix,iy,iz,it)

                        end
                    end
                    
                end
            end

            normalize!(u[mu])
            
        end

    end


    function heatbath!(u::Array{T,1},ranf,gparam,temps::Array{T_1d,1}) where {T <: SUNGaugeFields,T_1d <: SUNGaugeFields_1d}
        #println("Warning!!!!!!!!")
        #error("Heatbath update for SU(3) is not implemented")
        beta = gparam.β

    
        staple= temps[1]
        temp1= temps[2]
        temp2 = temps[3]
        temp3 = temps[4]



        NV = staple.NV
        ITERATION_MAX = 10000
        Wnew = zeros(ComplexF64,2,2)
        NX = u[1].NX
        NY = u[1].NY
        NZ = u[1].NZ
        NT = u[1].NT

        NC = u[1].NC
        V = zeros(ComplexF64,NC,NC)
        Vtemp = zeros(ComplexF64,NC,NC)

        a = zeros(Float64,4)


        for mu=1:4
            #make_staple_double!(staple,u,mu,temp1,temp2,temp3)
            if typeof(gparam) == GaugeActionParam_standard
                loops = make_plaq_staple(mu)
            end

            
            #for i=1:NV
            i = 0
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix = 1:NX
                            i += 1 

                            if typeof(gparam) == GaugeActionParam_standard
                                V .= 0
                                evaluate_wilson_loops!(V,loops,u,ix,iy,iz,it)
                            elseif typeof(gparam) == GaugeActionParam_autogenerator
                                V .= 0
                                Vtemp .= 0
                                for iloop = 1:gparam.numactions
                                    evaluate_wilson_loops!(Vtemp,gparam.staples[iloop][mu],u,ix,iy,iz,it)
                                    @. V += (gparam.βs[iloop]/beta)*Vtemp
                                end
                            end

                            #println("#Heatbath for one SU(3) link started")
                            for l=1:NC
                            #for l=1:2NC
                            

                                UV = u[mu][:,:,ix,iy,iz,it]*V

                                
                                n = rand(1:NC)#l
                                m = rand(1:NC)
                                while(n==m)
                                    m = rand(1:NC)
                                end
                                
                                #=
                                if l < NC
                                    n = l
                                    m = l+1
                                else
                                    n = rand(1:NC)#l
                                    m = rand(1:NC)
                                    while(n==m)
                                        m = rand(1:NC)
                                    end
                                end
                                =#



                                S = make_submatrix(UV,n,m)
                                #gramschmidt_special!(S)
                                project_onto_SU2!(S)

                                K = SU2update_KP(S,beta,NC,ITERATION_MAX)


                                A = make_largematrix(K,n,m,NC)

                                AU = A*u[mu][:,:,ix,iy,iz,it]

                                u[mu][:,:,ix,iy,iz,it] = AU
                                #println("det U ",det(AU))

                            end

                            AU = u[mu][:,:,ix,iy,iz,it]
                            normalizeN!(AU)
                            u[mu][:,:,ix,iy,iz,it] = AU

                            #exit()
                            #set_wing!(u[mu])
                            set_wing!(u[mu],ix,iy,iz,it)

                        end

                    end
                    #set_wing_y!(u[mu],iz,it)
                    
                end
                #set_wing_z!(u[mu],it)
            end
            #set_wing_t!(u[mu])

            normalize!(u[mu])
            
        end

    end

    function gramschmidt!(v)
        n = size(v)[1]
        for i=1:n
            for j=1:i-1
                v[:,i] = v[:,i] - v[:,j]'*v[:,i]*v[:,j]
            end
            v[:,i] = v[:,i]/norm(v[:,i])
        end
    end


    function gramschmidt_special!(v)
        n = size(v)[1]
        #vdet = det(v)


        vnorm1 = norm(v[:,1])
        for i=1:n
            #vnorm[i] = norm(v[:,i])
            for j=1:i-1
                v[:,i] = v[:,i] - v[:,j]'*v[:,i]*v[:,j]
            end
            v[:,i] = v[:,i]/norm(v[:,i])
        end
        for i=1:n
            #v[:,i] = v[:,i]*vnorm[i]
            v[:,i] = v[:,i]*vnorm1
        end
    end

    function project_onto_SU2!(S) # This project onto SU(2) up to normalization.
        #S2 = zeros(ComplexF64,2,2)
        α = S[1,1]*0.5 + conj(S[2,2])*0.5
        β = S[2,1]*0.5 - conj(S[1,2])*0.5
        S[1,1] = α
        S[2,1] = β
        S[1,2] = -conj(β)
        S[2,2] = conj(α)
        #return S2
    end

    function make_submatrix(UV,i,j)
        S = zeros(ComplexF64,2,2)
        S[1,1] = UV[i,i]
        S[1,2] = UV[i,j]
        S[2,1] = UV[j,i]
        S[2,2] = UV[j,j]
        return S
    end


    function make_largematrix(K,i,j,NC)
        A = zeros(ComplexF64,NC,NC)
        for n=1:NC
            A[n,n] = 1
        end
        #K = project_onto_su2(K)
        A[i,i] = K[1,1]
        A[i,j] = K[1,2] 
        A[j,i] = K[2,1]
        A[j,j] = K[2,2]  
        return A
    end

    const nhit = 6
    const rwidth = 0.4


    """
-------------------------------------------------c
     su2-submatrix(c) in su3 matrix(x)
            su2            su3
     k=1         <-    1-2 elements
     k=2         <-    2-3 elements
     k=3         <-    1-3 elements
     k=4          ->   1-2 elements
     k=5          ->   2-3 elements
     k=6          ->   1-3 elements
-------------------------------------------------c
    """
    function submat!(x,c,n,k,id)

        if k==1
            for i=1:n
                c[1,i] = real(x[1,1,i]+x[2,2,i])*0.5
                c[2,i] = imag(x[1,2,i]+x[2,1,i])*0.5
                c[3,i] = real(x[1,2,i]-x[2,1,i])*0.5
                c[4,i] = imag(x[1,1,i]-x[2,2,i])*0.5
            end
        elseif k==2
            for i=1:n
                c[1,i] = real(x[2,2,i]+x[3,3,i])*0.5
                c[2,i] = imag(x[3,2,i]+x[2,3,i])*0.5
                c[3,i] = real(x[3,2,i]-x[2,3,i])*0.5
                c[4,i] = imag(x[2,2,i]-x[3,3,i])*0.5
            end

        elseif k==3
            for i=1:n
                c[1,i] = real(x[1,1,i]+x[3,3,i])*0.5
                c[2,i] = imag(x[3,1,i]+x[1,3,i])*0.5
                c[3,i] = real(x[1,3,i]-x[3,1,i])*0.5
                c[4,i] = imag(x[1,1,i]-x[3,3,i])*0.5
            end
        elseif k==4

            for i=1:n
                #println("i = $i")
                #println(c[:,i])
                if id[i] == 1
                    x[1,1,i] = c[1,i] + im*c[4,i]
                    x[1,2,i] = c[3,i] + im*c[2,i]
                    x[1,3,i] = 0
                    x[2,1,i] = -c[3,i] + im*c[2,i]
                    x[2,2,i] = c[1,i] - im*c[4,i]
                    x[2,3,i] = 0
                    x[3,1,i] = 0
                    x[3,2,i] = 0
                    x[3,3,i] = 1

                elseif id[i] == 0
                    x[1,1,i] = 1
                    x[1,2,i] = 0
                    x[1,3,i] = 0
                    x[2,1,i] = 0
                    x[2,2,i] = 1
                    x[2,3,i] = 0
                    x[3,1,i] = 0
                    x[3,2,i] = 0
                    x[3,3,i] = 1
                end 
            end
        elseif k==5
            for i=1:n
                if id[i] == 1
                    x[1,1,i] = 1
                    x[1,2,i] = 0
                    x[1,3,i] = 0
                    x[2,1,i] = 0
                    x[2,2,i] = c[1,i] + im*c[4,i]
                    x[2,3,i] = -c[3,i] + im*c[2,i]
                    x[3,1,i] = 0
                    x[3,2,i] = c[3,i] + im*c[2,i]
                    x[3,3,i] = c[1,i] -im*c[4,i]

                elseif id[i] == 0
                    x[1,1,i] = 1
                    x[1,2,i] = 0
                    x[1,3,i] = 0
                    x[2,1,i] = 0
                    x[2,2,i] = 1
                    x[2,3,i] = 0
                    x[3,1,i] = 0
                    x[3,2,i] = 0
                    x[3,3,i] = 1
                end 
            end

        elseif k==6
            for i=1:n
                if id[i] == 1
                    x[1,1,i] = c[1,i] + im*c[4,i]
                    x[1,2,i] = 0
                    x[1,3,i] = c[3,i] + im*c[2,i]
                    x[2,1,i] = 0
                    x[2,2,i] = 1
                    x[2,3,i] = 0
                    x[3,1,i] = -c[3,i] + im*c[2,i]
                    x[3,2,i] = 0
                    x[3,3,i] = c[1,i] -im*c[4,i]

                elseif id[i] == 0
                    x[1,1,i] = 1
                    x[1,2,i] = 0
                    x[1,3,i] = 0
                    x[2,1,i] = 0
                    x[2,2,i] = 1
                    x[2,3,i] = 0
                    x[3,1,i] = 0
                    x[3,2,i] = 0
                    x[3,3,i] = 1
                end 
            end
        end
    end

    function rndprd!(ranf,n)
        rn = zeros(Float64,n)
        rndprd!(ranf,rn,n)
        return rn
    end

    function rndprd!(ranf,rn,n)
        for i=1:n
            rn[i] = ranf()
        end
        return rn
    end

    function rndprd2!(ranf,n)
        xrn = zeros(Float64,3,n)
        rndprd2!(ranf,xrn,n)
        return xrn
    end

    function rndprd2!(ranf,xrn,n)
        for j=1:n
            for i=1:3
                xrn[i,j] = ranf()
            end
        end
        return 
    end


end