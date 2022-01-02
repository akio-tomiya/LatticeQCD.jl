module Othermethods
    using Random,Distributions
    using LinearAlgebra

    using FastGaussQuadrature

    #using DoubleExponentialFormulas
    #using DifferentialEquations
    #using QuadGK

    import ..CGmethods:bicg,cg,shiftedcg,reducedshiftedcg
    import ..Diracoperators:DdagD_operator,DdagD_Staggered_operator,DdagDND_Staggered_operator
    import ..Fermionfields:gauss_distribution_fermi_Z2!
    import ..Fermionfields
    

    function Tn(n,x)
        return cos(n*acos(x))
    end

    function fitfunc(an,x)
        val = 0
        for n=1:2:length(an)-1
            val += an[n+1]*Tn(n,x)
        end 
        return val
    end

    function integcheb(an,a,b)
        wa = 0
        for k=1:length(an)-2
            wa += (Tn(k,b)-Tn(k,a))*(an[k-1+1]-an[k+1+1])/k
        end
        wa /= 2
        return wa
    end

    function calc_vecfunc_cg(func,vec_x,M,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1)
        DN = Adiag .+ ANmat

        numdata = length(vec_x)
        data = zeros(Float64,numdata)

        #Aϕ2 = zero(Aϕ)
        #ϕtemp2 = zero(ϕtemp1)


        #x = zero(ϕ)
        #Nx = zero(x)


        #=
        vec_Ax = Array{typeof(x),1}(undef,numdata)
        for l=1:numdata
            vec_Ax[l] = zero(x)
        end
        =#
        vec_tinv = Adiag[1,1]*(1 ./ vec_x  .- 1) #(D/t + A) = ((D(1/t - 1)I + D I + A)

        
        for k=1:M
            for i=1:N
                ϕ[i] = rand([-1,1])
                #ϕ[i] = randn()
                #ϕ[i] = rand([-1,1,-im,im])
                #r = rand()*2π
                #ϕ[i] = cos(r)+im*sin(r)
            end
            mul!(Aϕ,ANmat,ϕ)

            #cg(ϕtemp1,DtN,ϕ,eps = 1e-15,maxsteps= 3000)
            #shiftedcg(vec_x,vec_β,x,WdagW,univ.η,eps = univ.fparam.eps,maxsteps= univ.fparam.MaxCGstep)
            #@time shiftedcg(vec_Ax,vec_tinv,x,DN,ϕ,eps = 1e-15,maxsteps= 3000)
            ϕtemp1 .= 0
            θ = reducedshiftedcg(Aϕ,vec_tinv,ϕtemp1,DN,ϕ,eps = 1e-15,maxsteps= 3000)
            θ ./= (M .* vec_x)
            data .+= real.(θ)

            #=

            for l=1:numdata
                t = vec_x[l]
                DtN = Adiag .+ t*ANmat
                cg(ϕtemp2,DtN,ϕ,eps = 1e-15,maxsteps= 3000)
                mul!(Aϕ2,ANmat,ϕtemp2)
                val = real(ϕ ⋅ Aϕ2 /M)
                println("$t data,val ",real(θ[l]),"\t",val)
            end
            =#
        
            #=
            for l=1:numdata
                #mul!(Aϕ,ANmat,vec_Ax[l])
                valshift = real(θ[l]) #real(Nx ⋅ vec_Ax[l] /(M*vec_x[l]))
                #println("val1 ",valshift,"\t theta ",real(θ[l]))

                #=
                mul!(Aϕ,ANmat,vec_Ax[l])
                valshift = real(ϕ ⋅ Aϕ /(M*vec_x[l]))
                println("val2 ",valshift)
                =#

                data[l] += valshift
                

                t = vec_x[l]
                DtN = Adiag .+ t*ANmat
                cg(ϕtemp1,DtN,ϕ,eps = 1e-15,maxsteps= 3000)
                mul!(Aϕ,ANmat,ϕtemp1)
                val = real(ϕ ⋅ Aϕ /M)
                println("$t data,val ",valshift,"\t",val)
                #val += ϕ ⋅ Aϕ /M
            end
            ϕtemp1 .= 0
            =#
            #=
            for l=1:numdata
                vec_Ax[l] .= 0
            end
            =#

            

        end
        #exit()

        return data

        exit()


        for (i,x) in enumerate(vec_x)
            #val1 = func(x)
            #println("1 ",val1)
            val1,val2 = calc_fa(x,M,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1)
            #println("2 ",val1)
            #exit()
            println("1 ",data[i])
            data[i] = val1
            println("2 ",data[i])
        end
        return data
    end

    function calc_vecfunc(func,vec_x)
        N = length(vec_x)
        data = zeros(Float64,N)
        for (i,x) in enumerate(vec_x)
            data[i] = func(x)
        end
        return data
    end

    
    function calc_fa(t,M,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1)
        DtN = Adiag .+ t*ANmat

        if t < 1e-15
            DtN2 = Adiag
        else
            DtN2 = Adiag ./t .+ ANmat
        end
        #count += 1
        #println(count)
        val = 0
        val2 = 0
        #println("random")
        for k=1:M
            #gauss_distribution_fermi_Z2!(ϕ)

            for i=1:N
                ϕ[i] = rand([-1,1])
            end
            #rand!(d,ϕ)
            #Dinvϕ = DtN \ ϕ
            cg(ϕtemp1,DtN,ϕ,eps = 1e-15,maxsteps= 3000)
            mul!(Aϕ,ANmat,ϕtemp1)
            val += ϕ ⋅ Aϕ /M

            if t < 1e-15
                ϕtemp1 = ϕ/Adiag[1,1]
                mul!(Aϕ,ANmat,ϕtemp1)
                val2 += ϕ ⋅ Aϕ /M
            else
                cg(ϕtemp1,DtN2,ϕ,eps = 1e-15,maxsteps= 3000)
                mul!(Aϕ,ANmat,ϕtemp1)
                val2 += ϕ ⋅ Aϕ /(M*t)
            end

            
            #exit()


            #val += ϕ'*Aϕ/M#ANmat*Dinvϕ/M
        end
        println("ts: $t ",real(val),"\t",real(val2))
        return real(val),real(val2)
    end



    #calc_f,m,nc,M,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1
    function chebyshevfit(func,m::T,nc,M,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1) where T <: Int
        Δt = (π/2)/2m
        vec_x = zeros(Float64,2m+1)

        for i=0:2m
            t = i*Δt
            x = cos(t)
            vec_x[i+1] = x
        end
        @time data = calc_vecfunc_cg(func,vec_x,M,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1)

        #fp = open("hiseki.dat","w")
        #close(fp)

        #=
        data1 = calc_vecfunc(func,vec_x)
        for i=1:length(data)
            println(data[i],"\t",data1[i])
        end
        exit()
        =#
        an = zeros(Float64,nc+1)

        for n=1:2:nc
            
            for i=1:m
                k = 2i-2
                k1 = 2i-1
                k2 = 2i
                t = k*Δt
                t1 = k1*Δt
                t2 = k2*Δt

                val = (2/π)*2*cos(n*t)*data[k+1]
                val1 = (2/π)*2*cos(n*t1)*data[k1+1]
                val2 = (2/π)*2*cos(n*t2)*data[k2+1]

                an[n+1] += Δt*(val + 4*val1 + val2)/3
            end
        end

        return an

    end

    function chebyshevfit(func,m,nc)
        Δt = (π/2)/2m
        vec_x = zeros(Float64,2m+1)
        
        for i=0:2m
            t = i*Δt
            x = cos(t)
            vec_x[i+1] = x
        end
        data = calc_vecfunc(vec_x,func)
        an = zeros(Float64,nc+1)

        for n=1:2:nc
            
            for i=1:m
                k = 2i-2
                k1 = 2i-1
                k2 = 2i
                t = k*Δt
                t1 = k1*Δt
                t2 = k2*Δt

                val = (2/π)*2*cos(n*t)*data[k+1]
                val1 = (2/π)*2*cos(n*t1)*data[k1+1]
                val2 = (2/π)*2*cos(n*t2)*data[k2+1]

                an[n+1] += Δt*(val + 4*val1 + val2)/3
            end
        end

        return an




        data = zeros(Float64,m)
        data1 = zeros(Float64,m)
        data2 = zeros(Float64,m)
        count = 0        
        for i=1:m
            t = (2i-2)*Δt
            x = cos(t)

            t1 = (2i-1)*Δt
            x1 = cos(t1)

            t2 = 2i*Δt
            x2 = cos(t2)

            val = func(x)
            val1 = func(x1)
            val2 = func(x2)

            data[i] = val
            data1[i] = val1
            data2[i] = val2
        end
        an = zeros(Float64,nc+1)

        for n=1:2:nc
            
            for i=1:m
                t = (2i-2)*Δt
                t1 = (2i-1)*Δt
                t2 = 2i*Δt

                val = (2/π)*2*cos(n*t)*data[i]
                val1 = (2/π)*2*cos(n*t1)*data1[i]
                val2 = (2/π)*2*cos(n*t2)*data2[i]

                an[n+1] += Δt*(val + 4*val1 + val2)/3
            end
        end

        return an

        

        fp = open("hiseki.dat","w")

        for i=1:m
            t = (2i-2)*Δt
            x = cos(t)

            t1 = (2i-1)*Δt
            x1 = cos(t1)

            t2 = 2i*Δt
            x2 = cos(t2)

            val = data[i]
            val1 = data1[i]
            val2 = data2[i]

            valx = 0
            for n=1:2:nc
                valx += an[n+1]*Tn(n,x)
            end 

            println("$x $val $valx")
            println(fp,"$x $val $valx")
        end
        close(fp)
        
        wa2 = 0
        Δx = 1/2m
        for i=1:m
            x = (2i-2)*Δx
            x1 = (2i-1)*Δx
            x2 = 2i*Δx

            val = fitfunc(an,x)
            val1 = fitfunc(an,x1)
            val2 = fitfunc(an,x2)

            wa2 += Δx*(val + 4*val1 + val2)/3
        end

        println(wa2)
        return wa2
        exit()


        


        return wa

    end

    

    function tdlogdet(A,M,m,tempvecs;filename=nothing,nc = 20,nonc = false)
        N,_ = size(A)
        Adiag = zero(A)
        d = Normal(0,1)
        ϕ = tempvecs[1]
        Aϕ = tempvecs[2]
        ϕtemp1 = tempvecs[3]

        #=
        for i=1:M
            rand!(d,ϕ)
            mul!(Aϕ,A,ϕ)
            for k=1:N
                Adiag[k,k] += real(conj(ϕ[k])*Aϕ[k]/M)
            end
        end
        =#

        for k=1:N
            Adiag[k,k] = A[1,1]
        end

        
        #for k=1:N
        #    println("$k $(A[k,k]) $(Adiag[k,k])")
        #end
        
        ANmat = A .- Adiag 
        Δ0 = 0
        for k=1:N
            Δ0 += real(log(Adiag[k,k]))
        end
        Δx = 1/2m
        
       

        #@time ldet = real(logdet(Matrix(A)))
        #println("log det A(exact): $ldet")

        count = 0

        function calc_f(t)
            DtN = Adiag .+ t*ANmat

            if t < 1e-15
                DtN2 = Adiag
            else
                DtN2 = Adiag ./t .+ ANmat
            end
            #count += 1
            #println(count)
            val = 0
            val2 = 0
            #println("random")
            for k=1:M
                #gauss_distribution_fermi_Z2!(ϕ)

                for i=1:N
                    ϕ[i] = rand([-1,1])
                end
                #rand!(d,ϕ)
                #Dinvϕ = DtN \ ϕ
                cg(ϕtemp1,DtN,ϕ,eps = 1e-15,maxsteps= 3000)
                mul!(Aϕ,ANmat,ϕtemp1)
                val += ϕ ⋅ Aϕ /M

                if t < 1e-15
                    ϕtemp1 = ϕ/Adiag[1,1]
                    mul!(Aϕ,ANmat,ϕtemp1)
                    val2 += ϕ ⋅ Aϕ /M
                else
                    cg(ϕtemp1,DtN2,ϕ,eps = 1e-15,maxsteps= 3000)
                    mul!(Aϕ,ANmat,ϕtemp1)
                    val2 += ϕ ⋅ Aϕ /(M*t)
                end

                
                #exit()


                #val += ϕ'*Aϕ/M#ANmat*Dinvϕ/M
            end
            println("t: $t ",real(val),"\t",real(val2))
            return real(val)
        end

        function calc_finv(t)
            if t == 0
                DtN = Adiag
            else
                DtN = Adiag ./t .+ ANmat
            end
            #count += 1
            #println(count)
            val = 0
            #println("random")
            for k=1:M
                #gauss_distribution_fermi_Z2!(ϕ)

                for i=1:N
                    ϕ[i] = rand([-1,1])
                end
                #rand!(d,ϕ)
                #Dinvϕ = DtN \ ϕ
                if t == 0
                    ϕtemp1 = ϕ/Adiag[1,1]
                    mul!(Aϕ,ANmat,ϕtemp1)
                    val += ϕ ⋅ Aϕ /M
                else
                    cg(ϕtemp1,DtN,ϕ,eps = 1e-15,maxsteps= 3000)
                    mul!(Aϕ,ANmat,ϕtemp1)
                    val += ϕ ⋅ Aϕ /(M*t)
                end
                

                
                #exit()


                #val += ϕ'*Aϕ/M#ANmat*Dinvϕ/M
            end
            println("t: $t ",val)
            return real(val)
        end

        Δx = 1/2m
        Δ = Δ0
        Δc = Δ0

        #I,E = quadgk(calc_f, 0,1,maxevals=3m)
        #println(I)
        #println(E)
        #exit()

        nodes, weights = gausslegendre( 1000 )
        vec_x = nodes ./ 2 .+ 1/2

        @time data = calc_vecfunc_cg(calc_f,vec_x,M,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1)
        Δ += dot(weights,data)/2

        return Δ


        #nc = 20
        
        if nonc == false
            #@time ldet = real(logdet(Matrix(A)))
            ldet = 0
            println("log det A(exact): $ldet")

            nodes, weights = gausslegendre( 1000 )
            vec_x = nodes ./ 2 .+ 1/2
            #=
            @time data = calc_vecfunc_cg(calc_f,vec_x,M,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1)
            Δ += dot(weights,data)/2 
            println(" Δ ",Δ)
            exit()
            =#

            #=
            fpk = open("kdep_L8_nc10_complex.dat","w")
            kmax=32
            for k=1:kmax
                M2 = 10
                nc2 = 10
                m2 = 1000
                @time an = chebyshevfit(calc_f,m2,nc2,M2,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1)
                Δ = integcheb(an,0,1) + Δ0

                #@time data = calc_vecfunc_cg(calc_f,vec_x,M2,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1)
                #Δ = dot(weights,data)/2 + Δ0
                println("$k $Δ $ldet")
                println(fpk,"$k $Δ $ldet")
            end
            close(fpk)
            exit()
            =#

            fp4 = open("MGL10_complex_k4.dat","w")
            kmax=4
            for M2 in [16,32,64,128,256,512,1024]
                Δ = Δ0
                for k=1:kmax
                    @time data = calc_vecfunc_cg(calc_f,vec_x,M2,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1)
                    Δ += (dot(weights,data)/2)/kmax
                end
                #Δ = dot(weights,data)/2 + Δ0
                println("M = $M2, Δ = $Δ, exact: $ldet ")
                println(fp4,"$M2 $Δ $ldet ")
            end
            close(fp4)
            exit()

            #=
            vec_x = zeros(Float64,2m+1)
        
            for i=0:2m
                x = i*Δx
                vec_x[i+1] = x
            end
            @time data = calc_vecfunc_cg(calc_f,vec_x,M,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1)

            for i=1:m
                k = 2i-2
                k1 = 2i-1
                k2 = 2i
                x = k*Δx
                x1 = k1*Δx
                x2 = k2*Δx

                val = data[k+1]
                val1 = data[k1+1]
                val2 = data[k2+1]

                Δ += Δx*(val + 4*val1 + val2)/3
            end
            println(" Δ ",Δ)
            =#


            fp4 = open("Mnc.dat","w")
            an = chebyshevfit(calc_f,m,nc,M,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1)
            for M2 in [16,32,64,128,256]
                for nc2 in [10,20,50,100,150,200]
                    an = chebyshevfit(calc_f,m,nc2,M2,N,ANmat,Adiag,ϕ,Aϕ,ϕtemp1)
                    Δ3 = integcheb(an,0,1) + Δ0
                    println("M = $M2, nc = $nc2, Δ = $Δ3, exact: $ldet ")
                    println(fp4,"$M2 $nc2 $Δ3 $ldet ")
                end
            end
            close(fp4)
            exit()
            Δ3 = integcheb(an,0,1) + Δ0
            println("an: ",an," Δ ",Δ3)

            #aninv = chebyshevfit(calc_finv,m,nc)
            #println("aninv: ",aninv)
            #Δ3 = integcheb(aninv,0,1) + Δ0
            #println("aninv: ",aninv," Δ ",Δ3)
            exit()
        end
        
    

        if filename != nothing
            fp = open(filename,"w")
        end
        for i=1:m
            x = (2i-2)*Δx
            x1 = (2i-1)*Δx
            x2 = 2i*Δx

            if nonc
    
                val = calc_f(x)
                val1 = calc_f(x1)
                val2 = calc_f(x2)
            else
                valc = fitfunc(an,x)
                val1c = fitfunc(an,x1)
                val2c = fitfunc(an,x2)
            end

            
            if filename != nothing
                if nonc
                
                    println(fp,x,"\t",val)
                    println(fp,x1,"\t",val1)
                    println(fp,x2,"\t",val2)
                else
                
                    println(fp,x,"\t",valc)
                    println(fp,x1,"\t",val1c)
                    println(fp,x2,"\t",val2c)
                end
                
                flush(fp)
            end

            if nonc
                Δ += Δx*(val + 4*val1 + val2)/3
            else
                Δc += Δx*(valc + 4*val1c + val2c)/3
            end
            #println("i=$i ",Δ,"\t",ldet)

        end

        if nonc
            return real(Δ)
        else
            return real(Δ3)
        end

        println(Δ)
        println("cheb ",Δc)
        println("cheb2 ", Δ3)
        


        @time ldet = real(logdet(Matrix(A)))
        println("log det A(exact): $ldet")
        exit()

        if filename != nothing
            close(fp)
        end
        #exit()

        return real(Δ)



    end

    function tdlogdet(A::T,M,m;filename=nothing,nc = 20,nonc = false) where T <: DdagD_operator
        ic,ix,iy,iz,it,α = 1,1,1,1,1,1
        temps = A.dirac._temporal_fermi
        xi = temps[7]
        x0 = temps[8]
        Fermionfields.clear!(x0)
        Fermionfields.clear!(xi)
        x0[ic,ix,iy,iz,it,α] = 1
        mul!(xi,A,x0,(ix,iy,iz,it,α))
        diagvalue = xi[ic,ix,iy,iz,it,α]
        #println(diagvalue)
        

        N,_ = size(A)
        Δ0 = 0
        Δ0 = N*real(log(diagvalue))
        Δx = 1/2m

        ϕ = temps[7]
        ϕtemp1 = temps[8]
        Aϕ = temps[9]

        ANmat = DdagDND_Staggered_operator(A) #non-diagonal part of A
        
    
        #@time ldet = real(logdet(Matrix(A)))
        #println("log det A(exact): $ldet")

        count = 0

        function calc_f(t)
            DtN = DdagD_Staggered_operator(A,t)
            #DtN = Adiag .+ t*ANmat
            #count += 1
            #println(count)
            val = 0
            #println("random")
            for k=1:M
                gauss_distribution_fermi_Z2!(ϕ)

                cg(ϕtemp1,DtN,ϕ,eps = 1e-10,maxsteps= 3000)
                mul!(Aϕ,ANmat,ϕtemp1)
                val += ϕ ⋅ Aϕ /M
                println("t: $t ",val)
                exit()
                #val += ϕ'*Aϕ/M#ANmat*Dinvϕ/M
            end
            return real(val)
        end

        Δx = 1/2m
        Δ = Δ0
        Δc = Δ0

        #I,E = quadgk(calc_f, 0,1,maxevals=3m)
        #println(I)
        #println(E)
        #exit()


        #nc = 20
        
        if nonc == false
            an = chebyshevfit(calc_f,m,nc)
            Δ3 = integcheb(an,0,1) + Δ0
        end
        
    

        if filename != nothing
            fp = open(filename,"w")
        end
        for i=1:m
            x = (2i-2)*Δx
            x1 = (2i-1)*Δx
            x2 = 2i*Δx

            if nonc
    
                val = calc_f(x)
                val1 = calc_f(x1)
                val2 = calc_f(x2)
            else
                valc = fitfunc(an,x)
                val1c = fitfunc(an,x1)
                val2c = fitfunc(an,x2)
            end

            
            if filename != nothing
                if nonc
                
                    println(fp,x,"\t",val)
                    println(fp,x1,"\t",val1)
                    println(fp,x2,"\t",val2)
                else
                
                    println(fp,x,"\t",valc)
                    println(fp,x1,"\t",val1c)
                    println(fp,x2,"\t",val2c)
                end
                
                flush(fp)
            end

            if nonc
                Δ += Δx*(val + 4*val1 + val2)/3
            else
                Δc += Δx*(valc + 4*val1c + val2c)/3
            end
            #println("i=$i ",Δ,"\t",ldet)

        end

        if nonc
            return real(Δ)
        else
            return real(Δ3)
        end



    end
end