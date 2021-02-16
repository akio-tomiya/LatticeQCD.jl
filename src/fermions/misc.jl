module Othermethods
    using Random,Distributions
    using LinearAlgebra
    #using DoubleExponentialFormulas
    #using DifferentialEquations
    #using QuadGK

    import ..CGmethods:bicg,cg,shiftedcg

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

    function chebyshevfit(func,m,nc)
        Δt = (π/2)/2m
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
            #count += 1
            #println(count)
            val = 0
            #println("random")
            for k=1:M
                for i=1:N
                    ϕ[i] = rand([-1,1])
                end
                #rand!(d,ϕ)
                #Dinvϕ = DtN \ ϕ
                cg(ϕtemp1,DtN,ϕ,eps = 1e-10,maxsteps= 3000)
                mul!(Aϕ,ANmat,ϕtemp1)
                val += ϕ ⋅ Aϕ /M
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
end