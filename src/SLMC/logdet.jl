module Logdet
    using LinearAlgebra
    using Distributions
    using SparseArrays
    export trlog

    function trlog(A)
        itemax = 80
        num = 80
        return trlog(A,itemax,num)
    end

    function calclog(T)
        e,v = eigen(T)
        valuei = 0.0
        kmax = size(T)[1]
        for l=1:kmax
            valuei +=  log(e[l])*abs(v[1,l])^2
        end
        return valuei
    end

    function trlog(A,itemax,num)
        eps = 1e-4
        return trlog(A,itemax,num,eps)
    end

    function trlog(A,itemax,num,eps)
        println("Estimation of log(det(A))=Tr(log(A)) start------------------------------------------")
        println("Num. of Lanczos steps: $itemax")
        println("Num. of random vectors: $num")
        
        n = size(A)[1]
        value = 0.0
        oldvalue = 0.0
        eps2= 1e-6
        g = zeros(Float64,2,2)
        y = zeros(Float64,2)
        #value1 = 0
        #println("trace(A) = ",tr(A))
        #trA = tr(A)        
        #trvalue = 0
        #trvalue1 = 0
        d = Normal(0,1/sqrt(2))
        #d = Normal(0,1)
        #ene,ve = eigen(A)
        #println(ene)
        #logdet = sum(log.(ene))
#        logdet = tr(log(A))
        #meanvalue = 0.0
        icount = 0
        #fp = open("hi2.dat","w")
        for i=1:num
            #phii = (rand(n)+im*rand(n))
            #phii = (rand(n) .- 0.5).*2 .+ im*(rand(n) .- 0.5).*2 
            phii = rand(d,n)+im*rand(d,n)
            #nphi = norm(phii)
            #dotd = dot(phii,phii)
            #println("dot0 ", dotd)
            #phii = phii ./norm(phii)
            #println("dot ", dot(phii,phii))
            #phii = zeros(Float64,n)
            #phii[i] = 1
            #phii = rand([-1.0,1.0],n)
#            phii = rand([-1.0,1.0,im,-im],n)
            nphi = norm(phii)
            phii = phii/norm(phii)
            valuei = Lanczos(A,phii,itemax,eps)
            #T = Lanczos(A,phii,itemax)
            #e,v = eigen(T)
            
            #valuei = 0
            #valuei1 = 0
            #println(e)
            #exit()
            #for l=1:itemax
                #println(e[l],"\t",abs(v[1,l])^2)
            #    valuei +=  log(e[l])*abs(v[1,l])^2
            #    valuei1 +=  e[l]*abs(v[1,l])^2
            #end
            #println(valuei1,"\t",dot(phii,A*phii))
            #exit()
            #trvalue += dot(phii,A*phii)
            #trvalue1 += dot(phii,phii)
            value += valuei*nphi^2
            #=
            if i > itemax-20
                icount += 1
                meanvalue += value/i
            end
            =#
            if i % div(num,10) == 0
                println("$i-th ", value/i)
            end
            #=
            g[1,1] += 1
            g[1,2] += 1/i
            g[2,1] += 1/i
            g[2,2] += 1/i^2
            y[1] += value/i
            y[2] += (1/i)*value/i
            if i > 20
                ceff = g \ y
                #println(ceff)
                #println(ceff[1]+ceff[2]/i)
            end
            =#
        
            
            if i >1
                hi = abs(oldvalue-value/i)/abs(oldvalue)
                if abs(oldvalue-value/i) < eps2
                    println("In the $i-th step, this is converged. Absolute error might be ",abs(oldvalue-value/i))
                    println("Estimation of log(det(A))=Tr(log(A)) end------------------------------------------")
                    return value/i
                end
                ##println(i,"\t",hi,"\t",value/i,"\t",oldvalue-value/i)#,"\t",meanvalue/icount)
                #if i > 20
                    #println(fp,i,"\t",hi,"\t",value/i,"\t",oldvalue-value/i)#,"\t",ceff[1]+ceff[2]/i)
                #end
            end
            

            oldvalue = value/i
            #value1 += valuei1*nphi^2
            #println("$i-th value = ",value,"\t",logdet,"\t",abs(value-logdet)/abs(logdet))
            #println("$i-th trvalue = ",trvalue,"\t",trA,"\t",abs(trvalue-trA)/abs(trA),"\t",value1)
            #println("$i-th value = ",value/i,"\t",logdet,"\t",abs(value/i-logdet)/abs(logdet))
            #println("$i-th trvalue = ",trvalue/i,"\t",trA,"\t",abs(trvalue/i-trA)/abs(trA),"\t",value1/i)
            #println("$i-th trvalue1 = ",trvalue1/i)
            #println("eigenvalues = ",real.(e))
        end
        #close(fp)
        println("Estimation of log(det(A))=Tr(log(A)) end------------------------------------------")
        #exit()

        return value/num
        
    end

    function Lanczos(A,v0,itemax)
        eps = 1e-8
        return Lanczos(A,v0,itemax,eps)
    end

    function Lanczos(A,v0,itemax,eps)
        ipt = Int64[1,2]
        n = length(v0)
        ftype = eltype(v0)
        β = zeros(Float64,itemax+1)
        α = zeros(Float64,itemax+1)
        w = zeros(ftype,n,2)
        x = similar(v0)
        w[:,2] = v0[:]/norm(v0)
        T = zeros(ftype,itemax+1,itemax+1)
        trlog = 0.0
        e = 0
        V = 0
        oldvaluei = 0.0
        valuei = 0.0
        hi = 0.0

        @inbounds for k in 1:itemax
            #vj = w[:,ipt[2]]
            #vj-1 = w[:,ipt[1]]
            mul!(x,A,w[:,ipt[2]])
            #x = A*w[:,ipt[2]] #Avj = A*v
            #mul!(w[:,ipt[1]],A,w[:,ipt[2]],1,- β[k]) #mul!(C,A,B,alpha,beta) -> C = A*B*alpha + C*beta
            @. w[:,ipt[1]] = x[:] - β[k]*w[:,ipt[1]] #wj = Avj - betaj*vj-1
            α[k] = real(dot(w[:,ipt[1]],w[:,ipt[2]])) # alphaj = (wj,vj)
            @. w[:,ipt[1]] = w[:,ipt[1]]-α[k]*w[:,ipt[2]] #wj = wj - alphaj*vj
            β[k+1] = norm(w[:,ipt[1]]) #betaj+1 = norm(w[:,ipt[1]])
            @. w[:,ipt[1]] /= β[k+1] #vj = wj/betaj+1
            ii = ipt[1]
            ipt[1] = ipt[2]
            ipt[2] = ii
            T[k,k] = α[k]
            T[k,k+1] =  β[k+1]
            T[k+1,k] =  β[k+1] 
            valuei = calclog(T[1:k,1:k])
            #println("$k-th value: ",valuei)
            
            if k>1
                hi = abs(oldvaluei-valuei)/abs(valuei)
                #println("$k-th ",hi)
                if hi < eps
                    return  valuei
                end
            end
            oldvaluei = valuei 
            
            #println("$k-th alpha: ",abs(α[k]-α[k-1])/abs(α[k]))
            #println("$k-th beta: ",abs(β[k+1]-β[k])/abs(β[k+1]))
            #println("$k-th beta: ",β[k+1])
            
            #ista = ifelse(k+1-4 > 0,k+1-4,k+1)
            #trlog = 0.0
            #for i=1:k+1
            #    trlog += log(e[i])*abs(V[i,1])^2
            #end
            #println(e[ista:end])            
            #println(trlog)
            #e,V =  eigen(T[1:itemax,1:itemax])
            #println(e)
        end
        
        #T =  eigen(T[1:itemax,1:itemax])
        #println(e)
        #println(V*v0)
        #exit()
        println("Warning!: Lanczos in logdet.jl is not converged. error: ",hi)
        return valuei

        #return T[1:itemax,1:itemax]
    end
end

#=

using LinearAlgebra
using .Logdet

using Random
using Distributions
using SparseArrays
Random.seed!(2234)
d = Normal()

n = 2000
#A = sprand(Float64,n,n,3/n)
A = zeros(Float64,n,n)
t = -1
for i=1:n
    j = i+1
    j += ifelse(j > n,-n,0)
    A[i,j] = -t

    j = i-1
    j += ifelse(j <1,n,0)
    A[i,j] = -t
    A[i,i] = rand() 
end

#A = rand(d,n,n) ./n

AdagA = A'*A  + Diagonal(ones(n))


#AdagA = Matrix(AdagA)
e,v = eigen(AdagA)
println(e[end-4:end])

itemax = 50*2
v0 = rand(n)
#Logdet.Lanczos(AdagA,v0,itemax)

trlog = Logdet.trlog(AdagA,itemax,40)

#logdet = log(det(AdagA))
trlog_exact = tr(log(AdagA))
println("trlog_exact AdagA = ", trlog_exact)
println("trlog AdagA = ", trlog)

=#