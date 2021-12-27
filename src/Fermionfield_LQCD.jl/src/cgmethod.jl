module CGmethods
    using LinearAlgebra
    using InteractiveUtils
    import ..Gaugefield:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3
    #import ..Verbose_print:Verbose_level,Verbose_3,Verbose_2,Verbose_1,println_verbose3

    function add!(b,Y,a,X) #b*Y + a*X -> Y
        axpby!(a,X,b,Y) #X*a + Y*b -> Y
    end

    function add!(Y,a,X) #Y + a*X -> Y
        axpby!(a,X,1,Y) #X*a + Y -> Y
    end

    function bicg(x,A,b;eps=1e-10,maxsteps = 1000,verbose = Verbose_2()) #Ax=b
        println_verbose3(verbose,"--------------------------------------")
        println_verbose3(verbose,"bicg method")
        res = deepcopy(b)
        temp1 = similar(x)
        p = similar(x)
        q = similar(x)
        s = similar(x)

        mul!(temp1,A,x)
        add!(res,-1,temp1)

        rnorm = real(res⋅res)

        if rnorm < eps
            return
        end
        #println(rnorm)

        mul!(p,A',res)
        c1 = p ⋅ p

        for i=1:maxsteps
            mul!(q,A,p)
            #! ...  c2 = < q | q >
            c2 = q ⋅ q
            
            alpha = c1 / c2
            #! ...  x   = x   + alpha * p  
            add!(x,alpha,p)
            #...  res = res - alpha * q 
            add!(res,-alpha,q)
            rnorm = real(res ⋅ res) 
            println_verbose3(verbose,"$i-th eps: $rnorm")

            if rnorm < eps
                println_verbose3(verbose,"Converged at $i-th step. eps: $rnorm")
                println_verbose3(verbose,"--------------------------------------")
                return
            end

            mul!(s,A',res)

            #c3 = s * s
            c3 = s ⋅ s

            beta = c3 / c1
            c1 = c3

            add!(beta,p,1,s) #p = beta*p + s

        end

        
        error("""
        The BICG is not converged! with maxsteps = $(maxsteps)
        residual is $rnorm
        maxsteps should be larger.""")


    end

    function cg(x,A,b;eps=1e-10,maxsteps = 1000,verbose = Verbose_2()) #Ax=b
        println_verbose3(verbose,"--------------------------------------")
        println_verbose3(verbose,"cg method")
        res = deepcopy(b)
        temp1 = similar(x)
        
        q = similar(x)

        mul!(temp1,A,x)
        add!(res,-1,temp1)
        p = deepcopy(res)

        rnorm = real(res⋅res)
        #println(rnorm)

        if rnorm < eps
            return
        end

        c1 = p ⋅ p

        for i=1:maxsteps
            mul!(q,A,p)
            
            c2 = p ⋅ q
            
            α = c1 / c2
            #! ...  x   = x   + alpha * p  
            add!(x,α,p)
            #...  res = res - alpha * q 
            add!(res,-α,q)
            c3 = res ⋅ res
            rnorm = real(c3) 
            println_verbose3(verbose,"$i-th eps: $rnorm")
            

            if rnorm < eps
                println_verbose3(verbose,"Converged at $i-th step. eps: $rnorm")
                println_verbose3(verbose,"--------------------------------------")
                return
            end

            β = c3 / c1
            c1 = c3

            add!(β,p,1,res) #p = beta*p + s

        end

        
        error("""
        The CG is not converged! with maxsteps = $(maxsteps)
        residual is $rnorm
        maxsteps should be larger.""")


    end

    function Base.:*(x::Array{T,1},y::Array{T,1}) where T <: Number
        return x'*y
    end



    function shiftedcg(vec_x,vec_β,x,A,b;eps=1e-10,maxsteps = 1000,verbose = Verbose_2()) #Ax=b
        
        println_verbose3(verbose,"--------------------------------------")
        println_verbose3(verbose,"shifted cg method")
        N = length(vec_β)
        temp1 = similar(b)
        r = deepcopy(b)
        p = deepcopy(b)
        q = similar(b)

        vec_r = Array{typeof(r),1}(undef,N)
        vec_p = Array{typeof(p),1}(undef,N)
        for j=1:N
            vec_r[j] = deepcopy(b)
            vec_p[j] = deepcopy(b)
        end

        αm = 1.0
        βm = 0.0

        ρm = ones(ComplexF64,N)
        ρ0 = ones(ComplexF64,N)
        ρp = ones(ComplexF64,N)
        residual = 0


        for i=1:maxsteps
            mul!(q,A,p)

            pAp = p ⋅ q

            rr = dot(r,r)
            αk = rr / pAp


            #! ...  x   = x   + alpha * p   
            add!(x,αk,p)

            #...  r = r - alpha * q 
            add!(r,-αk,q)

            βk = dot(r,r)/ rr
            add!(βk,p,1,r) #p = beta*p + r

            for j=1:N
                ρkj = ρ0[j]
                if abs(ρkj) < eps
                    continue
                end
                ρkmj =ρm[j]
                ρp[j] = ρkj*ρkmj*αm/(ρkmj*αm*(1.0+αk*vec_β[j])+αk*βm*(ρkmj-ρkj))
                αkj = (ρp[j]/ρkj)*αk

                add!(vec_x[j],αkj,vec_p[j])
                βkj = (ρp[j]/ρkj)^2*βk
                add!(βkj,vec_p[j],ρp[j],r) 

            end

            ρm[:] = ρ0[:]
            ρ0[:] = ρp[:]
            αm = αk
            βm = βk


            ρMAX = maximum(abs.(ρp))^2
            residual = abs(rr*ρMAX)
            println_verbose3(verbose,"$i-th eps: $residual")

            if abs(residual) < eps
                println_verbose3(verbose,"Converged at $i-th step. eps: $residual")
                println_verbose3(verbose,"--------------------------------------")
                return
            end


        end

        
        error("""
        The shifted CG is not converged! with maxsteps = $(maxsteps)
        residual is $residual
        maxsteps should be larger.""")


    end

    function reducedshiftedcg(leftvec,vec_β,x,A,b;eps=1e-10,maxsteps = 1000,verbose = Verbose_2()) #Ax=b
        println_verbose3(verbose,"--------------------------------------")
        println_verbose3(verbose,"shifted cg method")
        N = length(vec_β)
        temp1 = similar(b)
        r = deepcopy(b)
        p = deepcopy(b)
        q = similar(b)

        #=
        vec_r = Array{typeof(r),1}(undef,N)
        vec_p = Array{typeof(p),1}(undef,N)
        for j=1:N
            vec_r[j] = deepcopy(b)
            vec_p[j] = deepcopy(b)
        end
        =#

        Σ = leftvec ⋅ b

        θ = zeros(ComplexF64,N)
        Π = ones(ComplexF64,N) .* Σ


        αm = 1.0
        βm = 0.0

        ρm = ones(ComplexF64,N)
        ρ0 = ones(ComplexF64,N)
        ρp = ones(ComplexF64,N)
        residual = 0


        for i=1:maxsteps
            mul!(q,A,p)

            pAp = p ⋅ q
            rr = r*r
            αk = rr / pAp

            #! ...  x   = x   + alpha * p   
            add!(x,αk,p)

            #...  r = r - alpha * q 
            add!(r,-αk,q)

            βk = (r*r)/ rr
            add!(βk,p,1,r) #p = beta*p + r

            Σ = leftvec ⋅ r

            for j=1:N
                ρkj = ρ0[j]
                if abs(ρkj) < eps
                    continue
                end
                ρkmj =ρm[j]
                ρp[j] = ρkj*ρkmj*αm/(ρkmj*αm*(1.0+αk*vec_β[j])+αk*βm*(ρkmj-ρkj))
                αkj = (ρp[j]/ρkj)*αk
                θ[j] += αkj*Π[j]

                #add!(vec_x[j],αkj,vec_p[j])
                βkj = (ρp[j]/ρkj)^2*βk
                Π[j] = βkj*Π[j] + ρp[j]*Σ

                #mul!(Π,ρp[j],Σ,1,βkj)
                #Π[:] = 
                #add!(βkj,vec_p[j],ρp[j],r) 

            end

            ρm[:] = ρ0[:]
            ρ0[:] = ρp[:]
            αm = αk
            βm = βk


            ρMAX = maximum(abs.(ρp))^2
            residual = abs(rr*ρMAX)
            println_verbose3(verbose,"$i-th eps: $residual")

            if abs(residual) < eps
                println_verbose3(verbose,"Converged at $i-th step. eps: $residual")
                println_verbose3(verbose,"--------------------------------------")
                return θ
            end


        end

        
        error("""
        The shifted CG is not converged! with maxsteps = $(maxsteps)
        residual is $residual
        maxsteps should be larger.""")


    end

end