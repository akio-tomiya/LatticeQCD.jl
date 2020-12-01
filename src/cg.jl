module tesmo
end

module CGfermion
    #using Base.Threads
    #using Distributed
    import ..Fermionfields:FermionFields,WilsonFermion,
                Wx!,Wdagx!,add!,WdagWx!
    import ..Actions:FermiActionParam,FermiActionParam_WilsonClover
    import ..Clover:Make_CloverFμν!
    import ..LTK_universe:Universe


    

"""
c-----------------------------------------------------------------------
c     this routine  solves
c                  W(u) * x = b      for iflag=1
c                  W_adj(u) * x = b  for iflag=2
c     by conjugate gradient 
c
c     eps and imax controle when we stop CG iteration, i.e.,
c     <res|res> < eps or (# of CG iteration)>imax, the iteration is
c     terminated.
c----------------------------------------------------------------------
"""
    function cg0!(x,b,iflag,univ::Universe)
        cg0!(x,b,iflag,univ.U,univ._temporal_gauge,univ._temporal_fermi,univ.fparam)
    end

    function cg0!(x,b,iflag,u,temps_g,temps,fparam::T) where T <: FermiActionParam
        #id = threadid()
        #println("id = $id")
        icheck = 1

        temp1 = similar(x)
        res = similar(x)
        p = similar(x)
        q = similar(x)
        s = similar(x)

        if typeof(fparam) == FermiActionParam_WilsonClover
            Make_CloverFμν!(fparam,u,temps_g)

        end

        #println(b*b)
        #exit()

        #the initial condition  ( for i=0 )
        #! ...  res = b - W*x        
        if iflag == 1
            Wx!(temp1,u,x,temps,fparam)

            add!(res,1,b,-1,temp1)


        elseif iflag == 2

            Wdagx!(temp1,u,x,temps,fparam)



            add!(res,1,b,-1,temp1)

            

        end 
        #display(res)

        rnorm = real(res * res)
        #println(rnorm)

        if rnorm < x.eps
            return
        end

        #println(res*res)
        # ...  p = W_adj * res
        if iflag == 1
            Wdagx!(p,u,res,temps,fparam)
        elseif iflag == 2
            #println("p ",p*p)
            Wx!(p,u,res,temps,fparam)

        end 


        # ...  c1 = < p | p >
        c1 = p * p

        #...  the iteration starts
        for i=1:fparam.MaxCGstep

            #...  q = W * p
            if iflag == 1
                #Field_f.wxvect!(q,u,p,temps,2)
                Wx!(q,u,p,temps,fparam)
            elseif iflag == 2
                Wdagx!(q,u,p,temps,fparam)
                #Field_f.wxvect!(q,u,p,temps,3)
            end 

            #! ...  c2 = < q | q >
            c2 = q * q 
            alpha = c1 / c2

            #! ...  x   = x   + alpha * p   
            add!(x,alpha,p)


            #...  res = res - alpha * q 
            add!(res,-alpha,q)

        

            #.....   check of the convergence   ...............
            rnorm = real(res * res)

            #println(rnorm)

            #if icheck == 1
                #println("$i rnorm : ",rnorm)
            #end

            if rnorm < fparam.eps
                #=
                if iflag == 1
                    Wx!(temp1,u,x,temps,fparam)
                    #Field_f.wxvect!(temp1,u,x,temps,2)
                    add!(res,1,b,-1,temp1)
                    #Field_f.add!(res,1,b,-1,temp1)
                elseif iflag == 2
                    Wdagx!(temp1,u,x,temps,fparam)
                    #Field_f.wxvect!(temp1,u,x,temps,3)
                    add!(res,1,b,-1,temp1)
                    #Field_f.add!(res,1,b,-1,temp1)
                end
                =#
                #println("converged! at $i : ",res*res,"\t x \t ",x*x)
                #println("\t")
                
                return
            end

            #! ...  s = W_adj * res 
            if iflag == 1
                Wdagx!(s,u,res,temps,fparam)
                #Field_f.wxvect!(s,u,res,temps,3)
            elseif iflag == 2
                Wx!(s,u,res,temps,fparam)
                #Field_f.wxvect!(s,u,res,temps,2)
            end 

        
            c3 = s * s

            beta = c3 / c1
            c1 = c3


            add!(beta,p,1,s) #p = beta*p + s


        
            if i==fparam.MaxCGstep
                error("""
                The CG is not converged! with fparam.MaxCGstep = $(fparam.MaxCGstep)
                residual is $rnorm
                fparam.MaxCGstep should be larger.""")

            end
        end


    end

    
"""
c-----------------------------------------------------------------------
c     this routine  solves
c                  (W(u)^+W(u)) * x = b    
c     by conjugate gradient 
c
c     eps and imax controle when we stop CG iteration, i.e.,
c     <res|res> < eps or (# of CG iteration)>imax, the iteration is
c     terminated.
c----------------------------------------------------------------------
"""
    function cg0_WdagW!(x,b,u,temps_g,temps,fparam::T) where T <: FermiActionParam
        temp1 = similar(x)
        res = similar(x)
        p = similar(x)
        q = similar(x)
        s = similar(x)

        WdagWx!(temp1,u,x,temps,fparam)
        add!(res,1,b,-1,temp1) # res = b - temp1
        add!(p,1,res,0,temp1) #p = res

        rnorm = real(res * res)
        #println("0 rnorm : ",rnorm)

        if rnorm < x.eps
            return
        end

        # ...  c1 = < p | p >
        c1 = p * p



        for i=1:fparam.MaxCGstep
            
            WdagWx!(q,u,p,temps,fparam)
            #println(q*q)
            c2 = p * q

            α = c1 / c2

            #! ...  x   = x   + alpha * p   
            add!(x,α,p)
            #...  res = res - alpha * q 
            add!(res,-α,q)

            c3 = res*res

            #.....   check of the convergence   ...............
            rnorm = real(c3)
            #println("$i rnorm : ",rnorm)

            if rnorm < fparam.eps
                #=
                WdagWx!(temp1,u,x,temps,fparam)
                add!(res,1,b,-1,temp1)

                println("converged! at $i : ",res*res,"\t x \t ",x*x)
                println("\t")
                =#
                return
            end

            β = c3 / c1
            c1 = c3

            add!(β,p,1,res) #p = beta*p + res
        
            if i==fparam.MaxCGstep
                error("""
                The CG is not converged! with fparam.MaxCGstep = $(fparam.MaxCGstep)
                residual is $rnorm
                fparam.MaxCGstep should be larger.""")

            end


        end


    end

    
    function shiftedcg0_WdagW!(vec_x,x,b,vec_β,u,temps_g,temps,fparam::T) where T <: FermiActionParam
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

        for i=1:fparam.MaxCGstep
            WdagWx!(q,u,p,temps,fparam)
            pAp = p * q
            rr = r*r
            αk = rr / pAp
            #! ...  x   = x   + alpha * p   
            add!(x,αk,p)

            #...  r = r - alpha * q 
            add!(r,-αk,q)

            

            βk = (r*r)/ rr
            add!(βk,p,1,r) #p = beta*p + r

            for j=1:N
                ρkj = ρ0[j]
                if abs(ρkj) < fparam.eps
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

            if abs(residual) < fparam.eps

                return
            end

            if i==fparam.MaxCGstep
                error("""
                The CG is not converged! with fparam.MaxCGstep = $(fparam.MaxCGstep)
                residual is $rnorm
                fparam.MaxCGstep should be larger.""")

            end

        end


        return

    end

end