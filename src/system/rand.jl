module Rand
    export Crand

    mutable struct Crand
        iw::Array{Int64,1}
        jr::Int64
        kr::Int64
    end

    mutable struct Random_LCGs
        crand_public::Crand
        Random_LCGs(ndelay) = new(cint3(ndelay))
    end

    function (g::Random_LCGs)()
        #return rand()
        return ranf!(g.crand_public)
    end


    const IP = 521
    const IQ = 32
    const MACRM = 40
    const MACRI = 1
    const FNORM=Float32(0.465661e-9)

    const intmax32 = Int32(2147483647)


    function cint3(ndelay)
        println("ndelay : $ndelay")
        #ccall((:readseed,"./randfortran.so"),Nothing,(),)
        crand_public = init3()
        delay3!(crand_public,ndelay)
        #display(crand_public.iw)
        
        return crand_public
    end

    function make_value_Int32_again(ix)
        while ix < -intmax32 - 1
            ix += 2*intmax32 + 1
        end
        while ix > intmax32
            ix -= 2*(intmax32 +1)
        end
        return ix
    end

    function init3()
        
        
        
        
        
        #println(2147483647+1-2*(intmax32+1))
        ib = zeros(Int64,IP) #INTEGER ib(0:IP-1)
        iw = zeros(Int64,IP) #INTEGER iw(0:IP-1)

        ix = Int32(MACRI)
        for i=0:IP-1
            ix = ix*69069
            ix = make_value_Int32_again(ix)

            #while ix < -intmax32-1 || intmax32 < ix
            #    ix -= Int32(Int32(2)*(intmax32+Int32(1)))
            #end
            #println(Int32(ix) >>> -(-31))
            ib[i+1] = Int32(ix) >>> -(-31)
            #println("$i ",ix," ",typeof(ix))
            #println(ix,"\t",ib[i+1])
            #ib[i+1] = ishift(ix,-31)
        end
        

        jr = 0
        kr = IP-IQ
        for j=0:IP-1
            iwork = 0
            for i=0:31
                iwork=iwork*2+ib[jr+1]
                #while iwork < -intmax32-1 || intmax32 < iwork
                #    iwork -= Int32(Int32(2)*(intmax32+Int32(1)))
                #end
                iwork = make_value_Int32_again(iwork)

                ib[jr+1] = Int32(ib[jr+1]) ⊻ Int32(ib[kr+1])
                jr += 1
                if jr == IP
                    jr = 0
                end
                kr += 1
                if kr == IP
                    kr = 0
                end
            end
            iw[j+1] = Int32(iwork) >>> -(-1)
            #println(j,"\t",iw[j+1])
        end

        return Crand(iw,jr,kr)

#        exit()

    end

    function delay3!(crand_public::Crand,lambda)
        
        

        iwk = zeros(Int32,2*IP-2+1)
        c = zeros(Int32,2*IP-1+1)
        ib = zeros(Int32,IP+31+1)

        mu = MACRM
        for i=0:IP-1
            iwk[i+1] = crand_public.iw[i+1]
        end

        for i=IP:2*IP-2
            iwk[i+1] = Int32(iwk[i-IP+1]) ⊻ Int32(iwk[i-IQ+1])
        end

        for i=0:mu-1
            ib[i+1] = 0
        end

        m = lambda
        nb = mu-1

        #count = 0
        while m>IP-1
            nb += 1
            ib[nb+1] = m % 2
            m = div(m,2)
            #count += 1
            #if count > 100000000
            #    error("something is wrong in rand.jl!")
            #end
        end

        for i=0:IP-1
            c[i+1] = 0
        end

        c[m+1] = 1


        for j=nb:-1:0
            for i=IP-1:-1:0
                c[2*i+ib[j+1]+1] = c[i+1]
                c[2*i+1-ib[j+1]+1] = 0
            end

            for i=2*IP-1:-1:IP
                c[i-IP+1] = Int32(c[i-IP+1]) ⊻ Int32(c[i+1])
                c[i-IQ+1] = Int32(c[i-IQ+1]) ⊻ Int32(c[i+1])
            end
        end

        for j=0:IP-1
            iwork = 0
            for i=0:IP-1
                #while iwork < -intmax32-1 || intmax32 < iwork
                #    iwork -= Int32(Int32(2)*(intmax32+Int32(1)))
                #end
                iwork = make_value_Int32_again(iwork)
                iwork = Int32(iwork) ⊻ Int32(c[i+1]*iwk[j+i+1])
            end
            crand_public.iw[j+1] = iwork
        end

    end

    function ranf!(crand_public)
        #ranfs = zeros(Float64,1)
        #ccall((:call_random,"./randfortran.so"),Nothing,(Ref{Float64},),ranfs)
        #return ranfs[1]


        ir = crand_public.iw
        j = crand_public.jr
        k = crand_public.kr
        crand_public.iw[j+1] = Int32(crand_public.iw[j+1]) ⊻ Int32(crand_public.iw[k+1])
        irnd = crand_public.iw[j+1]
        #println("irnd = $irnd ", FNORM,"\t",irnd*FNORM)
        #ranf = Float32(irnd*FNORM)
        #granf = Float64(irnd)*FNORM
        frnd = Float32(irnd*FNORM)
        #println(frnd,"\t",irnd*FNORM,"\t",Float64(irnd),"\t",granf)
        #exit()
        ranf = Float32(frnd)

        crand_public.jr = crand_public.jr + 1
        if crand_public.jr >= IP
            crand_public.jr = crand_public.jr - IP
        end
        crand_public.kr = crand_public.kr + 1
        if crand_public.kr >= IP
            crand_public.kr = crand_public.kr - IP
        end
        return ranf
    end

    function rgauss1!(crand_public,nup)
        xr = zeros(Float64,nup)
        
        nup2 = div(nup,2)
        if nup != 2*nup2
            error(
                """
                !!!!  Attention nup in rgauss should be even nup = $nup
                """
            )
        end

        for i=1:nup2
            i2 = 2*i
            i1 = i2 - 1
            v1 = sqrt(-log(ranf!(crand_public)+1e-10)*2)
            v2 = 2pi*ranf!(crand_public)

            xr[i1] = v1*cos(v2)
            xr[i2] = v1*sin(v2)
            
        end

        return xr
    end
end