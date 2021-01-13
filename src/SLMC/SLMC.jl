module SLMC
    import ..Actions:GaugeActionParam_standard,
                    GaugeActionParam,
                    GaugeActionParam_autogenerator

    using LinearAlgebra
    mutable struct SLMC_data
        KdagK::Array{Float64,2}
        KdagY::Array{Float64,1}
        factor::Float64


        function SLMC_data(n,NC)
            KdagK = zeros(Float64,n+1,n+1)
            KdagY= zeros(Float64,n+1)
            factor = -1/NC
            return new(KdagK,KdagY,factor)
        end
    end
    


    function update_slmcdata!(s::SLMC_data,inputdata,outputdata)
        n = length(s.KdagY)
        KdagK = zero(s.KdagK)
        KdagY = zero(s.KdagY)
        S = outputdata
        
        K = zeros(Float64,n)
        
        K[1] = 1
        for i=2:n
            K[i] = s.factor*inputdata[i-1]
        end
        println("K = ",K)

        s.KdagY += K[:]*S #KdagK 
        s.KdagK += K*K' #KdagY 

        #numdata = s.KdagK[1,1]
        #s.KdagY ./= numdata
        #s.KdagK ./= numdata
        return
    end

    function make_beta(s::SLMC_data,num)    
        #if num == 1
            #detKdagK = det(s.KdagK[1:num+1,1:num+1])
            #println("det",detKdagK)
            #invKdagK = [s.KdagK[2,2] -s.KdagK[1,2]
            #s.KdagK[2,1] s.KdagK[1,1] 
            #]/detKdagK
            #println("invKdagK ",invKdagK)
            #println("KdagY ",s.KdagY[1:num+1] )
            
        #    slmc_betas = invKdagK*s.KdagY[1:num+1]            
        #else       
            slmc_betas = s.KdagK[1:num+1,1:num+1] \ s.KdagY[1:num+1]
        #end
        return slmc_betas
    end

    function show_effbeta(s::SLMC_data,effectiveterms)
        #effectiveterms = ["plaq"]
        n = length(s.KdagY)
        be1=zeros(n)
        IsSucs=true

        try
            for i=1:n-1
                print("#Estimation$i learning ")
                print("$(effectiveterms[1:i]) only: ")
                be1[1:i+1] = make_beta(s,i)
                #println("slmc_betas ",be1 )
                for k = 1:length(be1)-1
                    print("$(be1[k+1]) \t")
                end
                println("constant is ",be1[1])
            end
        catch
            println("Failed to calculate on-fly SLMC parameter.")
            println(s.KdagK,"\t",s.KdagY)            
            
            IsSucs=false
        end
        return be1[2:end],be1[1],IsSucs
    end

    function show_effbeta(s::SLMC_data)
        effectiveterms = ["plaq"]
        return show_effbeta(s,effectiveterms)
    end

    function show_effbeta(s::SLMC_data,gparam::GaugeActionParam)
        
        if typeof(gparam) == GaugeActionParam_autogenerator
            effectiveterms = gparam.couplinglist
            #println(gparam.couplinglist )
        elseif typeof(gparam) == GaugeActionParam_standard
            effectiveterms = ["plaq"]
        else
            error("error!")
        end
        return show_effbeta(s,effectiveterms)
    end


end