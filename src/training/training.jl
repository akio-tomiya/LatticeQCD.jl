module Training
    abstract type TrainableWeights
    end

    ```
    Adam: A Method for Stochastic Optimization
    Diederik P. Kingma, Jimmy Ba, arXiv:1412.6980
    ```
    mutable struct TrainableWeights_ADAM <: TrainableWeights
        ρs::Array{Array{Float64,1},1}
        ms::Array{Array{Float64,1},1}
        vs::Array{Array{Float64,1},1}
        t::Int128
        η::Float64
        ε::Float64
        β1::Float64
        β2::Float64

        function TrainableWeights_ADAM(ρs::Array{Array{Float64,1},1};η = 0.0005,ε=10^(-8),β1=0.9,β2=0.999)
            ms = deepcopy(ρs)
            vs = deepcopy(ρs)
            for i=1:length(ms)
                ms[i] .= 0
                vs[i] .= 0
            end
            t = 1
            return new(ρs,ms,vs,t,η,ε,β1,β2)
        end
    end

    mutable struct TrainableWeights_SGD <: TrainableWeights
        ρs::Array{Array{Float64,1},1}
        η::Float64

        function TrainableWeights_SGD(ρs::Array{Array{Float64,1},1};η = 0.0005)
            return new(ρs,η)
        end
    end

    function train!(tw::T,g) where T <: TrainableWeights_SGD
        for i=1:length(tw.ρs)
            for j=1:length(tw.ρs[i])
                tw.ρs[i][j] = tw.ρs[i][j] - tw.η*g[i][j]
            end
        end
    end

    ```
    training with ADAM.
    Adam: A Method for Stochastic Optimization
    Diederik P. Kingma, Jimmy Ba, arXiv:1412.6980
    ```
    function train!(tw::T,g) where T <: TrainableWeights_ADAM
        t = tw.t
        for i=1:length(tw.ms)
            for j=1:length(tw.ms[i])
                tw.ms[i][j] = tw.β1*tw.ms[i][j] + (1-tw.β1)*g[i][j]
            end
        end
        #tw.ms = β1*θ.m + (1-β1)*g
        for i=1:length(tw.vs)
            for j=1:length(tw.vs[i])
                tw.vs[i][j] = tw.β2*tw.vs[i][j] + (1-tw.β2)*(g[i][j]^2)
            end
        end
        #θ.v = β2*θ.v + (1-β2)*(g^2)
        # correction
        for i=1:length(tw.ms)
            for j=1:length(tw.ms[i])
                mhat = tw.ms[i][j]/(1-tw.β1^tw.t) 
                vhat = tw.vs[i][j]/(1-tw.β2^tw.t)
                #println(tw.η*mhat/(sqrt(vhat)+tw.ε))
                #println("before ",tw.ρs[i][j])
                tw.ρs[i][j] = tw.ρs[i][j] - tw.η*mhat/(sqrt(vhat)+tw.ε)
                #println("after ",tw.ρs[i][j])
            end
        end
        #
        tw.t+=1
    end
end