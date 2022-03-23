struct Heatbathupdate{Dim,T} <: AbstractUpdate
    heatbath::Union{Heatbath_update{Dim,T},Heatbath{T}}
    isevenodd::Bool
    numOR::Int64
    useOR::Bool
end

function Heatbathupdate(U::Array{T,1},gauge_action,quench;
        useOR = false,numOR=0,isevenodd=false,β=2.3,ITERATION_MAX=10^5) where T
    @assert quench "Heatbath update is only for quanch case!"
    Dim = length(U)
    if isevenodd 
        println_verbose_level1(U[1],"Heatbath update with even-odd method. Only plaquette action can be treated. Now β = $β")
        
        hb = Heatbath(U,β,ITERATION_MAX=ITERATION_MAX) 
    else
        hb = Heatbath_update(U,gauge_action,ITERATION_MAX=ITERATION_MAX)
    end

    return Heatbathupdate{Dim,T}(hb,isevenodd,numOR,useOR)
#error("in StandardHMC!!")
end

function update!(updatemethod::T,U) where T <: Heatbathupdate
    heatbath!(U,updatemethod.heatbath)
    if updatemethod.useOR
        for i=1:updatemethod.numOR
            overrelaxation!(U,updatemethod.heatbath)
        end
    end
    #error("updatemethod type $(updatemethod) is not supported!!")
end
