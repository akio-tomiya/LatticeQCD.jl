mutable struct Gradientflow{TA,T} <: Abstractsmearing
    Nflow::Int64
    eps::Float64
    _temporal_TA_field::Array{TA,1} 
    _temporal_G_field::Array{T,1}
    _temporal_U_field::Array{Array{T,1},1}

    function Gradientflow(U::Array{T,1},Nflow = 1,eps = 0.01) where T <: AbstractGaugefields 
        F0 = initialize_TA_Gaugefields(U)
        Ftemps = Array{typeof(F0),1}(undef,4)
        Ftemps[1] = F0
        for i=2:4
            Ftemps[i] = initialize_TA_Gaugefields(U)
        end

        Utemps = Array{Array{T,1},1}(undef,2)
        for i=1:2
            Utemps[i] = similar(U)
        end

        tempG = Array{T,1}(undef,3)
        for i=1:3
            tempG[i] = similar(U[1])
        end

        return new{typeof(F0),T}(Nflow,eps,Ftemps,tempG,Utemps)
    end
    
end

function get_tempG(x::T) where T <: Gradientflow
    return x._temporal_G_field
end

function get_eps(x::T) where T <: Gradientflow
    return x.eps
end

function flow!(U,g::T) where T<: Gradientflow
    #@assert Dim == 4 "Dimension should be 4. But Dim = $Dim "
    Ftemps = g._temporal_TA_field
    Utemps = g._temporal_U_field
    temps = g._temporal_G_field
    
    F0 = Ftemps[1]
    F1 = Ftemps[2]
    F2 = Ftemps[3]
    Ftmp = Ftemps[4]

    #Ftmp = similar(U)
    W1 = Utemps[1]
    W2 = Utemps[2]
    temp1 = temps[1]
    temp2 = temps[2]
    temp3 = temps[3]
    eps = g.eps

    for istep=1:g.Nflow #RK4 integrator -> RK3?
        clear_U!(F0)
        add_force!(F0,U,temps,plaqonly = true)

        #add_force!(F0,U,[temp1,temp2,temp3],gparam)
        
        exp_aF_U!(W1,-eps*(1/4),F0,U,[temp1,temp2,temp3]) #exp(a*F)*U

        #println("W1 ",W1[1][1,1,1,1,1,1])
        #
        clear_U!(F1)
        add_force!(F1,W1,[temp1,temp2,temp3],plaqonly = true)
        #add_force!(F1,W1,[temp1,temp2,temp3],gparam) #F
        #println("F1 ",F1[1][1,1,1,1,1,1])
        clear_U!(Ftmp)
        add_U!(Ftmp,-(8/9*eps),F1)
        #println("Ftmp ",Ftmp[1][1,1,1,1,1,1])
        add_U!(Ftmp,(17/36*eps),F0)
        #println("Ftmp1 ",Ftmp[1][1,1,1,1,1,1])
        exp_aF_U!(W2,1,Ftmp,W1,[temp1,temp2,temp3]) #exp(a*F)*U
        #exp_aF_U!(W2,1,Ftmp,U,[temp1,temp2,temp3]) #exp(a*F)*U
        #println("W2 ",W2[1][1,1,1,1,1,1])
        
        #
        clear_U!(F2)
        add_force!(F2,W2,[temp1,temp2,temp3],plaqonly = true)
        #add_force!(F2,W2,[temp1,temp2,temp3],gparam) #F
        #calc_gaugeforce!(F2,W2,univ) #F
        clear_U!(Ftmp)
        add_U!(Ftmp,-(3/4*eps),F2)
        add_U!(Ftmp,(8/9*eps),F1)
        add_U!(Ftmp,-(17/36*eps),F0)
        #exp_aF_U!(W1,1,Ftmp,U,[temp1,temp2,temp3]) #exp(a*F)*U  
        exp_aF_U!(U,1,Ftmp,W2,[temp1,temp2,temp3]) #exp(a*F)*U  
        
        #println(typeof(U[1]))
        #println(U[1][1,1,1,1,1,1])
       
        #error("U")
    end
    
end