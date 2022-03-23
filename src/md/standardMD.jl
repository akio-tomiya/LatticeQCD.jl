struct StandardMD{Dim,TG,TA,quench,T_FA,TF} <: AbstractMD{Dim,TG}
    gauge_action::GaugeAction{Dim,TG}
    quench::Bool
    Δτ::Float64
    MDsteps::Float64
    p::Vector{TA}
    QPQ::Bool
    fermi_action::T_FA
    η::TF
    ξ::TF
    SextonWeingargten::Bool
    Nsw::Int64

    function StandardMD(U,gauge_action::GaugeAction{Dim,TG},quench,Δτ,MDsteps,fermi_action = nothing;QPQ=true,SextonWeingargten=false,Nsw=2) where {Dim,TG}
        p = initialize_TA_Gaugefields(U) #This is a traceless-antihermitian gauge fields. This has NC^2-1 real coefficients. 
        TA = eltype(p)
        T_FA = typeof(fermi_action)

        if quench 
            η = nothing
            ξ = nothing
            if SextonWeingargten 
                error("The quench update does not need the SextonWeingargten method. Put SextonWeingargten = false")
            end
        else
            if fermi_action == nothing
                η = nothing
                ξ = nothing
            else
                η = similar(fermi_action._temporary_fermionfields[1])
                ξ = similar(η)
            end
        end
        TF = typeof(η)

        @assert Nsw % 2 == 0 "Nsw should be even number! now Nsw = $Nsw"

        return new{Dim,TG,TA,quench,T_FA,TF}(gauge_action,quench,Δτ,MDsteps,p,QPQ,fermi_action,η,ξ,SextonWeingargten,Nsw)
    end
end

function initialize_MD!(U,md::StandardMD{Dim,TG,TA,quench,T_FA}) where {Dim,TG,TA,quench,T_FA}
    gauss_distribution!(md.p) #initial momentum

    if quench == false
        gauss_sampling_in_action!(md.ξ,U,md.fermi_action)
        sample_pseudofermions!(md.η,U,md.fermi_action,md.ξ)
        #error("not supported yet")
    end
end

function runMD!(U,md::StandardMD{Dim,TG,TA,quench}) where {Dim,TG,TA,quench}
    #p = md.p

    if md.QPQ
        if md.SextonWeingargten
            runMD_QPQ_sw!(U,md)
        else
            runMD_QPQ!(U,md)
        end
    else
        if md.SextonWeingargten
            error("PQP update with SextonWeingargten is not supported")
        else
            runMD_PQP!(U,md)
        end
    end

    #error("type $(typeof(md)) is not supported")
end

function runMD_QPQ!(U,md::StandardMD{Dim,TG,TA,quench}) where {Dim,TG,TA,quench}
    p = md.p

    for itrj=1:md.MDsteps
        U_update!(U,p,0.5,md)
        P_update!(U,p,1.0,md)
        if quench == false
            P_update_fermion!(U,p,1.0,md)
        end
        U_update!(U,p,0.5,md)
    end

    #error("type $(typeof(md)) is not supported")
end

function runMD_QPQ_sw!(U,md::StandardMD{Dim,TG,TA,quench}) where {Dim,TG,TA,quench}
    p = md.p

    for itrj=1:md.MDsteps
        for isw=1:div(md.Nsw,2)
            U_update!(U,p,0.5/md.Nsw,md)
            P_update!(U,p,1.0/md.Nsw,md)
            U_update!(U,p,0.5/md.Nsw,md)
        end
        if quench == false
            P_update_fermion!(U,p,1.0,md)
        end
        for isw=1:div(md.Nsw,2)
            U_update!(U,p,0.5/md.Nsw,md)
            P_update!(U,p,1.0/md.Nsw,md)
            U_update!(U,p,0.5/md.Nsw,md)
        end
    end

    #error("type $(typeof(md)) is not supported")
end


function runMD_PQP!(U,md::StandardMD{Dim,TG,TA,quench}) where {Dim,TG,TA,quench}
    p = md.p

    for itrj=1:md.MDsteps
        P_update!(U,p,0.5,md)
        if quench == false
            P_update_fermion!(U,p,0.5,md)
        end

        U_update!(U,p,1.0,md)
        P_update!(U,p,0.5,md)
        if quench == false
            P_update_fermion!(U,p,0.5,md)
        end
    end

    #error("type $(typeof(md)) is not supported")
end