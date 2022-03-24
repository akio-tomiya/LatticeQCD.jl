mutable struct Topological_charge_measurement{Dim,TG,S} <: AbstractMeasurement
    filename::String
    _temporary_gaugefields::Vector{TG}
    temp_UμνTA::Matrix{TG}
    Dim::Int8
    #factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    Usmr::Array{TG,1}
    smearing::S
    smearing_type::String
    numflow::Int64

    function Topological_charge_measurement(U::Vector{T};
            Nflowsteps = 1,
            eps_flow = 0.01,
            smearing_type = "gradient_flow",
            numflow = 10,
            filename = nothing,
            verbose_level = 2,
            printvalues = true) where T
        myrank = get_myrank(U)

        Usmr = deepcopy(U)


        if printvalues
            verbose_print = Verbose_print(verbose_level,myid = myrank,filename=filename)
        else
            verbose_print = nothing
        end
        Dim = length(U)

        temp_UμνTA = Array{T,2}(undef,Dim,Dim)

        for μ=1:Dim
            for ν=1:Dim
                temp_UμνTA[ν,μ] = similar(U[1])
            end
        end

        if smearing_type == "gradient_flow"
            smearing = Gradientflow(U,Nflow = Nflowsteps,eps = eps_flow)
            #smearing = Gradientflow(U,Nflow = Nflowsteps,eps = eps_flow)
            S = typeof(smearing)
        else
            error("smearing_type = $smearing_type is not supported")
        end

    

        numg = 3
        _temporary_gaugefields = Vector{T}(undef,numg)
        _temporary_gaugefields[1] = similar(U[1])
        for i=2:numg
            _temporary_gaugefields[i] = similar(U[1])
        end

        return new{Dim,T,S}(filename,_temporary_gaugefields,temp_UμνTA,
                        Dim,verbose_print,printvalues,
                        Usmr,
                        smearing,smearing_type,
                        numflow)

    end



end

function measure(m::M,itrj,U::Array{<: AbstractGaugefields{NC,Dim},1}) where {M <: Topological_charge_measurement,NC,Dim}
    temps = get_temporary_gaugefields(m)
    temp1 = temps[1]
    temp2 = temps[2]
    eps_flow = m.smearing.eps
    #Nflowsteps = m.smearing.Nflow

    if m.printvalues
        println_verbose_level2(U[1],"-----------------")
        println_verbose_level2(U[1],"# epsilon for the Wilson flow is $(eps_flow)")
    end

    substitute_U!(m.Usmr,U)  
    
    τ = 0.0
       
    if Dim == 4
        comb = 6 #4*3/2
    elseif Dim == 3
        comb = 3
    elseif Dim == 2
        comb = 1
    else
        error("dimension $Dim is not supported")
    end
    factor = 1/(comb*U[1].NV*NC)

    plaq =  calculate_Plaquette(U,temp1,temp2)*factor
    Qplaq = calculate_topological_charge_plaq(m.Usmr,m.temp_UμνTA,temps)
    #println("Qplaq = ",Qplaq)
    Qclover= calculate_topological_charge_clover(m.Usmr,m.temp_UμνTA,temps)
    Qimproved= calculate_topological_charge_improved(m.Usmr,m.temp_UμνTA,Qclover,temps)
    #println("Qimproved = ",Qimproved)
    clov = calculate_energy_density(m.Usmr,m.temp_UμνTA,temps)

    if m.printvalues
        println_verbose_level2(U[1],"$itrj $τ $plaq $clov $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #flow itrj flowtime plaq E Qplaq Qclov Qimproved")
        #println(m.fp,"$itrj $τ $plaq $clov $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #flow itrj flowtime plaq E Qplaq Qclov Qimproved")
        flush(stdout)
    end


    if m.smearing_type == "gradient_flow"
        for iflow = 1:m.numflow#5000 # eps=0.01: t_flow = 50
            flow!(m.Usmr,m.smearing)
            plaq = calculate_Plaquette(m.Usmr,temp1,temp2)*factor
            Qplaq = calculate_topological_charge_plaq(m.Usmr,m.temp_UμνTA,temps)
            Qclover= calculate_topological_charge_clover(m.Usmr,m.temp_UμνTA,temps)
            Qimproved= calculate_topological_charge_improved(m.Usmr,m.temp_UμνTA,Qclover,temps)
            clov = calculate_energy_density(m.Usmr,m.temp_UμνTA,temps)
            #@time Q = calc_topological_charge(Usmr)
            τ = iflow*eps_flow*m.smearing.Nflow
            if m.printvalues
                println_verbose_level2(m.verbose_print,"$itrj $(round(τ, digits=3)) $plaq $clov $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #flow itrj flowtime plaq E Qplaq Qclov Qimproved")
                #println(m.fp,"$itrj $(round(τ, digits=3)) $plaq $clov $(real(Qplaq)) $(real(Qclover)) $(real(Qimproved)) #flow itrj flowtime plaq E Qplaq Qclov Qimproved")
                #if iflow%10 == 0
                flush(stdout)
            end
        
        end
    else
        error("$(m.smearig_type) is not suppoorted")
    end


    return Qclover
end

function calculate_topological_charge_plaq(U::Array{T,1},temp_UμνTA,temps) where T
    UμνTA = temp_UμνTA
    numofloops = calc_UμνTA!(UμνTA,"plaq",U,temps)
    Q = calc_Q(UμνTA,numofloops,U)
    return Q
end

function calculate_topological_charge_clover(U::Array{T,1},temp_UμνTA,temps) where T 
    UμνTA = temp_UμνTA
    numofloops = calc_UμνTA!(UμνTA,"clover",U,temps)
    Q = calc_Q(UμνTA,numofloops,U)
    return Q
end

function calculate_topological_charge_improved(U::Array{T,1},temp_UμνTA,Qclover,temps) where T 
    UμνTA = temp_UμνTA
    #numofloops = calc_UμνTA!(UμνTA,"clover",U)
    #Qclover = calc_Q(UμνTA,numofloops,U)

    numofloops = calc_UμνTA!(UμνTA,"rect",U,temps)
    Qrect = 2*calc_Q(UμνTA,numofloops,U)
    c1 = -1/12
    c0 = 5/3
    Q = c0*Qclover + c1*Qrect
    return Q
end

function calc_UμνTA!(temp_UμνTA,name::String,U::Array{<: AbstractGaugefields{NC,Dim},1},temps) where {NC,Dim}
    loops_μν,numofloops = calc_loopset_μν_name(name,Dim)
    calc_UμνTA!(temp_UμνTA,loops_μν,U,temps)
    return numofloops
end


function calc_UμνTA!(temp_UμνTA,loops_μν,U::Array{<: AbstractGaugefields{NC,Dim},1},temps) where {NC,Dim}
    UμνTA = temp_UμνTA
    for μ=1:Dim
        for ν=1:Dim
            if ν == μ
                continue
            end

            evaluate_gaugelinks!(temps[1],loops_μν[μ,ν],U,temps[2:3])
            Traceless_antihermitian!(UμνTA[μ,ν],temps[1])
            #loopset = Loops(U,loops_μν[μ,ν])
            #UμνTA[μ,ν] = evaluate_loops(loopset,U)

            #UμνTA[μ,ν] = Traceless_antihermitian(UμνTA[μ,ν])
        end
    end
    return 
end


        #=
    implementation of topological charge is based on
    https://arxiv.org/abs/1509.04259
    =#
function calc_Q(UμνTA,numofloops,U::Array{<: AbstractGaugefields{NC,Dim},1}) where {NC,Dim}
    Q = 0.0
    if Dim == 4
        ε(μ,ν,ρ,σ) = epsilon_tensor(μ,ν,ρ,σ)  
    else
        error("Dimension $Dim is not supported")
    end
    for μ=1:Dim
        for ν=1:Dim
            if ν == μ
                continue
            end
            Uμν = UμνTA[μ,ν]                 
            for ρ =1:Dim
                for σ=1:Dim
                    if ρ == σ
                        continue
                    end
                    Uρσ = UμνTA[ρ,σ]
                    s = tr(Uμν,Uρσ)
                    Q += ε(μ,ν,ρ,σ)*s/numofloops^2
                end
            end
        end
    end

    return -Q/(32*(π^2))
end




#topological charge
function epsilon_tensor(mu::Int,nu::Int,rho::Int,sigma::Int) 
    sign=1 # (3) 1710.09474 extended epsilon tensor
    if mu < 0
        sign*=-1
        mu=-mu
    end
    if nu < 0
        sign*=-1
        nu=-nu
    end
    if rho < 0
        sign*=-1
        rh=-rho
    end
    if sigma < 0
        sign*=-1
        sigma=-sigma
    end
    epsilon = zeros(Int,4,4,4,4)
    epsilon[ 1, 2, 3, 4 ] = 1
    epsilon[ 1, 2, 4, 3 ] = -1
    epsilon[ 1, 3, 2, 4 ] = -1
    epsilon[ 1, 3, 4, 2 ] = 1
    epsilon[ 1, 4, 2, 3 ] = 1
    epsilon[ 1, 4, 3, 2 ] = -1
    epsilon[ 2, 1, 3, 4 ] = -1
    epsilon[ 2, 1, 4, 3 ] = 1
    epsilon[ 2, 3, 1, 4 ] = 1
    epsilon[ 2, 3, 4, 1 ] = -1
    epsilon[ 2, 4, 1, 3 ] = -1
    epsilon[ 2, 4, 3, 1 ] = 1
    epsilon[ 3, 1, 2, 4 ] = 1
    epsilon[ 3, 1, 4, 2 ] = -1
    epsilon[ 3, 2, 1, 4 ] = -1
    epsilon[ 3, 2, 4, 1 ] = 1
    epsilon[ 3, 4, 1, 2 ] = 1
    epsilon[ 3, 4, 2, 1 ] = -1
    epsilon[ 4, 1, 2, 3 ] = -1
    epsilon[ 4, 1, 3, 2 ] = 1
    epsilon[ 4, 2, 1, 3 ] = 1
    epsilon[ 4, 2, 3, 1 ] = -1
    epsilon[ 4, 3, 1, 2 ] = -1
    epsilon[ 4, 3, 2, 1 ] = 1
    return epsilon[mu,nu,rho,sigma]*sign
end


       # = = = calc energy density = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
function calculate_energy_density(U::Array{T,1}, Wmat,temps) where T <: AbstractGaugefields
# Making a ( Ls × Lt) Wilson loop operator for potential calculations
    WL = 0.0+0.0im
    NV = U[1].NV
    NC = U[1].NC
    #Wmat = Array{T,2}(undef,4,4)
    #
    make_energy_density!(Wmat,U,temps) # make wilon loop operator and evaluate as a field, not traced.
    WL =  make_energy_density_core(Wmat,U,NV) # tracing over color and average over spacetime and x,y,z.
    NDir = 4.0*3.0/2 # choice of 2 axis from 4.
    return real(WL)/NV/NDir/NC/8
end

function make_energy_density!(Wmat,U::Array{<: AbstractGaugefields{NC,Dim},1},temps) where {NC,Dim}
    W_operator,numofloops = calc_loopset_μν_name("clover",Dim)#make_Wilson_loop(Lt,Ls)
    calc_large_wilson_loop!(Wmat,W_operator,U,temps)
    return 
end

function  make_energy_density_core(Wmat::Array{<: AbstractGaugefields{NC,Dim},2}, U::Array{T,1} ,NV) where {T <: AbstractGaugefields,NC,Dim}
    @assert Dim == 4

    W = 0.0 + 0.0im
    for μ=1:Dim # all directions
        for ν=1:Dim
            if μ == ν
                continue
            end
            W += tr(Wmat[μ,ν],Wmat[μ,ν])/4
        end
    end
    return W
end

function calc_loopset_μν_name(name,Dim)
    loops_μν= Array{Vector{Wilsonline{Dim}},2}(undef,Dim,Dim)
    if name == "plaq"
        numofloops = 1
        for μ=1:Dim
            for ν=1:Dim
                loops_μν[μ,ν] = Wilsonline{Dim}[]
                if ν == μ
                    continue
                end
                plaq = make_plaq(μ,ν,Dim=Dim)
                push!(loops_μν[μ,ν],plaq)
            end
        end
    elseif name == "clover"
        numofloops = 4
        for μ=1:Dim
            for ν=1:Dim
                loops_μν[μ,ν] = Wilsonline{Dim}[]
                if ν == μ
                    continue
                end
                loops_μν[μ,ν] = make_cloverloops(μ,ν,Dim=Dim)
            end
        end
    elseif name == "rect"
        numofloops = 8
        for μ=1:4
            for ν=1:4
                if ν == μ
                    continue
                end
                loops = Wilsonline{Dim}[]
                loop_righttop = Wilsonline([(μ,2),(ν,1),(μ,-2),(ν,-1)])
                loop_lefttop = Wilsonline([(ν,1),(μ,-2),(ν,-1),(μ,2)])
                loop_rightbottom = Wilsonline([(ν,-1),(μ,2),(ν,1),(μ,-2)])
                loop_leftbottom= Wilsonline([(μ,-2),(ν,-1),(μ,2),(ν,1)])
                push!(loops,loop_righttop)
                push!(loops,loop_lefttop)
                push!(loops,loop_rightbottom)
                push!(loops,loop_leftbottom)

                loop_righttop = Wilsonline([(μ,1),(ν,2),(μ,-1),(ν,-2)])
                loop_lefttop = Wilsonline([(ν,2),(μ,-1),(ν,-2),(μ,1)])
                loop_rightbottom = Wilsonline([(ν,-2),(μ,1),(ν,2),(μ,-1)])
                loop_leftbottom= Wilsonline([(μ,-1),(ν,-2),(μ,1),(ν,2)])
                push!(loops,loop_righttop)
                push!(loops,loop_lefttop)
                push!(loops,loop_rightbottom)
                push!(loops,loop_leftbottom)

                loops_μν[μ,ν] = loops
            end
        end
    else
        error("$name is not supported")
    end
    return loops_μν,numofloops
end


function make_cloverloops(μ,ν;Dim=4)
    loops = Wilsonline{Dim}[]
    loop_righttop = Wilsonline([(μ,1),(ν,1),(μ,-1),(ν,-1)])
    loop_lefttop = Wilsonline([(ν,1),(μ,-1),(ν,-1),(μ,1)])
    loop_rightbottom = Wilsonline([(ν,-1),(μ,1),(ν,1),(μ,-1)])
    loop_leftbottom= Wilsonline([(μ,-1),(ν,-1),(μ,1),(ν,1)])
    push!(loops,loop_righttop)
    push!(loops,loop_lefttop)
    push!(loops,loop_rightbottom)
    push!(loops,loop_leftbottom)
    return loops
end

function calc_large_wilson_loop!(temp_Wmat::Array{<: AbstractGaugefields{NC,Dim},2},W_operator,U::Array{T,1},temps) where {T <: AbstractGaugefields,NC,Dim}
    W = temp_Wmat
    for μ=1:Dim
        for ν=1:Dim
            if μ == ν
                continue
            end
            #println(typeof(μ)," ",typeof(ν))
            #exit()
            #loopset = Loops(U,W_operator[μ,ν])
            evaluate_gaugelinks!(W[μ,ν],W_operator[μ,ν],U,temps)
            #W[μ,ν] = evaluate_loops(loopset,U)
        end
    end
    return 
end