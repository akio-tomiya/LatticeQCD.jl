mutable struct Topological_charge_measurement{Dim,TG} <: AbstractMeasurement
    filename::String
    _temporary_gaugefields::Vector{TG}
    temp_UμνTA::Matrix{TG}
    Dim::Int8
    #factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool
    TC_methods::Vector{String}

    function Topological_charge_measurement(
        U::Vector{T};
        filename = nothing,
        verbose_level = 2,
        printvalues = true,
        TC_methods = ["plaquette"],
    ) where {T}
        myrank = get_myrank(U)

        if printvalues
            verbose_print = Verbose_print(verbose_level, myid = myrank, filename = filename)
        else
            verbose_print = nothing
        end
        Dim = length(U)

        temp_UμνTA = Array{T,2}(undef, Dim, Dim)

        for μ = 1:Dim
            for ν = 1:Dim
                temp_UμνTA[ν, μ] = similar(U[1])
            end
        end


        numg = 3
        _temporary_gaugefields = Vector{T}(undef, numg)
        _temporary_gaugefields[1] = similar(U[1])
        for i = 2:numg
            _temporary_gaugefields[i] = similar(U[1])
        end



        return new{Dim,T}(
            filename,
            _temporary_gaugefields,
            temp_UμνTA,
            Dim,
            verbose_print,
            printvalues,
            TC_methods,
        )

    end



end

function Topological_charge_measurement(
    U::Vector{T},
    params::TopologicalCharge_parameters,
    filename,
) where {T}

    return Topological_charge_measurement(
        U,
        filename = filename,
        verbose_level = params.verbose_level,
        printvalues = params.printvalues,
        TC_methods = params.kinds_of_topological_charge, #["plaquette"]
    )

end

function measure(
    m::M,
    itrj,
    U::Array{<:AbstractGaugefields{NC,Dim},1};
    additional_string = "",
) where {M<:Topological_charge_measurement,NC,Dim}
    temps = get_temporary_gaugefields(m)
    temp1 = temps[1]
    temp2 = temps[2]
    measurestring = ""

    nummethod = length(m.TC_methods)
    values = Float64[]
    printstring = "$itrj " * additional_string
    for i = 1:nummethod
        methodname = m.TC_methods[i]
        if methodname == "plaquette"
            Qplaq = calculate_topological_charge_plaq(U, m.temp_UμνTA, temps)
            push!(values, real(Qplaq))
        elseif methodname == "clover"
            Qclover = calculate_topological_charge_clover(U, m.temp_UμνTA, temps)
            push!(values, real(Qclover))
            Qimproved =
                calculate_topological_charge_improved(U, m.temp_UμνTA, Qclover, temps)
            push!(values, real(Qimproved))
        else
            error("method $methodname is not supported in topological charge measurement")
        end
        #printstring *= "$(values[i]) "
    end
    for value in values
        printstring *= "$(value) "
    end
    printstring *= "# itrj "

    for i = 1:nummethod
        methodname = m.TC_methods[i]
        if methodname == "plaquette"
            printstring *= "Qplaq "
        elseif methodname == "clover"
            printstring *= "Qclover Qimproved "
        else
            error("method $methodname is not supported in topological charge measurement")
        end
    end

    if m.printvalues
        #println_verbose_level2(U[1],"-----------------")
        measurestring = printstring
        println_verbose_level2(m.verbose_print, printstring)
        #println_verbose_level2(U[1],"-----------------")
    end

    return values, measurestring
end

function calculate_topological_charge_plaq(U::Array{T,1}, temp_UμνTA, temps) where {T}
    UμνTA = temp_UμνTA
    numofloops = calc_UμνTA!(UμνTA, "plaq", U, temps)
    Q = calc_Q(UμνTA, numofloops, U)
    return Q
end

function calculate_topological_charge_clover(U::Array{T,1}, temp_UμνTA, temps) where {T}
    UμνTA = temp_UμνTA
    numofloops = calc_UμνTA!(UμνTA, "clover", U, temps)
    Q = calc_Q(UμνTA, numofloops, U)
    return Q
end

function calculate_topological_charge_improved(
    U::Array{T,1},
    temp_UμνTA,
    Qclover,
    temps,
) where {T}
    UμνTA = temp_UμνTA
    #numofloops = calc_UμνTA!(UμνTA,"clover",U)
    #Qclover = calc_Q(UμνTA,numofloops,U)

    numofloops = calc_UμνTA!(UμνTA, "rect", U, temps)
    Qrect = 2 * calc_Q(UμνTA, numofloops, U)
    c1 = -1 / 12
    c0 = 5 / 3
    Q = c0 * Qclover + c1 * Qrect
    return Q
end

function calc_UμνTA!(
    temp_UμνTA,
    name::String,
    U::Array{<:AbstractGaugefields{NC,Dim},1},
    temps,
) where {NC,Dim}
    loops_μν, numofloops = calc_loopset_μν_name(name, Dim)
    calc_UμνTA!(temp_UμνTA, loops_μν, U, temps)
    return numofloops
end


function calc_UμνTA!(
    temp_UμνTA,
    loops_μν,
    U::Array{<:AbstractGaugefields{NC,Dim},1},
    temps,
) where {NC,Dim}
    UμνTA = temp_UμνTA
    for μ = 1:Dim
        for ν = 1:Dim
            if ν == μ
                continue
            end

            evaluate_gaugelinks!(temps[1], loops_μν[μ, ν], U, temps[2:3])
            Traceless_antihermitian!(UμνTA[μ, ν], temps[1])
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
function calc_Q(UμνTA, numofloops, U::Array{<:AbstractGaugefields{NC,Dim},1}) where {NC,Dim}
    Q = 0.0
    if Dim == 4
        ε(μ, ν, ρ, σ) = epsilon_tensor(μ, ν, ρ, σ)
    else
        error("Dimension $Dim is not supported")
    end
    for μ = 1:Dim
        for ν = 1:Dim
            if ν == μ
                continue
            end
            Uμν = UμνTA[μ, ν]
            for ρ = 1:Dim
                for σ = 1:Dim
                    if ρ == σ
                        continue
                    end
                    Uρσ = UμνTA[ρ, σ]
                    s = tr(Uμν, Uρσ)
                    Q += ε(μ, ν, ρ, σ) * s / numofloops^2
                end
            end
        end
    end

    return -Q / (32 * (π^2))
end




#topological charge
function epsilon_tensor(mu::Int, nu::Int, rho::Int, sigma::Int)
    sign = 1 # (3) 1710.09474 extended epsilon tensor
    if mu < 0
        sign *= -1
        mu = -mu
    end
    if nu < 0
        sign *= -1
        nu = -nu
    end
    if rho < 0
        sign *= -1
        rh = -rho
    end
    if sigma < 0
        sign *= -1
        sigma = -sigma
    end
    epsilon = zeros(Int, 4, 4, 4, 4)
    epsilon[1, 2, 3, 4] = 1
    epsilon[1, 2, 4, 3] = -1
    epsilon[1, 3, 2, 4] = -1
    epsilon[1, 3, 4, 2] = 1
    epsilon[1, 4, 2, 3] = 1
    epsilon[1, 4, 3, 2] = -1
    epsilon[2, 1, 3, 4] = -1
    epsilon[2, 1, 4, 3] = 1
    epsilon[2, 3, 1, 4] = 1
    epsilon[2, 3, 4, 1] = -1
    epsilon[2, 4, 1, 3] = -1
    epsilon[2, 4, 3, 1] = 1
    epsilon[3, 1, 2, 4] = 1
    epsilon[3, 1, 4, 2] = -1
    epsilon[3, 2, 1, 4] = -1
    epsilon[3, 2, 4, 1] = 1
    epsilon[3, 4, 1, 2] = 1
    epsilon[3, 4, 2, 1] = -1
    epsilon[4, 1, 2, 3] = -1
    epsilon[4, 1, 3, 2] = 1
    epsilon[4, 2, 1, 3] = 1
    epsilon[4, 2, 3, 1] = -1
    epsilon[4, 3, 1, 2] = -1
    epsilon[4, 3, 2, 1] = 1
    return epsilon[mu, nu, rho, sigma] * sign
end



function calc_loopset_μν_name(name, Dim)
    loops_μν = Array{Vector{Wilsonline{Dim}},2}(undef, Dim, Dim)
    if name == "plaq"
        numofloops = 1
        for μ = 1:Dim
            for ν = 1:Dim
                loops_μν[μ, ν] = Wilsonline{Dim}[]
                if ν == μ
                    continue
                end
                plaq = make_plaq(μ, ν, Dim = Dim)
                push!(loops_μν[μ, ν], plaq)
            end
        end
    elseif name == "clover"
        numofloops = 4
        for μ = 1:Dim
            for ν = 1:Dim
                loops_μν[μ, ν] = Wilsonline{Dim}[]
                if ν == μ
                    continue
                end
                loops_μν[μ, ν] = make_cloverloops_topo(μ, ν, Dim = Dim)
            end
        end
    elseif name == "rect"
        numofloops = 8
        for μ = 1:4
            for ν = 1:4
                if ν == μ
                    continue
                end
                loops = Wilsonline{Dim}[]
                loop_righttop = Wilsonline([(μ, 2), (ν, 1), (μ, -2), (ν, -1)])
                loop_lefttop = Wilsonline([(ν, 1), (μ, -2), (ν, -1), (μ, 2)])
                loop_rightbottom = Wilsonline([(ν, -1), (μ, 2), (ν, 1), (μ, -2)])
                loop_leftbottom = Wilsonline([(μ, -2), (ν, -1), (μ, 2), (ν, 1)])
                push!(loops, loop_righttop)
                push!(loops, loop_lefttop)
                push!(loops, loop_rightbottom)
                push!(loops, loop_leftbottom)

                loop_righttop = Wilsonline([(μ, 1), (ν, 2), (μ, -1), (ν, -2)])
                loop_lefttop = Wilsonline([(ν, 2), (μ, -1), (ν, -2), (μ, 1)])
                loop_rightbottom = Wilsonline([(ν, -2), (μ, 1), (ν, 2), (μ, -1)])
                loop_leftbottom = Wilsonline([(μ, -1), (ν, -2), (μ, 1), (ν, 2)])
                push!(loops, loop_righttop)
                push!(loops, loop_lefttop)
                push!(loops, loop_rightbottom)
                push!(loops, loop_leftbottom)

                loops_μν[μ, ν] = loops
            end
        end
    else
        error("$name is not supported")
    end
    return loops_μν, numofloops
end


function make_cloverloops_topo(μ, ν; Dim = 4)
    loops = Wilsonline{Dim}[]
    loop_righttop = Wilsonline([(μ, 1), (ν, 1), (μ, -1), (ν, -1)])
    loop_lefttop = Wilsonline([(ν, 1), (μ, -1), (ν, -1), (μ, 1)])
    loop_rightbottom = Wilsonline([(ν, -1), (μ, 1), (ν, 1), (μ, -1)])
    loop_leftbottom = Wilsonline([(μ, -1), (ν, -1), (μ, 1), (ν, 1)])
    push!(loops, loop_righttop)
    push!(loops, loop_lefttop)
    push!(loops, loop_rightbottom)
    push!(loops, loop_leftbottom)
    return loops
end
