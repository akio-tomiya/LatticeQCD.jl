mutable struct Energy_density_measurement{Dim,TG} <: AbstractMeasurement
    filename::String
    _temporary_gaugefields::Vector{TG}
    temp_UμνTA::Matrix{TG}
    Dim::Int8
    #factor::Float64
    verbose_print::Union{Verbose_print,Nothing}
    printvalues::Bool

    function Energy_density_measurement(U::Vector{T};
            filename = nothing,
            verbose_level = 2,
            printvalues = true) where T
        myrank = get_myrank(U)

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
    

        numg = 3
        _temporary_gaugefields = Vector{T}(undef,numg)
        _temporary_gaugefields[1] = similar(U[1])
        for i=2:numg
            _temporary_gaugefields[i] = similar(U[1])
        end

        return new{Dim,T}(filename,_temporary_gaugefields,temp_UμνTA,
                        Dim,verbose_print,printvalues)

    end



end


function measure(m::M,itrj,U::Array{<: AbstractGaugefields{NC,Dim},1};additional_string = "") where {M <: Energy_density_measurement,NC,Dim}
    temps = get_temporary_gaugefields(m)
    value = calculate_energy_density(U,m.temp_UμνTA,temps)

    if m.printvalues
        println_verbose_level2(U[1],"-----------------")
        println_verbose_level2(m.verbose_print,"$itrj $additional_string $value # energydensity")
        println_verbose_level2(U[1],"-----------------")
    end

    return values
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