abstract type TA_Gaugefields_4D{NC} <: TA_Gaugefields{NC,4}
end

include("./TA_gaugefields_4D_serial.jl")

function TA_Gaugefields(NC,NX,NY,NZ,NT;mpi=false)
    if mpi
        error("mpi = $mpi is not supoorted")
    else
        return TA_Gaugefields_4D_serial(NC,NX,NY,NZ,NT)
    end 
end

function clear_U!(U::Array{T,1}) where T <: TA_Gaugefields_4D
    for μ=1:4
        clear_U!(U[μ])
    end
end

