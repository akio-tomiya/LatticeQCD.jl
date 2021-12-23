#abstract type StaggeredFermion{NC,Dim} <: AbstractFermionfields{NC,Dim} 
#end 

include("./StaggeredFermion_4D_wing.jl")

function StaggeredFermion(params,NC,NN...)
    Dim = length(NN)
    if Dim == 4
        fermion = StaggeredFermion_4D_wing(params,NC,NN...)
    else
        error("Dimension $Dim is not supported")
    end
    return fermion
end

function StaggeredFermion(fparam::FermiActionParam,NC,BoundaryCondition,NN...)
    Dim = length(NN)
    if Dim == 4
        fermion = StaggeredFermion_4D_wing(NC,NN...,fparam.mass,fparam.eps,
            fparam.MaxCGstep,BoundaryCondition)
        #fermion = StaggeredFermion_4D_wing(params,NC,NN...)
    else
        error("Dimension $Dim is not supported")
    end
    return fermion
end