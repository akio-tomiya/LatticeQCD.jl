module StaggeredFermion_module
using LinearAlgebra

import ..Actions:
    FermiActionParam,
    FermiActionParam_Wilson,
    FermiActionParam_WilsonClover,
    FermiActionParam_Staggered
import ..Gaugefields:
    GaugeFields,
    GaugeFields_1d,
    SU3GaugeFields,
    SU2GaugeFields,
    SU3GaugeFields_1d,
    SU2GaugeFields_1d,
    staggered_phase,
    SUn,
    SU2,
    SU3,
    SUNGaugeFields,
    SUNGaugeFields_1d,
    SU


#using ..Fermion
import ..AbstractFermion:
    FermionFields,
    Wx!,
    Wdagx!,
    clear!,
    substitute_fermion!,
    Dx!,
    Dxplus!,
    set_wing_fermi!,
    fermion_shift!,
    fermion_shiftB!,
    add!,
    WdagWx!,
    apply_periodicity


struct StaggeredFermion <: FermionFields
    NC::Int64
    NX::Int64
    NY::Int64
    NZ::Int64
    NT::Int64
    f::Array{ComplexF64,6}

    mass::Float64
    eps::Float64
    Dirac_operator::String
    MaxCGstep::Int64
    BoundaryCondition::Array{Int8,1}
end

function StaggeredFermion(NC, NX, NY, NZ, NT, fparam::FermiActionParam, BoundaryCondition)
    mass = fparam.mass
    eps = fparam.eps
    MaxCGstep = fparam.MaxCGstep
    return StaggeredFermion(NC, NX, NY, NZ, NT, mass, eps, MaxCGstep, BoundaryCondition)
end



function StaggeredFermion(NC, NX, NY, NZ, NT, mass, eps, MaxCGstep, BoundaryCondition)#r,hop,eps,MaxCGstep)
    Dirac_operator = "Staggered"
    return StaggeredFermion(
        NC,
        NX,
        NY,
        NZ,
        NT,
        zeros(ComplexF64, NC, NX + 2, NY + 2, NZ + 2, NT + 2, 1),
        mass,
        eps,
        Dirac_operator,
        MaxCGstep,
        BoundaryCondition,
    )
end



function StaggeredFermion(
    NC,
    NX,
    NY,
    NZ,
    NT,
    mass,
    eps,
    fermion,
    MaxCGstep,
    BoundaryCondition,
)
    return StaggeredFermion(
        NC,
        NX,
        NY,
        NZ,
        NT,
        zeros(ComplexF64, NC, NX + 2, NY + 2, NZ + 2, NT + 2, 1),
        mass,
        eps,
        fermion,
        MaxCGstep,
        BoundaryCondition,
    )
end

function Base.:*(a::StaggeredFermion, b::StaggeredFermion)
    c = 0.0im
    α = 1
    for it = 1:a.NT
        for iz = 1:a.NZ
            for iy = 1:a.NY
                for ix = 1:a.NX
                    @simd for ic = 1:a.NC
                        c += conj(a[ic, ix, iy, iz, it, α]) * b[ic, ix, iy, iz, it, α]
                    end
                end
            end
        end
    end

    return c
end

function Base.similar(x::StaggeredFermion)
    return StaggeredFermion(
        x.NC,
        x.NX,
        x.NY,
        x.NZ,
        x.NT,
        x.mass,
        x.eps,
        x.Dirac_operator,
        x.MaxCGstep,
        x.BoundaryCondition,
    )
end

function clear!(x::StaggeredFermion, evensite)
    ibush = ifelse(evensite, 0, 1)
    for it = 1:x.NT
        for iz = 1:x.NZ
            for iy = 1:x.NY
                xran = 1+(1+ibush+iy+iz+it)%2:2:x.NX
                for ix in xran
                    @simd for ic = 1:x.NC
                        x[ic, ix, iy, iz, it, 1] = 0
                    end
                end
            end
        end
    end
    return
end

function substitute_fermion!(H, j, x::StaggeredFermion)
    i = 0
    for ialpha = 1:1
        for it = 1:x.NT
            for iz = 1:x.NZ
                for iy = 1:x.NY
                    for ix = 1:x.NX
                        @simd for ic = 1:x.NC
                            i += 1
                            H[i, j] = x[ic, ix, iy, iz, it, ialpha]
                        end
                    end
                end
            end
        end
    end
end

function substitute_fermion!(y::Vector, x::StaggeredFermion)
    i = 0
    for ialpha = 1:1
        for it = 1:x.NT
            for iz = 1:x.NZ
                for iy = 1:x.NY
                    for ix = 1:x.NX
                        @simd for ic = 1:x.NC
                            i += 1
                            y[i] = x[ic, ix, iy, iz, it, ialpha]
                        end
                    end
                end
            end
        end
    end
end


"""
(D + m)*x
"""
function Wx!(
    xout::StaggeredFermion,
    U::Array{G,1},
    x::StaggeredFermion,
    temps::Array{T,1},
    fparam::FermiActionParam_Staggered,
) where {T<:FermionFields,G<:GaugeFields}

    temp = temps[4]
    Dx!(temp, U, x, temps, fparam)
    clear!(xout)
    add!(xout, fparam.mass, x, 1, temp)
    set_wing_fermi!(xout)
    return

end

function Wx!(
    xout::StaggeredFermion,
    U::Array{G,1},
    x::StaggeredFermion,
    temps::Array{T,1},
    fparam::FermiActionParam_Staggered,
    indices,
) where {T<:FermionFields,G<:GaugeFields}

    temp = temps[4]
    vec_indices1, vec_indices2 = Dx!(temp, U, x, temps, indices)
    clear!(xout)
    add!(xout, fparam.mass, x, 1, temp)
    set_wing_fermi!(xout)
    return vec_indices1, vec_indices2

end

function Wx!(
    xout::StaggeredFermion,
    evensite,
    U::Array{G,1},
    x::StaggeredFermion,
    temps::Array{T,1},
    fparam::FermiActionParam_Staggered,
) where {T<:FermionFields,G<:GaugeFields}


    temp = temps[4]
    Dx!(temp, evensite, U, x, temps, fparam)
    clear!(xout)
    add!(xout, fparam.mass, x, 1, temp)
    set_wing_fermi!(xout)
    return

end

function Dx!(
    xout::StaggeredFermion,
    U::Array{G,1},
    x::StaggeredFermion,
    temps::Array{T,1},
) where {T<:FermionFields,G<:GaugeFields}
    #temp = temps[4]
    temp1 = temps[1]
    temp2 = temps[2]

    #clear!(temp)
    set_wing_fermi!(x)
    clear!(xout)
    for ν = 1:4

        fermion_shift!(temp1, U, ν, x)

        fermion_shift!(temp2, U, -ν, x)

        add!(xout, 0.5, temp1, -0.5, temp2)

    end
    set_wing_fermi!(xout)

    return
end

function Dxplus!(
    xout::StaggeredFermion,
    ν::Int64,
    U::Array{G,1},
    x::StaggeredFermion,
    temps::Array{T,1},
) where {T<:FermionFields,G<:GaugeFields}
    #temp = temps[4]
    temp1 = temps[1]
    temp2 = temps[2]

    #clear!(temp)
    set_wing_fermi!(x)

    fermion_shift!(temp1, U, ν, x)

    fermion_shift!(temp2, U, -ν, x)

    add!(xout, 0.5, temp1, -0.5, temp2)

    set_wing_fermi!(xout)

    return
end


function Dx!(
    xout::StaggeredFermion,
    U::Array{G,1},
    x::StaggeredFermion,
    temps::Array{T,1},
    indices,
) where {T<:FermionFields,G<:GaugeFields}
    #temp = temps[4]
    temp1 = temps[1]
    temp2 = temps[2]

    #clear!(temp)
    set_wing_fermi!(x)
    clear!(xout)
    vec_indices1 = Array{Tuple,1}(undef, 4)
    vec_indices2 = Array{Tuple,1}(undef, 4)

    for ν = 1:4
        clear!(temp1)
        vec_indices1[ν] = fermion_shift!(temp1, U, ν, x, indices)

        #println(temp1*temp1)


        clear!(temp2)

        vec_indices2[ν] = fermion_shift!(temp2, U, -ν, x, indices)



        add!(xout, 0.5, temp1, -0.5, temp2)

    end
    set_wing_fermi!(xout)

    return vec_indices1, vec_indices2
end

function DdagDx!(
    xout::StaggeredFermion,
    U::Array{G,1},
    x::StaggeredFermion,
    temps::Array{T,1},
    indices,
) where {T<:FermionFields,G<:GaugeFields}
    #temp = temps[4]


    temp1 = temps[1]
    temp2 = temps[2]
    temp3 = temps[3]
    temp4 = temps[4]
    temp5 = temps[5]

    clear!(temp1)
    clear!(temp2)
    clear!(temp4)
    clear!(temp5)


    #clear!(temp)
    set_wing_fermi!(x)
    clear!(xout)
    vec_indices = Array{Tuple,1}(undef, 8)
    #vec_indices2 = Array{Tuple,1}(undef,4)

    clear!(temp3)
    for ν = 1:4
        vec_indices[2ν-1] = fermion_shift!(temp1, U, ν, x, [indices])[1]
        vec_indices[2ν] = fermion_shift!(temp2, U, -ν, x, [indices])[1]

    end
    add!(temp3, 0.5, temp1, -0.5, temp2)

    for ν = 1:4
        fermion_shift!(temp4, U, ν, temp3, vec_indices)
        fermion_shift!(temp5, U, -ν, temp3, vec_indices)

    end
    add!(xout, 0.5, temp4, -0.5, temp5)

    set_wing_fermi!(xout)
    return

end



function Dx!(
    xout::StaggeredFermion,
    U::Array{G,1},
    x::StaggeredFermion,
    temps::Array{T,1},
    fparam::FermiActionParam_Staggered,
) where {T<:FermionFields,G<:GaugeFields}
    Dx!(xout, U, x, temps)
    return
end

function Dx!(
    xout::StaggeredFermion,
    evensite,
    U::Array{G,1},
    x::StaggeredFermion,
    temps::Array{T,1},
    fparam::FermiActionParam_Staggered,
) where {T<:FermionFields,G<:GaugeFields}
    #temp = temps[4]
    temp1 = temps[1]
    temp2 = temps[2]

    #clear!(temp)
    set_wing_fermi!(x)
    clear!(xout)
    for ν = 1:4
        clear!(temp1)
        fermion_shift!(temp1, evensite, U, ν, x)

        #... Dirac multiplication
        #mul!(temp1,view(x.rminusγ,:,:,ν),temp1)

        #
        clear!(temp2)
        fermion_shift!(temp2, evensite, U, -ν, x)
        #mul!(temp2,view(x.rplusγ,:,:,ν),temp2)

        #add!(temp,1,temp1,-1,temp2)
        add!(xout, 0.5, temp1, -0.5, temp2)

    end
    set_wing_fermi!(xout)

    #clear!(xout)
    #add!(xout,fparam.mass,x,1,temp)

    #display(xout)
    #    exit()
    return
end


function Wdagx!(
    xout::StaggeredFermion,
    U::Array{G,1},
    x::StaggeredFermion,
    temps::Array{T,1},
    fparam::FermiActionParam_Staggered,
) where {T<:FermionFields,G<:GaugeFields}

    temp = temps[4]
    Dx!(temp, U, x, temps, fparam)
    clear!(xout)
    add!(xout, fparam.mass, x, -1, temp)
    return

end


function Wdagx!(
    xout::StaggeredFermion,
    U::Array{G,1},
    x::StaggeredFermion,
    temps::Array{T,1},
    fparam::FermiActionParam_Staggered,
    vec_indices1,
    vec_indices2,
) where {T<:FermionFields,G<:GaugeFields}

    temp = temps[4]
    Dx!(temp, U, x, temps, vec_indices1, vec_indices2)
    clear!(xout)
    add!(xout, fparam.mass, x, -1, temp)
    return

end


function WdagWx!(
    xout::T,
    U::Array{G,1},
    x::T,
    temps::Array{T,1},
    mass::Number,
    indices,
) where {T<:StaggeredFermion,G<:GaugeFields}

    temp = temps[6]
    #println("x = ", x*x)
    DdagDx!(temp, U, x, temps, indices)

    #println("temp ",temp*temp)
    clear!(xout)
    add!(xout, mass^2, x, -1, temp)
    #println("xout ",xout*xout)
    return
end

function WdagWx!(
    xout::T,
    U::Array{G,1},
    x::T,
    temps::Array{T,1},
    fparam::FermiActionParam_Staggered,
    indices,
) where {T<:StaggeredFermion,G<:GaugeFields}
    WdagWx!(xout, U, x, temps, fparam.mass, indices)
    return
end


function Ddagx!(
    xout::StaggeredFermion,
    U::Array{G,1},
    x::StaggeredFermion,
    temps::Array{T,1},
    fparam::FermiActionParam_Staggered,
) where {T<:FermionFields,G<:GaugeFields}
    #temp = temps[4]
    temp1 = temps[1]
    temp2 = temps[2]

    #clear!(temp)
    clear!(xout)
    set_wing_fermi!(x)
    for ν = 1:4
        fermion_shift!(temp1, U, ν, x)

        #... Dirac multiplication
        #mul!(temp1,view(x.rminusγ,:,:,ν),temp1)

        #
        fermion_shift!(temp2, U, -ν, x)
        #mul!(temp2,view(x.rplusγ,:,:,ν),temp2)

        #add!(temp,1,temp1,-1,temp2)

        add!(xout, -0.5, temp1, 0.5, temp2)
        #add!(xout,1,temp1,-1,temp2)

    end
    set_wing_fermi!(xout)

    #clear!(xout)
    #add!(xout,fparam.mass,x,1,temp)

    #display(xout)
    #    exit()
    return
end

function set_wing_fermi!(a::StaggeredFermion)
    NT = a.NT
    NZ = a.NZ
    NY = a.NY
    NX = a.NX
    NC = a.NC

    #!  X-direction
    for ialpha = 1:1
        for it = 1:NT
            for iz = 1:NZ
                for iy = 1:NY
                    @simd for k = 1:NC
                        a[k, 0, iy, iz, it, ialpha] =
                            a.BoundaryCondition[1] * a[k, NX, iy, iz, it, ialpha]
                    end
                end
            end
        end
    end

    for ialpha = 1:1
        for it = 1:NT
            for iz = 1:NZ
                for iy = 1:NY
                    @simd for k = 1:NC
                        a[k, NX+1, iy, iz, it, ialpha] =
                            a.BoundaryCondition[1] * a[k, 1, iy, iz, it, ialpha]
                    end
                end
            end
        end
    end

    #Y-direction
    for ialpha = 1:1
        for it = 1:NT
            for iz = 1:NZ
                for ix = 1:NX
                    @simd for k = 1:NC
                        a[k, ix, 0, iz, it, ialpha] =
                            a.BoundaryCondition[2] * a[k, ix, NY, iz, it, ialpha]
                    end
                end
            end
        end
    end

    for ialpha = 1:1
        for it = 1:NT
            for iz = 1:NZ
                for ix = 1:NX
                    @simd for k = 1:NC
                        a[k, ix, NY+1, iz, it, ialpha] =
                            a.BoundaryCondition[2] * a[k, ix, 1, iz, it, ialpha]
                    end
                end
            end
        end
    end


    for ialpha = 1:1
        # Z-direction
        for it = 1:NT
            for iy = 1:NY
                for ix = 1:NX
                    @simd for k = 1:NC
                        a[k, ix, iy, 0, it, ialpha] =
                            a.BoundaryCondition[3] * a[k, ix, iy, NZ, it, ialpha]
                        a[k, ix, iy, NZ+1, it, ialpha] =
                            a.BoundaryCondition[3] * a[k, ix, iy, 1, it, ialpha]

                    end
                end
            end
        end

        #T-direction
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    @simd for k = 1:NC
                        a[k, ix, iy, iz, 0, ialpha] =
                            a.BoundaryCondition[4] * a[k, ix, iy, iz, NT, ialpha]
                        a[k, ix, iy, iz, NT+1, ialpha] =
                            a.BoundaryCondition[4] * a[k, ix, iy, iz, 1, ialpha]
                    end
                end
            end
        end

    end





end

function fermion_shift!(
    b::StaggeredFermion,
    u::Array{T,1},
    μ::Int,
    a::StaggeredFermion,
) where {T<:SU3GaugeFields}
    if μ == 0
        substitute!(b, a)
        return
    end

    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    NC = a.NC

    #NTrange = get_looprange(NT)
    #println(NTrange)
    if μ > 0
        #idel = zeros(Int64,4)
        #idel[μ] = 1

        n6 = size(a.f)[6]
        for ialpha = 1:1
            #for it=NTrange
            for it = 1:NT
                it1 = it + ifelse(μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz + ifelse(μ == 3, 1, 0) #idel[3]
                    for iy = 1:NY
                        iy1 = iy + ifelse(μ == 2, 1, 0) #idel[2]
                        @simd for ix = 1:NX
                            ix1 = ix + ifelse(μ == 1, 1, 0) #idel[1]
                            η = staggered_phase(μ, ix, iy, iz, it, NX, NY, NZ, NT)

                            b[1, ix, iy, iz, it, ialpha] =
                                η * (
                                    u[μ][1, 1, ix, iy, iz, it] *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][1, 2, ix, iy, iz, it] *
                                    a[2, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][1, 3, ix, iy, iz, it] *
                                    a[3, ix1, iy1, iz1, it1, ialpha]
                                )

                            b[2, ix, iy, iz, it, ialpha] =
                                η * (
                                    u[μ][2, 1, ix, iy, iz, it] *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][2, 2, ix, iy, iz, it] *
                                    a[2, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][2, 3, ix, iy, iz, it] *
                                    a[3, ix1, iy1, iz1, it1, ialpha]
                                )

                            b[3, ix, iy, iz, it, ialpha] =
                                η * (
                                    u[μ][3, 1, ix, iy, iz, it] *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][3, 2, ix, iy, iz, it] *
                                    a[2, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][3, 3, ix, iy, iz, it] *
                                    a[3, ix1, iy1, iz1, it1, ialpha]
                                )

                            #=
                            for k1=1:3
                                b[k1,ix,iy,iz,it,ialpha] = 0
                                for k2=1:3
                                    b[k1,ix,iy,iz,it,ialpha] += u[μ][k1,k2,ix,iy,iz,it]*a[k2,ix1,iy1,iz1,it1,ialpha]
                                end
                            end
                            =#
                        end
                    end
                end
            end
        end

    elseif μ < 0
        #idel = zeros(Int64,4)
        #idel[-μ] = 1
        #n6 = size(b.f)[6]
        for ialpha = 1:1
            for it = 1:NT
                it1 = it - ifelse(-μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz - ifelse(-μ == 3, 1, 0) #idel[3]
                    for iy = 1:NY
                        iy1 = iy - ifelse(-μ == 2, 1, 0)  #idel[2]
                        @simd for ix = 1:NX
                            ix1 = ix - ifelse(-μ == 1, 1, 0) #idel[1]

                            η = staggered_phase(-μ, ix1, iy1, iz1, it1, NX, NY, NZ, NT)

                            b[1, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(u[-μ][1, 1, ix1, iy1, iz1, it1]) *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][2, 1, ix1, iy1, iz1, it1]) *
                                    a[2, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][3, 1, ix1, iy1, iz1, it1]) *
                                    a[3, ix1, iy1, iz1, it1, ialpha]
                                )

                            b[2, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(u[-μ][1, 2, ix1, iy1, iz1, it1]) *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][2, 2, ix1, iy1, iz1, it1]) *
                                    a[2, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][3, 2, ix1, iy1, iz1, it1]) *
                                    a[3, ix1, iy1, iz1, it1, ialpha]
                                )

                            b[3, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(u[-μ][1, 3, ix1, iy1, iz1, it1]) *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][2, 3, ix1, iy1, iz1, it1]) *
                                    a[2, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][3, 3, ix1, iy1, iz1, it1]) *
                                    a[3, ix1, iy1, iz1, it1, ialpha]
                                )
                            #=
                            for k1=1:3
                                b[k1,ix,iy,iz,it,ialpha] = 0
                                for k2=1:3
                                    b[k1,ix,iy,iz,it,ialpha] += conj(u[-μ][k2,k1,ix1,iy1,iz1,it1])*a[k2,ix1,iy1,iz1,it1,ialpha]
                                end
                            end
                            =#
                        end
                    end
                end
            end
        end
    end

end

function fermion_shift!(
    b::StaggeredFermion,
    u::Array{GaugeFields{SU{NC}},1},
    μ::Int,
    a::StaggeredFermion,
) where {NC}
    if μ == 0
        substitute!(b, a)
        return
    end

    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    #NC = a.NC

    #NTrange = get_looprange(NT)
    #println(NTrange)
    if μ > 0
        #idel = zeros(Int64,4)
        #idel[μ] = 1

        n6 = size(a.f)[6]
        for ialpha = 1:1
            #for it=NTrange
            for it = 1:NT
                it1 = it + ifelse(μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz + ifelse(μ == 3, 1, 0) #idel[3]
                    for iy = 1:NY
                        iy1 = iy + ifelse(μ == 2, 1, 0) #idel[2]
                        for ix = 1:NX
                            ix1 = ix + ifelse(μ == 1, 1, 0) #idel[1]
                            η = staggered_phase(μ, ix, iy, iz, it, NX, NY, NZ, NT)

                            for k1 = 1:NC
                                b[k1, ix, iy, iz, it, ialpha] = 0
                                @simd for k2 = 1:NC
                                    b[k1, ix, iy, iz, it, ialpha] +=
                                        η *
                                        u[μ][k1, k2, ix, iy, iz, it] *
                                        a[k2, ix1, iy1, iz1, it1, ialpha]
                                end
                            end

                        end
                    end
                end
            end
        end

    elseif μ < 0
        #idel = zeros(Int64,4)
        #idel[-μ] = 1
        #n6 = size(b.f)[6]
        for ialpha = 1:1
            for it = 1:NT
                it1 = it - ifelse(-μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz - ifelse(-μ == 3, 1, 0) #idel[3]
                    for iy = 1:NY
                        iy1 = iy - ifelse(-μ == 2, 1, 0)  #idel[2]
                        for ix = 1:NX
                            ix1 = ix - ifelse(-μ == 1, 1, 0) #idel[1]

                            η = staggered_phase(-μ, ix1, iy1, iz1, it1, NX, NY, NZ, NT)


                            for k1 = 1:NC
                                b[k1, ix, iy, iz, it, ialpha] = 0
                                @simd for k2 = 1:NC
                                    b[k1, ix, iy, iz, it, ialpha] +=
                                        η *
                                        conj(u[-μ][k2, k1, ix1, iy1, iz1, it1]) *
                                        a[k2, ix1, iy1, iz1, it1, ialpha]
                                end
                            end

                        end
                    end
                end
            end
        end
    end

end

function fermion_shift!(
    b::StaggeredFermion,
    evensite,
    u::Array{T,1},
    μ::Int,
    a::StaggeredFermion,
) where {T<:SU3GaugeFields}
    if μ == 0
        substitute!(b, a)
        return
    end

    ibush = ifelse(evensite, 0, 1)

    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    NC = a.NC

    #NTrange = get_looprange(NT)
    #println(NTrange)
    if μ > 0
        #idel = zeros(Int64,4)
        #idel[μ] = 1

        n6 = size(a.f)[6]
        for ialpha = 1:1
            #for it=NTrange
            for it = 1:NT
                it1 = it + ifelse(μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz + ifelse(μ == 3, 1, 0) #idel[3]
                    for iy = 1:NY
                        iy1 = iy + ifelse(μ == 2, 1, 0) #idel[2]
                        xran = 1+(1+ibush+iy+iz+it)%2:2:NX
                        @simd for ix in xran
                            #for ix=1:NX
                            ix1 = ix + ifelse(μ == 1, 1, 0) #idel[1]
                            η = staggered_phase(μ, ix, iy, iz, it, NX, NY, NZ, NT)

                            b[1, ix, iy, iz, it, ialpha] =
                                η * (
                                    u[μ][1, 1, ix, iy, iz, it] *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][1, 2, ix, iy, iz, it] *
                                    a[2, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][1, 3, ix, iy, iz, it] *
                                    a[3, ix1, iy1, iz1, it1, ialpha]
                                )

                            b[2, ix, iy, iz, it, ialpha] =
                                η * (
                                    u[μ][2, 1, ix, iy, iz, it] *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][2, 2, ix, iy, iz, it] *
                                    a[2, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][2, 3, ix, iy, iz, it] *
                                    a[3, ix1, iy1, iz1, it1, ialpha]
                                )

                            b[3, ix, iy, iz, it, ialpha] =
                                η * (
                                    u[μ][3, 1, ix, iy, iz, it] *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][3, 2, ix, iy, iz, it] *
                                    a[2, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][3, 3, ix, iy, iz, it] *
                                    a[3, ix1, iy1, iz1, it1, ialpha]
                                )

                            #=
                            for k1=1:3
                                b[k1,ix,iy,iz,it,ialpha] = 0
                                for k2=1:3
                                    b[k1,ix,iy,iz,it,ialpha] += u[μ][k1,k2,ix,iy,iz,it]*a[k2,ix1,iy1,iz1,it1,ialpha]
                                end
                            end
                            =#
                        end
                    end
                end
            end
        end

    elseif μ < 0
        #idel = zeros(Int64,4)
        #idel[-μ] = 1
        #n6 = size(b.f)[6]
        for ialpha = 1:1
            for it = 1:NT
                it1 = it - ifelse(-μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz - ifelse(-μ == 3, 1, 0) #idel[3]
                    for iy = 1:NY
                        iy1 = iy - ifelse(-μ == 2, 1, 0)  #idel[2]
                        xran = 1+(1+ibush+iy+iz+it)%2:2:NX
                        @simd for ix in xran
                            ix1 = ix - ifelse(-μ == 1, 1, 0) #idel[1]

                            η = staggered_phase(-μ, ix1, iy1, iz1, it1, NX, NY, NZ, NT)

                            b[1, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(u[-μ][1, 1, ix1, iy1, iz1, it1]) *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][2, 1, ix1, iy1, iz1, it1]) *
                                    a[2, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][3, 1, ix1, iy1, iz1, it1]) *
                                    a[3, ix1, iy1, iz1, it1, ialpha]
                                )

                            b[2, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(u[-μ][1, 2, ix1, iy1, iz1, it1]) *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][2, 2, ix1, iy1, iz1, it1]) *
                                    a[2, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][3, 2, ix1, iy1, iz1, it1]) *
                                    a[3, ix1, iy1, iz1, it1, ialpha]
                                )

                            b[3, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(u[-μ][1, 3, ix1, iy1, iz1, it1]) *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][2, 3, ix1, iy1, iz1, it1]) *
                                    a[2, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][3, 3, ix1, iy1, iz1, it1]) *
                                    a[3, ix1, iy1, iz1, it1, ialpha]
                                )
                            #=
                            for k1=1:3
                                b[k1,ix,iy,iz,it,ialpha] = 0
                                for k2=1:3
                                    b[k1,ix,iy,iz,it,ialpha] += conj(u[-μ][k2,k1,ix1,iy1,iz1,it1])*a[k2,ix1,iy1,iz1,it1,ialpha]
                                end
                            end
                            =#
                        end
                    end
                end
            end
        end
    end

end

function fermion_shift!(
    b::StaggeredFermion,
    u::Array{T,1},
    μ::Int,
    a::StaggeredFermion,
) where {T<:SU2GaugeFields}
    if μ == 0
        substitute!(b, a)
        return
    end

    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    NC = a.NC

    #NTrange = get_looprange(NT)
    #println(NTrange)
    if μ > 0
        #idel = zeros(Int64,4)
        #idel[μ] = 1

        n6 = size(a.f)[6]
        for ialpha = 1:1
            #for it=NTrange
            for it = 1:NT
                it1 = it + ifelse(μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz + ifelse(μ == 3, 1, 0) #idel[3]
                    for iy = 1:NY
                        iy1 = iy + ifelse(μ == 2, 1, 0) #idel[2]
                        for ix = 1:NX
                            ix1 = ix + ifelse(μ == 1, 1, 0) #idel[1]
                            η = staggered_phase(μ, ix, iy, iz, it, NX, NY, NZ, NT)

                            b[1, ix, iy, iz, it, ialpha] =
                                η * (
                                    u[μ][1, 1, ix, iy, iz, it] *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][1, 2, ix, iy, iz, it] *
                                    a[2, ix1, iy1, iz1, it1, ialpha]
                                )

                            b[2, ix, iy, iz, it, ialpha] =
                                η * (
                                    u[μ][2, 1, ix, iy, iz, it] *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][2, 2, ix, iy, iz, it] *
                                    a[2, ix1, iy1, iz1, it1, ialpha]
                                )


                            #=
                            for k1=1:3
                                b[k1,ix,iy,iz,it,ialpha] = 0
                                for k2=1:3
                                    b[k1,ix,iy,iz,it,ialpha] += u[μ][k1,k2,ix,iy,iz,it]*a[k2,ix1,iy1,iz1,it1,ialpha]
                                end
                            end
                            =#
                        end
                    end
                end
            end
        end

    elseif μ < 0
        #idel = zeros(Int64,4)
        #idel[-μ] = 1
        #n6 = size(b.f)[6]
        for ialpha = 1:1
            for it = 1:NT
                it1 = it - ifelse(-μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz - ifelse(-μ == 3, 1, 0) #idel[3]
                    for iy = 1:NY
                        iy1 = iy - ifelse(-μ == 2, 1, 0)  #idel[2]
                        for ix = 1:NX
                            ix1 = ix - ifelse(-μ == 1, 1, 0) #idel[1]
                            η = staggered_phase(-μ, ix1, iy1, iz1, it1, NX, NY, NZ, NT)

                            b[1, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(u[-μ][1, 1, ix1, iy1, iz1, it1]) *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][2, 1, ix1, iy1, iz1, it1]) *
                                    a[2, ix1, iy1, iz1, it1, ialpha]
                                )

                            b[2, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(u[-μ][1, 2, ix1, iy1, iz1, it1]) *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][2, 2, ix1, iy1, iz1, it1]) *
                                    a[2, ix1, iy1, iz1, it1, ialpha]
                                )

                        end
                    end
                end
            end
        end
    end

end



function fermion_shift!(
    b::StaggeredFermion,
    u::Array{T,1},
    μ::Int,
    a::StaggeredFermion,
    vec_indices,
) where {T<:SU2GaugeFields}
    if μ == 0
        substitute!(b, a)
        return
    end

    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    NC = a.NC

    #NTrange = get_looprange(NT)
    #println(NTrange)
    outindices = Array{Tuple,1}(undef, length(vec_indices))
    #clear!(b)
    if μ > 0
        for (k, indices) in enumerate(vec_indices)

            ix1 = indices[1]
            iy1 = indices[2]
            iz1 = indices[3]
            it1 = indices[4]
            ialpha = indices[5]

            ix = ix1 - ifelse(μ == 1, 1, 0)
            iy = iy1 - ifelse(μ == 2, 1, 0)
            iz = iz1 - ifelse(μ == 3, 1, 0)
            it = it1 - ifelse(μ == 4, 1, 0)

            η = staggered_phase(μ, ix, iy, iz, it, NX, NY, NZ, NT)

            ix, iy, iz, it, sign =
                apply_periodicity(ix, iy, iz, it, NX, NY, NZ, NT, b.BoundaryCondition)


            outindices[k] = (ix, iy, iz, it, ialpha)



            b[1, ix, iy, iz, it, ialpha] +=
                sign *
                η *
                (
                    u[μ][1, 1, ix, iy, iz, it] * a[1, ix1, iy1, iz1, it1, ialpha] +
                    u[μ][1, 2, ix, iy, iz, it] * a[2, ix1, iy1, iz1, it1, ialpha]
                )

            b[2, ix, iy, iz, it, ialpha] +=
                sign *
                η *
                (
                    u[μ][2, 1, ix, iy, iz, it] * a[1, ix1, iy1, iz1, it1, ialpha] +
                    u[μ][2, 2, ix, iy, iz, it] * a[2, ix1, iy1, iz1, it1, ialpha]
                )
        end

    elseif μ < 0
        for (k, indices) in enumerate(vec_indices)

            ix1 = indices[1]
            iy1 = indices[2]
            iz1 = indices[3]
            it1 = indices[4]
            ialpha = indices[5]

            ix = ix1 + ifelse(-μ == 1, 1, 0)
            iy = iy1 + ifelse(-μ == 2, 1, 0)
            iz = iz1 + ifelse(-μ == 3, 1, 0)
            it = it1 + ifelse(-μ == 4, 1, 0)

            η = staggered_phase(-μ, ix1, iy1, iz1, it1, NX, NY, NZ, NT)

            ix, iy, iz, it, sign =
                apply_periodicity(ix, iy, iz, it, NX, NY, NZ, NT, b.BoundaryCondition)

            outindices[k] = (ix, iy, iz, it, ialpha)


            b[1, ix, iy, iz, it, ialpha] +=
                sign *
                η *
                (
                    conj(u[-μ][1, 1, ix1, iy1, iz1, it1]) *
                    a[1, ix1, iy1, iz1, it1, ialpha] +
                    conj(u[-μ][2, 1, ix1, iy1, iz1, it1]) *
                    a[2, ix1, iy1, iz1, it1, ialpha]
                )

            b[2, ix, iy, iz, it, ialpha] +=
                sign *
                η *
                (
                    conj(u[-μ][1, 2, ix1, iy1, iz1, it1]) *
                    a[1, ix1, iy1, iz1, it1, ialpha] +
                    conj(u[-μ][2, 2, ix1, iy1, iz1, it1]) *
                    a[2, ix1, iy1, iz1, it1, ialpha]
                )
        end

    end

    return outindices

end

function fermion_shift!(
    b::StaggeredFermion,
    u::Array{T,1},
    μ::Int,
    a::StaggeredFermion,
    vec_indices,
) where {T<:SU3GaugeFields}
    if μ == 0
        substitute!(b, a)
        return
    end

    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    NC = a.NC

    #NTrange = get_looprange(NT)
    #println(NTrange)
    outindices = Array{Tuple,1}(undef, length(vec_indices))
    #clear!(b)
    if μ > 0
        for (k, indices) in enumerate(vec_indices)

            ix1 = indices[1]
            iy1 = indices[2]
            iz1 = indices[3]
            it1 = indices[4]
            ialpha = indices[5]

            ix = ix1 - ifelse(μ == 1, 1, 0)
            iy = iy1 - ifelse(μ == 2, 1, 0)
            iz = iz1 - ifelse(μ == 3, 1, 0)
            it = it1 - ifelse(μ == 4, 1, 0)

            η = staggered_phase(μ, ix, iy, iz, it, NX, NY, NZ, NT)

            ix, iy, iz, it, sign =
                apply_periodicity(ix, iy, iz, it, NX, NY, NZ, NT, b.BoundaryCondition)


            outindices[k] = (ix, iy, iz, it, ialpha)



            b[1, ix, iy, iz, it, ialpha] +=
                sign *
                η *
                (
                    u[μ][1, 1, ix, iy, iz, it] * a[1, ix1, iy1, iz1, it1, ialpha] +
                    u[μ][1, 2, ix, iy, iz, it] * a[2, ix1, iy1, iz1, it1, ialpha] +
                    u[μ][1, 3, ix, iy, iz, it] * a[3, ix1, iy1, iz1, it1, ialpha]
                )

            b[2, ix, iy, iz, it, ialpha] +=
                sign *
                η *
                (
                    u[μ][2, 1, ix, iy, iz, it] * a[1, ix1, iy1, iz1, it1, ialpha] +
                    u[μ][2, 2, ix, iy, iz, it] * a[2, ix1, iy1, iz1, it1, ialpha] +
                    u[μ][2, 3, ix, iy, iz, it] * a[3, ix1, iy1, iz1, it1, ialpha]
                )
            b[3, ix, iy, iz, it, ialpha] +=
                sign *
                η *
                (
                    u[μ][3, 1, ix, iy, iz, it] * a[1, ix1, iy1, iz1, it1, ialpha] +
                    u[μ][3, 2, ix, iy, iz, it] * a[2, ix1, iy1, iz1, it1, ialpha] +
                    u[μ][3, 3, ix, iy, iz, it] * a[3, ix1, iy1, iz1, it1, ialpha]
                )
        end

    elseif μ < 0
        for (k, indices) in enumerate(vec_indices)

            ix1 = indices[1]
            iy1 = indices[2]
            iz1 = indices[3]
            it1 = indices[4]
            ialpha = indices[5]

            ix = ix1 + ifelse(-μ == 1, 1, 0)
            iy = iy1 + ifelse(-μ == 2, 1, 0)
            iz = iz1 + ifelse(-μ == 3, 1, 0)
            it = it1 + ifelse(-μ == 4, 1, 0)

            η = staggered_phase(-μ, ix1, iy1, iz1, it1, NX, NY, NZ, NT)

            ix, iy, iz, it, sign =
                apply_periodicity(ix, iy, iz, it, NX, NY, NZ, NT, b.BoundaryCondition)

            outindices[k] = (ix, iy, iz, it, ialpha)


            b[1, ix, iy, iz, it, ialpha] +=
                sign *
                η *
                (
                    conj(u[-μ][1, 1, ix1, iy1, iz1, it1]) *
                    a[1, ix1, iy1, iz1, it1, ialpha] +
                    conj(u[-μ][2, 1, ix1, iy1, iz1, it1]) *
                    a[2, ix1, iy1, iz1, it1, ialpha] +
                    conj(u[-μ][3, 1, ix1, iy1, iz1, it1]) * a[3, ix1, iy1, iz1, it1, ialpha]
                )

            b[2, ix, iy, iz, it, ialpha] +=
                sign *
                η *
                (
                    conj(u[-μ][1, 2, ix1, iy1, iz1, it1]) *
                    a[1, ix1, iy1, iz1, it1, ialpha] +
                    conj(u[-μ][2, 2, ix1, iy1, iz1, it1]) *
                    a[2, ix1, iy1, iz1, it1, ialpha] +
                    conj(u[-μ][3, 2, ix1, iy1, iz1, it1]) * a[3, ix1, iy1, iz1, it1, ialpha]
                )
            b[3, ix, iy, iz, it, ialpha] +=
                sign *
                η *
                (
                    conj(u[-μ][1, 3, ix1, iy1, iz1, it1]) *
                    a[1, ix1, iy1, iz1, it1, ialpha] +
                    conj(u[-μ][2, 3, ix1, iy1, iz1, it1]) *
                    a[2, ix1, iy1, iz1, it1, ialpha] +
                    conj(u[-μ][3, 3, ix1, iy1, iz1, it1]) * a[3, ix1, iy1, iz1, it1, ialpha]
                )
        end

    end

    return outindices

end

function fermion_shift!(
    b::StaggeredFermion,
    u::Array{GaugeFields{SU{NC}},1},
    μ::Int,
    a::StaggeredFermion,
    vec_indices,
) where {NC}
    if μ == 0
        substitute!(b, a)
        return
    end

    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    #NC = a.NC

    #NTrange = get_looprange(NT)
    #println(NTrange)
    outindices = Array{Tuple,1}(undef, length(vec_indices))
    #clear!(b)
    if μ > 0
        for (k, indices) in enumerate(vec_indices)

            ix1 = indices[1]
            iy1 = indices[2]
            iz1 = indices[3]
            it1 = indices[4]
            ialpha = indices[5]

            ix = ix1 - ifelse(μ == 1, 1, 0)
            iy = iy1 - ifelse(μ == 2, 1, 0)
            iz = iz1 - ifelse(μ == 3, 1, 0)
            it = it1 - ifelse(μ == 4, 1, 0)

            η = staggered_phase(μ, ix, iy, iz, it, NX, NY, NZ, NT)

            ix, iy, iz, it, sign =
                apply_periodicity(ix, iy, iz, it, NX, NY, NZ, NT, b.BoundaryCondition)


            outindices[k] = (ix, iy, iz, it, ialpha)

            for k1 = 1:NC
                b[k1, ix, iy, iz, it, ialpha] = 0
                for k2 = 1:NC
                    b[k1, ix, iy, iz, it, ialpha] +=
                        sign *
                        η *
                        u[μ][k1, k2, ix, iy, iz, it] *
                        a[k2, ix1, iy1, iz1, it1, ialpha]
                end
            end


        end

    elseif μ < 0
        for (k, indices) in enumerate(vec_indices)

            ix1 = indices[1]
            iy1 = indices[2]
            iz1 = indices[3]
            it1 = indices[4]
            ialpha = indices[5]

            ix = ix1 + ifelse(-μ == 1, 1, 0)
            iy = iy1 + ifelse(-μ == 2, 1, 0)
            iz = iz1 + ifelse(-μ == 3, 1, 0)
            it = it1 + ifelse(-μ == 4, 1, 0)

            η = staggered_phase(-μ, ix1, iy1, iz1, it1, NX, NY, NZ, NT)

            ix, iy, iz, it, sign =
                apply_periodicity(ix, iy, iz, it, NX, NY, NZ, NT, b.BoundaryCondition)

            outindices[k] = (ix, iy, iz, it, ialpha)

            for k1 = 1:NC
                b[k1, ix, iy, iz, it, ialpha] = 0
                for k2 = 1:NC
                    b[k1, ix, iy, iz, it, ialpha] +=
                        sign *
                        η *
                        conj(u[-μ][k2, k1, ix1, iy1, iz1, it1]) *
                        a[k2, ix1, iy1, iz1, it1, ialpha]
                end
            end


        end

    end

    return outindices

end


function fermion_shift!(
    b::StaggeredFermion,
    evensite,
    u::Array{T,1},
    μ::Int,
    a::StaggeredFermion,
) where {T<:SU2GaugeFields}
    if μ == 0
        substitute!(b, a)
        return
    end

    ibush = ifelse(evensite, 0, 1)

    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    NC = a.NC

    #NTrange = get_looprange(NT)
    #println(NTrange)
    if μ > 0
        #idel = zeros(Int64,4)
        #idel[μ] = 1

        n6 = size(a.f)[6]
        for ialpha = 1:1
            #for it=NTrange
            for it = 1:NT
                it1 = it + ifelse(μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz + ifelse(μ == 3, 1, 0) #idel[3]
                    for iy = 1:NY
                        iy1 = iy + ifelse(μ == 2, 1, 0) #idel[2]
                        xran = 1+(1+ibush+iy+iz+it)%2:2:NX
                        for ix in xran
                            ix1 = ix + ifelse(μ == 1, 1, 0) #idel[1]
                            η = staggered_phase(μ, ix, iy, iz, it, NX, NY, NZ, NT)

                            b[1, ix, iy, iz, it, ialpha] =
                                η * (
                                    u[μ][1, 1, ix, iy, iz, it] *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][1, 2, ix, iy, iz, it] *
                                    a[2, ix1, iy1, iz1, it1, ialpha]
                                )

                            b[2, ix, iy, iz, it, ialpha] =
                                η * (
                                    u[μ][2, 1, ix, iy, iz, it] *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    u[μ][2, 2, ix, iy, iz, it] *
                                    a[2, ix1, iy1, iz1, it1, ialpha]
                                )


                            #=
                            for k1=1:3
                                b[k1,ix,iy,iz,it,ialpha] = 0
                                for k2=1:3
                                    b[k1,ix,iy,iz,it,ialpha] += u[μ][k1,k2,ix,iy,iz,it]*a[k2,ix1,iy1,iz1,it1,ialpha]
                                end
                            end
                            =#
                        end
                    end
                end
            end
        end

    elseif μ < 0
        #idel = zeros(Int64,4)
        #idel[-μ] = 1
        #n6 = size(b.f)[6]
        for ialpha = 1:1
            for it = 1:NT
                it1 = it - ifelse(-μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz - ifelse(-μ == 3, 1, 0) #idel[3]
                    for iy = 1:NY
                        iy1 = iy - ifelse(-μ == 2, 1, 0)  #idel[2]
                        xran = 1+(1+ibush+iy+iz+it)%2:2:NX
                        for ix in xran
                            ix1 = ix - ifelse(-μ == 1, 1, 0) #idel[1]
                            η = staggered_phase(-μ, ix1, iy1, iz1, it1, NX, NY, NZ, NT)

                            b[1, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(u[-μ][1, 1, ix1, iy1, iz1, it1]) *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][2, 1, ix1, iy1, iz1, it1]) *
                                    a[2, ix1, iy1, iz1, it1, ialpha]
                                )

                            b[2, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(u[-μ][1, 2, ix1, iy1, iz1, it1]) *
                                    a[1, ix1, iy1, iz1, it1, ialpha] +
                                    conj(u[-μ][2, 2, ix1, iy1, iz1, it1]) *
                                    a[2, ix1, iy1, iz1, it1, ialpha]
                                )

                        end
                    end
                end
            end
        end
    end

end

function fermion_shiftB!(
    b::StaggeredFermion,
    u::Array{GaugeFields{SU{NC}},1},
    μ,
    a::StaggeredFermion,
) where {NC}
    if μ == 0
        substitute!(b, a)
        return
    end


    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    #NC = a.NC


    if μ > 0
        error("""
        Sorry this case if not yet ready
        mu = $mu
        """)

    elseif μ < 0
        #idel = zeros(Int64,4)
        #idel[-μ] = 1
        #n6 = size(b.f)[6]
        for ialpha = 1:1
            for it = 1:NT
                it1 = it + ifelse(-μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz + ifelse(-μ == 3, 1, 0)  #idel[3]
                    for iy = 1:NY
                        iy1 = iy + ifelse(-μ == 2, 1, 0) #idel[2]
                        for ix = 1:NX
                            ix1 = ix + ifelse(-μ == 1, 1, 0)  #idel[1]
                            η = staggered_phase(-μ, ix, iy, iz, it, NX, NY, NZ, NT)

                            for k1 = 1:NC
                                b[k1, ix, iy, iz, it, ialpha] = 0
                                for k2 = 1:NC
                                    b[k1, ix, iy, iz, it, ialpha] +=
                                        η *
                                        conj(a[k2, ix1, iy1, iz1, it1, ialpha]) *
                                        conj(u[-μ][k1, k2, ix, iy, iz, it])
                                end
                            end


                        end
                    end
                end
            end
        end
    end

end

function fermion_shiftB!(
    b::StaggeredFermion,
    evensite,
    u::Array{GaugeFields{SU{NC}},1},
    μ,
    a::StaggeredFermion,
) where {NC} #T <: Union{SUNGaugeFields,U1GaugeFields}
    if μ == 0
        substitute!(b, a)
        return
    end

    ibush = ifelse(evensite, 0, 1)

    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    #NC = a.NC


    if μ > 0
        error("""
        Sorry this case if not yet ready
        mu = $mu
        """)

    elseif μ < 0
        #idel = zeros(Int64,4)
        #idel[-μ] = 1
        #n6 = size(b.f)[6]
        for ialpha = 1:1
            for it = 1:NT
                it1 = it + ifelse(-μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz + ifelse(-μ == 3, 1, 0)  #idel[3]
                    for iy = 1:NY
                        iy1 = iy + ifelse(-μ == 2, 1, 0) #idel[2]
                        xran = 1+(1+ibush+iy+iz+it)%2:2:NX
                        for ix in xran
                            #for ix=1:NX
                            ix1 = ix + ifelse(-μ == 1, 1, 0)  #idel[1]
                            η = staggered_phase(-μ, ix, iy, iz, it, NX, NY, NZ, NT)

                            for k1 = 1:NC
                                b[k1, ix, iy, iz, it, ialpha] = 0
                                for k2 = 1:NC
                                    b[k1, ix, iy, iz, it, ialpha] +=
                                        η *
                                        conj(a[k2, ix1, iy1, iz1, it1, ialpha]) *
                                        conj(u[-μ][k1, k2, ix, iy, iz, it])
                                end
                            end


                        end
                    end
                end
            end
        end
    end

end




function fermion_shiftB!(
    b::StaggeredFermion,
    u::Array{SU3GaugeFields,1},
    μ,
    a::StaggeredFermion,
)
    if μ == 0
        substitute!(b, a)
        return
    end


    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    NC = a.NC


    if μ > 0
        error("""
        Sorry this case if not yet ready
        mu = $mu
        """)

    elseif μ < 0
        #idel = zeros(Int64,4)
        #idel[-μ] = 1
        #n6 = size(b.f)[6]
        for ialpha = 1:1
            for it = 1:NT
                it1 = it + ifelse(-μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz + ifelse(-μ == 3, 1, 0)  #idel[3]
                    for iy = 1:NY
                        iy1 = iy + ifelse(-μ == 2, 1, 0) #idel[2]
                        for ix = 1:NX
                            ix1 = ix + ifelse(-μ == 1, 1, 0)  #idel[1]
                            η = staggered_phase(-μ, ix, iy, iz, it, NX, NY, NZ, NT)

                            b[1, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(a[1, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][1, 1, ix, iy, iz, it]) +
                                    conj(a[2, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][1, 2, ix, iy, iz, it]) +
                                    conj(a[3, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][1, 3, ix, iy, iz, it])
                                )

                            b[2, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(a[1, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][2, 1, ix, iy, iz, it]) +
                                    conj(a[2, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][2, 2, ix, iy, iz, it]) +
                                    conj(a[3, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][2, 3, ix, iy, iz, it])
                                )

                            b[3, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(a[1, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][3, 1, ix, iy, iz, it]) +
                                    conj(a[2, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][3, 2, ix, iy, iz, it]) +
                                    conj(a[3, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][3, 3, ix, iy, iz, it])
                                )


                            #=
                            for k1=1:NC
                                b[k1,ix,iy,iz,it,ialpha] = 0
                                for k2=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] += conj(a[k2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][k1,k2,ix,iy,iz,it])
                                end
                            end
                            =#

                        end
                    end
                end
            end
        end
    end

end

function fermion_shiftB!(
    b::StaggeredFermion,
    evensite,
    u::Array{SU3GaugeFields,1},
    μ,
    a::StaggeredFermion,
)
    if μ == 0
        substitute!(b, a)
        return
    end

    ibush = ifelse(evensite, 0, 1)

    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    NC = a.NC


    if μ > 0
        error("""
        Sorry this case if not yet ready
        mu = $mu
        """)

    elseif μ < 0
        #idel = zeros(Int64,4)
        #idel[-μ] = 1
        #n6 = size(b.f)[6]
        for ialpha = 1:1
            for it = 1:NT
                it1 = it + ifelse(-μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz + ifelse(-μ == 3, 1, 0)  #idel[3]
                    for iy = 1:NY
                        iy1 = iy + ifelse(-μ == 2, 1, 0) #idel[2]
                        xran = 1+(1+ibush+iy+iz+it)%2:2:NX
                        for ix in xran
                            #for ix=1:NX
                            ix1 = ix + ifelse(-μ == 1, 1, 0)  #idel[1]
                            η = staggered_phase(-μ, ix, iy, iz, it, NX, NY, NZ, NT)

                            b[1, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(a[1, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][1, 1, ix, iy, iz, it]) +
                                    conj(a[2, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][1, 2, ix, iy, iz, it]) +
                                    conj(a[3, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][1, 3, ix, iy, iz, it])
                                )

                            b[2, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(a[1, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][2, 1, ix, iy, iz, it]) +
                                    conj(a[2, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][2, 2, ix, iy, iz, it]) +
                                    conj(a[3, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][2, 3, ix, iy, iz, it])
                                )

                            b[3, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(a[1, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][3, 1, ix, iy, iz, it]) +
                                    conj(a[2, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][3, 2, ix, iy, iz, it]) +
                                    conj(a[3, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][3, 3, ix, iy, iz, it])
                                )


                            #=
                            for k1=1:NC
                                b[k1,ix,iy,iz,it,ialpha] = 0
                                for k2=1:NC
                                    b[k1,ix,iy,iz,it,ialpha] += conj(a[k2,ix1,iy1,iz1,it1,ialpha])*conj(u[-μ][k1,k2,ix,iy,iz,it])
                                end
                            end
                            =#

                        end
                    end
                end
            end
        end
    end

end


function fermion_shiftB!(
    b::StaggeredFermion,
    u::Array{SU2GaugeFields,1},
    μ,
    a::StaggeredFermion,
)
    if μ == 0
        substitute!(b, a)
        return
    end


    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    NC = a.NC


    if μ > 0
        error("""
        Sorry this case if not yet ready
        mu = $mu
        """)

    elseif μ < 0
        #idel = zeros(Int64,4)
        #idel[-μ] = 1
        #n6 = size(b.f)[6]
        for ialpha = 1:1
            for it = 1:NT
                it1 = it + ifelse(-μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz + ifelse(-μ == 3, 1, 0)  #idel[3]
                    for iy = 1:NY
                        iy1 = iy + ifelse(-μ == 2, 1, 0) #idel[2]
                        for ix = 1:NX
                            ix1 = ix + ifelse(-μ == 1, 1, 0)  #idel[1]
                            η = staggered_phase(-μ, ix, iy, iz, it, NX, NY, NZ, NT)


                            b[1, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(a[1, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][1, 1, ix, iy, iz, it]) +
                                    conj(a[2, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][1, 2, ix, iy, iz, it])
                                )

                            b[2, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(a[1, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][2, 1, ix, iy, iz, it]) +
                                    conj(a[2, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][2, 2, ix, iy, iz, it])
                                )


                        end
                    end
                end
            end
        end
    end

end

function fermion_shiftB!(
    b::StaggeredFermion,
    u::Array{SU2GaugeFields,1},
    μ,
    a::StaggeredFermion,
    indices,
)
    if μ == 0
        substitute!(b, a)
        return
    end


    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    NC = a.NC


    if μ > 0
        error("""
        Sorry this case if not yet ready
        mu = $mu
        """)

    elseif μ < 0

        ix1 = indices[1]
        iy1 = indices[2]
        iz1 = indices[3]
        it1 = indices[4]
        ialpha = indices[5]

        ix = ix1 - ifelse(-μ == 1, 1, 0)
        iy = iy1 - ifelse(-μ == 2, 1, 0)
        iz = iz1 - ifelse(-μ == 3, 1, 0)
        it = it1 - ifelse(-μ == 4, 1, 0)

        ix, iy, iz, it, sign =
            apply_periodicity(ix, iy, iz, it, NX, NY, NZ, NT, b.BoundaryCondition)

        η = staggered_phase(-μ, ix, iy, iz, it, NX, NY, NZ, NT)

        b[1, ix, iy, iz, it, ialpha] =
            sign *
            η *
            (
                conj(a[1, ix1, iy1, iz1, it1, ialpha]) * conj(u[-μ][1, 1, ix, iy, iz, it]) +
                conj(a[2, ix1, iy1, iz1, it1, ialpha]) * conj(u[-μ][1, 2, ix, iy, iz, it])
            )

        b[2, ix, iy, iz, it, ialpha] =
            sign *
            η *
            (
                conj(a[1, ix1, iy1, iz1, it1, ialpha]) * conj(u[-μ][2, 1, ix, iy, iz, it]) +
                conj(a[2, ix1, iy1, iz1, it1, ialpha]) * conj(u[-μ][2, 2, ix, iy, iz, it])
            )

    end

    return (ix, iy, iz, it, ialpha)

end

function fermion_shiftB!(
    b::StaggeredFermion,
    evensite,
    u::Array{SU2GaugeFields,1},
    μ,
    a::StaggeredFermion,
)
    if μ == 0
        substitute!(b, a)
        return
    end

    ibush = ifelse(evensite, 0, 1)


    NX = a.NX
    NY = a.NY
    NZ = a.NZ
    NT = a.NT
    NC = a.NC


    if μ > 0
        error("""
        Sorry this case if not yet ready
        mu = $mu
        """)

    elseif μ < 0
        #idel = zeros(Int64,4)
        #idel[-μ] = 1
        #n6 = size(b.f)[6]
        for ialpha = 1:1
            for it = 1:NT
                it1 = it + ifelse(-μ == 4, 1, 0) #idel[4]
                for iz = 1:NZ
                    iz1 = iz + ifelse(-μ == 3, 1, 0)  #idel[3]
                    for iy = 1:NY
                        iy1 = iy + ifelse(-μ == 2, 1, 0) #idel[2]
                        xran = 1+(1+ibush+iy+iz+it)%2:2:NX
                        for ix in xran
                            ix1 = ix + ifelse(-μ == 1, 1, 0)  #idel[1]
                            η = staggered_phase(-μ, ix, iy, iz, it, NX, NY, NZ, NT)


                            b[1, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(a[1, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][1, 1, ix, iy, iz, it]) +
                                    conj(a[2, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][1, 2, ix, iy, iz, it])
                                )

                            b[2, ix, iy, iz, it, ialpha] =
                                η * (
                                    conj(a[1, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][2, 1, ix, iy, iz, it]) +
                                    conj(a[2, ix1, iy1, iz1, it1, ialpha]) *
                                    conj(u[-μ][2, 2, ix, iy, iz, it])
                                )


                        end
                    end
                end
            end
        end
    end

end


end
