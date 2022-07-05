module Diracoperators
import ..Gaugefields: GaugeFields, GaugeFields_1d
import ..Fermionfields:
    FermionFields,
    WilsonFermion,
    Wx!,
    Wdagx!,
    WdagWx!,
    StaggeredFermion,
    set_wing_fermi!,
    Dx!,
    add!,
    clear!,
    Dxplus!,
    DomainwallFermion,
    mul_γ5x!,
    D5DWx!,
    D5DWdagx!,
    substitute_fermion!
import ..Clover: Make_CloverFμν
import ..Actions:
    FermiActionParam_Wilson,
    FermiActionParam_WilsonClover,
    FermiActionParam_Staggered,
    FermiActionParam_Domainwall,
    SmearingParam_nosmearing
import ..CGmethods: cg, bicg
using LinearAlgebra

import ..Verbose_print: Verbose_level, Verbose_3, Verbose_2, Verbose_1, println_verbose3

abstract type Operator end

abstract type Dirac_operator <: Operator end

abstract type Adjoint_Dirac_operator <: Operator end

abstract type DdagD_operator <: Operator end

abstract type γ5D_operator <: Operator end


function get_U(A::Dirac_operator)
    return A.U
end

function get_U(A::Adjoint_Dirac_operator)
    return A.parent.U
end

function get_U(A::DdagD_operator)
    return A.dirac.U
end

function get_U(A::Operator)
    error("get_U is not implemented for operator $(typeof(A))")
end

function get_temporal_fermi(A::Dirac_operator)
    return A._temporal_fermi
end

function get_temporal_fermi(A::Adjoint_Dirac_operator)
    return A.parent._temporal_fermi
end

function get_temporal_fermi(A::DdagD_operator)
    return A.dirac._temporal_fermi
end

function get_temporal_fermi(A::Operator)
    error("get_temporal_fermi is not implemented for operator $(typeof(A))")
end




function Dirac_operator(U::Array{T,1}, x, fparam) where {T<:GaugeFields}
    if typeof(fparam) == FermiActionParam_Wilson
        W = Wilson_operator(U, x, fparam)
    elseif typeof(fparam) == FermiActionParam_WilsonClover
        W = WilsonClover_operator(U, x, fparam)
    elseif typeof(fparam) == FermiActionParam_Staggered
        W = Staggered_operator(U, x, fparam)
    elseif typeof(fparam) == FermiActionParam_Domainwall
        W = Domainwall_operator(U, x, fparam)
    else
        error("unknown fparam: fparam is ", typeof(fparam))
    end
end

function γ5D_operator(U::Array{T,1}, x, fparam) where {T<:GaugeFields}
    if typeof(fparam) == FermiActionParam_Wilson
        W = γ5D_Wilson_operator(Wilson_operator(U, x, fparam))
    elseif typeof(fparam) == FermiActionParam_WilsonClover
        W = γ5D_WilsonClover_operator(U, x, fparam)
    elseif typeof(fparam) == FermiActionParam_Staggered
        W = γ5D_Staggered_operator(U, x, fparam)
    elseif typeof(fparam) == FermiActionParam_Domainwall
        W = γ5D_Domainwall_operator(U, x, fparam)
    else
        error("unknown fparam: fparam is ", typeof(fparam))
    end
end




struct Wilson_operator{T} <: Dirac_operator where {T<:GaugeFields}
    U::Array{T,1}
    _temporal_fermi::Array{WilsonFermion,1}

    function Wilson_operator(U::Array{T,1}, x, fparam) where {T<:GaugeFields}
        num = 5
        _temporal_fermi = Array{WilsonFermion,1}(undef, num)
        for i = 1:num
            _temporal_fermi[i] = similar(x)
        end
        return new{eltype(U)}(U, _temporal_fermi)
    end
end

struct DdagD_Wilson_operator <: DdagD_operator
    dirac::Wilson_operator
    function DdagD_Wilson_operator(U::Array{T,1}, x, fparam) where {T<:GaugeFields}
        return new(Wilson_operator(U, x, fparam))
    end
end

struct WilsonClover_operator{T} <: Dirac_operator where {T<:GaugeFields}
    U::Array{T,1}
    _temporal_fermi::Array{WilsonFermion,1}
    CloverFμν::Array{ComplexF64,4}


    function WilsonClover_operator(U::Array{T,1}, x, fparam) where {T<:GaugeFields}
        num = 5
        _temporal_fermi = Array{WilsonFermion,1}(undef, num)
        for i = 1:num
            _temporal_fermi[i] = similar(x)
        end

        NV = U[1].NV
        NX = U[1].NX
        NY = U[1].NY
        NZ = U[1].NZ
        NT = U[1].NT
        NC = U[1].NC

        temp1 = GaugeFields_1d(NC, NX, NY, NZ, NT)
        temps = Array{typeof(temp1),1}(undef, num)
        for i = 1:num
            temps[i] = deepcopy(temp1)
        end


        CloverFμν = Make_CloverFμν(fparam, U, temps)


        return new{eltype(U)}(U, _temporal_fermi, CloverFμν)
    end
end


struct DdagD_WilsonClover_operator <: DdagD_operator
    dirac::WilsonClover_operator
    function DdagD_WilsonClover_operator(U::Array{T,1}, x, fparam) where {T<:GaugeFields}
        return new(WilsonClover_operator(U, x, fparam))
    end
end

struct Staggered_operator{T} <: Dirac_operator where {T<:GaugeFields}
    U::Array{T,1}
    _temporal_fermi::Array{StaggeredFermion,1}
    mass::Float64

    function Staggered_operator(U::Array{T,1}, x, fparam) where {T<:GaugeFields}
        num = 9
        _temporal_fermi = Array{StaggeredFermion,1}(undef, num)
        for i = 1:num
            _temporal_fermi[i] = similar(x)
        end
        return new{eltype(U)}(U, _temporal_fermi, fparam.mass)
    end
end


struct DdagD_Staggered_operator{perturb} <: DdagD_operator
    dirac::Staggered_operator
    t::Float64
    function DdagD_Staggered_operator(
        U::Array{T,1},
        x,
        fparam;
        t = 1,
    ) where {T<:GaugeFields}
        if t == 1
            perturb = false
        else
            perturb = true
        end

        return new{perturb}(Staggered_operator(U, x, fparam), t)
    end

    function DdagD_Staggered_operator(A::DdagD_Staggered_operator, t)
        if t == 1
            perturb = false
        else
            perturb = true
        end
        return new{perturb}(A.dirac, t)
    end
end


struct DdagDND_Staggered_operator <: DdagD_operator #nondiagonal part of DdagD
    dirac::Staggered_operator

    function DdagDND_Staggered_operator(A::DdagD_Staggered_operator)
        return new(A.dirac)
    end
end


struct D5DW_Domainwall_operator{T} <: Dirac_operator where {T<:GaugeFields}
    U::Array{T,1}
    wilsonoperator::Wilson_operator{T}
    m::Float64
    _temporal_fermi::Array{DomainwallFermion,1}

    function D5DW_Domainwall_operator(U::Array{T,1}, x, fparam, m) where {T<:GaugeFields}
        r = fparam.r
        M = fparam.M
        hop = 1 / (8r + 2M)

        if typeof(fparam.smearing) == SmearingParam_nosmearing
            fparam_wilson = FermiActionParam_Wilson(
                hop,
                r,
                fparam.eps,
                fparam.Dirac_operator,
                fparam.MaxCGstep,
                fparam.quench,
            )
        else
            fparam_wilson = FermiActionParam_Wilson(
                hop,
                r,
                fparam.eps,
                fparam.Dirac_operator,
                fparam.MaxCGstep,
                fparam.quench,
                smearingparameters = "stout",
                loops_list = fparam.loops_list,
                coefficients = fparam.coefficients,
                numlayers = fparam.numlayers,
                L = fparam.L,
            )
        end
        num = 1
        _temporal_fermi = Array{DomainwallFermion,1}(undef, num)
        for i = 1:num
            _temporal_fermi[i] = similar(x)
        end
        wilsonoperator = Wilson_operator(U, x.f[1], fparam_wilson)

        return new{eltype(U)}(U, wilsonoperator, m, _temporal_fermi)
    end

    function D5DW_Domainwall_operator(U::Array{T,1}, x, fparam) where {T<:GaugeFields}
        return D5DW_Domainwall_operator(U, x, fparam, fparam.m)
    end

end

struct D5DWdagD5DW_Wilson_operator <: DdagD_operator
    dirac::D5DW_Domainwall_operator
    function D5DWdagD5DW_Wilson_operator(U::Array{T,1}, x, fparam) where {T<:GaugeFields}
        return new(D5DW_Domainwall_operator(U, x, fparam))
    end
end






struct Domainwall_operator{T} <: Dirac_operator where {T<:GaugeFields}
    D5DW::D5DW_Domainwall_operator{T}
    D5DW_PV::D5DW_Domainwall_operator{T}

    function Domainwall_operator(U::Array{T,1}, x, fparam) where {T<:GaugeFields}
        D5DW = D5DW_Domainwall_operator(U, x, fparam)
        D5DW_PV = D5DW_Domainwall_operator(U, x, fparam, 1)
        return new{eltype(U)}(D5DW, D5DW_PV)
    end
end

function get_U(A::Domainwall_operator)
    return get_U(A.D5DW)
end

function get_temporal_fermi(A::Domainwall_operator)
    return get_temporal_fermi(A.D5DW)
end

struct DdagD_Domainwall_operator <: DdagD_operator
    dirac::Domainwall_operator
    DdagD::D5DWdagD5DW_Wilson_operator

    function DdagD_Domainwall_operator(U::Array{T,1}, x, fparam) where {T<:GaugeFields}
        return new(
            Domainwall_operator(U, x, fparam),
            D5DWdagD5DW_Wilson_operator(U, x, fparam),
        )
    end
end

function get_U(A::DdagD_Domainwall_operator)
    return get_U(A.dirac)
end


function Base.size(A::Staggered_operator)
    NX = A.U[1].NX
    NY = A.U[1].NY
    NZ = A.U[1].NZ
    NT = A.U[1].NT
    NC = A.U[1].NC
    return NX * NY * NZ * NT * NC, NX * NY * NZ * NT * NC
end

function Base.size(A::DdagD_operator)
    return size(A.dirac)
end



function DdagD_operator(U::Array{T,1}, x, fparam) where {T<:GaugeFields}
    if typeof(fparam) == FermiActionParam_Wilson
        W = DdagD_Wilson_operator(U, x, fparam)
    elseif typeof(fparam) == FermiActionParam_WilsonClover
        W = DdagD_WilsonClover_operator(U, x, fparam)
    elseif typeof(fparam) == FermiActionParam_Staggered
        W = DdagD_Staggered_operator(U, x, fparam)
    elseif typeof(fparam) == FermiActionParam_Domainwall
        W = DdagD_Domainwall_operator(U, x, fparam)
    else
        error("unknown fparam: fparam is ", typeof(fparam))
    end
end

struct Adjoint_Wilson_operator <: Adjoint_Dirac_operator
    parent::Wilson_operator
end

struct Adjoint_WilsonClover_operator <: Adjoint_Dirac_operator
    parent::WilsonClover_operator
end

struct Adjoint_Staggered_operator <: Adjoint_Dirac_operator
    parent::Staggered_operator
end

struct Adjoint_D5DW_Domainwall_operator <: Adjoint_Dirac_operator
    parent::D5DW_Domainwall_operator
end

function Base.adjoint(A::Wilson_operator)
    Adjoint_Wilson_operator(A)
end

function Base.adjoint(A::WilsonClover_operator)
    Adjoint_WilsonClover_operator(A)
end

function Base.adjoint(A::Staggered_operator)
    Adjoint_Staggered_operator(A)
end

function Base.adjoint(A::D5DW_Domainwall_operator)
    Adjoint_D5DW_Domainwall_operator(A)
end

struct γ5D_Wilson_operator <: γ5D_operator
    parent::Wilson_operator
end

struct γ5D_WilsonClover_operator <: γ5D_operator
    parent::WilsonClover_operator
end

struct γ5D_Staggered_operator <: γ5D_operator
    parent::Staggered_operator
end

struct γ5D_D5DW_Domainwall_operator <: γ5D_operator
    parent::D5DW_Domainwall_operator
end

Base.adjoint(A::Adjoint_Dirac_operator) = A.parent

function LinearAlgebra.mul!(y::WilsonFermion, A::Wilson_operator, x::WilsonFermion) #y = A*x
    Wx!(y, A.U, x, A._temporal_fermi)
    return
end

function LinearAlgebra.mul!(y::WilsonFermion, A::γ5D_Wilson_operator, x::WilsonFermion) #y = A*x
    temp = A.parent._temporal_fermi[5]
    mul!(temp, A.parent, x)
    mul_γ5x!(y, temp)
    return
end

function LinearAlgebra.mul!(y::WilsonFermion, A::WilsonClover_operator, x::WilsonFermion) #y = A*x
    Wx!(y, A.U, x, A._temporal_fermi, A.CloverFμν)
    return
end

function LinearAlgebra.mul!(
    y::DomainwallFermion,
    A::D5DW_Domainwall_operator,
    x::DomainwallFermion,
) #y = A*x
    D5DWx!(y, A.U, x, A.m, A.wilsonoperator._temporal_fermi)
    return
end

function LinearAlgebra.mul!(
    y::DomainwallFermion,
    A::Domainwall_operator,
    x::DomainwallFermion,
) #y = A*x

    #A = D5DW(m)*D5DW(m=1))^(-1)
    #y = A*x = D5DW(m)*D5DW(m=1))^(-1)*x
    bicg(A.D5DW_PV._temporal_fermi[1], A.D5DW_PV, x)
    mul!(y, A.D5DW, A.D5DW_PV._temporal_fermi[1])

    #error("Do not use Domainwall_operator directory. Use D5DW_Domainwall_operator M = D5DW(m)*D5DW(-1)^{-1}")
    #D5DWx!(y,A.U,x,A.m,A.wilsonoperator._temporal_fermi) 
    return
end

function LinearAlgebra.mul!(
    y::FermionFields,
    A::T,
    x::FermionFields,
) where {T<:DdagD_operator} #y = A*x
    temp = A.dirac._temporal_fermi[5]
    mul!(temp, A.dirac, x)
    mul!(y, A.dirac', temp)

    return
end

function LinearAlgebra.mul!(
    y::FermionFields,
    A::T,
    x::FermionFields,
) where {T<:D5DWdagD5DW_Wilson_operator}#y = A*x
    temp = A.dirac._temporal_fermi[1]
    mul!(temp, A.dirac, x)
    mul!(y, A.dirac', temp)

    return
end

function LinearAlgebra.mul!(
    y::FermionFields,
    A::DdagD_Staggered_operator{perturb},
    x::FermionFields,
) where {perturb}
    temps = A.dirac._temporal_fermi
    temp = temps[5]
    temp2 = temps[6]
    Dx!(temp, A.dirac.U, x, temps)
    Dx!(temp2, A.dirac.U, temp, temps)
    clear!(y)


    if perturb
        add!(y, A.dirac.mass^2, x, -A.t, temp2)


        for ν = 1:4
            Dxplus!(temp, ν, A.dirac.U, x, temps)
            Dxplus!(temp2, ν, A.dirac.U, temp, temps)

            add!(y, 1, x, -(1 - A.t), temp2)
        end

        set_wing_fermi!(y)
    else

        add!(y, A.dirac.mass^2, x, -1, temp2)
        set_wing_fermi!(y)
    end

    return
end

function LinearAlgebra.mul!(
    y::FermionFields,
    A::DdagDND_Staggered_operator,
    x::FermionFields,
) where {perturb}
    temps = A.dirac._temporal_fermi
    temp = temps[5]
    temp2 = temps[6]
    Dx!(temp, A.dirac.U, x, temps)
    Dx!(temp2, A.dirac.U, temp, temps)
    clear!(y)
    add!(y, 1, x, -1, temp2)

    for ν = 1:4
        Dxplus!(temp, ν, A.dirac.U, x, temps)
        Dxplus!(temp2, ν, A.dirac.U, temp, temps)

        add!(y, 1, x, 1, temp2)
    end
    set_wing_fermi!(y)

    return
end

function LinearAlgebra.mul!(y::FermionFields, A::DdagD_operator, x::FermionFields, indices) #y = A*x        
    mul!(y, A, x)
    return
end

function LinearAlgebra.mul!(
    y::FermionFields,
    A::DdagD_Staggered_operator,
    x::FermionFields,
    indices,
) #y = A*x        
    WdagWx!(y, A.dirac.U, x, A.dirac._temporal_fermi, A.dirac.mass, indices)
    return
end


function LinearAlgebra.mul!(y::StaggeredFermion, A::Staggered_operator, x::StaggeredFermion) #y = A*x
    temps = A._temporal_fermi
    temp = temps[4]
    Dx!(temp, A.U, x, temps)
    clear!(y)
    add!(y, A.mass, x, 1, temp)
    set_wing_fermi!(y)
    return
end

function Base.:*(A::Dirac_operator, x::FermionFields)
    y = similar(x)
    mul!(y, A, x)
    return y
end

function Base.:*(A::DdagD_operator, x::FermionFields)
    y = similar(x)
    mul!(y, A, x)
    return y
end

function LinearAlgebra.mul!(y::WilsonFermion, A::Adjoint_Wilson_operator, x::WilsonFermion) #y = A*x
    Wdagx!(y, A.parent.U, x, A.parent._temporal_fermi)
end

function LinearAlgebra.mul!(
    y::WilsonFermion,
    A::Adjoint_WilsonClover_operator,
    x::WilsonFermion,
) #y = A*x
    Wdagx!(y, A.parent.U, x, A.parent._temporal_fermi, A.parent.CloverFμν)
end

function LinearAlgebra.mul!(
    y::StaggeredFermion,
    A::Adjoint_Staggered_operator,
    x::StaggeredFermion,
) #y = A*x
    temps = A.parent._temporal_fermi
    temp = temps[4]
    Dx!(temp, A.parent.U, x, temps)
    clear!(y)
    add!(y, A.parent.mass, x, -1, temp)
    set_wing_fermi!(y)
    return
end

function LinearAlgebra.mul!(
    y::DomainwallFermion,
    A::Adjoint_D5DW_Domainwall_operator,
    x::DomainwallFermion,
) #y = A*x
    D5DWdagx!(y, A.parent.U, x, A.parent.m, A.parent.wilsonoperator._temporal_fermi)
end



function Base.:*(A::Adjoint_Dirac_operator, x::WilsonFermion)
    y = similar(x)
    mul!(y, A, x)
    return y
end

function bicg(
    x,
    A::Domainwall_operator,
    b;
    eps = 1e-10,
    maxsteps = 1000,
    verbose = Verbose_2(),
) #A*x = b -> x = A^-1*b
    #A = D5DW(m)*D5DW(m=1))^(-1)
    #A^-1 = D5DW(m=1)*DsDW(m)^-1
    #x = A^-1*b = D5DW(m=1)*DsDW(m)^-1*b
    bicg(
        A.D5DW._temporal_fermi[1],
        A.D5DW,
        b;
        eps = eps,
        maxsteps = maxsteps,
        verbose = verbose,
    )
    mul!(x, A.D5DW_PV, A.D5DW._temporal_fermi[1])
end

function cg(
    x,
    A::DdagD_Domainwall_operator,
    b;
    eps = 1e-10,
    maxsteps = 1000,
    verbose = Verbose_2(),
)
    #=
    A^-1 = ( (D5DW(m)*D5DW(m=1))^(-1))^+ D5DW(m)*D5DW(m=1))^(-1) )^-1
      = ( D5DW(m=1)^+)^(-1) D5DW(m)^+ D5DW(m)*D5DW(m=1))^(-1) )^-1
      = D5DW(m=1) (  D5DW(m)^+ D5DW(m) )^(-1) )^-1 D5DW(m=1)^+
    x = A^-1*b = D5DW(m=1) (  D5DW(m)^+ D5DW(m) )^(-1)  D5DW(m=1)^+*b
    =#
    mul!(A.dirac.D5DW_PV._temporal_fermi[1], A.dirac.DSDW_PV', b) #D5DW(m=1)^+*b

    temp = A.dirac.D5DW_PV._temporal_fermi[1]
    cg(
        A.dirac.D5DW._temporal_fermi[1],
        A.DdagD,
        temp;
        eps = eps,
        maxsteps = maxsteps,
        verbose = verbose,
    ) #(  D5DW(m)^+ D5DW(m) )^(-1)  D5DW(m=1)^+*b
    mul!(x, A.dirac.D5DW_PV, A.dirac.D5DW._temporal_fermi[1])
end


function make_densematrix(A::T) where {T<:Operator}
    x = get_temporal_fermi(A)[1]
    xi = similar(x)
    clear!(x)
    j = 0
    NV = length(x)
    A_dense = zeros(ComplexF64, NV, NV)
    println("Matrix size: ", NV)
    for i = 1:NV
        clear!(x)
        j += 1
        x[i] = 1
        set_wing_fermi!(x)

        #println("x*x ",x*x)

        mul!(xi, A, x)
        #for j=1:NV
        #    println("$i $j ",xi[j])
        #end


        #error("i = $i")
        println("i = $i ", xi * xi)
        substitute_fermion!(A_dense, j, xi)
        #=
        for j=1:NV
            if abs(A_dense[i,j]) != 0
                println("$i $j $(A_dense[i,j])")
            end
        end
        =#
    end
    return A_dense

end
end
