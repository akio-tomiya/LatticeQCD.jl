module LieAlgebrafields
using Random
using LinearAlgebra
import ..Gaugefields:
    SU3GaugeFields,
    SU3GaugeFields_1d,
    SU2GaugeFields,
    SU2GaugeFields_1d,
    GaugeFields,
    GaugeFields_1d,
    make_staple_double!,
    substitute!,
    projlink!,
    set_wing!,
    evaluate_wilson_loops!,
    muladd!,
    SUNGaugeFields_1d,
    SUNGaugeFields,
    U1GaugeFields_1d,
    U1GaugeFields,
    evaluate_tensor_lines,
    SU
import ..SUN_generator: Generator, lie2matrix!, matrix2lie!
import ..Gaugefields
import ..Actions:
    GaugeActionParam_autogenerator, SmearingParam_single, SmearingParam_multi, SmearingParam
import ..Wilsonloops: make_plaq_staple_prime, make_plaq_staple

abstract type LieAlgebraFields end

struct SU3AlgebraFields <: LieAlgebraFields
    a::Array{Float64,5}
    NX::Int64
    NY::Int64
    NZ::Int64
    NT::Int64
    NC::Int64
    NumofBasis::Int64

    function SU3AlgebraFields(NC, NX, NY, NZ, NT)
        @assert NC == 3
        NumofBasis = 8
        return new(
            zeros(Float64, NumofBasis, NX, NY, NZ, NT),
            NX,
            NY,
            NZ,
            NT,
            NC,
            NumofBasis,
        )
    end
end

struct SU2AlgebraFields <: LieAlgebraFields
    a::Array{Float64,5}
    NX::Int64
    NY::Int64
    NZ::Int64
    NT::Int64
    NC::Int64
    NumofBasis::Int64

    function SU2AlgebraFields(NC, NX, NY, NZ, NT)
        @assert NC == 2
        NumofBasis = 3
        return new(
            zeros(Float64, NumofBasis, NX, NY, NZ, NT),
            NX,
            NY,
            NZ,
            NT,
            NC,
            NumofBasis,
        )
    end
end

struct SUNAlgebraFields <: LieAlgebraFields
    a::Array{Float64,5}
    NX::Int64
    NY::Int64
    NZ::Int64
    NT::Int64
    NC::Int64
    NumofBasis::Int64
    generators::Generator

    function SUNAlgebraFields(NC, NX, NY, NZ, NT)
        NumofBasis = NC^2 - 1
        generators = Generator(NC)
        return new(
            zeros(Float64, NumofBasis, NX, NY, NZ, NT),
            NX,
            NY,
            NZ,
            NT,
            NC,
            NumofBasis,
            generators,
        )
    end
end

struct U1AlgebraFields <: LieAlgebraFields
    a::Array{Float64,5}
    NX::Int64
    NY::Int64
    NZ::Int64
    NT::Int64
    NC::Int64
    NumofBasis::Int64
    #generators::Generator

    function U1AlgebraFields(NC, NX, NY, NZ, NT)
        NumofBasis = 1
        #generators = Generator(NC)
        return new(
            zeros(Float64, NumofBasis, NX, NY, NZ, NT),
            NX,
            NY,
            NZ,
            NT,
            NC,
            NumofBasis,
        )
    end
end


function Base.setindex!(x::LieAlgebraFields, v, i...)
    x.a[i...] = v
end

function Base.getindex(x::LieAlgebraFields, i...)
    return x.a[i...]
end

function LieAlgebraFields(NC, NX, NY, NZ, NT)
    if NC == 3
        return SU3AlgebraFields(NC, NX, NY, NZ, NT)
    elseif NC == 2
        return SU2AlgebraFields(NC, NX, NY, NZ, NT)
    elseif NC ≥ 4
        return SUNAlgebraFields(NC, NX, NY, NZ, NT)
    elseif NC == 1
        return U1AlgebraFields(NC, NX, NY, NZ, NT)
    end
end

function Base.similar(x::LieAlgebraFields)
    return LieAlgebraFields(x.NC, x.NX, x.NY, x.NZ, x.NT)
end

function Gaugefields.substitute!(x::LieAlgebraFields, pwork)
    n1, n2, n3, n4, n5 = size(x.a)
    #println(size(pwork))
    x.a[:, :, :, :, :] = reshape(pwork, (n1, n2, n3, n4, n5))
end

function Base.display(a::LieAlgebraFields)
    n1, n2, n3, n4, n5 = size(a.a)
    for i5 = 1:n5
        for i4 = 1:n4
            for i3 = 1:n3
                for i2 = 1:n2
                    for i1 = 1:n1
                        println("$i1 $i2 $i3 $i4 $i5 ", a.a[i1, i2, i3, i4, i5])# = α*a.a[i1,i2,i3,i4,i5]
                    end
                end
            end
        end
    end
    return
end

function clear!(a::Array{T,1}) where {T<:LieAlgebraFields}
    for μ = 1:4
        clear!(a[μ])
    end
    return
end

function clear!(a::LieAlgebraFields)
    n1, n2, n3, n4, n5 = size(a.a)
    for i5 = 1:n5
        for i4 = 1:n4
            for i3 = 1:n3
                for i2 = 1:n2
                    for i1 = 1:n1
                        a.a[i1, i2, i3, i4, i5] = 0
                    end
                end
            end
        end
    end
    return
end

function LinearAlgebra.mul!(c::LieAlgebraFields, α::Number, a::LieAlgebraFields)
    #n1,n2,n3,n4,n5 = size(a.a)
    @. c.a = α * a.a

    return
    n1, n2, n3, n4, n5 = size(a.a)

    for i5 = 1:n5
        for i4 = 1:n4
            for i3 = 1:n3
                for i2 = 1:n2
                    for i1 = 1:n1
                        c.a[i1, i2, i3, i4, i5] = α * a.a[i1, i2, i3, i4, i5]
                    end
                end
            end
        end
    end
    return
end

function Base.:*(x::Array{T,1}, y::Array{T,1}) where {T<:LieAlgebraFields}
    s = 0
    for μ = 1:4
        s += x[μ] * y[μ]
    end

    return s
end


function Base.:*(x::T, y::T) where {T<:LieAlgebraFields}
    n1, n2, n3, n4, n5 = size(x.a)
    s = 0
    for i5 = 1:n5
        for i4 = 1:n4
            for i3 = 1:n3
                for i2 = 1:n2
                    for i1 = 1:n1
                        s += x.a[i1, i2, i3, i4, i5] * y.a[i1, i2, i3, i4, i5]
                    end
                end
            end
        end
    end
    return s
end


function Base.:*(x::SU3AlgebraFields, y::SU3AlgebraFields)
    n1, n2, n3, n4, n5 = size(x.a)
    s = 0
    for i5 = 1:n5
        for i4 = 1:n4
            for i3 = 1:n3
                for i2 = 1:n2
                    for i1 = 1:n1
                        s += x.a[i1, i2, i3, i4, i5] * y.a[i1, i2, i3, i4, i5]
                    end
                end
            end
        end
    end
    return s
end

function Base.:*(x::SU2AlgebraFields, y::SU2AlgebraFields)
    n1, n2, n3, n4, n5 = size(x.a)
    s = 0
    for i5 = 1:n5
        for i4 = 1:n4
            for i3 = 1:n3
                for i2 = 1:n2
                    for i1 = 1:n1
                        s += x.a[i1, i2, i3, i4, i5] * y.a[i1, i2, i3, i4, i5]
                    end
                end
            end
        end
    end
    return s
end


function add!(a::LieAlgebraFields, α, b::LieAlgebraFields)
    n1, n2, n3, n4, n5 = size(a.a)
    for i5 = 1:n5
        for i4 = 1:n4
            for i3 = 1:n3
                for i2 = 1:n2
                    for i1 = 1:n1
                        a.a[i1, i2, i3, i4, i5] += α * b.a[i1, i2, i3, i4, i5]
                    end
                end
            end
        end
    end
    return
end

function gauss_distribution_lie!(p::Array{LieAlgebraFields,1})
    for i = 1:length(p)
        gauss_distribution_lie!(p[i])
    end
end

function gauss_distribution_lie!(p::LieAlgebraFields)
    n1, n2, n3, n4, n5 = size(p.a)

    for i5 = 1:n5
        for i4 = 1:n4
            for i3 = 1:n3
                for i2 = 1:n2
                    for i1 = 1:n1
                        p.a[i1, i2, i3, i4, i5] = randn()
                    end
                end
            end
        end
    end

    return
end

const sr3 = sqrt(3)
const sr3i = 1 / sr3
const sr3i2 = 2 * sr3i

function Gauge2Lie!(c::SU3AlgebraFields, x::SU3GaugeFields)
    NX = x.NX
    NY = x.NY
    NZ = x.NZ
    NT = x.NT

    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    x11 = x[1, 1, ix, iy, iz, it]
                    x12 = x[1, 2, ix, iy, iz, it]
                    x13 = x[1, 3, ix, iy, iz, it]
                    x21 = x[2, 1, ix, iy, iz, it]
                    x22 = x[2, 2, ix, iy, iz, it]
                    x23 = x[2, 3, ix, iy, iz, it]
                    x31 = x[3, 1, ix, iy, iz, it]
                    x32 = x[3, 2, ix, iy, iz, it]
                    x33 = x[3, 3, ix, iy, iz, it]

                    c[1, ix, iy, iz, it] = (imag(x12) + imag(x21))
                    c[2, ix, iy, iz, it] = (real(x12) - real(x21))
                    c[3, ix, iy, iz, it] = (imag(x11) - imag(x22))
                    c[4, ix, iy, iz, it] = (imag(x13) + imag(x31))
                    c[5, ix, iy, iz, it] = (real(x13) - real(x31))

                    c[6, ix, iy, iz, it] = (imag(x23) + imag(x32))
                    c[7, ix, iy, iz, it] = (real(x23) - real(x32))
                    c[8, ix, iy, iz, it] = sr3i * (imag(x11) + imag(x22) - 2 * imag(x33))
                end
            end

        end

    end


end

function Gauge2Lie!(c::SU2AlgebraFields, x::SU2GaugeFields)
    NX = x.NX
    NY = x.NY
    NZ = x.NZ
    NT = x.NT

    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX


                    x11 = x[1, 1, ix, iy, iz, it]
                    x12 = x[1, 2, ix, iy, iz, it]
                    x21 = x[2, 1, ix, iy, iz, it]
                    x22 = x[2, 2, ix, iy, iz, it]

                    c[1, ix, iy, iz, it] = (imag(x12) + imag(x21))
                    c[2, ix, iy, iz, it] = (real(x12) - real(x21))
                    c[3, ix, iy, iz, it] = (imag(x11) - imag(x22))

                end
            end

        end

    end


end

function Gauge2Lie!(c::SUNAlgebraFields, x::SUNGaugeFields)
    NX = x.NX
    NY = x.NY
    NZ = x.NZ
    NT = x.NT
    g = c.generators
    NC = c.NC
    matrix = zeros(ComplexF64, NC, NC)
    a = zeros(ComplexF64, length(g))

    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    for k2 = 1:NC
                        for k1 = 1:NC
                            matrix[k1, k2] = x[k1, k2, ix, iy, iz, it]
                        end
                    end

                    matrix2lie!(a, g, matrix)
                    for k = 1:length(g)
                        c[k, ix, iy, iz, it] = 2 * imag(a[k])
                    end

                end
            end

        end

    end


end


function Gauge2Lie!(c::SU3AlgebraFields, x::SU3GaugeFields_1d)
    NX = c.NX
    NY = c.NY
    NZ = c.NZ
    NT = c.NT

    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    i = (((it - 1) * NX + iz - 1) * NY + iy - 1) * NX + ix

                    x11 = x[1, 1, i]
                    x12 = x[1, 2, i]
                    x13 = x[1, 3, i]
                    x21 = x[2, 1, i]
                    x22 = x[2, 2, i]
                    x23 = x[2, 3, i]
                    x31 = x[3, 1, i]
                    x32 = x[3, 2, i]
                    x33 = x[3, 3, i]

                    c[1, ix, iy, iz, it] = (imag(x12) + imag(x21))
                    c[2, ix, iy, iz, it] = (real(x12) - real(x21))
                    c[3, ix, iy, iz, it] = (imag(x11) - imag(x22))
                    c[4, ix, iy, iz, it] = (imag(x13) + imag(x31))
                    c[5, ix, iy, iz, it] = (real(x13) - real(x31))

                    c[6, ix, iy, iz, it] = (imag(x23) + imag(x32))
                    c[7, ix, iy, iz, it] = (real(x23) - real(x32))
                    c[8, ix, iy, iz, it] = sr3i * (imag(x11) + imag(x22) - 2 * imag(x33))


                end
            end

        end

    end


end

function Gauge2Lie!(c::SU2AlgebraFields, x::SU2GaugeFields_1d)
    NX = c.NX
    NY = c.NY
    NZ = c.NZ
    NT = c.NT

    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    i = (((it - 1) * NX + iz - 1) * NY + iy - 1) * NX + ix

                    x11 = x[1, 1, i]
                    x12 = x[1, 2, i]
                    #x13 = x[1,3,i]
                    x21 = x[2, 1, i]
                    x22 = x[2, 2, i]


                    c[1, ix, iy, iz, it] = (imag(x12) + imag(x21))
                    c[2, ix, iy, iz, it] = (real(x12) - real(x21))
                    c[3, ix, iy, iz, it] = (imag(x11) - imag(x22))

                end
            end

        end

    end


end

function Gauge2Lie!(c::SUNAlgebraFields, x::SUNGaugeFields_1d)
    NX = c.NX
    NY = c.NY
    NZ = c.NZ
    NT = c.NT

    g = c.generators
    NC = c.NC
    matrix = zeros(ComplexF64, NC, NC)
    a = zeros(ComplexF64, length(g))

    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    i = (((it - 1) * NX + iz - 1) * NY + iy - 1) * NX + ix

                    for k2 = 1:NC
                        for k1 = 1:NC
                            matrix[k1, k2] = x[k1, k2, i]
                        end
                    end
                    matrix2lie!(a, g, matrix)
                    for k = 1:length(g)
                        c[k, ix, iy, iz, it] = 2 * imag(a[k])
                    end

                end
            end

        end

    end


end

function Gauge2Lie!(c::U1AlgebraFields, x::U1GaugeFields_1d)
    NX = c.NX
    NY = c.NY
    NZ = c.NZ
    NT = c.NT

    NC = c.NC


    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    i = (((it - 1) * NX + iz - 1) * NY + iy - 1) * NX + ix
                    c[1, ix, iy, iz, it] = real(x[1, 1, i])
                end
            end

        end

    end


end




function expF_U!(
    U::Array{T,1},
    F::Array{N,1},
    Δτ,
    temps::Array{T_1d,1},
    temp_a::N,
) where {T<:GaugeFields,N<:LieAlgebraFields,T_1d<:GaugeFields_1d}
    @assert Δτ != 0 "Δτ should not be zero in expF_U! function!"
    temp1 = temps[1]
    temp2 = temps[2]
    temp3 = temps[3]
    temp4 = temps[4]

    c = temp_a

    for μ = 1:4
        substitute!(temp1, U[μ])

        mul!(c, Δτ, F[μ]) #c = Δτ*F

        expA!(temp2, c, temp3, temp4) #temp2 = exp(c)

        mul!(temp3, temp2, temp1) #temp3 = temp2*temp1

        substitute!(U[μ], temp3)
        set_wing!(U[μ])

    end
end


function add_gaugeforce!(
    F::Array{N,1},
    U::Array{T,1},
    temps::Array{T_1d,1},
    temp_a::N;
    fac = 1,
) where {T<:GaugeFields,N<:LieAlgebraFields,T_1d<:GaugeFields_1d}
    temp1 = temps[1]
    temp2 = temps[2]
    temp3 = temps[3]
    staple = temps[4]

    c = temp_a

    for μ = 1:4
        #println(1)
        make_staple_double!(staple, U, μ, temp1, temp2, temp3)

        substitute!(temp1, U[μ])

        #=
        c
        c                   ---
        c                  |   |           
        c            tmp2 = ---    +   ---
        c                             |   |
        c                              ---
        c
        =#

        mul!(temp2, temp1, staple')




        #.....   Projection onto Lie Algebra   .....

        projlink!(temp3, temp2)

        #display(temp3)
        #exit()


        Gauge2Lie!(c, temp3)

        #display(c)
        #exit()

        #c     .....   p(new) = p(old) + fac * c  .....
        add!(F[μ], fac, c)

    end
    return
end

function add_gaugeforce!(
    F::Array{N,1},
    U::Array{T,1},
    temps::Array{T_1d,1},
    temp_a::N,
    gparam::GaugeActionParam_autogenerator;
    fac = 1,
) where {T<:GaugeFields,N<:LieAlgebraFields,T_1d<:GaugeFields_1d}
    temp1 = temps[1]
    temp2 = temps[2]
    temp3 = temps[3]
    staple = temps[4]
    Utemp = temps[5]
    temp4 = temps[6]

    c = temp_a

    #=
    NX = U[1].NX
    NY = U[1].NY
    NZ = U[1].NZ
    NT = U[1].NT
    NC = U[1].NC
    V = zeros(ComplexF64,NC,NC)
    Vtemp = zeros(ComplexF64,NC,NC)
    =#


    for μ = 1:4
        substitute!(Utemp, U[μ])
        Gaugefields.clear!(staple)
        for i = 1:gparam.numactions
            evaluate_wilson_loops!(temp4, gparam.staples[i][μ], U, [temp1, temp2, temp3])
            muladd!(staple, gparam.βs[i] / gparam.β, temp4)


            #ix,iy,iz,it = 1,1,1,1
            #icum = (((it-1)*NZ+iz-1)*NY+iy-1)*NX+ix
            #evaluate_wilson_loops!(Vtemp,gparam.staples[i][μ],U,ix,iy,iz,it)
            #println("μ = $μ, i = $i")
            #println("all")
            #display(temp4[:,:,icum])
            #println("\t")
            #println("partial")
            #display(Vtemp)
            #println("\t")
        end
        #exit()

        #substitute!(temp1,U[μ])
        #mul!(temp2,temp1,staple)

        mul!(temp2, Utemp, staple)

        projlink!(temp3, temp2)
        Gauge2Lie!(c, temp3)
        add!(F[μ], fac, c)

    end
    return
end


const pi23 = 2pi / 3
const tinyvalue = 1e-100

function expA!(v::SU2GaugeFields_1d, u::SU2AlgebraFields, temp1, temp2)
    NX = u.NX
    NY = u.NY
    NZ = u.NZ
    NT = u.NT
    NC = u.NC
    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    icum = (((it - 1) * NX + iz - 1) * NY + iy - 1) * NX + ix
                    u1 = u[1, ix, iy, iz, it] / 2
                    u2 = u[2, ix, iy, iz, it] / 2
                    u3 = u[3, ix, iy, iz, it] / 2
                    R = sqrt(u1^2 + u2^2 + u3^2) + tinyvalue
                    sR = sin(R) / R
                    #sR = ifelse(R == 0,1,sR)
                    a0 = cos(R)
                    a1 = u1 * sR
                    a2 = u2 * sR
                    a3 = u3 * sR

                    v[1, 1, icum] = cos(R) + im * a3
                    v[1, 2, icum] = im * a1 + a2
                    v[2, 1, icum] = im * a1 - a2
                    v[2, 2, icum] = cos(R) - im * a3
                end
            end
        end
    end

end


function expA!(v::SUNGaugeFields_1d, u::SUNAlgebraFields, temp1, temp2)
    g = u.generators


    NX = u.NX
    NY = u.NY
    NZ = u.NZ
    NT = u.NT
    NC = u.NC
    u0 = zeros(ComplexF64, NC, NC)
    a = zeros(Float64, length(g))

    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    icum = (((it - 1) * NX + iz - 1) * NY + iy - 1) * NX + ix
                    for k = 1:length(a)
                        a[k] = u[k, ix, iy, iz, it]
                    end

                    lie2matrix!(u0, g, a)
                    v[:, :, icum] = exp((im / 2) * u0)
                end
            end
        end
    end

end

function expA!(v::U1GaugeFields_1d, u::U1AlgebraFields, temp1, temp2)


    NX = u.NX
    NY = u.NY
    NZ = u.NZ
    NT = u.NT
    NC = u.NC

    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    icum = (((it - 1) * NX + iz - 1) * NY + iy - 1) * NX + ix
                    v[1, 1, icum] = exp(im * u[1, ix, iy, iz, it])
                end
            end
        end
    end

end



"""
    The descrition is based on the LTK in fortran.
c----------------------------------------------------------------c
c     Project Lie algebra u to Lie group element g.
c     U = u(1)*lambda(1)/2 + u(2)*lambda(2)/2) + ....
c     g = exp(i*U) 
c     ( by exact diagonalization                           )
c     (  X = Adj(V) * D * V                                )
c     (  where D = diag(d1,d2,d3)i                         )
c     (  g = Adj(V) * diag(exp(i*d1),exp(i*d2),exp(i*d3))) )
c----------------------------------------------------------------c
c for i = 1 to n, v(i) is the SU(3) projection of Lie algebra element u(i)
c     !!!!!   This is for SU(3)   !!!!!
c----------------------------------------------------------------c
    """

function expA!(v::SU3GaugeFields_1d, u::SU3AlgebraFields, temp1, temp2)
    ww = temp1
    w = temp2

    #ww = similar(v)
    #w =similar(v)

    algemat!(v, u)


    for i = 1:v.NV
        v1 = real(v[1, 1, i])
        v2 = imag(v[1, 1, i])
        v3 = real(v[1, 2, i])
        v4 = imag(v[1, 2, i])
        v5 = real(v[1, 3, i])
        v6 = imag(v[1, 3, i])
        v7 = real(v[2, 1, i])
        v8 = imag(v[2, 1, i])
        v9 = real(v[2, 2, i])
        v10 = imag(v[2, 2, i])
        v11 = real(v[2, 3, i])
        v12 = imag(v[2, 3, i])
        v13 = real(v[3, 1, i])
        v14 = imag(v[3, 1, i])
        v15 = real(v[3, 2, i])
        v16 = imag(v[3, 2, i])
        v17 = real(v[3, 3, i])
        v18 = imag(v[3, 3, i])

        #c find eigenvalues of v
        trv3 = (v1 + v9 + v17) / 3.0
        cofac = v1 * v9 - v3^2 - v4^2 + v1 * v17 - v5^2 - v6^2 + v9 * v17 - v11^2 - v12^2
        det =
            v1 * v9 * v17 - v1 * (v11^2 + v12^2) - v9 * (v5^2 + v6^2) -
            v17 * (v3^2 + v4^2) +
            (v5 * (v3 * v11 - v4 * v12) + v6 * (v3 * v12 + v4 * v11)) * 2.0
        p3 = cofac / 3.0 - trv3^2
        q = trv3 * cofac - det - 2.0 * trv3^3
        x = sqrt(-4.0 * p3) + tinyvalue
        arg = q / (x * p3)

        arg = min(1, max(-1, arg))
        theta = acos(arg) / 3.0
        e1 = x * cos(theta) + trv3
        theta = theta + pi23
        e2 = x * cos(theta) + trv3
        #       theta = theta + pi23
        #       e3 = x * cos(theta) + trv3
        e3 = 3.0 * trv3 - e1 - e2

        # solve for eigenvectors

        w1 = v5 * (v9 - e1) - v3 * v11 + v4 * v12
        w2 = -v6 * (v9 - e1) + v4 * v11 + v3 * v12
        w3 = (v1 - e1) * v11 - v3 * v5 - v4 * v6
        w4 = -(v1 - e1) * v12 - v4 * v5 + v3 * v6
        w5 = -(v1 - e1) * (v9 - e1) + v3^2 + v4^2
        w6 = 0.0

        coeff = 1.0 / sqrt(w1^2 + w2^2 + w3^2 + w4^2 + w5^2)

        w1 = w1 * coeff
        w2 = w2 * coeff
        w3 = w3 * coeff
        w4 = w4 * coeff
        w5 = w5 * coeff

        w7 = v5 * (v9 - e2) - v3 * v11 + v4 * v12
        w8 = -v6 * (v9 - e2) + v4 * v11 + v3 * v12
        w9 = (v1 - e2) * v11 - v3 * v5 - v4 * v6
        w10 = -(v1 - e2) * v12 - v4 * v5 + v3 * v6
        w11 = -(v1 - e2) * (v9 - e2) + v3^2 + v4^2
        w12 = 0.0

        coeff = 1.0 / sqrt(w7^2 + w8^2 + w9^2 + w10^2 + w11^2)
        w7 = w7 * coeff
        w8 = w8 * coeff
        w9 = w9 * coeff
        w10 = w10 * coeff
        w11 = w11 * coeff

        w13 = v5 * (v9 - e3) - v3 * v11 + v4 * v12
        w14 = -v6 * (v9 - e3) + v4 * v11 + v3 * v12
        w15 = (v1 - e3) * v11 - v3 * v5 - v4 * v6
        w16 = -(v1 - e3) * v12 - v4 * v5 + v3 * v6
        w17 = -(v1 - e3) * (v9 - e3) + v3^2 + v4^2
        w18 = 0.0

        coeff = 1.0 / sqrt(w13^2 + w14^2 + w15^2 + w16^2 + w17^2)
        w13 = w13 * coeff
        w14 = w14 * coeff
        w15 = w15 * coeff
        w16 = w16 * coeff
        w17 = w17 * coeff

        # construct the projection v
        c1 = cos(e1)
        s1 = sin(e1)
        ww1 = w1 * c1 - w2 * s1
        ww2 = w2 * c1 + w1 * s1
        ww3 = w3 * c1 - w4 * s1
        ww4 = w4 * c1 + w3 * s1
        ww5 = w5 * c1 - w6 * s1
        ww6 = w6 * c1 + w5 * s1

        c2 = cos(e2)
        s2 = sin(e2)
        ww7 = w7 * c2 - w8 * s2
        ww8 = w8 * c2 + w7 * s2
        ww9 = w9 * c2 - w10 * s2
        ww10 = w10 * c2 + w9 * s2
        ww11 = w11 * c2 - w12 * s2
        ww12 = w12 * c2 + w11 * s2

        c3 = cos(e3)
        s3 = sin(e3)
        ww13 = w13 * c3 - w14 * s3
        ww14 = w14 * c3 + w13 * s3
        ww15 = w15 * c3 - w16 * s3
        ww16 = w16 * c3 + w15 * s3
        ww17 = w17 * c3 - w18 * s3
        ww18 = w18 * c3 + w17 * s3

        w[1, 1, i] = w1 + im * w2
        w[1, 2, i] = w3 + im * w4
        w[1, 3, i] = w5 + im * w6
        w[2, 1, i] = w7 + im * w8
        w[2, 2, i] = w9 + im * w10
        w[2, 3, i] = w11 + im * w12
        w[3, 1, i] = w13 + im * w14
        w[3, 2, i] = w15 + im * w16
        w[3, 3, i] = w17 + im * w18

        ww[1, 1, i] = ww1 + im * ww2
        ww[1, 2, i] = ww3 + im * ww4
        ww[1, 3, i] = ww5 + im * ww6
        ww[2, 1, i] = ww7 + im * ww8
        ww[2, 2, i] = ww9 + im * ww10
        ww[2, 3, i] = ww11 + im * ww12
        ww[3, 1, i] = ww13 + im * ww14
        ww[3, 2, i] = ww15 + im * ww16
        ww[3, 3, i] = ww17 + im * ww18
    end

    mul!(v, w', ww)
    return



end




const sr3 = sqrt(3)
const sr3i = 1 / sr3
const sr3i2 = 2 * sr3i

function algemat!(x::SU3GaugeFields_1d, t::SU3AlgebraFields)
    NX = t.NX
    NY = t.NY
    NZ = t.NZ
    NT = t.NT
    NC = t.NC
    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    icum = (((it - 1) * NX + iz - 1) * NY + iy - 1) * NX + ix



                    c1 = t[1, ix, iy, iz, it] * 0.5
                    c2 = t[2, ix, iy, iz, it] * 0.5
                    c3 = t[3, ix, iy, iz, it] * 0.5
                    c4 = t[4, ix, iy, iz, it] * 0.5
                    c5 = t[5, ix, iy, iz, it] * 0.5
                    c6 = t[6, ix, iy, iz, it] * 0.5
                    c7 = t[7, ix, iy, iz, it] * 0.5
                    c8 = t[8, ix, iy, iz, it] * 0.5

                    x[1, 1, icum] = c3 + sr3i * c8 + im * (0.0)
                    x[1, 2, icum] = c1 + im * (-c2)
                    x[1, 3, icum] = c4 + im * (-c5)

                    x[2, 1, icum] = c1 + im * (c2)
                    x[2, 2, icum] = -c3 + sr3i * c8 + im * (0.0)
                    x[2, 3, icum] = c6 + im * (-c7)

                    x[3, 1, icum] = c4 + im * (c5)
                    x[3, 2, icum] = c6 + im * (c7)
                    x[3, 3, icum] = -sr3i2 * c8 + im * (0.0)
                end

            end

        end
    end
end


function stoutfource(dSdU, U, smearing::T) where {T<:SmearingParam_single}
    dSdUnew = deepcopy(dSdU)
    dSdρs = stoutfource!(dSdUnew, dSdU, U, smearing, smearing.ρs)
    return dSdUnew, dSdρs
end

#evaluate_tensor_lines!(V,nu,dSdU::Array{T_1d,1},gparam::GaugeActionParam_autogenerator,U::Array{T,1},umu::Tuple{I,I},ix,iy,iz,it) where {T <: GaugeFields,I <: Int,T_1d <: GaugeFields_1d}
function stoutfource!(dSdUnew, dSdU, U, smearing::T, ρs) where {T<:SmearingParam}
    #nu = 1
    #evaluate_tensor_lines(nu,dSdU,smearing,U,(nu,1),ρs)


    #println("stout force!!")
    #dSdUnew = deepcopy(dSdU)

    dSdρs = zero(ρs)
    NX = U[1].NX
    NY = U[1].NY
    NZ = U[1].NZ
    NT = U[1].NT
    NC = U[1].NC
    V = zeros(ComplexF64, NC, NC)
    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    icum = (((it - 1) * NZ + iz - 1) * NY + iy - 1) * NX + ix
                    for nu = 1:4

                        #println(ρs)
                        #V,dSdρs_site = evaluate_tensor_lines_2(nu,dSdU,smearing,U,(nu,1),ix,iy,iz,it,ρs)
                        V, dSdρs_site = evaluate_tensor_lines(
                            nu,
                            dSdU,
                            smearing,
                            U,
                            (nu, 1),
                            ix,
                            iy,
                            iz,
                            it,
                            ρs,
                        )
                        #V = evaluate_tensor_lines_2(nu,dSdU,gparam,U,(nu,1),ρ,ix,iy,iz,it) 
                        #V = evaluate_SU2_force(nu,dSdU,gparam,U,(nu,1),ρ,ix,iy,iz,it) 
                        #V = evaluate_tensor_lines(nu,dSdU,gparam,U,(nu,1),ρ,ix,iy,iz,it) 
                        #println("dSdU\t")
                        #println("VV", V)
                        for i = 1:NC
                            for j = 1:NC
                                dSdUnew[nu][j, i, icum] = V[j, i]
                            end
                        end

                        @. dSdρs += dSdρs_site
                        #dSdUnew[nu][:,:,icum] = V2[:,:]

                        #=
                        if nu == 1
                           if (ix,iy,iz,it) == (1,1,1,1) || (ix,iy,iz,it) ==(1,2,2,2) || (ix,iy,iz,it) ==(2,2,2,2) || (ix,iy,iz,it) ==(4,4,4,4)
                               println("x ",(ix,iy,iz,it))
                           # #=
                               println("U ", U[nu][:,:,ix,iy,iz,it])

                               println("dS[Ufat]/dUfat ")
                               display(dSdU[nu][:,:,icum])
                               println("\t")
                               println("dS[Ufat]/dU ")
                               display(dSdUnew[nu][:,:,icum])
                               println("\t")
                               println("new2")
                               #display(V2)
                               println("\t")
                               println("old")
                               display(V)
                               println("\t")

                               #=
                               V = evaluate_tensor_lines(nu,dSdU,gparam,U,(nu,1),-ρ,ix,iy,iz,it) 
                               #V = evaluate_tensor_lines(nu,dSdU,gparam,U,(nu,1),-ρ,ix,iy,iz,it) 
                               println("new")
                               display(V)
                               println("\t")

                               V2 = evaluate_tensor_lines_2(nu,dSdU,gparam,U,(nu,1),-ρ,ix,iy,iz,it) 
                               println("new2")
                               display(V2)
                               println("\t")
                               =#
                               exit()

                            # =#
                           end
                        end
                        =#




                    end


                end
            end
        end
    end
    #exit()

    return dSdρs
    #exit()
    #exit()
    #return dSdUnew
end

function stoutfource(dSdU, U_multi, U, smearing::T) where {T<:SmearingParam_multi}
    dSdUnew = deepcopy(dSdU)
    dSdρs = stoutfource!(dSdUnew, dSdU, U_multi, U, smearing)
    return dSdUnew, dSdρs
end


function stoutfource!(
    dSdUnew,
    dSdU,
    U_multi,
    U::Array{GaugeFields{SU{NC}},1},
    smearing::T,
) where {NC,T<:SmearingParam_multi}
    num = length(U_multi)
    dSdU_m = deepcopy(dSdU)
    dSdU_m2 = deepcopy(dSdU)
    dSdρs = deepcopy(smearing.ρs)
    if num == 1
        dSdρs[1] = stoutfource!(dSdUnew, dSdU, U, smearing, smearing.ρs[1])
    else
        for i = num:-1:2
            #println(i)
            #println(smearing.ρs[i])
            dSdρs[i] = stoutfource!(dSdU_m2, dSdU_m, U_multi[i-1], smearing, smearing.ρs[i])
            dSdU_m, dSdU_m2 = dSdU_m2, dSdU_m
        end
        dSdρs[1] = stoutfource!(dSdUnew, dSdU_m, U, smearing, smearing.ρs[1])
    end
    return dSdρs

end

end
