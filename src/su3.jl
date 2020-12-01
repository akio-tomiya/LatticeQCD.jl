
struct SU3GaugeFields <: GaugeFields
    g::Array{ComplexF64,6}
    NX::Int64
    NY::Int64
    NZ::Int64
    NT::Int64
    NC::Int64
    NDW::Int64
    NV::Int64

    function SU3GaugeFields(NC,NDW,NX,NY,NZ,NT)
        @assert NC == 3
        NV = NX*NY*NZ*NT
        g = zeros(ComplexF64,NC,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW)
        return new(g,NX,NY,NZ,NT,NC,NDW,NV)
    end
end

struct SU3GaugeFields_1d <: GaugeFields_1d
    g::Array{ComplexF64,3}
    NC::Int64
    NV::Int64

    function SU3GaugeFields_1d(NC,NX,NY,NZ,NT)
        @assert NC == 3
        NV = NX*NY*NZ*NT
        g = zeros(ComplexF64,NC,NC,NV)
        return new(g,NC,NV)
    end

    function SU3GaugeFields_1d(NC,NV) 
        g = zeros(ComplexF64,NC,NC,NV)
        return new(g,NC,NV)
    end
end

function Base.similar(x::SU3GaugeFields_1d)
    return SU3GaugeFields_1d(x.NC,x.NV) 

end

function Base.setindex!(x::SU3GaugeFields,v,i1,i2,i3,i4,i5,i6) 
    x.g[i1,i2,i3 + x.NDW,i4 + x.NDW,i5 + x.NDW,i6 + x.NDW] = v
end

function Base.getindex(x::SU3GaugeFields,i1,i2,i3,i4,i5,i6)
    return x.g[i1,i2,i3 .+ x.NDW,i4 .+ x.NDW,i5 .+ x.NDW,i6 .+ x.NDW]
end

function Base.setindex!(x::SU3GaugeFields_1d,v,i1,i2,i3) 
    x.g[i1,i2,i3] = v
end

function Base.getindex(x::SU3GaugeFields_1d,i1,i2,i3)
    return x.g[i1,i2,i3]
end

struct Adjoint_SU3GaugeFields_1d <: Adjoint_GaugeFields_1d
    parent::SU3GaugeFields_1d
end

struct Adjoint_SU3GaugeFields <: Adjoint_GaugeFields
    parent::SU3GaugeFields
end

function Base.adjoint(x::SU3GaugeFields_1d)
    Adjoint_SU3GaugeFields_1d(x)
end

function Base.adjoint(x::SU3GaugeFields)
    Adjoint_SU3GaugeFields(x)
end

function IdentitySU3Gauges(NX,NY,NZ,NT,NDW)
    NC = 3
    U = SU3GaugeFields(NC,NDW,NX,NY,NZ,NT)
    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    for ic=1:3
                        U[ic,ic,ix,iy,iz,it] = 1 
                    end
                end
            end
        end
    end
    return U
end

function RandomSU3Gauges(NX,NY,NZ,NT,NDW)
    NC = 3
    U = SU3GaugeFields(NC,NDW,NX,NY,NZ,NT)
    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    U[1,1,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
                    U[1,2,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
                    U[1,3,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
                    U[2,1,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
                    U[2,2,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
                    U[2,3,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
                end
            end
        end
    end

    normalize!(U)
    return U
end

"""
normalize!(u)
c----------------------------------------------------------------------c
c     normalizes links                                                 c
c     input  x : arbitrary 3*3 complex matrix                          c
c     output x : SU(3) matrix 
"""
function normalize!(u::SU3GaugeFields)
    NX = u.NX
    NY = u.NY
    NZ = u.NZ
    NT = u.NT

    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    w1 = 0
                    w2 = 0
                    for ic=1:3
                        w1 += u[2,ic,ix,iy,iz,it]*conj(u[1,ic,ix,iy,iz,it])
                        w2 += u[1,ic,ix,iy,iz,it]*conj(u[1,ic,ix,iy,iz,it])
                    end
                    zerock2 = w2
                    if zerock2 == 0 
                        println("w2 is zero  !!  (in normlz)")
                        println("u[1,1),u[1,2),u[1,3) : ",u[1,1,ix,iy,iz,it], "\t",u[1,2,ix,iy,iz,it],"\t", u[1,3,ix,iy,iz,it])
                    end

                    w1 = -w1/w2

                    x4 = (u[2,1,ix,iy,iz,it]) + w1*u[1,1,ix,iy,iz,it]
                    x5 = (u[2,2,ix,iy,iz,it]) + w1*u[1,2,ix,iy,iz,it]
                    x6 = (u[2,3,ix,iy,iz,it]) + w1*u[1,3,ix,iy,iz,it]

                    w3 = x4*conj(x4) + x5*conj(x5) + x6*conj(x6)

                    zerock3 = w3
                    if zerock3 == 0
                        println("w3 is zero  !!  (in normlz)")
                        println("x4, x5, x6 : $x4, $x5, $x6")
                        exit()
                    end

                    u[2,1,ix,iy,iz,it] = x4
                    u[2,2,ix,iy,iz,it] = x5
                    u[2,3,ix,iy,iz,it] = x6

                    w3 = 1/sqrt(w3)
                    w2 = 1/sqrt(w2)

                    u[1,1,ix,iy,iz,it] = u[1,1,ix,iy,iz,it]*w2
                    u[1,2,ix,iy,iz,it] = u[1,2,ix,iy,iz,it]*w2
                    u[1,3,ix,iy,iz,it] = u[1,3,ix,iy,iz,it]*w2
                    u[2,1,ix,iy,iz,it] = u[2,1,ix,iy,iz,it]*w3
                    u[2,2,ix,iy,iz,it] = u[2,2,ix,iy,iz,it]*w3
                    u[2,3,ix,iy,iz,it] = u[2,3,ix,iy,iz,it]*w3

                    if zerock2*zerock3 == 0 
                        println("!! devided by zero !! (in normalize)")
                        println("w2 or w3 in normlz is zero !!")
                        println("w2, w3 : $w2, $w3   ")
                        exit()
                    end
                end
            end
        end
    end
end

function m3complv!(a::SU3GaugeFields)
    aa = zeros(Float64,18)
    NX = u.NX
    NY = u.NY
    NZ = u.NZ
    NT = u.NT

    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX

                    aa[ 1] = real( a[1,1,ix,iy,iz,it])
                    aa[ 2] = imag(a[1,1,ix,iy,iz,it])
                    aa[ 3] = real( a[1,2,ix,iy,iz,it])
                    aa[ 4] = imag(a[1,2,ix,iy,iz,it])
                    aa[ 5] = real( a[1,3,ix,iy,iz,it])
                    aa[ 6] = imag(a[1,3,ix,iy,iz,it])
                    aa[ 7] = real( a[2,1,ix,iy,iz,it])
                    aa[ 8] = imag(a[2,1,ix,iy,iz,it])
                    aa[ 9] = real( a[2,2,ix,iy,iz,it])
                    aa[10] = imag(a[2,2,ix,iy,iz,it])
                    aa[11] = real( a[2,3,ix,iy,iz,it])
                    aa[12] = imag(a[2,3,ix,iy,iz,it])

                    aa[13] = aa[ 3]*aa[11] - aa[ 4]*aa[12] -
                                aa[ 5]*aa[ 9] + aa[ 6]*aa[10]
                    aa[14] = aa[ 5]*aa[10] + aa[ 6]*aa[ 9] -
                                aa[ 3]*aa[12] - aa[ 4]*aa[11]
                    aa[15] = aa[ 5]*aa[ 7] - aa[ 6]*aa[ 8] -
                                aa[ 1]*aa[11] + aa[ 2]*aa[12]
                    aa[16] = aa[ 1]*aa[12] + aa[ 2]*aa[11] -
                                aa[ 5]*aa[ 8] - aa[ 6]*aa[ 7]
                    aa[17] = aa[ 1]*aa[ 9] - aa[ 2]*aa[10] -
                                aa[ 3]*aa[ 7] + aa[ 4]*aa[ 8]
                    aa[18] = aa[ 3]*aa[ 8] + aa[ 4]*aa[ 7] -
                                aa[ 1]*aa[10] - aa[ 2]*aa[ 9]

                    a[3,1,ix,iy,iz,it] = aa[13]+im*aa[14]
                    a[3,2,ix,iy,iz,it] = aa[15]+im*aa[16]
                    a[3,3,ix,iy,iz,it] = aa[17]+im*aa[18]
                end
            end
        end
    end
end

"""
    c-----------------------------------------------------c
c     !!!!!   vin and vout should be different vectors
c
c     Projectin of the etraceless antiermite part 
c     vout = x/2 - Tr(x)/6
c     wher   x = vin - Conjg(vin)      
c-----------------------------------------------------c
    """
    function projlink!(vout::SU3GaugeFields,vin::SU3GaugeFields)
        
        fac13 = 1/3
        NX = vin.NX
        NY = vin.NY
        NZ = vin.NZ
        NT = vin.NT

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX
                        v11 = vin[1,1,ix,iy,iz,it]
                        v22 = vin[2,2,ix,iy,iz,it]
                        v33 = vin[3,3,ix,iy,iz,it]

                        tri = fac13*(imag(v11)+imag(v22)+imag(v33))

                        vout[1,1,ix,iy,iz,it] = (imag(v11)-tri)*im
                        vout[2,2,ix,iy,iz,it] = (imag(v22)-tri)*im
                        vout[3,3,ix,iy,iz,it] = (imag(v33)-tri)*im
                    end
                end
            end
        end

        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for ix=1:NX

                        v12 = vin[1,2,ix,iy,iz,it]
                        v13 = vin[1,3,ix,iy,iz,it]
                        v21 = vin[2,1,ix,iy,iz,it]
                        v23 = vin[2,3,ix,iy,iz,it]
                        v31 = vin[3,1,ix,iy,iz,it]
                        v32 = vin[3,2,ix,iy,iz,it]

                        x12 = v12 - conj(v21)
                        x13 = v13 - conj(v31)
                        x23 = v23 - conj(v32)
                    
                        x21 = - conj(x12)
                        x31 = - conj(x13)
                        x32 = - conj(x23)

                        vout[1,2,ix,iy,iz,it] = 0.5  * x12
                        vout[1,3,ix,iy,iz,it] = 0.5  * x13
                        vout[2,1,ix,iy,iz,it] = 0.5  * x21
                        vout[2,3,ix,iy,iz,it] = 0.5  * x23
                        vout[3,1,ix,iy,iz,it] = 0.5  * x31
                        vout[3,2,ix,iy,iz,it] = 0.5  * x32
                    end
                end
            end
        end

    end

    const sr3 = sqrt(3)
    const sr3i = 1/sr3

    function projlink!(vout::SU3GaugeFields_1d,vin::SU3GaugeFields_1d)
        fac13 = 1/3
        nv = vin.NV

        for i=1:nv
            v11 = vin[1,1,i]
            v22 = vin[2,2,i]
            v33 = vin[3,3,i]

            tri = fac13*(imag(v11)+imag(v22)+imag(v33))

            vout[1,1,i] = (imag(v11)-tri)*im
            vout[2,2,i] = (imag(v22)-tri)*im
            vout[3,3,i] = (imag(v33)-tri)*im

        end

        for i=1:nv
            v12 = vin[1,2,i]
            v13 = vin[1,3,i]
            v21 = vin[2,1,i]
            v23 = vin[2,3,i]
            v31 = vin[3,1,i]
            v32 = vin[3,2,i]

            x12 = v12 - conj(v21)
            x13 = v13 - conj(v31)
            x23 = v23 - conj(v32)
        
            x21 = - conj(x12)
            x31 = - conj(x13)
            x32 = - conj(x23)

            vout[1,2,i] = 0.5  * x12
            vout[1,3,i] = 0.5  * x13
            vout[2,1,i] = 0.5  * x21
            vout[2,3,i] = 0.5  * x23
            vout[3,1,i] = 0.5  * x31
            vout[3,2,i] = 0.5  * x32
            
        end
    end






