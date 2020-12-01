
struct SU2GaugeFields <: GaugeFields
    g::Array{ComplexF64,6}
    NX::Int64
    NY::Int64
    NZ::Int64
    NT::Int64
    NC::Int64
    NDW::Int64
    NV::Int64

    function SU2GaugeFields(NC,NDW,NX,NY,NZ,NT)
        @assert NC == 2
        NV = NX*NY*NZ*NT
        g = zeros(ComplexF64,NC,NC,NX+2NDW,NY+2NDW,NZ+2NDW,NT+2NDW)
        return new(g,NX,NY,NZ,NT,NC,NDW,NV)
    end
end

struct SU2GaugeFields_1d <: GaugeFields_1d
    g::Array{ComplexF64,3}
    NC::Int64
    NV::Int64

    function SU2GaugeFields_1d(NC,NX,NY,NZ,NT)
        @assert NC == 2
        NV = NX*NY*NZ*NT
        g = zeros(ComplexF64,NC,NC,NV)
        return new(g,NC,NV)
    end

    function SU2GaugeFields_1d(NC,NV) 
        g = zeros(ComplexF64,NC,NC,NV)
        return new(g,NC,NV)
    end
end

function Base.similar(x::SU2GaugeFields_1d)
    return SU2GaugeFields_1d(x.NC,x.NV) 

end

function Base.setindex!(x::SU2GaugeFields,v,i1,i2,i3,i4,i5,i6) 
    x.g[i1,i2,i3 + x.NDW,i4 + x.NDW,i5 + x.NDW,i6 + x.NDW] = v
end

function Base.getindex(x::SU2GaugeFields,i1,i2,i3,i4,i5,i6)
    return x.g[i1,i2,i3 .+ x.NDW,i4 .+ x.NDW,i5 .+ x.NDW,i6 .+ x.NDW]
end

function Base.setindex!(x::SU2GaugeFields_1d,v,i1,i2,i3) 
    x.g[i1,i2,i3] = v
end

function Base.getindex(x::SU2GaugeFields_1d,i1,i2,i3)
    return x.g[i1,i2,i3]
end

struct Adjoint_SU2GaugeFields_1d <: Adjoint_GaugeFields_1d
    parent::SU2GaugeFields_1d
end

struct Adjoint_SU2GaugeFields <: Adjoint_GaugeFields
    parent::SU2GaugeFields
end

function Base.adjoint(x::SU2GaugeFields_1d)
    Adjoint_SU2GaugeFields_1d(x)
end

function Base.adjoint(x::SU2GaugeFields)
    Adjoint_SU2GaugeFields(x)
end

function IdentitySU2Gauges(NX,NY,NZ,NT,NDW)
    NC = 2
    U = SU2GaugeFields(NC,NDW,NX,NY,NZ,NT)
    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    for ic=1:2
                        U[ic,ic,ix,iy,iz,it] = 1 
                    end
                end
            end
        end
    end
    return U
end

function RandomSU2Gauges(NX,NY,NZ,NT,NDW)
    NC = 2
    U = SU2GaugeFields(NC,NDW,NX,NY,NZ,NT)
    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    U[1,1,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
                    U[1,2,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
                    U[2,1,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
                    U[2,2,ix,iy,iz,it] = rand()-0.5 + im*(rand()-0.5)
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
c     input  x : arbitrary 2*2 complex matrix                          c
c     output x : SU(2) matrix 
"""
function normalize!(u::SU2GaugeFields)
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
                    for ic=1:2
                        w1 += u[2,ic,ix,iy,iz,it]*conj(u[1,ic,ix,iy,iz,it])
                        w2 += u[1,ic,ix,iy,iz,it]*conj(u[1,ic,ix,iy,iz,it])
                    end
                    zerock2 = w2
                    if zerock2 == 0 
                        println("w2 is zero  !!  (in normlz)")
                        println("u[1,1),u[1,2) : ",u[1,1,ix,iy,iz,it], "\t",u[1,2,ix,iy,iz,it])
                    end

                    w1 = -w1/w2

                    x4 = (u[2,1,ix,iy,iz,it]) + w1*u[1,1,ix,iy,iz,it]
                    x5 = (u[2,2,ix,iy,iz,it]) + w1*u[1,2,ix,iy,iz,it]
                    #x6 = (u[2,3,ix,iy,iz,it]) + w1*u[1,3,ix,iy,iz,it]

                    w3 = x4*conj(x4) + x5*conj(x5) #+ x6*conj(x6)

                    zerock3 = w3
                    if zerock3 == 0
                        println("w3 is zero  !!  (in normlz)")
                        println("x4, x5: $x4, $x5")
                        exit()
                    end

                    u[2,1,ix,iy,iz,it] = x4
                    u[2,2,ix,iy,iz,it] = x5
                    #u[2,3,ix,iy,iz,it] = x6

                    w3 = 1/sqrt(w3)
                    w2 = 1/sqrt(w2)

                    u[1,1,ix,iy,iz,it] = u[1,1,ix,iy,iz,it]*w2
                    u[1,2,ix,iy,iz,it] = u[1,2,ix,iy,iz,it]*w2
                    #u[1,3,ix,iy,iz,it] = u[1,3,ix,iy,iz,it]*w2
                    u[2,1,ix,iy,iz,it] = u[2,1,ix,iy,iz,it]*w3
                    u[2,2,ix,iy,iz,it] = u[2,2,ix,iy,iz,it]*w3
                    #u[2,3,ix,iy,iz,it] = u[2,3,ix,iy,iz,it]*w3

                    if zerock2*zerock3 == 0 
                        println("!! devided by zero !! (in normalize)")
                        println("w2 or w3 in normlz is zero !!")
                        println("w2, w3 : $w2, $w3   ")
                        exit()
                    end
                    #println(u[:,:,ix,iy,iz,it]'*u[:,:,ix,iy,iz,it])
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
function projlink!(vout::SU2GaugeFields,vin::SU2GaugeFields)
    
    fac12 = 1/2
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

                    tri = fac12*(imag(v11)+imag(v22))

                    vout[1,1,ix,iy,iz,it] = (imag(v11)-tri)*im
                    vout[2,2,ix,iy,iz,it] = (imag(v22)-tri)*im
                end
            end
        end
    end

    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX

                    v12 = vin[1,2,ix,iy,iz,it]
                    #v13 = vin[1,3,ix,iy,iz,it]
                    v21 = vin[2,1,ix,iy,iz,it]

                    x12 = v12 - conj(v21)

                    x21 = - conj(x12)

                    vout[1,2,ix,iy,iz,it] = 0.5  * x12

                    vout[2,1,ix,iy,iz,it] = 0.5  * x21


                end
            end
        end
    end

end



function projlink!(vout::SU2GaugeFields_1d,vin::SU2GaugeFields_1d)
    fac12 = 1/2
    nv = vin.NV

    for i=1:nv
        v11 = vin[1,1,i]
        v22 = vin[2,2,i]

        tri = fac12*(imag(v11)+imag(v22))

        vout[1,1,i] = (imag(v11)-tri)*im
        vout[2,2,i] = (imag(v22)-tri)*im

    end

    for i=1:nv
        v12 = vin[1,2,i]
        v21 = vin[2,1,i]


        x12 = v12 - conj(v21)

    
        x21 = - conj(x12)

        vout[1,2,i] = 0.5  * x12
        vout[2,1,i] = 0.5  * x21

        
    end

end

function lambdamul(b::SU2GaugeFields_1d,a::SU2GaugeFields_1d,k)
    #=
    c----------------------------------------------------------------------c
    c     b = (lambda_k/2)*a
    C             lambda_k : GellMann matrices. k=1, 8 
    c----------------------------------------------------------------------c
    =#
    NV = a.NV


    if k==1
        for i=1:NV
            b[1,1,i] = -0.5*im* a[2,1,i]*im
            b[1,2,i] = -0.5*im * a[2,2,i]*im

            b[2,1,i] = -0.5*im * a[1,1,i]*im
            b[2,2,i] = -0.5*im * a[1,2,i]*im


        end
    elseif k==2
        for i=1:NV
            b[1,1,i] = -0.5 * a[2,1,i] *im
            b[1,2,i] = -0.5 * a[2,2,i]*im

            b[2,1,i] =  0.5 * a[1,1,i]*im*
            b[2,2,i] =  0.5 * a[1,2,i]*im

        end
    elseif k==3
        for i=1:NV
            b[1,1,i] =  -0.5*im * a[1,1,i] *im
            b[1,2,i] =  -0.5*im * a[1,2,i]*im

            b[2,1,i] = 0.5*im * a[2,1,i]*im
            b[2,2,i] = 0.5*im * a[2,2,i]*im

        end
    else
        error("k should be k <= 3 but k = $k")
    end

    return
end







