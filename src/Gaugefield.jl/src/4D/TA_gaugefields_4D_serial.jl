struct TA_Gaugefields_4D_serial{NC,NumofBasis} <: TA_Gaugefields_4D{NC}
    a::Array{Float64,5}
    NX::Int64
    NY::Int64
    NZ::Int64
    NT::Int64
    NC::Int64
    NumofBasis::Int64
    generators::Union{Nothing,Generator}

    function TA_Gaugefields_4D_serial(NC,NX,NY,NZ,NT)
        NumofBasis = ifelse(NC == 1,1,NC^2-1)
        if NC <= 3
            generators = nothing
        else
            generators = Generator(NC)
        end
        
        return new{NC,NumofBasis}(zeros(Float64,NumofBasis,NX,NY,NZ,NT),NX,NY,NZ,NT,NC,NumofBasis,generators)
    end
end

function Base.setindex!(x::T,v,i...)  where T<: TA_Gaugefields_4D_serial
    x.a[i...] = v 
end

function Base.getindex(x::T,i...) where T<: TA_Gaugefields_4D_serial
    return x.a[i...]
end


function Base.similar(u::TA_Gaugefields_4D_serial{NC,NumofBasis}) where {NC,NumofBasis}
    return TA_Gaugefields_4D_serial(NC,u.NX,u.NY,u.NZ,u.NT)
    #error("similar! is not implemented in type $(typeof(U)) ")
end

function substitute_U!(Uμ::TA_Gaugefields_4D_serial{NC,NumofBasis},pwork) where {NC,NumofBasis}
    NT = Uμ.NT
    NZ = Uμ.NZ
    NY = Uμ.NY
    NX = Uμ.NX
    #NumofBasis = Uμ.NumofBasis
    icount = 0
    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    for k=1:NumofBasis 
                        icount += 1
                        Uμ[k,ix,iy,iz,it] = pwork[icount]
                    end
                end
            end
        end
    end


    #n1,n2,n3,n4,n5 = size(x.a)
    #println(size(pwork))
    #x.a[:,:,:,:,:] = reshape(pwork,(n1,n2,n3,n4,n5))
end


    
function Base.:*(x::TA_Gaugefields_4D_serial{NC,NumofBasis},y::TA_Gaugefields_4D_serial{NC,NumofBasis}) where {NC,NumofBasis}
    NT = x.NT
    NZ = x.NZ
    NY = x.NY
    NX = x.NX
    #NumofBasis = Uμ.NumofBasis
    s = 0.0
    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    for k=1:NumofBasis 
                        s += x[k,ix,iy,iz,it]*y[k,ix,iy,iz,it]
                    end
                end
            end
        end
    end

    return s
end

function clear_U!(Uμ::TA_Gaugefields_4D_serial{NC,NumofBasis}) where {NC,NumofBasis}
    NT = Uμ.NT
    NZ = Uμ.NZ
    NY = Uμ.NY
    NX = Uμ.NX
    #NumofBasis = Uμ.NumofBasis
    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    for k=1:NumofBasis 
                        Uμ[k,ix,iy,iz,it] = 0
                    end
                end
            end
        end
    end
end

function add_U!(c::TA_Gaugefields_4D_serial{NC,NumofBasis},α::N,a::TA_Gaugefields_4D_serial{NC,NumofBasis}) where {NC, N<:Number,NumofBasis}
    NT = c.NT
    NZ = c.NZ
    NY = c.NY
    NX = c.NX
    #NumofBasis = c.NumofBasis
    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    for k=1:NumofBasis 
                        c[k,ix,iy,iz,it] = c[k,ix,iy,iz,it] + α*a[k,ix,iy,iz,it]
                    end
                end
            end
        end
    end
    #error("add_U! is not implemented in type $(typeof(c)) ")
end



const sr3 = sqrt(3)
const sr3i = 1/sr3
const sr3i2 = 2*sr3i

function Traceless_antihermitian_add!(c::TA_Gaugefields_4D_serial{3,NumofBasis},factor,vin::Gaugefields_4D_wing{3}) where NumofBasis
    #error("Traceless_antihermitian! is not implemented in type $(typeof(vout)) ")
    fac13 = 1/3
    NX = vin.NX
    NY = vin.NY
    NZ = vin.NZ
    NT = vin.NT

    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                @simd for ix=1:NX
                    v11 = vin[1,1,ix,iy,iz,it]
                    v22 = vin[2,2,ix,iy,iz,it] 
                    v33 = vin[3,3,ix,iy,iz,it]

                    tri = fac13*(imag(v11)+imag(v22)+imag(v33))

                    #=
                    vout[1,1,ix,iy,iz,it] = (imag(v11)-tri)*im
                    vout[2,2,ix,iy,iz,it] = (imag(v22)-tri)*im
                    vout[3,3,ix,iy,iz,it] = (imag(v33)-tri)*im
                    =#
                    y11 = (imag(v11)-tri)*im
                    y22 = (imag(v22)-tri)*im
                    y33 = (imag(v33)-tri)*im

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

                    #=
                    vout[1,2,ix,iy,iz,it] = 0.5  * x12
                    vout[1,3,ix,iy,iz,it] = 0.5  * x13
                    vout[2,1,ix,iy,iz,it] = 0.5  * x21
                    vout[2,3,ix,iy,iz,it] = 0.5  * x23
                    vout[3,1,ix,iy,iz,it] = 0.5  * x31
                    vout[3,2,ix,iy,iz,it] = 0.5  * x32
                    =#
                    y12 = 0.5  * x12
                    y13 = 0.5  * x13
                    y21 = 0.5  * x21
                    y23 = 0.5  * x23
                    y31 = 0.5  * x31
                    y32 =  0.5  * x32


                    c[1,ix,iy,iz,it] = ( imag(y12) + imag(y21) )*factor + c[1,ix,iy,iz,it] 
                    c[2,ix,iy,iz,it] = ( real(y12) - real(y21) )*factor + c[2,ix,iy,iz,it] 
                    c[3,ix,iy,iz,it] = ( imag(y11) - imag(y22) )*factor + c[3,ix,iy,iz,it] 
                    c[4,ix,iy,iz,it] = ( imag(y13) + imag(y31) )*factor  + c[4,ix,iy,iz,it] 
                    c[5,ix,iy,iz,it] = ( real(y13) - real(y31) )*factor  + c[5,ix,iy,iz,it] 
                    
                    c[6,ix,iy,iz,it] = ( imag(y23) + imag(y32) )*factor  + c[6,ix,iy,iz,it] 
                    c[7,ix,iy,iz,it] = ( real(y23) - real(y32) )*factor + c[7,ix,iy,iz,it] 
                    c[8,ix,iy,iz,it] = sr3i *
                            ( imag(y11) + imag(y22) -
                                    2*imag(y33) )*factor  + c[8,ix,iy,iz,it] 
                end
            end
        end
    end


end

function Traceless_antihermitian_add!(c::TA_Gaugefields_4D_serial{2,NumofBasis},factor,vin::Gaugefields_4D_wing{2}) where NumofBasis
    #error("Traceless_antihermitian! is not implemented in type $(typeof(vout)) ")
    fac12 = 1/2
    NX = vin.NX
    NY = vin.NY
    NZ = vin.NZ
    NT = vin.NT

    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                @simd for ix=1:NX
                    v11 = vin[1,1,ix,iy,iz,it]
                    v22 = vin[2,2,ix,iy,iz,it]

                    tri = fac12*(imag(v11)+imag(v22))

                    

                    v12 = vin[1,2,ix,iy,iz,it]
                    #v13 = vin[1,3,ix,iy,iz,it]
                    v21 = vin[2,1,ix,iy,iz,it]

                    x12 = v12 - conj(v21)

                    x21 = - conj(x12)

                    y11 = (imag(v11)-tri)*im
                    y12 = 0.5  * x12
                    y21 = 0.5  * x21
                    y22 = (imag(v22)-tri)*im

                    c[1,ix,iy,iz,it] = (imag(y12)+imag(y21))*factor  + c[1,ix,iy,iz,it]
                    c[2,ix,iy,iz,it] = (real(y12)-real(y21))*factor  + c[2,ix,iy,iz,it]
                    c[3,ix,iy,iz,it] = (imag(y11)-imag(y22))*factor  + c[3,ix,iy,iz,it]

                end
            end
        end
    end


end

"""
-----------------------------------------------------c
     !!!!!   vin and vout should be different vectors

     Projectin of the etraceless antiermite part 
     vout = x/2 - Tr(x)/6
     wher   x = vin - Conjg(vin)      
-----------------------------------------------------c
    """
function Traceless_antihermitian!(c::TA_Gaugefields_4D_serial{3,NumofBasis},vin::Gaugefields_4D_wing{3}) where NumofBasis
    #error("Traceless_antihermitian! is not implemented in type $(typeof(vout)) ")
    fac13 = 1/3
    NX = vin.NX
    NY = vin.NY
    NZ = vin.NZ
    NT = vin.NT

    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                @simd for ix=1:NX
                    v11 = vin[1,1,ix,iy,iz,it]
                    v22 = vin[2,2,ix,iy,iz,it]
                    v33 = vin[3,3,ix,iy,iz,it]

                    tri = fac13*(imag(v11)+imag(v22)+imag(v33))

                    #=
                    vout[1,1,ix,iy,iz,it] = (imag(v11)-tri)*im
                    vout[2,2,ix,iy,iz,it] = (imag(v22)-tri)*im
                    vout[3,3,ix,iy,iz,it] = (imag(v33)-tri)*im
                    =#
                    y11 = (imag(v11)-tri)*im
                    y22 = (imag(v22)-tri)*im
                    y33 = (imag(v33)-tri)*im

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

                    #=
                    vout[1,2,ix,iy,iz,it] = 0.5  * x12
                    vout[1,3,ix,iy,iz,it] = 0.5  * x13
                    vout[2,1,ix,iy,iz,it] = 0.5  * x21
                    vout[2,3,ix,iy,iz,it] = 0.5  * x23
                    vout[3,1,ix,iy,iz,it] = 0.5  * x31
                    vout[3,2,ix,iy,iz,it] = 0.5  * x32
                    =#
                    y12 = 0.5  * x12
                    y13 = 0.5  * x13
                    y21 = 0.5  * x21
                    y23 = 0.5  * x23
                    y31 = 0.5  * x31
                    y32 =  0.5  * x32


                    c[1,ix,iy,iz,it] = ( imag(y12) + imag(y21) )
                    c[2,ix,iy,iz,it] = ( real(y12) - real(y21) )
                    c[3,ix,iy,iz,it] = ( imag(y11) - imag(y22) )
                    c[4,ix,iy,iz,it] = ( imag(y13) + imag(y31) )
                    c[5,ix,iy,iz,it] = ( real(y13) - real(y31) )
                    
                    c[6,ix,iy,iz,it] = ( imag(y23) + imag(y32) )
                    c[7,ix,iy,iz,it] = ( real(y23) - real(y32) )
                    c[8,ix,iy,iz,it] = sr3i *
                            ( imag(y11) + imag(y22) -
                                    2*imag(y33) )
                end
            end
        end
    end


end

function Traceless_antihermitian!(c::TA_Gaugefields_4D_serial{NC,NumofBasis},vin::Gaugefields_4D_wing{NC}) where {NC,NumofBasis}
    @assert NC != 3 && NC != 2 
    #NC = vout.NC
    fac1N = 1/NC
    nv = vin.NV

    NX = vin.NX
    NY = vin.NY
    NZ = vin.NZ
    NT = vin.NT
    g = c.generators
    matrix = zeros(ComplexF64,NC,NC)
    a = zeros(ComplexF64,length(g))

    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                @simd for ix=1:NX
                    tri = 0.0
                    @simd for k=1:NC
                        tri += imag(vin[k,k,ix,iy,iz,it])
                    end
                    tri *= fac1N
                    @simd for k=1:NC
                        #vout[k,k,ix,iy,iz,it] = (imag(vin[k,k,ix,iy,iz,it])-tri)*im
                        matrix[k,k] = (imag(vin[k,k,ix,iy,iz,it])-tri)*im
                    end

                    @simd for k2=k1+1:NC
                        vv = 0.5*(vin[k1,k2,ix,iy,iz,it] - conj(vin[k2,k1,ix,iy,iz,it]))
                        #vout[k1,k2,ix,iy,iz,it] = vv
                        #vout[k2,k1,ix,iy,iz,it] = -conj(vv)
                        matrix[k1,k2] = vv
                        matrix[k2,k1] = -conj(vv)
                    end

                    matrix2lie!(a,g,matrix)
                    for k = 1:length(g)
                        c[k,ix,iy,iz,it] = 2*imag(a[k])
                    end

                end
            end
        end
    end

        
    
end

function Traceless_antihermitian_add!(c::TA_Gaugefields_4D_serial{NC,NumofBasis},factor,vin::Gaugefields_4D_wing{NC}) where {NC,NumofBasis}
    @assert NC != 3 && NC != 2 "NC should be NC >4! in this function"
    #NC = vout.NC
    fac1N = 1/NC
    nv = vin.NV

    NX = vin.NX
    NY = vin.NY
    NZ = vin.NZ
    NT = vin.NT
    g = c.generators
    matrix = zeros(ComplexF64,NC,NC)
    a = zeros(ComplexF64,length(g))

    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                @simd for ix=1:NX
                    tri = 0.0
                    @simd for k=1:NC
                        tri += imag(vin[k,k,ix,iy,iz,it])
                    end
                    tri *= fac1N
                    @simd for k=1:NC
                        #vout[k,k,ix,iy,iz,it] = (imag(vin[k,k,ix,iy,iz,it])-tri)*im
                        matrix[k,k] = (imag(vin[k,k,ix,iy,iz,it])-tri)*im
                    end

                    for k1=1:NC
                        @simd for k2=k1+1:NC
                            vv = 0.5*(vin[k1,k2,ix,iy,iz,it] - conj(vin[k2,k1,ix,iy,iz,it]))
                            #vout[k1,k2,ix,iy,iz,it] = vv
                            #vout[k2,k1,ix,iy,iz,it] = -conj(vv)
                            matrix[k1,k2] = vv
                            matrix[k2,k1] = -conj(vv)
                        end
                    end

                    matrix2lie!(a,g,matrix)
                    for k = 1:length(g)
                        c[k,ix,iy,iz,it] = 2*imag(a[k])*factor + c[k,ix,iy,iz,it]
                    end

                end
            end
        end
    end

        
    
end




function exptU!(uout::T,t::N,u::TA_Gaugefields_4D_serial{NC,NumofBasis},temps::Array{T,1}) where {N <: Number, T <: Gaugefields_4D_wing, NC,NumofBasis} #uout = exp(t*u)
    @assert NC != 3 && NC != 2 "This function is for NC != 2,3"
    g = u.generators
    NT = u.NT
    NZ = u.NZ
    NY = u.NY
    NX = u.NX

    u0 = zeros(ComplexF64,NC,NC)
    a = zeros(Float64,length(g))
    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    for k=1:length(a)
                        a[k] = u[k,ix,iy,iz,it]
                    end

                    lie2matrix!(u0,g,a)
                    uout[:,:,ix,iy,iz,it] = exp(t*(im/2)*u0)

                end
            end
        end
    end
    #error("exptU! is not implemented in type $(typeof(u)) ")
end

function exptU!(uout::T,t::N,u::TA_Gaugefields_4D_serial{3,NumofBasis},temps::Array{T,1}) where {N <: Number, T <: Gaugefields_4D_wing,NumofBasis} #uout = exp(t*u)     
    ww = temps[1]
    w = temps[2]
    NT = u.NT
    NZ = u.NZ
    NY = u.NY
    NX = u.NX
    


    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    c1 = t*u[1,ix,iy,iz,it] * 0.5 
                    c2 = t*u[2,ix,iy,iz,it] * 0.5
                    c3 = t*u[3,ix,iy,iz,it] * 0.5
                    c4 = t*u[4,ix,iy,iz,it] * 0.5
                    c5 = t*u[5,ix,iy,iz,it] * 0.5
                    c6 = t*u[6,ix,iy,iz,it] * 0.5
                    c7 = t*u[7,ix,iy,iz,it] * 0.5
                    c8 = t*u[8,ix,iy,iz,it] * 0.5
                    csum = c1+c2+c3+c4+c5+c6+c7+c8
                    if csum == 0
                        w[1,1,ix,iy,iz,it]=   1 
                        w[1,2,ix,iy,iz,it]=   0
                        w[1,3,ix,iy,iz,it]=   0 
                        w[2,1,ix,iy,iz,it]=  0
                        w[2,2,ix,iy,iz,it]=   1 
                        w[2,3,ix,iy,iz,it]=   0
                        w[3,1,ix,iy,iz,it]=   0 
                        w[3,2,ix,iy,iz,it]=   0  
                        w[3,3,ix,iy,iz,it]=   1  
                
                        ww[1,1,ix,iy,iz,it]=   1
                        ww[1,2,ix,iy,iz,it]=   0
                        ww[1,3,ix,iy,iz,it]=   0
                        ww[2,1,ix,iy,iz,it]=   0
                        ww[2,2,ix,iy,iz,it]=   1 
                        ww[2,3,ix,iy,iz,it]=   0 
                        ww[3,1,ix,iy,iz,it]=   0
                        ww[3,2,ix,iy,iz,it]=   0  
                        ww[3,3,ix,iy,iz,it]=   1
                        continue
                    end

                
                    #x[1,1,icum] =  c3+sr3i*c8 +im*(  0.0 )
                    v1 = c3+sr3i*c8
                    v2 = 0.0
                    #x[1,2,icum] =  c1         +im*( -c2   )
                    v3 = c1
                    v4 = -c2
                    #x[1,3,icum] =  c4         +im*(-c5   )
                    v5 = c4
                    v6 = -c5
                
                    #x[2,1,icum] =  c1         +im*(  c2   )
                    v7 = c1
                    v8 = c2

                    #x[2,2,icum] =  -c3+sr3i*c8+im*(  0.0 )
                    v9 =-c3+sr3i*c8
                    v10 = 0.0

                    #x[2,3,icum] =  c6         +im*( -c7   )
                    v11 = c6
                    v12 = -c7
                
                    #x[3,1,icum] =  c4         +im*(  c5   )
                    v13 = c4
                    v14 = c5

                    #x[3,2,icum] =  c6         +im*(  c7   )
                    v15 = c6 
                    v16 = c7
                    #x[3,3,icum] =  -sr3i2*c8  +im*(  0.0 )
                    v17 = -sr3i2*c8
                    v18 = 0.0


                    #c find eigenvalues of v
                    trv3 = (v1 + v9 + v17) / 3.0
                    cofac = v1 *  v9  - v3^2  - v4^2 + 
                            v1 * v17 -  v5^2  - v6^2 + 
                            v9 * v17 - v11^2 - v12^2
                    det = v1 * v9 * v17 -
                            v1 * (v11^2 + v12^2) -
                            v9 * (v5^2 + v6^2) -
                            v17 * (v3^2 + v4^2) +
                            (v5 * (v3 * v11 - v4 * v12) +
                            v6 * (v3 * v12 + v4 * v11)) * 2.0
                    
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

                    w1 =   v5 * (v9 - e1) - v3 * v11 + v4 * v12 
                    w2 = - v6 * (v9 - e1) + v4 * v11 + v3 * v12
                    w3 =   (v1 - e1) * v11 - v3 * v5 - v4 * v6
                    w4 = - (v1 - e1) * v12 - v4 * v5 + v3 * v6
                    w5 = - (v1 - e1) * (v9 - e1) +  v3^2 + v4^2
                    w6 = 0.0
                    #println("1c $w1 $w2 $w3 $w4 $w5 $w6 ",)
                    #coeffv = sqrt(w1^2 + w2^2 + w3^2 + w4^2 + w5^2)
                    
                    #coeff = ifelse(coeffv == zero(coeffv),0,coeffv)
                    coeff = 1.0 / sqrt(w1^2 + w2^2 + w3^2 + w4^2 + w5^2)
                    #println("1 ",coeff)

                    w1 = w1 * coeff
                    w2 = w2 * coeff
                    w3 = w3 * coeff
                    w4 = w4 * coeff
                    w5 = w5 * coeff
            
                    w7 =   v5 * (v9 - e2) - v3 * v11 + v4 * v12
                    w8 = - v6 * (v9 - e2) + v4 * v11 + v3 * v12
                    w9 = (v1 - e2) * v11 - v3 * v5 - v4 * v6
                    w10 = - (v1 - e2) * v12 - v4 * v5 + v3 * v6
                    w11 = - (v1 - e2) * (v9 - e2) +  v3^2 + v4^2
                    w12 = 0.0
            
                    coeff = 1.0 / sqrt(w7^2  + w8^2 + w9^2+ w10^2 + w11^2)
                    
                    w7  = w7  * coeff
                    w8  = w8  * coeff
                    w9  = w9  * coeff
                    w10 = w10 * coeff
                    w11 = w11 * coeff
            
                    w13 =   v5 * (v9 - e3) - v3 * v11 + v4 * v12
                    w14 = - v6 * (v9 - e3) + v4 * v11 + v3 * v12
                    w15 =   (v1 - e3) * v11 - v3 * v5 - v4 * v6
                    w16 = - (v1 - e3) * v12 - v4 * v5 + v3 * v6
                    w17 = - (v1 - e3) * (v9 - e3) +  v3^2 + v4^2
                    w18 = 0.0

                    coeff = 1.0 / sqrt(w13^2 + w14^2 + w15^2+ w16^2 + w17^2)
                    w13 = w13 * coeff
                    w14 = w14 * coeff
                    w15 = w15 * coeff
                    w16 = w16 * coeff
                    w17 = w17 * coeff
            
            # construct the projection v
                    c1 = cos(e1)
                    s1 = sin(e1)
                    ww1  = w1  * c1 - w2 * s1
                    ww2  = w2  * c1 + w1 * s1
                    ww3  = w3  * c1 - w4 * s1
                    ww4  = w4  * c1 + w3 * s1
                    ww5  = w5  * c1 - w6 * s1
                    ww6  = w6  * c1 + w5 * s1

                    c2 = cos(e2)
                    s2 = sin(e2)
                    ww7  = w7  * c2 - w8 * s2
                    ww8  = w8  * c2 + w7 * s2
                    ww9  = w9  * c2 - w10 * s2
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

                    w[1,1,ix,iy,iz,it]=   w1+im*w2 
                    w[1,2,ix,iy,iz,it]=   w3+im*w4 
                    w[1,3,ix,iy,iz,it]=   w5+im*w6 
                    w[2,1,ix,iy,iz,it]=   w7+im*w8 
                    w[2,2,ix,iy,iz,it]=   w9+im*w10  
                    w[2,3,ix,iy,iz,it]=   w11+im*w12  
                    w[3,1,ix,iy,iz,it]=   w13+im*w14  
                    w[3,2,ix,iy,iz,it]=   w15+im*w16  
                    w[3,3,ix,iy,iz,it]=   w17+im*w18  
            
                    ww[1,1,ix,iy,iz,it]=   ww1+im*ww2 
                    ww[1,2,ix,iy,iz,it]=   ww3+im*ww4 
                    ww[1,3,ix,iy,iz,it]=   ww5+im*ww6 
                    ww[2,1,ix,iy,iz,it]=   ww7+im*ww8 
                    ww[2,2,ix,iy,iz,it]=   ww9+im*ww10  
                    ww[2,3,ix,iy,iz,it]=   ww11+im*ww12  
                    ww[3,1,ix,iy,iz,it]=   ww13+im*ww14  
                    ww[3,2,ix,iy,iz,it]=   ww15+im*ww16  
                    ww[3,3,ix,iy,iz,it]=   ww17+im*ww18  

                end
            end
        end
    end

    mul!(uout,w',ww)


end
const tinyvalue =1e-100


function exptU!(uout::T,t::N,u::TA_Gaugefields_4D_serial{2,NumofBasis},temps::Array{T,1}) where {N <: Number, T <: Gaugefields_4D_wing,NumofBasis} #uout = exp(t*u)     
    NT = u.NT
    NZ = u.NZ
    NY = u.NY
    NX = u.NX


    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    #icum = (((it-1)*NX+iz-1)*NY+iy-1)*NX+ix  
                    u1 = t*u[1,ix,iy,iz,it]/2
                    u2 = t*u[2,ix,iy,iz,it]/2
                    u3 = t*u[3,ix,iy,iz,it]/2
                    R = sqrt(u1^2+u2^2+u3^2) +  tinyvalue
                    sR = sin(R)/R
                    #sR = ifelse(R == 0,1,sR)
                    a0 = cos(R)
                    a1 = u1*sR
                    a2 = u2*sR
                    a3 = u3*sR

                    uout[1,1,ix,iy,iz,it] = cos(R) + im*a3
                    uout[1,2,ix,iy,iz,it] = im*a1 + a2
                    uout[2,1,ix,iy,iz,it] = im*a1 - a2
                    uout[2,2,ix,iy,iz,it]= cos(R) - im*a3

                end
            end
        end
    end




end
