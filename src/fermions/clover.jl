module Clover
using LinearAlgebra
import ..Actions: FermiActionParam_WilsonClover
import ..LieAlgebrafields: clear!
import ..Gaugefields:
    GaugeFields_1d,
    GaugeFields,
    set_wing!,
    substitute!,
    gauge_shift!,
    lambdamul,
    SU3GaugeFields_1d,
    SU2GaugeFields_1d,
    SUNGaugeFields_1d,
    Loops,
    evaluate_loops!,
    evaluate_loops,
    SU
import ..Wilsonloops: Wilson_loop, Wilson_loop_set

function Make_CloverFμν(
    fparam::FermiActionParam_WilsonClover,
    U,
    temps::Array{T,1},
) where {T<:GaugeFields_1d}
    NV = temps[1].NV
    NC = temps[1].NC
    CloverFμν = zeros(ComplexF64, NC, NC, NV, 6)
    Make_CloverFμν!(CloverFμν, fparam, U, temps)
    return CloverFμν
end

function Make_CloverFμν!(
    CloverFμν,
    fparam::FermiActionParam_WilsonClover,
    U,
    temps::Array{T,1},
) where {T<:GaugeFields_1d}
    work1 = temps[1]
    work2 = temps[2]
    work3 = temps[3]
    work4 = temps[4]

    #CloverFμν .= 0
    set_wing!(U)

    NV = work1.NV
    NC = work1.NC

    # ... Calculation of 4 leaves under the counter clock order.
    μν = 0
    for μ = 1:3
        for ν = μ+1:4
            μν += 1
            if μν > 6
                error("μν > 6 ?")
            end

            loopset = Loops(U, fparam._cloverloops[μ, ν], [work2, work3, work4])
            evaluate_loops!(work1, loopset, U)
            setFμν!(CloverFμν, μν, work1)

        end
    end

    for μν = 1:6
        cimaglink!(CloverFμν, μν, NV, NC)
    end

    coe = im * 0.125 * fparam.hop * fparam.Clover_coefficient
    CloverFμν .*= coe

    #println("clover ",sum(abs.(fparam.CloverFμν)))


    return

end


function Make_CloverFμν!(
    fparam::FermiActionParam_WilsonClover,
    U,
    temps::Array{T,1},
) where {T<:GaugeFields_1d}
    Make_CloverFμν!(fparam.CloverFμν, fparam, U, temps)
    return
end

"""
 y = x - x_aj
"""
function cimaglink!(Fμν, μν, NV, NC)
    if NC == 3
        for i = 1:NV
            z11 = Fμν[1, 1, i, μν] - conj(Fμν[1, 1, i, μν])
            z12 = Fμν[1, 2, i, μν] - conj(Fμν[2, 1, i, μν])
            z13 = Fμν[1, 3, i, μν] - conj(Fμν[3, 1, i, μν])

            z22 = Fμν[2, 2, i, μν] - conj(Fμν[2, 2, i, μν])
            z23 = Fμν[2, 3, i, μν] - conj(Fμν[3, 2, i, μν])

            z33 = Fμν[3, 3, i, μν] - conj(Fμν[3, 3, i, μν])

            Fμν[1, 1, i, μν] = z11
            Fμν[1, 2, i, μν] = z12
            Fμν[1, 3, i, μν] = z13

            Fμν[2, 1, i, μν] = -conj(z12)
            Fμν[2, 2, i, μν] = z22
            Fμν[2, 3, i, μν] = z23

            Fμν[3, 1, i, μν] = -conj(z13)
            Fμν[3, 2, i, μν] = -conj(z23)
            Fμν[3, 3, i, μν] = z33
        end
    elseif NC == 2
        for i = 1:NV
            z11 = Fμν[1, 1, i, μν] - conj(Fμν[1, 1, i, μν])
            z12 = Fμν[1, 2, i, μν] - conj(Fμν[2, 1, i, μν])

            z22 = Fμν[2, 2, i, μν] - conj(Fμν[2, 2, i, μν])


            Fμν[1, 1, i, μν] = z11
            Fμν[1, 2, i, μν] = z12


            Fμν[2, 1, i, μν] = -conj(z12)
            Fμν[2, 2, i, μν] = z22

        end
    else
        for i = 1:NV
            for k2 = 1:NC
                for k1 = k2:NC
                    z = Fμν[k1, k2, i, μν] - conj(Fμν[k2, k1, i, μν])
                    Fμν[k1, k2, i, μν] = z
                    if k1 != k2
                        Fμν[k2, k1, i, μν] = -conj(z)
                    end
                end
            end
        end
        #error("error!")
    end
    return
end

function setFμν!(Fμν, μν, v::GaugeFields_1d{SU{NC}}) where {NC}
    NV = v.NV
    for i = 1:NV
        for k2 = 1:NC
            for k1 = 1:NC
                Fμν[k1, k2, i, μν] = v[k1, k2, i]
            end
        end
    end
    return
end


function addFμν!(Fμν, μν, v::SU3GaugeFields_1d)
    NV = v.NV

    for i = 1:NV
        Fμν[1, 1, i, μν] += v[1, 1, i]
        Fμν[1, 2, i, μν] += v[1, 2, i]
        Fμν[1, 3, i, μν] += v[1, 3, i]
        Fμν[2, 1, i, μν] += v[2, 1, i]
        Fμν[2, 2, i, μν] += v[2, 2, i]
        Fμν[2, 3, i, μν] += v[2, 3, i]
        Fμν[3, 1, i, μν] += v[3, 1, i]
        Fμν[3, 2, i, μν] += v[3, 2, i]
        Fμν[3, 3, i, μν] += v[3, 3, i]
    end
    return
end

function addFμν!(Fμν, μν, v::SU2GaugeFields_1d)
    NV = v.NV

    for i = 1:NV
        Fμν[1, 1, i, μν] += v[1, 1, i]
        Fμν[1, 2, i, μν] += v[1, 2, i]
        Fμν[2, 1, i, μν] += v[2, 1, i]
        Fμν[2, 2, i, μν] += v[2, 2, i]

    end
    return
end

function addFμν!(Fμν, μν, v::SUNGaugeFields_1d)
    NV = v.NV
    NC = v.NC

    for i = 1:NV
        for k2 = 1:NC
            for k1 = 1:NC
                Fμν[k1, k2, i, μν] += v[k1, k2, i]
            end
        end
    end
    return
end

jf(I, L) = (I + 1) % L
jb(I, L) = (I - 1 + L) % L
ic(ix, iy, iz, it, N1, N2, N3) = 1 + ix + N1 * (iy + N2 * (iz + N3 * it))


function fprep3!(NX, NY, NZ, NT, inn)
    for it = 0:NT-1
        for iz = 0:NZ-1
            for iy = 0:NY-1
                for ix = 0:NX-1
                    i = ic(ix, iy, iz, it, NX, NY, NZ)
                    inn[i, 1, 1] = ic(jf(ix, NX), iy, iz, it, NX, NY, NZ)
                    inn[i, 2, 1] = ic(ix, jf(iy, NY), iz, it, NX, NY, NZ)
                    inn[i, 3, 1] = ic(ix, iy, jf(iz, NZ), it, NX, NY, NZ)
                    inn[i, 4, 1] = ic(ix, iy, iz, jf(it, NT), NX, NY, NZ)
                    inn[i, 1, 2] = ic(jb(ix, NX), iy, iz, it, NX, NY, NZ)
                    inn[i, 2, 2] = ic(ix, jb(iy, NY), iz, it, NX, NY, NZ)
                    inn[i, 3, 2] = ic(ix, iy, jb(iz, NZ), it, NX, NY, NZ)
                    inn[i, 4, 2] = ic(ix, iy, iz, jb(it, NT), NX, NY, NZ)
                end
            end
        end
    end
    return
end


"""
Calculate   dS_clover/dA_mu(x)
"""
function dSclover!(
    z,
    μ,
    X,
    Y,
    U,
    fparam::FermiActionParam_WilsonClover,
    temps::Array{T,1},
) where {T<:GaugeFields_1d}
    NX = X.NX
    NY = X.NY
    NZ = X.NZ
    NT = X.NT
    NC = X.NC
    NV = NX * NY * NZ * NT
    numbasis = NC^2 - 1
    #numbasis = ifelse(NC==3,8,3)


    work1 = temps[1]
    work2 = temps[2]
    work3 = temps[3]
    work4 = temps[4]

    gtmp1 = temps[5]
    gtmp2 = temps[6]
    gtmp3 = temps[7]
    gtmp4 = temps[8]

    dF1 = temps[8+1:8+numbasis]
    #dF1 = temps[9:16]
    dF2 = temps[8+numbasis+1:8+2numbasis]
    #dF2 = temps[17:24]




    if fparam.internal_flags[1] == false
        fprep3!(NX, NY, NZ, NT, fparam.inn_table)
        fparam.internal_flags[1] = true
    end

    clear!(z)
    veta = fparam._ftmp_vectors[1]
    vxi = fparam._ftmp_vectors[2]
    ftmp1 = fparam._ftmp_vectors[3]
    ftmp2 = fparam._ftmp_vectors[4]
    v1 = fparam._ftmp_vectors[5]
    v2 = fparam._ftmp_vectors[6]
    is1 = fparam._is1
    is2 = fparam._is2
    inn = fparam.inn_table



    for ialpha = 1:4
        for ic = 1:NC
            is = 0
            for it = 1:NT
                for iz = 1:NZ
                    for iy = 1:NY
                        for ix = 1:NX
                            is += 1
                            veta[ic, is, ialpha] = Y[ic, ix, iy, iz, it, ialpha]
                            vxi[ic, is, ialpha] = X[ic, ix, iy, iz, it, ialpha]
                        end
                    end
                end
            end
        end
    end



    for ν = 1:4
        if ν == μ
            continue
        end

        #=
         .....................  
            Case 1 and 3        
         .....................  
        =#
        iflag = 1
        cal_dFμν!(
            dF1,
            dF2,
            U,
            fparam,
            work1,
            work2,
            work3,
            work4,
            gtmp1,
            gtmp2,
            gtmp3,
            gtmp4,
            μ,
            ν,
            iflag,
        )

        # ... Case 1
        VxSigxV!(veta, vxi, dF1, z, ftmp1, ftmp2, μ, ν)
        #println("z1 = ",z*z)

        # ... Case 3
        for is = 1:NV
            is1[is] = inn[is, μ, 1]
        end

        for is = 1:NV
            is2[is] = inn[is1[is], ν, 1]
        end

        for ic = 1:NC
            for ialpha = 1:4
                for is = 1:NV
                    v1[ic, is, ialpha] = veta[ic, is2[is], ialpha]
                    v2[ic, is, ialpha] = vxi[ic, is2[is], ialpha]
                end
            end
        end

        VxSigxV!(v1, v2, dF2, z, ftmp1, ftmp2, μ, ν)
        #println("z2 = ",z*z)

        #=
         .....................  
            Case 2 and 4        
         .....................  
        =#

        iflag = 2

        cal_dFμν!(
            dF1,
            dF2,
            U,
            fparam,
            work1,
            work2,
            work3,
            work4,
            gtmp1,
            gtmp2,
            gtmp3,
            gtmp4,
            μ,
            ν,
            iflag,
        )

        # ... Case 2
        for is = 1:NV
            is1[is] = inn[is, μ, 1]
        end

        for ic = 1:NC
            for ialpha = 1:4
                for is = 1:NV
                    v1[ic, is, ialpha] = veta[ic, is1[is], ialpha]
                    v2[ic, is, ialpha] = vxi[ic, is1[is], ialpha]
                end
            end
        end

        VxSigxV!(v1, v2, dF1, z, ftmp1, ftmp2, μ, ν)
        #println("z3 = ",z*z)

        # ... Case 4
        for is = 1:NV
            is1[is] = inn[is, ν, 1]
        end

        for ic = 1:NC
            for ialpha = 1:4
                for is = 1:NV
                    v1[ic, is, ialpha] = veta[ic, is1[is], ialpha]
                    v2[ic, is, ialpha] = vxi[ic, is1[is], ialpha]
                end
            end
        end

        VxSigxV!(v1, v2, dF2, z, ftmp1, ftmp2, μ, ν)
        #println("z4 = ",z*z)

        #=
         .....................  
            Case 4' and 2'      
         .....................  
        =#

        iflag = 3

        cal_dFμν!(
            dF1,
            dF2,
            U,
            fparam,
            work1,
            work2,
            work3,
            work4,
            gtmp1,
            gtmp2,
            gtmp3,
            gtmp4,
            μ,
            ν,
            iflag,
        )

        # ... Case 4'
        VxSigxV!(veta, vxi, dF1, z, ftmp1, ftmp2, μ, ν)
        #println("z5 = ",z*z)

        # ... Case 2' 
        for is = 1:NV
            is1[is] = inn[is, μ, 1]
        end

        for is = 1:NV
            is2[is] = inn[is1[is], ν, 2]
        end

        for ic = 1:NC
            for ialpha = 1:4
                for is = 1:NV
                    v1[ic, is, ialpha] = veta[ic, is2[is], ialpha]
                    v2[ic, is, ialpha] = vxi[ic, is2[is], ialpha]
                end
            end
        end

        VxSigxV!(v1, v2, dF2, z, ftmp1, ftmp2, μ, ν)
        #println("z6 = ",z*z)

        #=
         .....................  
            Case 3' and 1'     
         .....................  
        =#
        iflag = 4

        cal_dFμν!(
            dF1,
            dF2,
            U,
            fparam,
            work1,
            work2,
            work3,
            work4,
            gtmp1,
            gtmp2,
            gtmp3,
            gtmp4,
            μ,
            ν,
            iflag,
        )




        # ... Case 3'
        for is = 1:NV
            is1[is] = inn[is, μ, 1]
        end

        for ic = 1:NC
            for ialpha = 1:4
                for is = 1:NV
                    v1[ic, is, ialpha] = veta[ic, is1[is], ialpha]
                    v2[ic, is, ialpha] = vxi[ic, is1[is], ialpha]
                end
            end
        end


        VxSigxV!(v1, v2, dF1, z, ftmp1, ftmp2, μ, ν)


        # ... Case 1' 
        for is = 1:NV
            is1[is] = inn[is, ν, 2]
        end

        for ic = 1:NC
            for ialpha = 1:4
                for is = 1:NV
                    v1[ic, is, ialpha] = veta[ic, is1[is], ialpha]
                    v2[ic, is, ialpha] = vxi[ic, is1[is], ialpha]
                end
            end
        end

        VxSigxV!(v1, v2, dF2, z, ftmp1, ftmp2, μ, ν)


    end


    for ia = 1:numbasis
        is = 0
        for it = 1:NT
            for iz = 1:NZ
                for iy = 1:NY
                    for ix = 1:NX
                        is += 1
                        z[ia, ix, iy, iz, it] = -2 * z[ia, ix, iy, iz, it]
                    end
                end
            end
        end

    end



end


function VxSigxV!(v1, v2, u, z, tmp1, tmp2, μ, ν)
    NC, NV, _ = size(v1)
    NX = z.NX
    NY = z.NY
    NZ = z.NZ
    NT = z.NT
    #numbasis = ifelse(NC==3,8,3)
    numbasis = NC^2 - 1


    if μ == ν
        error("""μ should not be equal to ν (VsSigV)
        μ,ν : $μ $ν
        """)
    end

    if μ < ν
        facμν = 1
        μ0 = μ
        ν0 = ν
    else
        facμν = -1
        μ0 = ν
        ν0 = μ
    end

    # ... σ_μν x v2 (in Dirac space)  .... 
    if μ0 == 1
        if ν0 == 2
            for ic = 1:NC
                for is = 1:NV
                    tmp1[ic, is, 1] = -v2[ic, is, 1]
                    tmp1[ic, is, 2] = +v2[ic, is, 2]
                    tmp1[ic, is, 3] = -v2[ic, is, 3]
                    tmp1[ic, is, 4] = +v2[ic, is, 4]
                end
            end

        elseif ν0 == 3
            for ic = 1:NC
                for is = 1:NV
                    tmp1[ic, is, 1] = -im * v2[ic, is, 2]
                    tmp1[ic, is, 2] = +im * v2[ic, is, 1]
                    tmp1[ic, is, 3] = -im * v2[ic, is, 4]
                    tmp1[ic, is, 4] = +im * v2[ic, is, 3]
                end
            end
        elseif ν0 == 4
            for ic = 1:NC
                for is = 1:NV
                    tmp1[ic, is, 1] = -v2[ic, is, 2]
                    tmp1[ic, is, 2] = -v2[ic, is, 1]
                    tmp1[ic, is, 3] = +v2[ic, is, 4]
                    tmp1[ic, is, 4] = +v2[ic, is, 3]
                end
            end
        else
            error("""something is wrong in VsSigV
            μ,ν : $μ $ν
            """)
        end
    elseif μ0 == 2
        if ν0 == 3
            for ic = 1:NC
                for is = 1:NV
                    tmp1[ic, is, 1] = -v2[ic, is, 2]
                    tmp1[ic, is, 2] = -v2[ic, is, 1]
                    tmp1[ic, is, 3] = -v2[ic, is, 4]
                    tmp1[ic, is, 4] = -v2[ic, is, 3]
                end
            end
        elseif ν0 == 4
            for ic = 1:NC
                for is = 1:NV
                    tmp1[ic, is, 1] = +im * v2[ic, is, 2]
                    tmp1[ic, is, 2] = -im * v2[ic, is, 1]
                    tmp1[ic, is, 3] = -im * v2[ic, is, 4]
                    tmp1[ic, is, 4] = +im * v2[ic, is, 3]
                end
            end
        else
            error("""something is wrong in VsSigV
            μ,ν : $μ $ν
            """)
        end
    elseif μ0 == 3
        if ν0 == 4
            for ic = 1:NC
                for is = 1:NV
                    tmp1[ic, is, 1] = -v2[ic, is, 1]
                    tmp1[ic, is, 2] = +v2[ic, is, 2]
                    tmp1[ic, is, 3] = +v2[ic, is, 3]
                    tmp1[ic, is, 4] = -v2[ic, is, 4]
                end
            end

        else
            error("""something is wrong in VsSigV
            μ,ν : $μ $ν
            """)
        end
    else
        error("""something is wrong in VsSigV
            μ,ν : $μ $ν
            """)
    end

    for ic = 1:NC
        for is = 1:NV
            tmp2[ic, is, 1] = facμν * tmp1[ic, is, 1]
            tmp2[ic, is, 2] = facμν * tmp1[ic, is, 2]
            tmp2[ic, is, 3] = facμν * tmp1[ic, is, 3]
            tmp2[ic, is, 4] = facμν * tmp1[ic, is, 4]
        end
    end

    if NC == 3

        # ... tmp1 = u x tmp2 (in color space) 
        for ia = 1:8
            for ialpha = 1:4
                for is = 1:NV
                    s1 =
                        +u[ia][1, 1, is] * tmp2[1, is, ialpha] +
                        u[ia][1, 2, is] * tmp2[2, is, ialpha] +
                        +u[ia][1, 3, is] * tmp2[3, is, ialpha]
                    s2 =
                        +u[ia][2, 1, is] * tmp2[1, is, ialpha] +
                        u[ia][2, 2, is] * tmp2[2, is, ialpha] +
                        u[ia][2, 3, is] * tmp2[3, is, ialpha]
                    s3 =
                        +u[ia][3, 1, is] * tmp2[1, is, ialpha] +
                        +u[ia][3, 2, is] * tmp2[2, is, ialpha] +
                        +u[ia][3, 3, is] * tmp2[3, is, ialpha]


                    tmp1[1, is, ialpha] = s1
                    tmp1[2, is, ialpha] = s2
                    tmp1[3, is, ialpha] = s3
                end
            end



            # ... v1^{\dagger} * tmp1
            is = 0
            for it = 1:NT

                for iz = 1:NZ
                    for iy = 1:NY
                        for ix = 1:NX
                            is += 1
                            s1 =
                                conj(v1[1, is, 1]) * tmp1[1, is, 1] +
                                conj(v1[1, is, 2]) * tmp1[1, is, 2] +
                                conj(v1[1, is, 3]) * tmp1[1, is, 3] +
                                conj(v1[1, is, 4]) * tmp1[1, is, 4]

                            s2 =
                                conj(v1[2, is, 1]) * tmp1[2, is, 1] +
                                conj(v1[2, is, 2]) * tmp1[2, is, 2] +
                                conj(v1[2, is, 3]) * tmp1[2, is, 3] +
                                conj(v1[2, is, 4]) * tmp1[2, is, 4]

                            s3 =
                                conj(v1[3, is, 1]) * tmp1[3, is, 1] +
                                conj(v1[3, is, 2]) * tmp1[3, is, 2] +
                                conj(v1[3, is, 3]) * tmp1[3, is, 3] +
                                conj(v1[3, is, 4]) * tmp1[3, is, 4]

                            z[ia, ix, iy, iz, it] += real(s1 + s2 + s3)
                            #println(z[ia,ix,iy,iz,it])
                        end
                    end
                end
            end


        end
    elseif NC == 2

        # ... tmp1 = u x tmp2 (in color space) 
        for ia = 1:3
            for ialpha = 1:4
                for is = 1:NV
                    s1 =
                        +u[ia][1, 1, is] * tmp2[1, is, ialpha] +
                        u[ia][1, 2, is] * tmp2[2, is, ialpha]
                    s2 =
                        +u[ia][2, 1, is] * tmp2[1, is, ialpha] +
                        u[ia][2, 2, is] * tmp2[2, is, ialpha]

                    tmp1[1, is, ialpha] = s1
                    tmp1[2, is, ialpha] = s2

                end
            end

            # ... v1^{\dagger} * tmp1
            is = 0
            for it = 1:NT
                for iz = 1:NZ
                    for iy = 1:NY
                        for ix = 1:NX
                            is += 1
                            s1 =
                                conj(v1[1, is, 1]) * tmp1[1, is, 1] +
                                conj(v1[1, is, 2]) * tmp1[1, is, 2] +
                                conj(v1[1, is, 3]) * tmp1[1, is, 3] +
                                conj(v1[1, is, 4]) * tmp1[1, is, 4]

                            s2 =
                                conj(v1[2, is, 1]) * tmp1[2, is, 1] +
                                conj(v1[2, is, 2]) * tmp1[2, is, 2] +
                                conj(v1[2, is, 3]) * tmp1[2, is, 3] +
                                conj(v1[2, is, 4]) * tmp1[2, is, 4]

                            z[ia, ix, iy, iz, it] += real(s1 + s2)
                        end
                    end
                end
            end

        end
    else
        # ... tmp1 = u x tmp2 (in color space) 
        for ia = 1:NC^2-1
            for ialpha = 1:4
                for is = 1:NV
                    for k1 = 1:NC
                        tmp1[k1, is, ialpha] = 0


                        for k2 = 1:NC
                            tmp1[k1, is, ialpha] += u[ia][k1, k2, is] * tmp2[k2, is, ialpha]
                        end

                    end
                end
            end



            # ... v1^{\dagger} * tmp1
            is = 0
            for it = 1:NT

                for iz = 1:NZ
                    for iy = 1:NY
                        for ix = 1:NX
                            is += 1

                            for k1 = 1:NC
                                z[ia, ix, iy, iz, it] += real(
                                    conj(v1[k1, is, 1]) * tmp1[k1, is, 1] +
                                    conj(v1[k1, is, 2]) * tmp1[k1, is, 2] +
                                    conj(v1[k1, is, 3]) * tmp1[k1, is, 3] +
                                    conj(v1[k1, is, 4]) * tmp1[k1, is, 4],
                                )
                            end

                            #println(z[ia,ix,iy,iz,it])
                        end
                    end
                end
            end


        end
    end

    #println("z = ",z*z)

end

function cal_dFμν!(
    dFμν1,
    dFμν2,
    U,
    fparam::FermiActionParam_WilsonClover,
    work1::T,
    work2,
    work3,
    work4,
    temp1,
    temp2,
    temp3,
    temp4,
    μ,
    ν,
    iflag,
) where {T<:GaugeFields_1d}
    coe = im * 0.125 * fparam.hop * fparam.Clover_coefficient
    NC = U[1].NC
    #numbasis = ifelse(NC == 3,8,3)
    numbasis = NC^2 - 1


    #=
        ...  A) the upper staple  ....................
                                   w3
           nu                 .-----------+
                              |           | 
           /|             w4  |           | w2
            |                 |           |
            |                 |           |
            +----> mu         .-----------.

    =#
    if iflag == 1 || iflag == 2
        substitute!(work1, U[μ])
        gauge_shift!(work2, μ, U[ν])
        gauge_shift!(work3, ν, U[μ])
        substitute!(work4, U[ν])
    end

    # ...  iflag = 1 (Case 1 and 3)
    if iflag == 1
        #=
                                       w3
                                  .-----------+
          Case 1                  |           |
               nu                 |           |
                |              w4 |           | w2
                |                 |           |
                +----> mu         o-----------. 
                                  x    w1 

                                       w3
                                  .<----------o
          Case 3                  |           |
               nu                 |           |
                |              w4 |           | w2
                |                 |           |
                +----> mu         .-----------. 

        =#
        mul!(temp1, work1, work2)
        mul!(temp2, work4, work3)

        for ia = 1:numbasis
            lambdamul(temp3, temp1, ia, fparam.SUNgenerator)
            mul!(temp4, temp3, temp2')
            mul!(temp4, im)
            cimaglnk!(temp4)
            substitute!(dFμν1[ia], coe, temp4)

            mul!(temp4, temp2', temp3)
            mul!(temp4, im)
            cimaglnk!(temp4)
            substitute!(dFμν2[ia], coe, temp4)

        end


    end
    # ...  iflag = 2 (Case 2 and 4)
    if iflag == 2
        #=
                                        w3
                                  .-----------+
          Case 2                  |           |
               nu                 |           |
                |             w4  |           | w2
                |                 |           |
                +----> mu         .---------->o 
                                  x    w1 

                                        w3
                                  o-----------+
          Case 4                  |           |
               nu                 |           |
                |              w4 |           | w2
                |                 |           |
                +----> mu         .-----------. 
                                  x    w1 
        =#
        mul!(temp1, work2, work3')

        for ia = 1:numbasis
            lambdamul(temp2, work1, ia, fparam.SUNgenerator) # temp2=(lambda_a/2)*work1
            mul!(temp3, work4', temp2)
            mul!(temp4, temp1, temp3)
            mul!(temp4, im)
            cimaglnk!(temp4)
            substitute!(dFμν1[ia], coe, temp4)

            mul!(temp4, temp3, temp1)
            mul!(temp4, im)
            cimaglnk!(temp4)
            substitute!(dFμν2[ia], coe, temp4)
        end

    end

    #=
        ...  B) the lower staple  ....................

           nu
           /|
            |
            |
            |                 x    w1
            +----> mu         .-----------+
                              |           |
                          w4  |           | w2
                              |           |
                              |           |
                              .-----------.

    =#
    if iflag == 3 || iflag == 4
        substitute!(work1, U[μ])
        idir2 = (μ, -ν)
        gauge_shift!(work2, idir2, U[ν])
        gauge_shift!(work3, -ν, U[μ])
        gauge_shift!(work4, -ν, U[ν])
        #println(work2[1,1,1])
        #exit()
    end
    #   ...  iflag = 3 (Case 4' and 2')
    if iflag == 3
        #=
          Case 4'
               nu
                |
                |                 x    w1
                +----> mu         o<----------. 
                                  |           |
                              w4  |           | w2
                                  |           |
                                  |           |
                                  .-----------.
                                       w3
          Case 2'
               nu
                |
                |                 x    w1
                +----> mu         .<----------. 
                                  |           |
                              w4  |           | w2
                                  |           |
                                  |           |
                                  .-----------o
                                       w3
        =#

        mul!(temp1, work4', work3)

        for ia = 1:numbasis
            lambdamul(temp2, work1, ia, fparam.SUNgenerator) # temp2=(lambda_a/2)*work1
            mul!(temp3, work2, temp2')

            mul!(temp4, temp1, temp3)
            mul!(temp4, -im)
            cimaglnk!(temp4)
            substitute!(dFμν1[ia], coe, temp4)

            mul!(temp4, temp3, temp1)
            mul!(temp4, -im)
            cimaglnk!(temp4)
            substitute!(dFμν2[ia], coe, temp4)
        end
    end

    # ...  iflag = 4 (Case 3' and 1')
    if iflag == 4
        #=
          Case 3' 
               nu
                |
                |                 x    w1
                +----> mu         .<----------o 
                                  |           |
                              w4  |           | w2
                                  |           |
                                  |           |
                                  .-----------.
                                       w3
          Case 1' 
               nu
                |
                |                 x    w1
                +----> mu         .<----------. 
                                  |           |
                              w4  |           | w2
                                  |           |
                                  |           |
                                  o-----------.

        =#

        mul!(temp1, work3, work2)

        for ia = 1:numbasis
            lambdamul(temp2, work1, ia, fparam.SUNgenerator) # temp2=(lambda_a/2)*work1
            mul!(temp3, work4, temp2)

            mul!(temp4, temp3', temp1)
            mul!(temp4, -im)
            cimaglnk!(temp4)
            #println(temp4.g)
            #exit()

            substitute!(dFμν1[ia], coe, temp4)

            mul!(temp4, temp1, temp3')
            mul!(temp4, -im)
            cimaglnk!(temp4)

            substitute!(dFμν2[ia], coe, temp4)
        end
    end
    return



end

"""
    x <= x - x_aj
"""
function cimaglnk!(x::SUNGaugeFields_1d)
    NV = x.NV
    NC = x.NC
    for i = 1:NV
        for k2 = 1:NC
            for k1 = k2:NC
                xdiff = x[k1, k2, i] - conj(x[k2, k1, i])
                x[k1, k2, i] = xdiff
                if k1 != k2
                    x[k2, k1, i] = -conj(xdiff)
                end
            end
        end
    end
end

function cimaglnk!(x::SU3GaugeFields_1d)
    NV = x.NV
    for i = 1:NV
        z11 = x[1, 1, i] - conj(x[1, 1, i])
        z12 = x[1, 2, i] - conj(x[2, 1, i])
        z13 = x[1, 3, i] - conj(x[3, 1, i])

        z22 = x[2, 2, i] - conj(x[2, 2, i])
        z23 = x[2, 3, i] - conj(x[3, 2, i])

        z33 = x[3, 3, i] - conj(x[3, 3, i])

        x[1, 1, i] = z11
        x[1, 2, i] = z12
        x[1, 3, i] = z13
        x[2, 1, i] = -conj(z12)
        x[2, 2, i] = z22
        x[2, 3, i] = z23
        x[3, 1, i] = -conj(z13)
        x[3, 2, i] = -conj(z23)
        x[3, 3, i] = z33
    end
end


function cimaglnk!(x::SU2GaugeFields_1d)
    NV = x.NV
    for i = 1:NV
        z11 = x[1, 1, i] - conj(x[1, 1, i])
        z12 = x[1, 2, i] - conj(x[2, 1, i])

        z22 = x[2, 2, i] - conj(x[2, 2, i])

        x[1, 1, i] = z11
        x[1, 2, i] = z12
        x[2, 1, i] = -conj(z12)
        x[2, 2, i] = z22

    end
end





end
