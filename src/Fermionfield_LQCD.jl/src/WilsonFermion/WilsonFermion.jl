include("./WilsonFermion_4D_wing.jl")

function WilsonFermion(params,NC,NN...)
    Dim = length(NN)
    if Dim == 4
        fermion = WilsonFermion_4D_wing(params,NC,NN...)
    else
        error("Dimension $Dim is not supported. NN = $NN")
    end
    return fermion
end


function WilsonFermion(fparam::FermiActionParam,NC,BoundaryCondition,NN...)
    Dim = length(NN)
    if Dim == 4
        fermion = WilsonFermion_4D_wing(NC,NN...,fparam.r,fparam.hop,fparam.eps,
            fparam.MaxCGstep,BoundaryCondition)
        #fermion = StaggeredFermion_4D_wing(params,NC,NN...)
    else
        error("Dimension $Dim is not supported. NN = $NN")
    end
    return fermion
end

"""
mk_gamma()
c----------------------------------------------------------------------c
c     Make gamma matrix
c----------------------------------------------------------------------c
C     THE CONVENTION OF THE GAMMA MATRIX HERE
C     ( EUCLIDEAN CHIRAL REPRESENTATION )
C
C               (       -i )              (       -1 )
C     GAMMA1 =  (     -i   )     GAMMA2 = (     +1   )
C               (   +i     )              (   +1     )
C               ( +i       )              ( -1       )
C
C               (     -i   )              (     -1   )
C     GAMMA3 =  (       +i )     GAMMA4 = (       -1 )
C               ( +i       )              ( -1       )
C               (   -i     )              (   -1     )
C
C               ( -1       )
C     GAMMA5 =  (   -1     )
C               (     +1   )
C               (       +1 )
C
C     ( GAMMA_MU, GAMMA_NU ) = 2*DEL_MU,NU   FOR MU,NU=1,2,3,4   
c----------------------------------------------------------------------c
"""
function mk_gamma(r)
    g0 = zeros(ComplexF64,4,4)
    g1 = zero(g0)
    g2 = zero(g1)
    g3 = zero(g1)
    g4 = zero(g1)
    g5 = zero(g1)
    gamma = zeros(ComplexF64,4,4,5)
    rpg = zero(gamma)
    rmg = zero(gamma)


    g0[1,1]=1.0; g0[1,2]=0.0; g0[1,3]=0.0; g0[1,4]=0.0
    g0[2,1]=0.0; g0[2,2]=1.0; g0[2,3]=0.0; g0[2,4]=0.0
    g0[3,1]=0.0; g0[3,2]=0.0; g0[3,3]=1.0; g0[3,4]=0.0
    g0[4,1]=0.0; g0[4,2]=0.0; g0[4,3]=0.0; g0[4,4]=1.0

    g1[1,1]=0.0; g1[1,2]=0.0; g1[1,3]=0.0; g1[1,4]=-im
    g1[2,1]=0.0; g1[2,2]=0.0; g1[2,3]=-im;  g1[2,4]=0.0
    g1[3,1]=0.0; g1[3,2]=+im;  g1[3,3]=0.0; g1[3,4]=0.0
    g1[4,1]=+im;  g1[4,2]=0.0; g1[4,3]=0.0; g1[4,4]=0.0

    g2[1,1]=0.0; g2[1,2]=0.0; g2[1,3]=0.0; g2[1,4]=-1.0
    g2[2,1]=0.0; g2[2,2]=0.0; g2[2,3]=1.0; g2[2,4]=0.0
    g2[3,1]=0.0; g2[3,2]=1.0; g2[3,3]=0.0; g2[3,4]=0.0
    g2[4,1]=-1.0;g2[4,2]=0.0; g2[4,3]=0.0; g2[4,4]=0.0

    g3[1,1]=0.0; g3[1,2]=0.0; g3[1,3]=-im;  g3[1,4]=0.0
    g3[2,1]=0.0; g3[2,2]=0.0; g3[2,3]=0.0; g3[2,4]=+im
    g3[3,1]=+im;  g3[3,2]=0.0; g3[3,3]=0.0; g3[3,4]=0.0
    g3[4,1]=0.0; g3[4,2]=-im;  g3[4,3]=0.0; g3[4,4]=0.0

    g4[1,1]=0.0; g4[1,2]=0.0; g4[1,3]=-1.0;g4[1,4]=0.0
    g4[2,1]=0.0; g4[2,2]=0.0; g4[2,3]=0.0; g4[2,4]=-1.0
    g4[3,1]=-1.0;g4[3,2]=0.0; g4[3,3]=0.0; g4[3,4]=0.0
    g4[4,1]=0.0; g4[4,2]=-1.0;g4[4,3]=0.0; g4[4,4]=0.0

    g5[1,1]=-1.0;g5[1,2]=0.0; g5[1,3]=0.0; g5[1,4]=0.0
    g5[2,1]=0.0; g5[2,2]=-1.0;g5[2,3]=0.0; g5[2,4]=0.0
    g5[3,1]=0.0; g5[3,2]=0.0; g5[3,3]=1.0; g5[3,4]=0.0
    g5[4,1]=0.0; g5[4,2]=0.0; g5[4,3]=0.0; g5[4,4]=1.0

    gamma[:,:,1] = g1[:,:]
    gamma[:,:,2] = g2[:,:]
    gamma[:,:,3] = g3[:,:]
    gamma[:,:,4] = g4[:,:]
    gamma[:,:,5] = g5[:,:]

    for mu=1:4
        for j=1:4
            for i=1:4
                rpg[i,j,mu] = r*g0[i,j] + gamma[i,j,mu]
                rmg[i,j,mu] = r*g0[i,j] - gamma[i,j,mu]
            end
        end
    end 

    return gamma,rpg,rmg


end