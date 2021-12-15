module Bridge_format
    import ..Gaugefields:GaugeFields,SU3GaugeFields,SU2GaugeFields,set_wing!
#=
Bridge++ Text file format
U(x,y,z,t;mu)_{ab} is like 
(outer) [mu][t][z][y][x][a][b][re/im] (inner)
re/im; 
  0: Real
  1: Imaginary
mu: (x,y,z,t)=(0,1,2,3)

Re of U(0,0,0,0;mu=0)_00  # mu=0, site (x,y,z,t)=(0,0,0,0)
Im of U(0,0,0,0;mu=)_00
Re of U(0,0,0,0;mu=0)_01
Im of U(0,0,0,0;mu=0)_01
...
Re of U(1,0,0,0;0)_00     # mu=0, site (x,y,z,t)=(1,0,0,0)
Im of U(1,0,0,0;0)_00
...
Re of U(0,0,0,0;mu=1)_00  # mu=1, site (x,y,z,t)=(0,0,0,0)
Im of U(0,0,0,0;mu=1)_00

=#
    function load_BridgeText!(initial,U,L,NC)
        NX = L[1]
        NY = L[2]
        NZ = L[3]
        NT = L[4]
        @assert U[1].NX == NX "NX mismatch"
        @assert U[1].NY == NY "NY mismatch"
        @assert U[1].NZ == NZ "NZ mismatch"
        @assert U[1].NT == NT "NT mismatch"
        @assert U[1].NC == NC "NC mismatch"
        fp = open(initial,"r")
        numdata = countlines(initial)
        @assert numdata == 4*NX*NY*NT*NZ*NC*NC*2 "data shape is wrong"


        for μ=1:4
            for it=1:NT
                for iz=1:NZ
                    for iy=1:NY
                        for ix=1:NX
                            for a=1:NC
                                for b=1:NC
                                    u = split(readline(fp))
                                    #println(u)
                                    rvalue = parse(Float64,u[1])
                                    u = split(readline(fp))
                                    ivalue = parse(Float64,u[1])
                                    U[μ][a,b,ix,iy,iz,it] = rvalue + im*ivalue
                                end
                            end
                            #u = U[μ][:,:,ix,iy,iz,it] 
                            #println(u'*u)
                        end
                    end
                end
            end
        end

        
        set_wing!(U)
 

        close(fp)


    end
end