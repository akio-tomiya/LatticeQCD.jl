

function set_wing_fermi!(a::FermionFields)
    NT = a.NT
    NZ = a.NZ
    NY = a.NY
    NX = a.NX
    NC = a.NC

    #!  X-direction
    for ialpha=1:4
        for it=1:NT
            for iz = 1:NZ
                for iy=1:NY
                    for k=1:NC
                        a[k,0,iy,iz,it,ialpha] = a.BoundaryCondition[1]*a[k,NX,iy,iz,it,ialpha]
                    end
                end
            end
        end
    end

    for ialpha=1:4
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for k=1:NC
                        a[k,NX+1,iy,iz,it,ialpha] = a.BoundaryCondition[1]*a[k,1,iy,iz,it,ialpha]
                    end
                end
            end
        end
    end

    #Y-direction
    for ialpha = 1:4
        for it=1:NT
            for iz=1:NZ
                for ix=1:NX
                    for k=1:NC
                        a[k,ix,0,iz,it,ialpha] = a.BoundaryCondition[2]*a[k,ix,NY,iz,it,ialpha]
                    end
                end
            end
        end
    end

    for ialpha=1:4
        for it=1:NT
            for iz=1:NZ
                for ix=1:NX
                    for k=1:NC
                        a[k,ix,NY+1,iz,it,ialpha] = a.BoundaryCondition[2]*a[k,ix,1,iz,it,ialpha]
                    end
                end
            end
        end
    end

    
    for ialpha=1:4
        # Z-direction
        for it=1:NT
            for iy=1:NY
                for ix=1:NX
                    for k=1:NC
                        a[k,ix,iy,0,it,ialpha] = a.BoundaryCondition[3]*a[k,ix,iy,NZ,it,ialpha]
                        a[k,ix,iy,NZ+1,it,ialpha] = a.BoundaryCondition[3]*a[k,ix,iy,1,it,ialpha]

                    end
                end
            end
        end

        #T-direction
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    for k=1:NC
                        a[k,ix,iy,iz,0,ialpha] = a.BoundaryCondition[4]*a[k,ix,iy,iz,NT,ialpha]
                        a[k,ix,iy,iz,NT+1,ialpha] = a.BoundaryCondition[4]*a[k,ix,iy,iz,1,ialpha]
                    end
                end
            end
        end

    end

    

    

end


function set_wing_fermi!(a::StaggeredFermion)
    NT = a.NT
    NZ = a.NZ
    NY = a.NY
    NX = a.NX
    NC = a.NC

    #!  X-direction
    for ialpha=1:1
        for it=1:NT
            for iz = 1:NZ
                for iy=1:NY
                    for k=1:NC
                        a[k,0,iy,iz,it,ialpha] = a.BoundaryCondition[1]*a[k,NX,iy,iz,it,ialpha]
                    end
                end
            end
        end
    end

    for ialpha=1:1
        for it=1:NT
            for iz=1:NZ
                for iy=1:NY
                    for k=1:NC
                        a[k,NX+1,iy,iz,it,ialpha] = a.BoundaryCondition[1]*a[k,1,iy,iz,it,ialpha]
                    end
                end
            end
        end
    end

    #Y-direction
    for ialpha = 1:1
        for it=1:NT
            for iz=1:NZ
                for ix=1:NX
                    for k=1:NC
                        a[k,ix,0,iz,it,ialpha] = a.BoundaryCondition[2]*a[k,ix,NY,iz,it,ialpha]
                    end
                end
            end
        end
    end

    for ialpha=1:1
        for it=1:NT
            for iz=1:NZ
                for ix=1:NX
                    for k=1:NC
                        a[k,ix,NY+1,iz,it,ialpha] = a.BoundaryCondition[2]*a[k,ix,1,iz,it,ialpha]
                    end
                end
            end
        end
    end

    
    for ialpha=1:1
        # Z-direction
        for it=1:NT
            for iy=1:NY
                for ix=1:NX
                    for k=1:NC
                        a[k,ix,iy,0,it,ialpha] = a.BoundaryCondition[3]*a[k,ix,iy,NZ,it,ialpha]
                        a[k,ix,iy,NZ+1,it,ialpha] = a.BoundaryCondition[3]*a[k,ix,iy,1,it,ialpha]

                    end
                end
            end
        end

        #T-direction
        for iz=1:NZ
            for iy=1:NY
                for ix=1:NX
                    for k=1:NC
                        a[k,ix,iy,iz,0,ialpha] = a.BoundaryCondition[4]*a[k,ix,iy,iz,NT,ialpha]
                        a[k,ix,iy,iz,NT+1,ialpha] = a.BoundaryCondition[4]*a[k,ix,iy,iz,1,ialpha]
                    end
                end
            end
        end

    end

    

    

end



