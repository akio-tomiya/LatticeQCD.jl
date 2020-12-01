function set_wing!(u::Array{T,1}) where T <: GaugeFields
    for μ=1:4
        set_wing!(u[μ])
    end
end

function set_wing!(u::GaugeFields)
    NT = u.NT
    NY = u.NY
    NZ = u.NZ
    NX = u.NX
    NDW = u.NDW
    NC = u.NC

    #X direction 
    #Now we send data

    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for id=1:NDW
                    for k2=1:NC
                        for k1=1:NC
                            u[k1,k2,-NDW+id,iy,iz,it] = u[k1,k2,NX+(id-NDW),iy,iz,it]
                        end
                    end
                end
            end
        end
    end

    for it=1:NT
        for iz=1:NZ
            for iy=1:NY
                for id=1:NDW
                    for k2=1:NC
                        for k1=1:NC
                            u[k1,k2,NX+id,iy,iz,it] = u[k1,k2,id,iy,iz,it]
                        end
                    end
                end
            end
        end
    end


    #Y direction 
    #Now we send data
    for it=1:NT
        for iz=1:NZ
            for ix=-NDW+1:NX+NDW
                for id=1:NDW
                    for k1=1:NC
                        for k2=1:NC
                            u[k1,k2,ix,-NDW+id,iz,it] = u[k1,k2,ix,NY+(id-NDW),iz,it]
                        end
                    end
                end
            end
        end
    end

    for it=1:NT
        for iz=1:NZ
            for ix=-NDW+1:NX+NDW
                for id=1:NDW
                    for k1=1:NC
                        for k2=1:NC
                            u[k1,k2,ix,NY+id,iz,it] = u[k1,k2,ix,id,iz,it]
                        end
                    end
                end
            end
        end
    end

    #Z direction 
    #Now we send data
    for id=1:NDW
        for it=1:NT
            for iy=-NDW+1:NY+NDW
                for ix=-NDW+1:NX+NDW
                    for k1=1:NC
                        for k2=1:NC
                            u[k1,k2,ix,iy,id-NDW,it] = u[k1,k2,ix,iy,NZ+(id-NDW),it]
                            u[k1,k2,ix,iy,NZ+id,it] = u[k1,k2,ix,iy,id,it]
                        end
                    end
                end
            end
        end
    end


    for id=1:NDW
        for iz=-NDW+1:NZ+NDW
            for iy=-NDW+1:NY+NDW
                for ix=-NDW+1:NX+NDW
                    for k1=1:NC
                        for k2=1:NC
                            u[k1,k2,ix,iy,iz,id-NDW] = u[k1,k2,ix,iy,iz,NT+(id-NDW)]
                            u[k1,k2,ix,iy,iz,NT+id] = u[k1,k2,ix,iy,iz,id]
                        end
                    end
                end
            end
        end
    end

    #display(u.g)
    #exit()

    return
end