module Print_config
import ..LTK_universe: Universe
function write_config(filename::String, univ::Universe)
    fp = open(filename, "w")
    NT = univ.NT
    NZ = univ.NZ
    NY = univ.NY
    NX = univ.NX
    NC = univ.NC

    for it = 1:NT
        for iz = 1:NZ
            for iy = 1:NY
                for ix = 1:NX
                    for μ = 1:4
                        for a = 1:NC
                            for b = 1:NC
                                #[:,:,ix,iy,iz,it]
                                write(fp, real(univ.U[μ][a, b, ix, iy, iz, it]))
                                write(fp, imag(univ.U[μ][a, b, ix, iy, iz, it]))
                            end
                        end

                    end
                end
            end
        end
    end
    close(fp)


end
end
