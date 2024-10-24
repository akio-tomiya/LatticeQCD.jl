using Gaugefields
function main(filename)
    NX = 4
    NY = 4
    NZ = 4
    NT = 4
    NC = 4
    Nwing = 1
    Dim = 4

    U = Initialize_Gaugefields(NC, Nwing, NX, NY, NZ, NT, condition="cold")

    ildg = ILDG(filename)
    i = 1
    L = [NX, NY, NZ, NT]
    load_gaugefield!(U, i, ildg, L, NC)

    save_textdata(U, filename * ".txt")
end
main(ARGS[1])