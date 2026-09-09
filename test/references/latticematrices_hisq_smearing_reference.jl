using LatticeQCD
import LatticeDiracOperators

const LM = LatticeDiracOperators.LatticeMatrices

const LATTICE = (4, 4, 4, 4)
const NC = 3
const NAIK_EPSILON = -0.083

function deterministic_links()
    links = [zeros(ComplexF64, NC, NC, LATTICE...) for _ in 1:4]
    for site in CartesianIndices(LATTICE)
        x = Tuple(site)
        coordinate = x[1] + 3x[2] + 5x[3] + 7x[4]
        for mu in 1:4, column in 1:NC, row in 1:NC
            re = 0.05 * 0.013 *
                 (2row - column + coordinate + 3mu)
            im = 0.05 * 0.017 *
                 (row + 2column - coordinate + mu)
            links[mu][row, column, x...] =
                complex(re + (row == column), im)
        end
    end
    return links
end

function fingerprint(field)
    values = vec(field)
    return (
        sum(real, values),
        sum(imag, values),
        sum(i * real(values[i]) for i in eachindex(values)),
        sum(i * imag(values[i]) for i in eachindex(values)),
        sum(abs2, values),
    )
end

thin = [
    LM.LatticeMatrix(link, 4, (1, 1, 1, 1);
        nw=3, comm0=LM.SerialCommunicator())
    for link in deterministic_links()
]
links = LM.hisq_links_from_thin(thin; naik_epsilon=NAIK_EPSILON)

for (stage, fields) in (("fat", links.fat_links), ("long", links.long_links))
    for mu in 1:4
        values = fingerprint(LM.gather_matrix(fields[mu]))
        println(
            "FINGERPRINT implementation=LatticeMatrices",
            " stage=", stage,
            " mu=", mu,
            " sum_re=", values[1],
            " sum_im=", values[2],
            " weighted_re=", values[3],
            " weighted_im=", values[4],
            " norm2=", values[5],
        )
    end
end
