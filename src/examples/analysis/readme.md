# Analysis code



Please  ``samples.zip`` decompress.



``analysis_sample.jl`` is the analysis code.

```
using LatticeQCD

betas = ["1.40","1.65"]
println("# beta plaquette")
for beta in betas
        dir = "samples/HMC_L08080808_su2nf8_beta$(beta)_Staggered_mass0.015"
        p = get_plaquette_average(dir)
        println("$beta $p #$dir")
end
println("# ")
println("# beta Polyakov loop")
for beta in betas
        dir = "samples/HMC_L08080808_su2nf8_beta$(beta)_Staggered_mass0.015"
        p = get_polyakov_average(dir)
        println("$beta $p #$dir")
end
```

