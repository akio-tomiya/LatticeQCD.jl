# LatticeQCD.jl
This is the lattice QCD package purely written in Julia language.
Lattice QCD is a well-established non-perturbative approach to solving the quantum chromodynamics (QCD) theory of quarks and gluons.

We confirmed that it works in Julia 1.5 or later. 

This code is inspired by the Lattice Tool Kit (LTK) written in [Fortran](https://nio-mon.riise.hiroshima-u.ac.jp/LTK/).
With the use of a modern programing language, it is easy to understand how the code works. 
The part of the codes is translated from the LTK. 


What LatticeQCD.jl can do is :

- Two flavor SU(3) Hybrid Monte Carlo with the Wilson Fermion (We want to add more functionalities.). 
- The speed is the faster than the original LTK code written in fortran. 

## How to do 
### simple version
You can try this code with 

```
julia ./src/run.jl
```

or 

```
julia ./src/run.jl params.jl
```
Here, in params.jl, we can see 

```julia
L = (4,4,4,4)
β = 6
#gparam = Setup_Gauge_action(β)

NTRACE = 3
gparam =  GaugeActionParam_standard(β,NTRACE)
#gparam = Setup_Gauge_action(β)

hop= 0.141139#Hopping parameter
r= 1#Wilson term
eps= 1e-19
Dirac_operator= "Wilson"
MaxCGstep= 3000
#fparam = Setup_Fermi_action()
fparam = FermiActionParam_Wilson(hop,r,eps,Dirac_operator,MaxCGstep)
```
You can change the parameters. 
If you set ```fparam=nothing``` in this file, you can do the quench HMC.

### more details

In the LatticeQCD.jl, the Universe type is an important type for simulations. 
At first, you have to generate your "universe". 

```julia
univ = Universe(file)
```
The file is like: 

```julia
L = (4,4,4,4)
β = 6
NTRACE = 3
#gparam = Setup_Gauge_action(β)
gparam =  GaugeActionParam_standard(β,NTRACE)

BoundaryCondition=[1,1,1,-1]
Nwing = 1
initial="cold"
NC =3


hop= 0.141139#Hopping parameter
r= 1#Wilson term
eps= 1e-19
Dirac_operator= "Wilson"
MaxCGstep= 3000

fparam = FermiActionParam_Wilson(hop,r,eps,Dirac_operator,MaxCGstep)

```
The parameters that you do not provide are set by the default values. 

Then, you can calculate physical obserbables. 
For example, if you want to calculate a plaquette, just do 

```julia
plaq = calc_plaquette(univ)
println("plaq = ",plaq)
```

If you want to do the HMC simulation, set the MD parameters: 

```julia
Δτ = 0.1
MDsteps = 10
βMD = β

mdparams =MD_parameters_standard(gparam,Δτ,MDsteps,βMD)
```

and do it like: 

```julia
for i=1:20
    Sold = md_initialize!(univ)
    Snew = md!(univ,mdparams)

    metropolis_update!(univ,Sold,Snew)
    plaq = calc_plaquette(univ)
    println("-------------------------------------")
    println("$i-th plaq = ",plaq)
    println("-------------------------------------")
end 
```



## Benchmarks
The speed is important for the Lattice QCD simulation. 
Remarkably, this Julia code is faster than s similar code in Fortran. 

For example, with the following parameters: 

```
6.0d0     6.0d0         beta, betamd
0.141139d0  1.d0         hop,  r (Hopping parametger, Wilson term)
.false.                   Clover term
(0.0d0,0.0d0)            cmu (Chemical potential)
1                        istart (1:Cold, 2:Hot, 3:File)
001       020            ntraj0, ntraj1 
1.d0                     gamma_G
10     0.1d0            nstep, dtau
.true.                   fermions
0                        flagMD
```
, where this is an input file of the LTK. 
We set eps=1e-19 as a convergence criteria in a CG solver. 

The elapsed time of the original LTK code on Mac mini (2018) with 3.2Ghz Intel Core i7 (6 cores) is 

```
Eold,Enew,Diff,accept:   0.1613911E+04    0.1614279E+04   -0.3678492E+00   T
 Plaq :   0.60513145439721250     
 Pol :         (1.1156431755476819,-3.20744714128515240E-002)
./a.out < input  227.40s user 0.07s system 99% cpu 3:47.51 total
```

On the other hand, the elapsed time of this LatticeQCD.jl is 

```
-------------------------------------
20-th plaq = 0.6051309989225465
-------------------------------------
julia --sysimage ~/sys_plots.so run.jl  180.41s user 0.25s system 100% cpu 3:00.62 total
```
The LatticeQCD.jl is faster than the Fortran-based code. 

We note that the plaquette value is consistently in single precision floating point numbers, 
since the random number generation is based on the original Fortran code and random numbers are in single precision floating point.



```@autodocs
Modules = [LatticeQCD.LTK_universe]
```

```@docs
Setup_Gauge_action
Setup_Fermi_action
```