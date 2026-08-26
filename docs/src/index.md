# LatticeQCD.jl
This is the lattice QCD package purely written in Julia language.
Lattice QCD is a well-established non-perturbative approach to solving the quantum chromodynamics (QCD) theory of quarks and gluons.

We confirmed that it works in Julia 1.5 or later. 

This code is inspired by the Lattice Tool Kit (LTK) written in [Fortran](https://nio-mon.riise.hiroshima-u.ac.jp/LTK/).
With the use of a modern programing language, it is easy to understand how the code works. 
The part of the codes is translated from the LTK. 


What LatticeQCD.jl can do includes:

- Hybrid Monte Carlo with Wilson and Wilson–clover fermions.
- HMC/RHMC with one-link staggered and SU(N) HISQ fermions.
- Standard and Möbius domain-wall fermions (experimental).
- Quenched HMC and heatbath updates for general gauge actions.

## Current HMC interfaces

The navigable Wizard returns a typed `SimulationSpec` without constructing
legacy `Params`:

```julia
using LatticeQCD

spec = run_wizard()
session = build_simulation(spec, GaugefieldsEnvironment())
summary = run!(session)
```

By default `run!` prints rank-zero trajectory progress, HMC acceptance
diagnostics, measurements, configuration saves, and a final summary. A GUI
can use its own event sink without terminal output:

```julia
summary = run!(session; verbose=false)
```

The Wizard also writes the selected TOML file. The conventional file runner
remains available as `run_LQCD("my_parameters.toml")`. For compatibility and
debug comparisons, `run_wizard_legacy()` runs the original Wizard and returns
the historical `Params` object. `run_wizardv2` is retained as an alias of the
new `run_wizard` and is also Param-free.

The typed `Simulation` interface assembles the same HMC update from explicit
configuration objects and runs it through the Gaugefields MD driver.  This is
the preferred interface for notebooks, MPI/GPU applications, and future GUIs:

```julia
using LatticeQCD
import Gaugefields

lattice = LatticeConfig((4, 4, 4, 4))
gauge = GaugeConfig(3, 1, HotStartConfig(1234))
target_action = GaugeActionConfig(
    GaugeActionTermConfig(:gauge_plaquette, "plaquette", 5.7),
)
integrator = LeapfrogConfig(QPQConfig(), ForceGroupConfig(:gauge))
md = MDConfig(0.01, 10, integrator)
momentum = GaussianMomentumConfig(
    1.0,
    RandomStreamConfig(5678, :momentum),
)
acceptance = RankZeroMetropolisConfig(
    RandomStreamConfig(9012, :metropolis),
)
hmc = HMCConfig(md, momentum, acceptance)
input = LQCDConfig(lattice, gauge, target_action, hmc)

environment = GaugefieldsEnvironment(
    communicator=Gaugefields.SerialCommunicator(),
    process_grid=(1, 1, 1, 1),
)
simulation = build_simulation(input, environment)
result = update!(simulation)
```

For self-learning HMC, the target action belongs to `LQCDConfig`, while the
approximate action used only for the MD proposal belongs to `SLHMCConfig`.
The initial and final target Hamiltonians are always used for the exact
Metropolis decision:

```julia
md_action = GaugeActionConfig(
    GaugeActionTermConfig(:gauge_plaquette, "plaquette", 5.2),
)
slhmc = SLHMCConfig(md, momentum, acceptance, md_action)
slhmc_input = LQCDConfig(lattice, gauge, target_action, slhmc)
slhmc_simulation = build_simulation(slhmc_input, environment)
slhmc_result = update!(slhmc_simulation)

# Exact target and approximate-MD diagnostics are both retained.
slhmc_result.delta_hamiltonian
slhmc_result.md_delta_hamiltonian
```

Runtime configurations use portable JLD2 gauge-field output by default.
Pseudofermions are trajectory workspaces and are not written:

```julia
save_configuration(
    "conf_00000001.jld2",
    simulation.configuration;
    format=:jld2,
)
load_configuration!(
    simulation.configuration,
    "conf_00000001.jld2";
    format=:jld2,
)
```

## Typed TOML input without `Params`

Existing Wizard TOML files can be read directly into a `SimulationSpec`.
This path does not construct `Params`, open a logfile, create output
directories, or retain an untyped dictionary in the simulation settings:

```julia
using LatticeQCD

spec = load_simulation_spec("my_parameters.toml")
show_config(spec)

environment = GaugefieldsEnvironment()
session = build_simulation(spec, environment)
summary = run!(session)
```

`load_simulation_spec` also reads the canonical, versioned TOML format. A
Wizard/legacy file can therefore be normalized once and read back without the
legacy section layout:

```julia
spec = load_simulation_spec("my_parameters.toml")
write_simulation_spec("simulation-spec.toml", spec)
restored = load_simulation_spec("simulation-spec.toml")
```

The canonical format represents the lattice, gauge initialization, all gauge
action terms, zero or more fermion actions, MD force groups, update method,
measurements, gradient flow, schedule, and configuration output separately.
It starts with the following format marker so that future schema migrations
are explicit:

```toml
format = "LatticeQCD.SimulationSpec"
schema_version = 1
```

The original `Params` route remains available as a compatibility and debugging
oracle. It deliberately retains the historical side effects of creating log
and measurement directories and opening `parameters.load_fp`:

```julia
using TOML

document = TOML.parsefile("my_parameters.toml")
parameters =
    LatticeQCD.Parameters_TOML.construct_Params_from_TOML(document)
legacy_spec = SimulationSpec(parameters)

try
    # Compare or run the legacy-compatible path here.
finally
    close(parameters.load_fp)
end
```

Both paths use the same typed conversion routines. Omitted legacy values are
taken from the existing `Print_*_parameters` defaults. The test suite compares
the complete `SimulationSpec` produced by both paths for every current Wizard
file and also checks a minimal input in which all optional settings are
omitted.

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
