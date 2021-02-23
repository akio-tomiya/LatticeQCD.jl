# LatticeQCD.jl (ver. 0.5)

[![Build Status](https://travis-ci.com/cometscome/LatticeQCD.jl.svg?branch=master)](https://travis-ci.com/cometscome/LatticeQCD.jl)
[![Coverage](https://codecov.io/gh/cometscome/LatticeQCD.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cometscome/LatticeQCD.jl)
[![Coverage](https://coveralls.io/repos/github/cometscome/LatticeQCD.jl/badge.svg?branch=master)](https://coveralls.io/github/cometscome/LatticeQCD.jl?branch=master)


![LatticeQCD.jl](logo.png)

This code enabales you to perform lattice QCD calculations! Native Julia code for Lattice QCD.

- [What is lattice QCD? (PDG)](https://pdg.lbl.gov/2019/reviews/rpp2018-rev-lattice-qcd.pdf) : Lattice regulated quantum chromo-dynamics used in high energy physics.

- [What is Julia?](https://julialang.org/) : An easy and fast scientific programming launguage with the JIT compilar


![demo](demo.gif)


# Quick start

You can start lattice QCD in 5 steps!



1.Download a Julia binary from [Julialang.org](https://julialang.org/downloads/). Set a path to the binary of Julia. 



2.In Julia REPL, push "]" key to enter the package mode and type

```
add https://github.com/akio-tomiya/LatticeQCD.jl
```
and "return" key. Press "backspace" key ( "delete" key for mac) to exit the package mode.

(All dependence will be solved automatically)



3.Include the package with

```
using LatticeQCD
```



4.Make a parameter file with wizard,

```
run_wizard()
```

Choose parameters as you want!



5.Start simulation with created your parameter file!

```
 run_LQCD("my_parameters.jl")
```

You'll get results!

Of cource, you can write/modify a parameter file by yourself.

Enjoy life with lattice QCD.

To see our demonstration above, execute,
```julia
using Plots
using LatticeQCD 
demo()
```
it takes time before showing up the window.

# What is supported?

We support 4 dimensional euclidean spacetime and following things:

- Gauge field
  - Optimized SU(2), SU(3)
  - General SU(N)
  - General gauge action = plaquette + rect + etc action
- Fermions
  - Wilson/Wilson-clover (2flavor)
  - Staggered fermion (4tastes)
- Configuration generation algorithms
  - Cold/Hot start for SU(2), SU(3) and SU(N). One instanton for SU(2)
  - Heatbath for SU(2), SU(3) and SU(N) & overrelaxation for a general gauge action
  - Quenched HMC with SU(2), SU(3), SU(N) for a general gauge action
  - HMC (2 flavor Wilson/Clover) with SU(2), SU(3), SU(N) with a general gauge action
  - HMC (4 flavor Staggered fermions) with SU(2), SU(3), SU(N) with a general gauge action
  - <s>RHMC (any flavor staggered) with SU(2), SU(3) and SU(N) for a general gauge action</s> (This will be supported)
  - Fermion integrated HMC with a general gauge action
  - Self-learning Monte-Carlo with plaquette action (simpler version of [this paper](https://arxiv.org/abs/2010.11900) )
  - Self-learning HMC with plaquette action
- Measurements
  - Plaquette
  - Polyakov loop
  - Chiral condensates (Wilson, Wilson-Clover/Staggered)
  - Momentum projected pion correlator (Wilson& Wilson-Clover fermion)
  - Topological charge(plaquette and clover definition)
  - Load & measurement mode (load and measure all configurations in a directory)
- Smearing
  - <s>APE, stout</s> (This will be public)
  - <s>Gradient flow for a genelic action</s>   (Not well tested)
- I/O for configurations
  - ILDG format
  - JLD format (default binary file for Julia, one of HDF5)

Many of smearing and improved fermion actions and parallelazation have not supported yet.


# USAGE/User interface

We support following two user interfaces

1. Julia REPL interface (For beginners, just after the lattice QCD textbook)
2. Genral run code (Experience with another code, for batch job, customised purpose)

Usage 1 was already explained. 

For Usage 2, in Julia REPL, push "]" key to enter the package mode and type

```
add https://github.com/akio-tomiya/LatticeQCD.jl
```
Then, LatticeQCD.jl is installed on your machine. 

You write the following code:
```julia:run.jl
using LatticeQCD
run_LQCD(ARGS[1])
```

Then, you can execute like

```
julia run.jl PARAMETER_FILE
```

then, you get results though standard I/O. 

The "PARAMETER_FILE" can be created through the wizard. To use the wizard on the shell, you write the following code: 

```julia:wizard.jl
using LatticeQCD
run_wizard()
```

Then, you can run the wizard:

```
julia wizard.jl
```

# Purpose of the code
We develop this code to achive following things:

1. Good portability (If one has Julia, this code is runnable. All dependence is under control.)
2. Easy to start/ pedagogical (start in 10 minutes)
3. Suite (configuration generation with and without fermions, and measurements)
4. Easy to modify (Good for prototyping)
5. Compatitive speed with Fortran 90 codes

This is the first open source Julia code for lattice QCD. High performance is out of our scope.

# How has it been tested?

We compared results to following papers/codes 

- SU(3) staggered HMC with https://inspirehep.net/literature/283285

# Reference

We refer "Lattice Tool Kit" https://nio-mon.riise.hiroshima-u.ac.jp/LTK/ written in Fortran 90.
