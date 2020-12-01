# LatticeQCD.jl

[![Build Status](https://travis-ci.com/cometscome/LatticeQCD.jl.svg?branch=master)](https://travis-ci.com/cometscome/LatticeQCD.jl)
[![Coverage](https://codecov.io/gh/cometscome/LatticeQCD.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cometscome/LatticeQCD.jl)
[![Coverage](https://coveralls.io/repos/github/cometscome/LatticeQCD.jl/badge.svg?branch=master)](https://coveralls.io/github/cometscome/LatticeQCD.jl?branch=master)


![LatticeQCD.jl](logo.png)

This code enabales you to perform lattice QCD calculations!

[What is lattice QCD? (PDG)](https://pdg.lbl.gov/2019/reviews/rpp2018-rev-lattice-qcd.pdf)



# Quick start

You can start lattice QCD in 5 steps!



1.Download a Julia binary from [Julialang.org](https://julialang.org/)



2.In Julia REPL, push "]" key to enter the package mode

```
add https://github.com/akio-tomiya/LatticeQCD.jl
```
press ESC key to exit the package mode.

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
 run_LQCD("parameters_used.jl")
```

You'll get results!

Of cource, you can write/modify a parameter file by yourself.

Enjoy life with lattice QCD.



# What is supported?

- Configuration generation
  - Cold/Hot start, One instanton for SU(2)
  - Heatbath for SU(2) and SU(3) for plaquette gauge action
  - <s>Quenched HMC with SU(2) or SU(3) for general gauge action </s> (This will be supported)
  - HMC (2 flavor Wilson/Clover) with SU(2) or SU(3) for plaquette gauge action
  - HMC (4 flavor Staggered fermions) with SU(2) or SU(3) for plaquette gauge action
  - <s>RHMC (any flavor staggered) with SU(2) or SU(3) for general gauge action</s> (This will be supported)
  - <s> General gauge action = plaquette+rect + etc action </s>  (This will be supported)
  - SLMC with plaquette action
  - Load & measurement mode (load and measure all configurations in a directory)
- Measurements
  - Plaquette
  - Polyakov loop
  - Chiral condensates (Wilson/Staggered)
  - Momentum projected pion correlator (Wilson fermion)
  - <s>Topological charge with the Wilson flow (Clover operaor definition)</s>   (Not well tested)
- I/O for configurations
  - <s>ILDG format</s> (This will be supported)
  - JLD format (Defoult binary file for Julia, one of HDF5)

Many of smearing and improved fermion actions have not supported yet.


# USAGE/User interface

We support following user interfaces

1. Julia REPL interface (For beginners, just after the lattice QCD textbook)
2. Genral run code (Experience with another code, for batch job, customised purpose)

Usage 1 was already explained. You can call wizard through Julia Repel.

For Usage 2, you can clone this repository with Git commands.
And in the directory, you can execute like,

```
julia run.jl PARAMETER_FILE
```

then, you get results though standard I/O.


# Purpose of the code
We develop this code to achive following thing:

1. Good portability (If one has Julia, this code is runnable)
2. Easy to start/ pedagogical
3. Suite
4. Easy to modify (Good for prototyping)
5. Compatitive speed with Fortran 90 codes

First open source Julia code for lattice QCD. High performance is out of our scope.

# How has it been tested?

We compared results to following papers/codes 

- SU(3) staggered HMC with https://inspirehep.net/literature/283285

# Reference

We refer "Lattice Tool Kit" https://nio-mon.riise.hiroshima-u.ac.jp/LTK/ written in Fortran 90.
