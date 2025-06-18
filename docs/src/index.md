# FewBodyToolkit.jl

[![Build Status](https://github.com/lhapp27/FewBodyToolkit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lhapp27/FewBodyToolkit.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lhapp27.github.io/FewBodyToolkit.jl/dev/)

A Julia package for solving quantum systems of 2 or 3 particles with general potentials in 1D, 2D, and 3D.

## Features and modules overview

The FewBodyToolkit.jl package currently provides three modules with the following features:

### `GEM2B` - Two-body solver (1D,2D,3D)
* Supports 1D, 2D, and 3D geometries
* Symmetric interaction potentials
* Real-ranged Gaussian basis functions
* Optional complex-ranged Gaussian basis functions for highly excited states
* Complex scaling method (CSM) for resonances
* Basis optimization routine
* Routine for scaling the interaction to yield a desired eigenenergy (inverse problem)
* Coupled-channel problems with additional derivative terms (up to 4th order)
* Output of the wave function

### `GEM3B1D` - Three-body solver (1D)
* Supports 1D three-body problems
* Symmetric two-body pair-interactions, no three-body forces
* Expansion in up to three Fadeev components (rearrangement channels)
* Automatic (anti-)symmetrization for identical bosons or fermions
* Product basis in Jacobi coordinates
* Complex scaling method (CSM) for resonances
* No separate output of the wave function

### `ISGL` - Three-body solver (3D)
* Supports 3D three-body problems
* Central two-body pair-interactions, no three-body forces
* Employs infinitesimally-shifted Gaussian basis functions to handle arbitrary high intrinsic angular momenta (computationally expensive for high values)
* Expansion in up to three Fadeev components (rearrangement channels)
* Automatic (anti-) symmetrization for identical bosons or fermions
* Product basis in Jacobi coordinates
* Complex scaling method (CSM) for resonances
* On-the-fly calculation of central observables and mean-square radii
* No separate output of the wave function


## Installation
To install FewBodyToolkit.jl you can use Julia's package manager
```
using Pkg
Pkg.add("FewBodyToolkit") # or ] add FewBodyToolkit
```

## Usage & Examples
All modules follow a similar usage pattern:

1.Casting a few-body system in physical parameters: (reduced) masses, interactions, angular momentum, ...
```
# 2-body system in 3D, with reduced mass mur and potential as interaction
using FewBodyToolkit
potential(r) = -10/(1+r)
phys_params = make_phys_params2B(;mur=2.0,vint_arr=[potential],dim=3)
```

2.Setting up numerical parameters: number and range of Gaussian basis functions, complex-rotation angle, ...
```
# 10 basis functions with minimum and maximum range parameters 0.5, and 30.0, respectively.
num_params = make_num_params2B(;gem_params=(nmax=10,r1=0.5,rnmax=30.0)) 
```

3.Solving the system using the solver provided by the corresponding module.
```
energies, vectors = GEM2B_solve(phys_params, num_params; wf_bool=1)
```

4.Post-processing via Plots, optional calculation of wave-functions or mean value of observables.
```
using Plots
rgrid = range(0.0,5.0,length=100)
firstexcited = GEM2B.wavefun_arr(rgrid,phys_params,num_params,vectors[:,2])
plot(rgrid,firstexcited)
```


Explicit examples showing this procedure for each of the modules can be found under [Examples](@ref) with fully runnable scripts in the /examples subfolder.

## Method and advanced options

More information on the underlying method and basis functions can be found under [Gaussian Basis Functions](@ref). An explanation of advanced options are listed in [Advanced options](@ref).


## API Reference

See the [API section](@ref API) for a list of all exported types and functions.


## Related packages:
- [Antique.jl](https://github.com/ohno/Antique.jl) - Analytical solutions to solvable quantum mechanical models
- [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl) - Approximating functions and operators using spectral methods
- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) - Comprehensive package for differential equations
- [FewBodyECG.jl](https://github.com/JuliaFewBody/FewBodyECG.jl) - For coulombic few-body systems in 3D, using explicitly correlated Gaussians
- [TwoBody.jl](https://github.com/ohno/TwoBody.jl) - Solutions to two-body systems using various methods, e.g. finite differences