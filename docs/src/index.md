# FewBodyToolkit.jl

[![Build Status](https://github.com/lhapp27/FewBodyToolkit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lhapp27/FewBodyToolkit.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://lhapp27.github.io/FewBodyToolkit.jl/dev/)

A Julia package for solving quantum few-body systems of 2 or 3 particles with general potentials in various dimensions.

## Features and modules overview

The FewBodyToolkit.jl package currently provides three modules with the following features:

### `GEM2B` - Two-body solver (1D,2D,3D)
* Two-body systems in 1D, 2D, or 3D
* Symmetric interaction potentials of arbitrary shape
* Real- and complex-ranged Gaussian basis functions
* Complex scaling method (CSM) for resonances
* Basis optimization
* Inverse problem: scaling the interaction to yield a desired eigenenergy
* Coupled-channel problems with additional derivative terms
* Output of the wave function

### `GEM3B1D` - Three-body solver (1D)
* Three-body systems in 1D
* Pair-interactions of arbitrary shape and symmetry (e.g. parity-violating)
* Expansion in up to three Fadeev components (rearrangement channels)
* Automatic handling of identical particles (symmetrization, reduction of Faddeev components)
* Complex scaling method (CSM) for resonances

### `ISGL` - Three-body solver (3D)
* Three-body systems in 3D
* Central-symmetric pair-interactions of arbitrary shape
* Expansion in up to three Fadeev components (rearrangement channels)
* Automatic handling of identical particles (symmetrization, reduction of Faddeev components)
* Complex scaling method (CSM) for resonances
* Arbitrary high intrinsic angular momenta (computationally expensive for high values) via infinitesimally-shifted Gaussian basis functions
* On-the-fly calculation of central observables and mean-square radii


## Installation
To install FewBodyToolkit.jl you can use Julia's package manager
```
using Pkg
Pkg.add("FewBodyToolkit") # or ] add FewBodyToolkit
```

## Usage & Examples
All modules follow a similar usage pattern:

1.Casting a few-body system in physical parameters: (reduced) masses, interactions, parity, etc. :
```
# 2-body system in 3D, with reduced mass mur and potential as interaction
using FewBodyToolkit
potential(r) = -10/(1+r)
phys_params = make_phys_params2B(;mur=2.0,vint_arr=[potential],dim=3)
```

2.Setting up numerical parameters: number and range of Gaussian basis functions, complex-rotation angle, etc. :
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


Explicit examples showing this procedure for each of the modules can be found under the Examples page with fully runnable scripts in the `/examples` subfolder of the repository.

## Method and advanced options

More information on the underlying method and basis functions can be found under [Gaussian Basis Functions](@ref). An explanation of advanced options are listed in [Advanced options](@ref).


## API Reference

See the [API section](@ref API) for a list of all exported types and functions.


## Related packages:
- [Antique.jl](https://github.com/ohno/Antique.jl) - Analytical solutions to solvable quantum mechanical models
- [FewBodyECG.jl](https://github.com/JuliaFewBody/FewBodyECG.jl) - For coulombic few-body systems in 3D, using explicitly correlated Gaussians
- [TwoBody.jl](https://github.com/ohno/TwoBody.jl) - Solutions to two-body systems using various methods, e.g. finite differences
- [QMsolve](https://github.com/quantum-visualizations/qmsolve) - Solving and visualizing the Schrödinger equation in Python
- [MOLSCAT](https://github.com/molscat/molscat) -  Atom-molecule scattering in Fortran
- [JPublicThreeBodySolver] (https://github.com/roudnev/JPublicThreeBodySolver) - Solver for Faddeev equations in Java