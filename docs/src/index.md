# FewBodyToolkit.jl

A Julia package for solving 2- and 3-body quantum systems in 1D, 2D, and 3D.


## Features

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
* Symmetric pairwise two-body interactions (no three-body forces)
* Expansion in up to three Fadeev components (rearrangement channels)
* Product basis in Jacobi coordinates
* Automatic (anti-)symmetrization for identical bosons or fermions
* Complex scaling method (CSM) for resonances
* No separate output of the wave function

### `ISGL` - Three-body solver (3D)
* Supports 3D three-body problems
* Central two-body pair-interactions, no three-body forces
* Employs infinitesimally-shifted Gaussian basis functions to handle arbitrary high intrinsic angular momenta (computationally expensive for high values)
* Expansion in up to three Fadeev components (rearrangement channels)
* Product basis in Jacobi coordinates
* Automatic (anti-) symmetrization for identical bosons or fermions
* Complex scaling method (CSM) for resonances
* On-the-fly calculation of central observables and mean-square radii
* No separate output of the wave function


## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/lhapp27/FewBodyToolkit.jl")
```

## Quick Start

```julia
using FewBodyToolkit

# 1) Build physical and numerical parameters for a 2-body problem: 3D Coulomb problem
phys = make_phys_params2B(vint_arr=[r->-1/r], dim=3)
num  = make_num_params2B(gem_params=(nmax=10, r1=0.3, rnmax=20.0))

# 2) Solve generalized eigenvalue problem:
E = GEM_solve(phys,num)

# 3) Print 4 lowest eigenenergies
println("E_1 = ", E[1]) # should be close to -0.50000
println("E_2 = ", E[2]) # should be close to -0.12500
println("E_3 = ", E[3]) # should be close to -0.05555
println("E_4 = ", E[4]) # should be close to -0.03125
```

## Modules

- [`GEM2B`](@ref) — two-body systems in 1D, 2D, and 3D  
- [`GEM3B1D`](@ref) — three-body systems in 1D  
- [`ISGL`](@ref) — three-body systems in 3D  


## Examples
Check out the [1D example](example1D.md) to for an example in 1D.

## API Reference

See the [full API documentation](@ref API) for a list of all exported types and functions.


### `GEM2B` Module

The `GEM2B` module provides tools for solving two-body problems using the Gaussian Expansion Method (GEM). The following functions and types are exported:

#### Parameter Setup

- [`make_phys_params2B`](): Create a named tuple with the physical parameters of the two-body system.
- [`make_num_params2B`](): Set up numerical parameters such as the number of basis functions.
- [`PreallocStruct2B`](): Preallocates and stores relevant arrays for two-body calculations as well as the outputs (energies, coefficient-eigenvector)

#### Solvers and Optimization

- [`GEM_solve`](): Solves the two-body Schrödinger equation for given physical and numerical parameters
- [`v0GEMOptim`](): Optimizes the interaction strength to match a target eigenenergy.
- [`GEM_Optim_2B`](): Optimizes Gaussian basis widths for improved convergence.

#### Wave Function Utilities

- [`wavefun_arr`]: Returns the full wave function as an array sampled on a grid.
- [`wavefun_point`]: Evaluates the wave function at a given point.
