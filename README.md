# FewBodyToolkit.jl
FewBodyToolkit.jl is a Julia package for solving quantum few-body systems. Currently it supports
- two-body systems in 1D,2D,3D
- three-body systems in 1D and 3D.

A more detailed list of features can be found in the [Documentation](https://lhapp27.github.io/FewBodyToolkit.jl/dev/).

Feel free to try it out! If you have questions, encounter issues, or want to extend the package, please open an issue or submit a pull request.

## Installation
To install FewBodyToolkit.jl you can use Julia's package manager
```
using Pkg
Pkg.add("FewBodyToolkit") # or ] add FewBodyToolkit
```

## Example usage
```
using FewBodyToolkit

# Define the interaction and masses of the particles:
v12(r) = -10/(1+r^4)
v23(r) = -8/(1+r^5)
masses = [1.0,1.0,2.0]

# Define physical and numerical parameters
pp = make_phys_params3B3D(;mass_arr = masses, vint_arr=[[v23],[v23],[v12]])
np = make_num_params3B3D(;gem_params=(nmax=10,r1=0.2,rnmax=20.0,Nmax=10,R1=0.2,RNmax=20.0))

# Solve a 3-body, 3D quantum system with your interaction and masses.
@time energies = ISGL_solve(pp,np) # ~0.5s on an average laptop
```

## Method & Documentation
This package implements an expansion into central Gaussian basis functions. For two-body systems in the center-of-mass frame this can be used directly. For three-body systems, a decomposition into Faddeev components is established and subsequently, each component is expanded into a a set of products of two Gaussian basis functions, one for each Jacobi coordinate.

More information on the underlying method and basis functions, as well as examples based on actual research articles can be found in the [Documentation](https://lhapp27.github.io/FewBodyToolkit.jl/dev/).


## Related packages:
- [Antique.jl](https://github.com/ohno/Antique.jl) - Analytical solutions to solvable quantum mechanical models
- [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl) - Approximating functions and operators using spectral methods
- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) - Comprehensive package for differential equations
- [FewBodyECG.jl](https://github.com/JuliaFewBody/FewBodyECG.jl) - For coulombic few-body systems in 3D, using explicitly correlated Gaussians
- [TwoBody.jl](https://github.com/ohno/TwoBody.jl) - Solutions to two-body systems using various methods, e.g. finite differences