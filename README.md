# FewBodyToolkit.jl
Julia package for solving 2- and 3-body quantum systems in 1D–3D with general potentials


# Example usage
```julia
using Pkg;Pkg.activate(".")
using FewBodyToolkit.ISGL

# Define default physical and numerical parameters:
pp = make_phys_params()
np = make_num_params()

# specify Gaussian potentials:
gaussopt=[[[1,-10.0,1.0]],[[1,-10.0,1.0]],[[1,-10.0,1.0]]]

# Solve a 3-body, 3D quantum system
energies = ISGL_solve(pp,np;gaussopt)

```

# Related packages:
- [Antique.jl](https://github.com/ohno/Antique.jl.git) - Analytical solutions to solvable quantum mechanical models
- [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl.git) - Approximating functions and operators using spectral methods
- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl.git) - Comprehensive package for differential equations
- [FewBodyPhysics.jl](https://github.com/MartinMikkelsen/FewBodyPhysics.jl.git) - For coulombic few-body systems in 3D, using explicitly correlated Gaussians
- [TwoBody.jl](https://github.com/ohno/TwoBody.jl.git) - Solutions to two-body systems using various methods, e.g. finite differences