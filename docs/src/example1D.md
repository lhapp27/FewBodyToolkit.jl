```@meta
EditURL = "../../examples/example1D.jl"
```

# 1D Example: Two particles with Pöschl–Teller interaction

This example demonstrates how to use the `FewBodyToolkit.jl` package to compute bound states for two particles in 1D. Here we use the Pöschl–Teller interaction, since it has analytic solutions. In relative coordinates, this system is equivalent to a single particle in a potential. It is governed by the following Schrödinger equation (``\hbar, \mu=1``)
\\[ -\frac{1}{2} \frac{d^2}{dr^2}\psi + V(r)\psi = E\psi \\]
with the Pöschl–Teller potential
\\[ V(r) = -\frac{\lambda(\lambda+1)}{2} \frac{1}{\cosh^2(r)}. \\]

## Setup

````@example example1D
using Printf, Interpolations, Antique, FewBodyToolkit
````

## Input parameters

#### Physical parameters

````@example example1D
mass_arr=[1.0,Inf] # this ensures a reduced mass of 1.0
mur = 1/(1/mass_arr[1]+1/mass_arr[2]) # reduced mass
lambda=8.0

function v_poschl(r)
    return -lambda*(lambda+1)/2/mur*1/cosh(r)^2
end;
nothing #hide
````

We define the physical parameters as a `NamedTuple` which carries the information about the Hamiltonian.

````@example example1D
phys_params = make_phys_params2B(;mur,vint_arr=[v_poschl],dim=1)
````

By leaving out the optional parameters, we use the defaults:
- `lmin = lmax = 0`: minimum and maximum angular momentum (in 1D this corresponds to even states)
- `hbar = 1.0`: when working in dimensionless units

#### Numerical parameters

````@example example1D
nmax=6 # number of Gaussian basis functions
r1=0.1;rnmax=10.0; # r1 and rnmax defining the widths of the basis functions
gem_params = (;nmax,r1,rnmax);
nothing #hide
````

We define the numerical parameters as a `NamedTuple`:

````@example example1D
num_params = make_num_params2B(;gem_params)
````

## 1. Numerical solution

We solve the two-body system by calling `GEM2B_solve`.

````@example example1D
energies_arr = GEM2B_solve(phys_params,num_params)
````

Determine the number of bound states

````@example example1D
simax = findlast(energies_arr.<0)
````

The Pöschl–Teller potential has `lambda = 8` eigenvalues. In this example we focus on the even states, hence there are only 4 bound states. Their energies can be found exactly:
\\[ E_n = -\frac{(\lambda-n)^2}{2\mu} \\]
where `` n = 1, 2, ... , \lambda-1 `` is the state index. The package [Antique.jl](https://github.com/ohno/Antique.jl) provides these exact energies in a convenient way.

````@example example1D
PT = Antique.PoschlTeller(λ=8)
ex_arr = [Antique.E(PT,n=i) for i=0:2:Int(floor(lambda-1))]

println("1. Numerical solution of the 1D problem:")
comparison(energies_arr, ex_arr, simax)
````

So far, the numerical solutions are not very accurate. This is because the basis parmeters are not optimal.

## 2. Optimization of basis parameters

We can optimize the basis parameters for a specific state indicated by `stateindex` using `GEM_Optim_2B`.

````@example example1D
stateindex = 4 # which state to optimize for
params_opt = GEM2B.GEM_Optim_2B(phys_params, num_params, stateindex)
gem_params_opt = (;nmax, r1 = params_opt[1], rnmax = params_opt[2])
num_params_opt = make_num_params2B(;gem_params=gem_params_opt)
e2_opt = GEM2B.GEM2B_solve(phys_params,num_params_opt)

println("\n2. Optimization of GEM parameters for E2[$stateindex]:")
@printf("%-15s %-15s %-15s %-15s %-15s\n", "r1", "rnmax", "E2[$(stateindex-1)]", "E2[$stateindex]", "E2[$(stateindex+1)]")

println("Before optimization:")
@printf("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", gem_params.r1, gem_params.rnmax, energies_arr[stateindex-1], energies_arr[stateindex], energies_arr[stateindex+1])

println("After optimization:")
@printf("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", gem_params_opt.r1, gem_params_opt.rnmax, e2_opt[stateindex-1], e2_opt[stateindex], e2_opt[stateindex+1])
````

With the optimized parameters, the exact energies are reproduced very well, using only 6 basis functions.

````@example example1D
comparison(e2_opt, ex_arr, simax; s1="Optimized")
````

## 3. Using an interpolated interaction
We can also create a potential from interpolated data.

````@example example1D
r_arr = -10.0:0.5:10.0
v_arr = v_poschl.(r_arr)
v_interpol = cubic_spline_interpolation(r_arr,v_arr,extrapolation_bc=Line())
v_int(r) = v_interpol(r); # we have to transform the interaction to an object of type "function"
nothing #hide
````

As input to the solver we need to define new physical parameters with the interpolated interaction. Moreover, we use the optimized numerical parameters from the previous step.

````@example example1D
phys_params_interpol = make_phys_params2B(;mur,vint_arr=[v_int],dim=1)

println("\n3. Numerical solution using an interpolated interaction:")
energies_interpol = GEM2B_solve(phys_params_interpol,num_params_opt)
comparison(energies_interpol, e2_opt, simax;s1="Interpolated", s2="Optimized")
````

## 4. Inverse problem: Tuning the potential strength

We can use `v0GEMOptim` to scale the interaction such that the state indicated by `stateindex` has a fixed energy `target_e2`. At the same time, the basis parameters are optimized for this state.

````@example example1D
stateindex = 3; target_e2 = -18.0;
println("\n4. Scaling the potential such that E2[$stateindex] = $target_e2:")
phys_params_scaled,num_params_scaled,vscale = GEM2B.v0GEMOptim(phys_params,num_params,stateindex,target_e2)
e2_v0 = GEM2B.GEM2B_solve(phys_params_scaled,num_params_scaled)

println("After scaling:")
@printf("%-15s %-15s %-15s %-15s %-15s\n", "r1", "rnmax", "E2[$(stateindex-1)]", "E2[$stateindex]", "E2[$(stateindex+1)]")
@printf("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", num_params_scaled.gem_params.r1, num_params_scaled.gem_params.rnmax, e2_v0[stateindex-1], e2_v0[stateindex], e2_v0[stateindex+1])
````

Here, we scale the potential such that the energy of the state with `stateindex = 3` is equal to `target_e2 = -18.0`. So far this was the energy of the state with index 2. For this special potential, this corresponds therefore to increasing the number of states and `lambda` by 2. Hence, the we expect the scaling factor to be approximately ``(\lambda+2)(\lambda+2+1)/(\lambda(\lambda+1))``

````@example example1D
println("vscale = $(round(vscale,digits=8)) should be approximately (λ+2)*(λ+2+1)/(λ*(λ+1)) = ", round((lambda+2)*(lambda+2+1)/(lambda*(lambda+1)),digits=8) )
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

