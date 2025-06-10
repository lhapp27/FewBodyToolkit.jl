# # 2D Example: Two particles with Harmonic oscillator interaction
#
# This example demonstrates how to use the `FewBodyToolkit.jl` package to compute bound states for two particles in 2D. Here we use the harmonic oscillator, since it has analytic solutions. In relative coordinates, this system is equivalent to a single particle in a potential. It is governed by the following Schrödinger equation (we set the magnetic quantum number ``m = 0``, and ``\hbar = 1``)
# \\[ -\frac{1}{2 \mu} \left[\frac{d^2}{dr^2} + \frac{1}{r} \frac{d}{dr} \right] \psi + V(r)\psi = E\psi \\]
# with the Harmonic oscillator potential
# \\[ V(r) = -\frac{1}{2} \mu \omega^2 r^2. \\]


# ## Setup

using Printf, Interpolations, FewBodyToolkit

# ## Input parameters

# #### Physical parameters

mass_arr=[1.0,10.0] # finite masses of the two particles
mur = 1/(1/mass_arr[1]+1/mass_arr[2]) # reduced mass
omega = 0.5

function v_ho(r)
    return 0.5*mur*omega^2*r^2
end;

# We define the physical parameters as a `NamedTuple` which carries the information about the Hamiltonian.
phys_params = make_phys_params2B(;mur,vint_arr=[v_ho],dim=2)
# By leaving out the optional parameters, we use the defaults:
# - `lmin = lmax = 0`: minimum and maximum angular momentum
# - `hbar = 1.0`: when working in dimensionless units

# #### Numerical parameters

nmax=14 # number of Gaussian basis functions
r1=0.5;rnmax=10.0;
gem_params = (;nmax,r1,rnmax);

# We define the numerical parameters as a `NamedTuple`:
num_params = make_num_params2B(;gem_params,threshold=10^-8)



# ## 1. Numerical solution

# We solve the two-body system by calling `GEM2B_solve`.
energies = GEM2B.GEM2B_solve(phys_params,num_params)


# The Harmonic Oscillator has infinitely many eigenvalues. For the radially symmetric states (m = 0) their energies are given by
# \\[ E_n = \omega (2n+1) \\]

energies_exact = ([2*i for i=0:15] .+ 1) .*omega

println("1. Numerical solution of the 2D problem:")
simax=10;
comparison(energies, energies_exact, simax)

# Especially for the first few states, the numerical energies are already reasonably close to the exact energies.

# ## 2. Optimization of basis parameters

# Still, we can try to improve the accuracy by optimizing the basis parameters for a given state indicated by `stateindex`.

stateindex = 6
params_opt = GEM2B.GEM_Optim_2B(phys_params, num_params, stateindex)
gem_params_opt = (;nmax, r1 = params_opt[1], rnmax = params_opt[2])
num_params_opt = make_num_params2B(;gem_params=gem_params_opt)
energies_opt = GEM2B.GEM2B_solve(phys_params,num_params_opt)

println("\n2. Optimization of GEM parameters for E2[$stateindex]:")
@printf("%-15s %-15s %-15s %-15s %-15s\n", "r1", "rnmax", "E2[$(stateindex-1)]", "E2[$stateindex]", "E2[$(stateindex+1)]")

println("Before optimization:")
@printf("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", gem_params.r1, gem_params.rnmax, energies[stateindex-1], energies[stateindex], energies[stateindex+1])

println("after optimization:")
@printf("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", params_opt[1], params_opt[2], energies_opt[stateindex-1], energies_opt[stateindex], energies_opt[stateindex+1])

# The optimized parameters yield not much of an iprovement, since the basis parameters were already quite good. We can improve the results by either using more basis functions, or by complex-ranged basis functions (see below)
comparison(energies_opt,energies_exact,simax; s1="Optimized")


# ## 3. Inverse problem: Tuning the potential strength

# We can use `v0GEMOptim` to scale the interaction such that the state indicated by `stateindex` has a fixed energy `target_e2`. At the same time, the basis parameters are optimized for this state.

stateindex = 2; target_e2 = 3.0;
println("\n3. Scaling the potential such that E2[$stateindex] = $target_e2:")
phys_params_scaled,num_params_scaled,vscale = GEM2B.v0GEMOptim(phys_params,num_params_opt,stateindex,target_e2)
energies_v0 = GEM2B.GEM2B_solve(phys_params_scaled,num_params_scaled)

@printf("%-15s %-15s %-15s %-15s %-15s\n", "r1", "rnmax", "E2[$(stateindex-1)]", "E2[$stateindex]", "E2[$(stateindex+1)]")
@printf("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", num_params_scaled.gem_params.r1, num_params_scaled.gem_params.rnmax, energies_v0[stateindex-1], energies_v0[stateindex], energies_v0[stateindex+1])

# Here, we scale the potential such that the energy of the state with `stateindex = 2` is equal to `target_e2 = 3.0`, i.e. twice its original value. Since the potential scales quadratically  with ``\omega``, we expect a scaling factor of 4.0.
println("vscale = $(round(vscale,digits=8)) should be approximately 4.0")


# ## 4. Using complex-ranged basis functions
# We can also use complex-ranged basis functions, which are useful for more oscillatory bound states, i.e. highly excited states. Note that `cr_bool=1` effectively employs twice the number of basis functions, hence for a fair comparison we choose `nmaxC = nmax/2 = 7`. Keep in mind that the optimal parameters for the complex-ranged basis functions are usually differnt. Hence, we optimize them separately.
nmaxC = 7
r1C = 0.5; rnmaxC = 10.0;
gem_paramsC = (;nmax=nmaxC,r1=r1C,rnmax=rnmaxC);
num_paramsC = make_num_params2B(;gem_params=gem_paramsC)

stateindex = 6
params_opt = GEM2B.GEM_Optim_2B(phys_params, num_paramsC, stateindex; cr_bool = 1)
gem_params_optC = (;nmax = nmaxC, r1 = params_opt[1], rnmax = params_opt[2])
num_params_optC = make_num_params2B(;gem_params=gem_params_optC)
energies_optC= GEM2B.GEM2B_solve(phys_params,num_params_optC; cr_bool = 1)

println("\n4. Using complex-ranged basis functions:")
comparison(energies_optC, energies_exact, 14; s1="Complex-ranged", s2="Exact")

# Using effectively 14 basis functions, we can reproduce the exact energies for the first 10 states with less than `1e-3` deviation, and even up to 14 states with reasonable accuracy.