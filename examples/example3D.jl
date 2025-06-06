# # 3D Example: Two particles with Coulomb interaction
#
# This example demonstrates how to use the `FewBodyToolkit.jl` package to compute bound states for two particles in 3D. Here we use the Coulomb interaction, since it has analytic solutions. In relative coordinates, this system is equivalent to a single particle in a potential. It is governed by the following Schrödinger equation:
# \\[ -\frac{1}{2} \frac{1}{r}\frac{d^2}{dr^2}\left( r\psi \right) + V(r)\psi = E\psi \\]
# with the Colomb potential
# \\[ V(r) = -\frac{Z}{r}. \\]

# ## Setup

using Printf, Interpolations, Plots, Antique, FewBodyToolkit#.GEM2B

# ## Input parameters

# #### Physical parameters

mass_arr = [1.0, Inf] # array of masses of particles [m1,m2]
mur = 1 / (1/mass_arr[1] + 1/mass_arr[2])
Z = 1.0

function v_coulomb(r)
    return -Z/r
end

# We define the physical parameters as a `NamedTuple` which carries the information about the Hamiltonian.
phys_params = make_phys_params2B(;vint_arr=[v_coulomb])
# By leaving out the optional parameters, we use the defaults:
# - `mur = 1.0`: reduced mass
# - `dim = 3`: dimension of the problem
# - `lmin = lmax = 0`: minimum and maximum angular momentum (in 1D this corresponds to even states)
# - `hbar = 1.0`: when working in dimensionless units

# #### Numerical parameters

nmax=10 # number of Gaussian basis functions
r1=0.1;rnmax=30.0; # r1 and rnmax defining the widths of the basis functions
gem_params = (;nmax,r1,rnmax);

# We define the numerical parameters as a `NamedTuple`:
num_params = make_num_params2B(;gem_params)



# ## 1. Numerical solution

# We solve the two-body system by calling `GEM2B_solve`.
energies = GEM2B.GEM2B_solve(phys_params,num_params)

# Number of bound states to consider:
simax = min(lastindex(energies),6); # max state index

# The Coulomb potential has infinitely many bound states, whose energies can be found exaclty. We can use the package Antique.jl to provide these energies.
CTB = Antique.CoulombTwoBody(m₁=mass_arr[1], m₂=mass_arr[2])
energies_exact = [Antique.E(CTB,n=i) for i=1:40]

println("1. Numerical solution of the 3D problem:")
comparison(energies,energies_exact,simax)

# The numerical solutions are good only for the few lowest state. Also, we only find 5 bound states.


# ## 2. Optimization of basis parameters

# We can optimize the basis parameters for a specific state indicated by `stateindex` using `GEM_Optim_2B`.

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

# Optimizing the parameters for the 6th excited states finds more bound states, while loosing some accuracy for the lower states. Here, only more basis functions would help.
comparison(energies_opt,energies_exact,simax; s1="Optimized")


# ## 3. Example with many basis functions

# Highly accurate results can indeed be obtained by using a larger basis. For a two-body system this comes only at a moderate computational cost. Here, we reproduce the table 2 of Ref. [hiyama2003](@cite) with the following numerical parameters:
println("\n3. Highly accurate solution using many basis functions:")
np = make_num_params2B(;gem_params=(;nmax=80,r1=0.015,rnmax=2000.0),omega_cr=1.5,threshold=10^-11)
@time energies_accurate = GEM2B.GEM2B_solve(phys_params,np;cr_bool=1) # ~3s on an average laptop
nlist=[1,2,3,4,5,10,14,18,22,26,30,32,34,36,38,40]; # states shown in the article
comparison(energies_accurate,energies_exact,lastindex(nlist); s1="Numerical", s2="Exact",indexlist=nlist)


# ## 4. Using an interpolated interaction
# We can also create a potential from interpolated data. Since the Coulomb potential diverges at the origin, a relatively fine grid is required.

r_arr = 0.001:0.01:50.001
v_arr = v_coulomb.(r_arr)
v_interpol = cubic_spline_interpolation(r_arr,v_arr,extrapolation_bc=Line())
v_int(r) = v_interpol(r); # we have to transform the interaction to an object of type "function"

# As input to the solver we need to define new physical parameters with the interpolated interaction. Moreover, we use the optimized numerical parameters from the previous step.
phys_params_i = make_phys_params2B(;mur,vint_arr=[v_int],dim=3)

println("\n4. Numerical solution using an interpolated interaction:")
energies_interpol = GEM2B.GEM2B_solve(phys_params_i,num_params_opt)
comparison(energies_interpol, energies_opt, simax;s1="Interpolated", s2="Optimized")


# ## 5. Coupled-channel problem

# The package also supports coupled-channel problems via `GEM2B_solveCC`. In this case the interaction is not provided via phys_params
phys_paramsCC = make_phys_params2B(;vint_arr=[r->0.0])

# but via extra arguments WCC: `wfun` on the diagonal; `wfun2` for off-diagonal couplings, and DCC: derivative terms of order `dor`, and radial prefactors `dfun` (diagonal) and `dfun2` (off-diagonal).

wfun(r) = v_coulomb(r);wfun2(r) = 0.05*exp(-r^2) 
dfun(r) = 0.0; dfun2(r) = 0.0
WCC = [wfun wfun2; wfun2 wfun]
dor = 1; #derivative-order
DCC = reshape([ [dor, dfun], [dor, dfun2], [dor, dfun2], [dor, dfun] ], 2, 2)

# As a test-case we use the coulomb interaction on the diagonal, and a weak repulsion on the off-diagonal. We don't consider any extra derivatives. Since the coupling is weak, we get approximately twofold degenerate eigenvalues, splitted around the original Coulomb results.
energiesCC = GEM2B.GEM2B_solveCC(phys_paramsCC, num_params_opt, WCC, DCC; diff_bool=0)

energies_exactCC = repeat(energies_exact, inner=(2,))

println("\n5. Coupled channel calculation:")
comparison(energiesCC, energies_exactCC, simax; s1="Coupled-Channel")




# ## 6. Calculation of the wave function

# Adding the optional argument `wf_bool=1` to `GEM2B_solve` also computes and returns a matrix of eigenvectors (in each column). These eigenvectors contain the weights of the basis functions.

energiesw,wfs = GEM2B.GEM2B_solve(phys_params,num_params_opt;wf_bool=1,cr_bool=0);

# We can use the functions `GEM2B.wavefun_arr` and `GEM2B.wavefun_point` to compute the wave function at a set of points or at a specific point, respectively. The information on the basis functions is provided via the optimized numerical parameters `num_params_opt`, the vector of Gaussian widths. We compare the numerical wave function with the exact solutions provided by Antique.jl and find very good agreement.

dr = 0.1
r_arr = 0.0:dr:50.0
redind = vcat(1:2:30,31:5:50,51:10:lastindex(r_arr)) # We evaluate the analytical solutions at a coarser grid to avoid overloading the plot.

wfA(r,n) = Antique.R(CTB, r; n, l=0) # Exact wave function for the n-th state

p = plot(xlabel="\$ r \$", ylabel="\$ r^2\\,|\\psi(r)|^2 \$", title="Two-body radial s-wave densities\n for a 3D Coulomb system", guidefont=18,legendfont=10)
density = zeros(length(r_arr),4)
density_exact = zeros(length(r_arr),4)
markers = [:circ, :square, :utriangle, :star]
for si = 1:4
    wf = wfs[:,si]
    psi_arr = GEM2B.wavefun_arr(r_arr,phys_params,num_params_opt,wf;cr_bool=0)
    
    density[:,si] .= abs2.(psi_arr).*r_arr.^2
    density_exact[:,si] = abs2.(wfA.(r_arr,si)).*r_arr.^2

    scatter!(r_arr[redind],density_exact[redind,si],label="n=$(si), exact", lw=2,color=si, marker=markers[si],markersize=3)
    plot!(r_arr, density[:,si], label="n=$(si), num", lw=2, color=si)
end

p

# The normalization of the wave function can be checked by integrating the density:
norms = density[:,1:4]'*fill(dr,lastindex(r_arr)) # A simple Riemann sum is sufficient here
println("\n6. Norms of the wave functions:")
comparison(norms, ones(4), 4; s1="Norm", s2="Exact")

