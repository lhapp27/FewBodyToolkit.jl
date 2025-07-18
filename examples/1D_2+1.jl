# # 2+1 system in 1D

# This example reproduces the results in the article [happ2019](@cite). It studies a one-dimensional 2+1 system of two identical particles interacting with a third particle via a central potential. Here, the interaction is taken to be a Gaussian potential, which supports a weakly-bound ground state. The two identical particles do not interact.
 
# ## Setup
using Printf, FewBodyToolkit


# This function returns the universal energy ratios of Table I in the article for a given mass ratio:
function exfun(mr)
    if mr == 2.2
        return [-2.1966, -1.0520]
    elseif mr == 12.4
        return [-2.5963, -1.4818, -1.1970, -1.0377, -1.0002]
    elseif mr == 22.2
        return [-2.7515, -1.6904, -1.3604, -1.1479, -1.0525, -1.0040]
    else
        error("Unknown mass ratio: $mr")
    end
end;


# ## Two-body inverse problem

# First we define the 2+1 system via the mass ratio. Here we choose 22.2 since it results in the most bound states. Feel free to change it to either 2.2 or 12.4.
massratio = 22.2 # other values are 2.2 and 12.4
mass_arr = [1.0, massratio, massratio]
mur = 1/(1/mass_arr[1]+1/mass_arr[2]) # reduced mass
println("Mass ratio: ", massratio, ", reduced mass: ", round(mur,digits=4))

# We need a potential whose ground state is weakly bound. Here, we use a Gaussian potential and set the target energy to -10^-3. Then, we find the required potential strength via the GEM2B code.
v0 = -1.0; mu_g = 1.0;
vg = GaussianPotential(v0,mu_g)

phys_params2B = make_phys_params2B(;mur,vint_arr=[vg],dim=1)
num_params2B = make_num_params2B(;gem_params=(;nmax=16, r1=1.0, rnmax=120.0))

stateindex = 1; target_e2 = -1e-3;
println("1. Two-body inverse problem")
pps,nps,vscale = GEM2B.v0GEMOptim(phys_params2B,num_params2B,stateindex,target_e2)
println("Found potential parameters: v0 = ", vscale*v0, ", mu_g = ", mu_g)

# We define the rescaled potential and corresponding physical parameters
vgscaled = GaussianPotential(v0*vscale,mu_g)
pps = make_phys_params2B(;mur,vint_arr=[vgscaled],dim=1)

#Check if the two-body system indeed has the desired binding energy:
e2s = GEM2B.GEM2B_solve(pps,nps)
println("Two-body binding energy: ", e2s[1], " (target: $(target_e2) )")


# ## Three-body problem

# ### Two identical bosons

# #### Inputs

# Having found the potential parameters, we can now set up the three-body problem with the scaled potential. For bosons we can use their symmetry with the argument `svals=["x","b","b"]`, where `x` is the different particle and `b` denotes two identical bosons.
vint_arr=[[],[vgscaled],[vgscaled]] #[[v23],[v31],[v12]]
phys_params3B = make_phys_params3B1D(;mass_arr=mass_arr,svals=["x","b","b"],vint_arr)

# For the numerical parameters `nmax`, `r1`, and `rnmax` we use the optimized ones, found by the two-body inverse problem. The parameters for the other Jacobi coordinates are set manually.
nmax = nps.gem_params.nmax;
r1 = nps.gem_params.r1; rnmax = nps.gem_params.rnmax;
num_params3B = make_num_params3B1D(;gem_params=(;nmax, r1, rnmax, Nmax=16, R1=1.5, RNmax=250.0))

# #### Solving the three-body problem

println("\n2. Solving the three-body problem and comparing to the article's results")
e3 = GEM3B1D.GEM3B1D_solve(phys_params3B,num_params3B);

# We compute the ratio of three-body to two-body binding energies.
println("Results for the bosonic case, mass ratio = $massratio")
epsilon = e3 /abs(e2s[1])

ex_arr = exfun(massratio)[1:2:end]
comparison(epsilon, ex_arr, min(length(epsilon),length(ex_arr)); s1="Gaussian", s2="Contact")



# ## Two identical fermions

# #### Inputs

# For fermions we can use the same potential. To account for their different statistics and parity, we use `svals=["x","f","f"]` and `parity=-1`. To allow for basis functions that obey these requirements, we need to set `lmax, Lmax=1`.
println("\n3. Results for the fermionic case, mass ratio = $massratio")
phys_params3B_F = make_phys_params3B1D(;mass_arr=mass_arr,svals=["x","f","f"],vint_arr,parity=-1)
num_params3B_F = make_num_params3B1D(;gem_params=(;nmax, r1, rnmax, Nmax=16, R1=1.5, RNmax=250.0), lmin=0, Lmin=0, lmax=1, Lmax=1)
e3_F = GEM3B1D.GEM3B1D_solve(phys_params3B_F,num_params3B_F);

# We compute the ratio of three-body to two-body binding energies.
epsilon_F = e3_F /abs(e2s[1])

ex_arr_F = exfun(massratio)[2:2:end]
comparison(epsilon_F, ex_arr_F, min(length(epsilon_F),length(ex_arr_F)); s1="Gaussian", s2="Contact")

# Overall, we can reproduce the article's results quite well for both bosonic and fermionic systems. Better results could be obtained with more basis functions and/or optimized basis parameters. Note, however, that perfect agreement cannot be reached, since here we used a finite-range interaction whereas the results taken from the article are for a contact interaction. Only in the limit of vanishing two-body binding energy, the two potentials should yield the same results.

# ## Page References

# ```@bibliography
# Pages = ["1D_2+1.md"]
# Canonical = false
# ```

# See also the [full bibliography](@ref References) for further references cited throughout this documentation.