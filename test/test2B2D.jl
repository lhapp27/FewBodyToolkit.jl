# Tests for the module GEM2B (2D scenario)
# Tests via harmonic oscillator potential.

## Setup:
# physical parameters:
mass_arr=[1.0,10.0] # finite masses of the two particles
mur = 1/(1/mass_arr[1]+1/mass_arr[2]) # reduced mass
omega = 0.5
v_ho(r) = 0.5*mur*omega^2*r^2
phys_params = make_phys_params2B(;mur,interactions=[v_ho],dim=2)

# numerical parameters:
gem_params = (;nmax=14,r1=0.82,rnmax=10.62) # gem_params
num_params = make_num_params2B(;gem_params)


## Tests:
# Exact results for the harmoic oscillator potential
exact_results = ([2*i for i=0:15] .+ 1) .*omega

# 1. Standard inputs: csm_bool = 0, cr_bool = 0
energies_arr = GEM2B.GEM2B_solve(phys_params,num_params)
@test all(isapprox.(energies_arr[1:4], exact_results[1:4]; atol=1e-3))

# 2. Using complex-ranged basis functions: complex_scaling = false, complex_ranged = true
num_paramsCR = make_num_params2B(;gem_params=(nmax=7,r1=1.3398861224184124,rnmax=4.4150781608810705))
energies_arr = GEM2B.GEM2B_solve(phys_params,num_paramsCR;complex_ranged=true)
@test all(isapprox.(energies_arr[1:10], exact_results[1:10]; atol=1e-3))

# 3. Complex scaling (angle 0°) should have no effect: complex_scaling = true
energies_arr = GEM2B.GEM2B_solve(phys_params,num_params;complex_scaling=true)
@test all(isapprox.(real.(energies_arr[1:4]), exact_results[1:4]; atol=1e-3))

# 4. Finite complex scaling angle (1°) should have very little effect on the bound states
num_paramsC = make_num_params2B(;gem_params,complex_scaling_angle=1.0)
energies_arr = GEM2B.GEM2B_solve(phys_params,num_paramsC;complex_scaling=true)
@test all(isapprox.(real.(energies_arr[1:4]), exact_results[1:4]; atol=1e-3))
