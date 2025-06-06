# Tests for the module GEM2B (1D scenario)
# Tests via Poschl-Teller potential.

## Setup:
# physical parameters:
lambda=8.0
v_poschl(r) = -lambda*(lambda+1)/2/mur*1/cosh(r)^2
mur = 1.0
phys_params = make_phys_params2B(;vint_arr=[v_poschl],dim=1)

# numerical parameters:
gem_params = (;nmax=10,r1=0.3,rnmax=3.0) # gem_params
num_params = make_num_params2B(;gem_params)


## Tests:
# Exact results for the Pöschl-Teller potential
exact_results = [-(lambda-i)^2/2/mur for i=0:2:Int(floor(lambda-1))] #

# 1. Standard inputs: csm_bool = 0, cr_bool = 0
energies_arr = GEM2B.GEM2B_solve(phys_params,num_params)
@test all(isapprox.(energies_arr[1:4], exact_results; atol=1e-3))

# 2. Using complex-ranged basis functions: csm_bool = 0, cr_bool = 1
energies_arr = GEM2B.GEM2B_solve(phys_params,num_params;cr_bool=1)
@test all(isapprox.(energies_arr[1:4], exact_results; atol=1e-3))

# 3. Complex scaling (angle 0°) should have no effect: csm_bool = 1
energies_arr = GEM2B.GEM2B_solve(phys_params,num_params;csm_bool=1)
@test all(isapprox.(real.(energies_arr[1:4]), exact_results; atol=1e-3))

# 4. Finite complex scaling angle (5°) should have very little effect on the bound states
num_paramsC = make_num_params2B(;gem_params,theta_csm=5.0)
energies_arr = GEM2B.GEM2B_solve(phys_params,num_paramsC;csm_bool=1)
@test all(isapprox.(real.(energies_arr[1:4]), exact_results; atol=1e-3))

# 5. Finite complex scaling angle (5°) , with simultaneous complex-ranged basis functions
num_paramsC = make_num_params2B(;gem_params,theta_csm=5.0)
energies_arr = GEM2B.GEM2B_solve(phys_params,num_paramsC;csm_bool=1,cr_bool=1)
@test all(isapprox.(real.(energies_arr[1:4]), exact_results; atol=1e-3))
