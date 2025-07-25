# Tests for the module GEM2B (3D scenario)
# Tests via Coulomb potential

## Setup:
# physical parameters:
Z=1.0
v_coulomb(r) = -Z/r
mur = 1.0
phys_params = make_phys_params2B(;vint_arr=[v_coulomb])

# numerical parameters:
gem_params = (;nmax=10,r1=0.1,rnmax=50.0) # gem_params
num_params = make_num_params2B(;gem_params)


## Tests:
# Exact results for the Coulomb potential
exact_results = [-1/(2*i^2) for i=1:4]

# 1. Standard inputs: csm_bool = 0, cr_bool = 0
energies_arr = GEM2B.GEM2B_solve(phys_params,num_params)
@test all(isapprox.(energies_arr[1:4], exact_results; atol=1e-3))

# 2. Test of preallocation
pa = GEM2B.PreallocStruct2B(num_params, 0, 0) #csm_bool = 0, cr_bool = 0
GEM2B.GEM2B_solve!(pa,phys_params,num_params,0,0,0,0,0,0.0)
energies_arr = pa.energies
@test all(isapprox.(energies_arr[1:4], exact_results; atol=1e-3))

# 3. Complex scaling (angle 0°) should have no effect: csm_bool = 1
energies_arr = GEM2B.GEM2B_solve(phys_params,num_params;csm_bool=1)
@test all(isapprox.(real.(energies_arr[1:4]), exact_results; atol=1e-3))

# 4. Finite complex scaling angle (5°) should have very little effect on the bound states
num_paramsC = make_num_params2B(;gem_params,theta_csm=5.0)
energies_arr = GEM2B.GEM2B_solve(phys_params,num_paramsC;csm_bool=1)
@test all(isapprox.(real.(energies_arr[1:4]), exact_results; atol=1e-3))

# 5. Coupled channels:
phys_paramsCC = make_phys_params2B(;vint_arr=[r->0.0])
wfun(r) = v_coulomb(r);wfun2(r) = 0.05*exp(-r^2) 
dfun(r) = 0.0; dfun2(r) = 0.0
WCC = [wfun wfun2; wfun2 wfun]
dor = 1; #derivative-order
DCC = reshape([ [dor, dfun], [dor, dfun2], [dor, dfun2], [dor, dfun] ], 2, 2)
energies_arr = GEM2B.GEM2B_solveCC(phys_paramsCC, num_params, WCC, DCC; diff_bool=0)
# Results of the current code:
exact_resultsCC = [-0.513423475743586, -0.486063239098383, -0.1260365599478753, -0.12340375414310918, -0.05573101564535911, -0.054935509952023065, -0.031248675303256424, -0.0309785242186285]
@test all(isapprox.(energies_arr[1:8], exact_resultsCC; atol=1e-5))
