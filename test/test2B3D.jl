# Tests for the module GEM2B (3D scenario)
# Tests via Coulomb potential

## Setup:
# physical parameters:
Z=1.0
v_coulomb(r) = -Z/r
mur = 1.0
phys_params = make_phys_params2B(;interactions=[v_coulomb])

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
pa = GEM2B.PreallocStruct2B(num_params, false, false) #complex_scaling = false, complex_ranged = false
GEM2B.GEM2B_solve!(pa,phys_params,num_params,false,false,false,false,false,0.0)
energies_arr = pa.energies
@test all(isapprox.(energies_arr[1:4], exact_results; atol=1e-3))

# 3. Complex scaling (angle 0°) should have no effect: complex_scaling = true
energies_arr = GEM2B.GEM2B_solve(phys_params,num_params;complex_scaling=true)
@test all(isapprox.(real.(energies_arr[1:4]), exact_results; atol=1e-3))

# 4. Finite complex scaling angle (5°) should have very little effect on the bound states
num_paramsC = make_num_params2B(;gem_params,complex_scaling_angle=5.0)
energies_arr = GEM2B.GEM2B_solve(phys_params,num_paramsC;complex_scaling=true)
@test all(isapprox.(real.(energies_arr[1:4]), exact_results; atol=1e-3))

# 5. Coupled channels:
phys_paramsCC = make_phys_params2B(;interactions=[r->0.0])
wfun(r) = v_coulomb(r);wfun2(r) = 0.05*exp(-r^2) 
dfun(r) = 0.0; dfun2(r) = 0.0
WCC = [wfun wfun2; wfun2 wfun]
dor = 1; #derivative-order
DCC = reshape([ [dor, dfun], [dor, dfun2], [dor, dfun2], [dor, dfun] ], 2, 2)
energies_arr = GEM2B.GEM2B_solveCC(phys_paramsCC, num_params, WCC, DCC; return_diff=false)
# Results of the current code:
exact_resultsCC = [-0.513423475743586, -0.486063239098383, -0.1260365599478753, -0.12340375414310918, -0.05573101564535911, -0.054935509952023065, -0.031248675303256424, -0.0309785242186285]
@test all(isapprox.(energies_arr[1:8], exact_resultsCC; atol=1e-5))

# 6. Debug return path should return the preallocation struct
debug_out = GEM2B.GEM2B_solve(phys_params, num_params; debug=true)
@test debug_out isa GEM2B.PreallocStruct2B

debug_out_wf = GEM2B.GEM2B_solve(phys_params, num_params; debug=true, return_wavefunctions=true)
@test debug_out_wf isa GEM2B.PreallocStruct2B

# 7. Inverse branch (both energies-only and energies+wavefunctions)
inv_vals = GEM2B.GEM2B_solve(phys_params, num_params; inverse=true, target_energy=-0.5)
@test all(isfinite.(inv_vals[1:4]))

inv_vals_wf, inv_wf = GEM2B.GEM2B_solve(phys_params, num_params; inverse=true, target_energy=-0.5, return_wavefunctions=true)
@test all(isfinite.(inv_vals_wf[1:4]))
@test size(inv_wf, 1) == num_params.gem_params.nmax

# 8. Deprecated keyword aliases should emit warnings and still execute
@test_logs (:warn, r"wf_bool is deprecated") GEM2B.GEM2B_solve(phys_params, num_params; wf_bool=true)
@test_logs (:warn, r"cr_bool is deprecated") GEM2B.GEM2B_solve(phys_params, num_params; cr_bool=false)
@test_logs (:warn, r"csm_bool is deprecated") GEM2B.GEM2B_solve(phys_params, num_params; csm_bool=false)
@test_logs (:warn, r"debug_bool is deprecated") GEM2B.GEM2B_solve(phys_params, num_params; debug_bool=false)
@test_logs (:warn, r"inverse_bool is deprecated") GEM2B.GEM2B_solve(phys_params, num_params; inverse_bool=false)
@test_logs (:warn, r"diff_bool is deprecated") GEM2B.GEM2B_solveCC(phys_paramsCC, num_params, WCC, DCC; diff_bool=false)

# 9. Incompatible coupled-channel options should throw
@test_throws ErrorException GEM2B.GEM2B_solveCC(phys_paramsCC, num_params, WCC, DCC; complex_ranged=true, complex_scaling=true)

# 10. Coupled channels with derivative contribution enabled
energies_arr_diff = GEM2B.GEM2B_solveCC(phys_paramsCC, num_params, WCC, DCC; return_diff=true)
@test length(energies_arr_diff) == length(energies_arr)
@test all(isfinite.(energies_arr_diff[1:8]))

# 11. GEM_Optim_2B: optimize GEM ranges for ground state of Coulomb potential
np_optim = make_num_params2B(;gem_params=(;nmax=8, r1=0.5, rnmax=20.0))
result_optim = GEM_Optim_2B(phys_params, np_optim, 1)
@test length(result_optim) == 3          # [r1_opt, rnmax_opt, energy]
@test isapprox(result_optim[3], -0.5; atol=1e-2)  # hydrogen ground state E = -0.5
