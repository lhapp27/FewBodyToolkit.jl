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
GEM2B.GEM2B_solve!(pa,phys_params,num_params,false,false,false)
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

# 7. Matrix export and inverse problem
# T is the kinetic energy alone, so T+V reproduces the forward energies
Tm,Vm,Sm = GEM2B_matrices(phys_params, num_params)
e_matrices = zeros(size(Tm,1))
FewBodyToolkit.eigen2step(e_matrices, Tm .+ Vm, Sm; threshold=num_params.threshold)
@test all(isapprox.(e_matrices[1:4], exact_results; atol=1e-3))

# reduce_basis: L orthonormalizes the basis
L = reduce_basis(Sm; threshold=num_params.threshold)
@test L' * Sm * L ≈ I

# inverse problem: the strengths are exact, scaling the interaction by inv_vals[1] puts the ground state at target_energy
inv_vals = inverse_solve(Tm, Vm, Sm; target_energy=-0.5, threshold=num_params.threshold)
@test isapprox(GEM2B.GEM2B_solve(make_phys_params2B(;interactions=[r -> inv_vals[1]*v_coulomb(r)]), num_params)[1], -0.5; atol=1e-10)

inv_vals_vec, inv_vecs = inverse_solve(Tm, Vm, Sm; target_energy=-0.5, threshold=num_params.threshold, return_vectors=true)
@test all(isapprox.(inv_vals_vec, inv_vals; rtol=1e-10)) # the vector path uses a different LAPACK driver
@test size(inv_vecs, 1) == num_params.gem_params.nmax
@test inv_vecs[:,1]' * Sm * inv_vecs[:,1] ≈ 1 # normalized as x'*S*x == 1

# purely attractive V: the S-reduced pencil reproduces the old -V-metric implementation
old_vals = zeros(size(Tm,1))
FewBodyToolkit.eigen2step(old_vals, Tm .+ 0.5.*Sm, -Vm; threshold=num_params.threshold)
@test all(isapprox.(inv_vals[1:4], old_vals[1:4]; rtol=1e-6))

# target_energy above the lowest eigenvalue of T makes T - target_energy*S indefinite
@test_throws ErrorException inverse_solve(Tm, Vm, Sm; target_energy=1.0, threshold=num_params.threshold)

# complex scaling leaves T and V complex-symmetric instead of hermitian and is rejected
Tc,Vc,Sc = GEM2B_matrices(phys_params, num_paramsC; complex_scaling=true)
@test_throws ErrorException inverse_solve(Tc, Vc, Sc; target_energy=-0.5, threshold=num_params.threshold)

# 8. Deprecated keyword aliases should emit warnings and still execute
@test_logs (:warn, r"wf_bool is deprecated") GEM2B.GEM2B_solve(phys_params, num_params; wf_bool=true)
@test_logs (:warn, r"cr_bool is deprecated") GEM2B.GEM2B_solve(phys_params, num_params; cr_bool=false)
@test_logs (:warn, r"csm_bool is deprecated") GEM2B.GEM2B_solve(phys_params, num_params; csm_bool=false)
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

# 12. PowerLawPotential (3D): analytic treatment of V(r) = v0*|r|^p
# 12a. Coulomb (p=-1) against the exact results and against the numerical path
v_coulomb_pow = PowerLawPotential(-Z,-1.0)
pp_pow = make_phys_params2B(;interactions=[v_coulomb_pow])
e_pow = GEM2B.GEM2B_solve(pp_pow,num_params)
@test all(isapprox.(e_pow[1:4], exact_results; atol=1e-3))
@test all(isapprox.(e_pow[1:4], GEM2B.GEM2B_solve(phys_params,num_params)[1:4]; rtol=1e-8))

# 12b. harmonic oscillator (p=2) against the exact 3D spectrum E=(2n+l+3/2)*omega
omega_ho = 0.7
v_ho_pow  = PowerLawPotential(0.5*1.0*omega_ho^2, 2.0)
v_ho_cent(r) = 0.5*1.0*omega_ho^2*r^2
gp_ho = (;nmax=24,r1=0.2,rnmax=12.0)
np_ho = make_num_params2B(;gem_params=gp_ho)
for l in [0,1,2]
    pp_ho_pow  = make_phys_params2B(;interactions=[v_ho_pow], lmax=l, lmin=l)
    pp_ho_cent = make_phys_params2B(;interactions=[v_ho_cent], lmax=l, lmin=l)
    e_ho_pow  = GEM2B.GEM2B_solve(pp_ho_pow, np_ho)
    e_ho_cent = GEM2B.GEM2B_solve(pp_ho_cent,np_ho)
    exact_ho = [(2*n+l+1.5)*omega_ho for n=0:3]
    @test all(isapprox.(e_ho_pow[1:4], exact_ho; atol=1e-2))
    @test all(isapprox.(e_ho_pow[1:4], e_ho_cent[1:4]; rtol=1e-8))
end

# 12c. non-integer exponent
v_pl_pow = PowerLawPotential(0.9,1.5)
v_pl_cent(r) = 0.9*abs(r)^1.5
@test all(isapprox.(GEM2B.GEM2B_solve(make_phys_params2B(;interactions=[v_pl_pow]),num_params)[1:4],
                    GEM2B.GEM2B_solve(make_phys_params2B(;interactions=[v_pl_cent]),num_params)[1:4]; rtol=1e-8))

# 12d. complex scaling at finite angle: analytic vs numerical
np_csm = make_num_params2B(;gem_params,theta_csm=8.0)
@test all(isapprox.(GEM2B.GEM2B_solve(pp_pow,np_csm;complex_scaling=true)[1:4],
                    GEM2B.GEM2B_solve(phys_params,np_csm;complex_scaling=true)[1:4]; atol=1e-6))

# 12e. validity check: in 3D with lmax=0 the bound is p > -3
@test_throws ErrorException GEM2B.GEM2B_solve(make_phys_params2B(;interactions=[PowerLawPotential(1.0,-3.0)]),num_params)
@test_throws ErrorException GEM2B.GEM2B_solve(make_phys_params2B(;interactions=[PowerLawPotential(1.0,-4.0)]),num_params)
# but for lmax=1 the bound is p > -5, so p=-4 is fine there
@test all(isfinite.(GEM2B.GEM2B_solve(make_phys_params2B(;interactions=[PowerLawPotential(1.0,-4.0)],lmax=1,lmin=1),num_params)[1:2]))
