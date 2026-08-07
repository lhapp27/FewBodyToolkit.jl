# Tests for the module GEM3B1D (three-body, 1D)
# Tests via Gaussian potential

## Setup:
# physical parameters:
v0 = -25.0; mu_g = 1.0
v_gauss = GaussianPotential(v0,mu_g)
v_gauss2(r) = v0*exp(-mu_g*r^2)
phys_params = make_phys_params3B1D(;species=[:b,:b,:b], interactions=[[v_gauss],[v_gauss],[v_gauss]])
phys_paramsN = make_phys_params3B1D(;species=[:b,:b,:b], interactions=[[v_gauss2],[v_gauss2],[v_gauss2]])

vg(r) = -5*exp(-r^2)
vga = GaussianPotential(-5.0, 1.0)

# numerical parameters:
gp = (;nmax=15,Nmax=15,r1=0.1,rnmax=8.0,R1=0.1,RNmax=7.0)
num_params = make_num_params3B1D(;gem_params=gp,theta_csm=5.0)


## Tests:
# 1. There are no exact results, but we can check consistency between analytical treatment (GaussianPotential) and numerical treatment (v_gauss2)
energies_arrA = GEM3B1D_solve(phys_params,num_params)
energies_arrN = GEM3B1D_solve(phys_paramsN,num_params)
@test all(isapprox.(energies_arrN[1:7], energies_arrA[1:7]; atol=1e-6)) # should be consistent

# 2. Complex scaling angle should almost not affect the bound states
energies_arr = GEM3B1D_solve(phys_params,num_params,complex_scaling=true)
@test all(isapprox.(energies_arr[1:4], energies_arrA[1:4]; atol=1e-3)) # 5th state already belongs to the continuum

# 3. complex scaling
# check if complex scaling with zero angle yields the same as without complex scaling
pp4 = make_phys_params3B1D(;interactions=[[vg],[vg],[vg]])
np4_nocsm = make_num_params3B1D(;lmax=0,Lmax=0,gem_params=gp, theta_csm = 0.0)
np4_00csm = make_num_params3B1D(;lmax=0,Lmax=0,gem_params=gp, theta_csm = 0.0)
e4_noscm = GEM3B1D_solve(pp4, np4_nocsm, complex_scaling=false)
e4_00csm = GEM3B1D_solve(pp4, np4_00csm, complex_scaling=true)
@test all(isapprox.(e4_noscm[1:5], e4_00csm[1:5]; atol=1e-3))

# check if complex scaling in basis functions and potential (only for gaussian!) gives the same result
pp4a = make_phys_params3B1D(;interactions=[[vga],[vga],[vga]])
np4_10csm = make_num_params3B1D(;lmax=0,Lmax=0,gem_params=gp, theta_csm = 10.0)

e4_csm_basisfkt = GEM3B1D_solve(pp4, np4_10csm, complex_scaling=true)
e4_csm_analytical = GEM3B1D_solve(pp4a, np4_10csm, complex_scaling=true)
@test all(isapprox.(e4_csm_basisfkt[1:5], e4_csm_analytical[1:5]; atol=1e-3))

# 4. Reconstruct two-body results with a spectator particle
# Use the analytical Gaussian potential here; the generic interpolated potential path
# does not support complex-ranged basis functions in this reduced-limit regression.
phys_params2B = make_phys_params2B(;mur=1/(1/1.0 + 1/20.0), interactions=[v_gauss], dim=1)
phys_params3B23 = make_phys_params3B1D(;masses=[1.0,20.0,20.0], interactions=[[],[],[v_gauss]])

gp23 = (;nmax=8, Nmax=1, r1=1.0, rnmax=10.0, R1=10000.0, RNmax=10000.0)
num_params2B = make_num_params2B(;gem_params=(;nmax=8, r1=1.0, rnmax=10.0))
num_params3B23 = make_num_params3B1D(;gem_params=gp23, theta_csm=0.0, complex_ranged_r=true, complex_ranged_R=false)

e2 = GEM2B.GEM2B_solve(phys_params2B, num_params2B)
e2cr = GEM2B.GEM2B_solve(phys_params2B, num_params2B; complex_ranged=true)
e2csmcr = GEM2B.GEM2B_solve(phys_params2B, num_params2B; complex_ranged=true, complex_scaling=true)
e3 = GEM3B1D_solve(phys_params3B23, num_params3B23)
e3cr = GEM3B1D_solve(phys_params3B23, num_params3B23; complex_ranged=true, complex_scaling=false)
e3csmcr = GEM3B1D_solve(phys_params3B23, num_params3B23; complex_ranged=true, complex_scaling=true)

simax23 = findlast(real.(e2) .< 0)
@test simax23 !== nothing
@test length(e3) >= simax23
@test length(e3cr) >= simax23
@test length(e3csmcr) >= simax23
@test all(isapprox.(e3[1:simax23], e2[1:simax23]; atol=1e-3))
@test all(isapprox.(e3cr[1:simax23], e2cr[1:simax23]; atol=2e-3))
@test all(isapprox.(e3csmcr[1:simax23], e2csmcr[1:simax23]; atol=2e-3))

# 5. return_wavefunctions=true
energies_wf, wfs = GEM3B1D_solve(phys_params, num_params; return_wavefunctions=true)
@test length(energies_wf) > 0
@test all(isfinite.(energies_wf[1:4]))
@test size(wfs, 1) == size(wfs, 2)

# 6. deprecated keyword aliases
@test_logs (:warn, r"wf_bool is deprecated") GEM3B1D_solve(phys_params, num_params; wf_bool=false)
@test_logs (:warn, r"csm_bool is deprecated") GEM3B1D_solve(phys_params, num_params; csm_bool=false)
@test_logs (:warn, r"debug_bool is deprecated") GEM3B1D_solve(phys_params, num_params; debug_bool=false)

# 7. sanity check failure paths (sanity_checks3B throws ErrorException)
pp_badsize1d  = (masses=[1.0,1.0], species=[:b,:b,:b], interactions=[[],[],[]], parity=1)
@test_throws ErrorException GEM3B1D_solve(pp_badsize1d, num_params)  # wrong masses size

pp_bad1ident1d = (masses=[1.0,1.0,1.0], species=[:b,:x,:y], interactions=[[],[],[]], parity=1)
@test_throws ErrorException GEM3B1D_solve(pp_bad1ident1d, num_params) # 1 boson: invalid

pp_badmasses1d = (masses=[1.0,2.0,1.0], species=[:b,:b,:y], interactions=[[],[],[]], parity=1)
@test_throws ErrorException GEM3B1D_solve(pp_badmasses1d, num_params) # bosons with unequal masses