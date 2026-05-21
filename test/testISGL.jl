# Tests for the module ISGL (three-body, 3D)
# Tests via harmonic oscillator potential

## Setup:
# physical parameters:
m=1/40.0;
omega = 15; a = 3.0;
vcent_ho(r) = 1/a*(m)/2*omega^2*r^2
phys_params = make_phys_params3B3D(;masses=[m,m,m], species=[:x,:y,:z], interactions=[[vcent_ho],[vcent_ho],[vcent_ho]])
vg(r) = -5*exp(-r^2)
vga = GaussianPotential(-5.0, 1.0)

# numerical parameters:
gp = (;nmax=10,Nmax=10,r1=0.5,rnmax=8.0,R1=0.5,RNmax=7.0)
num_params = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp)

## Tests:
# Exact results for harmonic oscillator potentials
exact_results = vcat([0],2*ones(3),4*ones(6),6*ones(10),8*ones(15),10*ones(21),12*ones(28))

# 1. lmax=Lmax=0
energies_arr = ISGL_solve(phys_params,num_params) /omega .- a;
@test all(isapprox.(energies_arr[1:12], exact_results[1:12]; atol=1e-2))

# 2. lmax=Lmax=2
num_params22 = make_num_params3B3D(;lmax=2,Lmax=2,gem_params=gp)
energies_arr = ISGL_solve(phys_params,num_params22) /omega .- a;
@test all(isapprox.(energies_arr[1:12], exact_results[1:12]; atol=1e-3)) # improved accuracy

# 3. svals = [:b,:b,:b]
exact_resultsBBB = vcat([0],2*ones(1),4*ones(2),6*ones(3),8*ones(4),10*ones(5)) # without degeneracy
phys_paramsBBB = make_phys_params3B3D(;masses=[m,m,m], species=[:b,:b,:b], interactions=[[vcent_ho],[vcent_ho],[vcent_ho]])
energies_arrBBB = ISGL_solve(phys_paramsBBB,num_params22) /omega .- a;
@test all(isapprox.(energies_arrBBB[1:5], exact_resultsBBB[1:5]; atol=1e-3)) # improved accuracy

# 4. complex scaling
# check if complex scaling with zero angle yields the same as without complex scaling
pp4 = make_phys_params3B3D(;interactions=[[vg],[vg],[vg]])
np4_nocsm = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp, theta_csm = 0.0)
np4_00csm = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp, theta_csm = 0.0)
e4_noscm = ISGL_solve(pp4, np4_nocsm, complex_scaling=false)
e4_00csm = ISGL_solve(pp4, np4_00csm, complex_scaling=true)
@test all(isapprox.(e4_noscm[1:5], e4_00csm[1:5]; atol=1e-3))

# check if complex scaling in basis functions and potential (only for gaussian!) gives the same result
pp4a = make_phys_params3B3D(;interactions=[[vga],[vga],[vga]])
np4_10csm = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp, theta_csm = 10.0)

e4_csm_basisfkt = ISGL_solve(pp4, np4_10csm, complex_scaling=true)
e4_csm_analytical = ISGL_solve(pp4a, np4_10csm, complex_scaling=true)
@test all(isapprox.(e4_csm_basisfkt[1:5], e4_csm_analytical[1:5]; atol=1e-3))

# 5. return_wavefunctions=true with R^2 observables
obs_params = (;stateindices=[1], centobs_arr=[Vector{PotentialFunction}() for _ in 1:3], R2_arr=[1,0,0])
energies_obs, wfs_obs, centobs_out, R2_out = ISGL_solve(phys_params, num_params; return_wavefunctions=true, observ_params=obs_params)
@test length(energies_obs) > 0
@test all(isfinite.(energies_obs[1:3]))
@test size(R2_out, 1) == 3   # 3 Jacobi sets
@test size(R2_out, 2) == 1   # 1 state requested
@test isfinite(R2_out[1,1])
@test R2_out[1,1] > 0        # <R^2> should be positive for HO

# 6. return_wavefunctions=true without observables (also exercises return path)
energies_wf, wfs_wf, co_wf, r2_wf = ISGL_solve(phys_params, num_params; return_wavefunctions=true)
@test size(wfs_wf, 1) == size(wfs_wf, 2)

# 7. deprecated keyword aliases
@test_logs (:warn, r"wf_bool is deprecated") ISGL_solve(phys_params, num_params; wf_bool=false)
@test_logs (:warn, r"csm_bool is deprecated") ISGL_solve(phys_params, num_params; csm_bool=false)
@test_logs (:warn, r"debug_bool is deprecated") ISGL_solve(phys_params, num_params; debug_bool=false)

# 8. observables error for complex_scaling=true
@test_throws ErrorException ISGL_solve(phys_params, num_params; return_wavefunctions=true, complex_scaling=true, observ_params=obs_params)

# 9. sanity check failure paths (ISGL returns nothing instead of throwing)
pp_badsize  = (masses=[m,m], species=[:x,:y,:z], interactions=[[],[],[]], J_tot=0, parity=1)
@test isnothing(ISGL_solve(pp_badsize, num_params))   # wrong masses/species size

pp_bad1ident = (masses=[m,m,m], species=[:b,:x,:y], interactions=[[],[],[]], J_tot=0, parity=1)
@test isnothing(ISGL_solve(pp_bad1ident, num_params)) # 1 boson: invalid identical count

pp_badmasses = (masses=[m,2m,m], species=[:b,:b,:y], interactions=[[],[],[]], J_tot=0, parity=1)
@test isnothing(ISGL_solve(pp_badmasses, num_params))  # bosons with unequal masses

