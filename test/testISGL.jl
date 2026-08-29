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

# 10. regression: mixed potential lists with wrapped functions should not fail,
# and equivalent combinations should produce the same low-lying energies.
vg_reg = GaussianPotential(-3.0, 0.5)
vgr_reg(r) = vg_reg(r)
vg0_reg(r) = -3.0 * exp(-0.5 * r^2)
vg2_reg(r) = -6.0 * exp(-0.5 * r^2)

gp_reg = (;nmax=6,Nmax=6,r1=0.1,rnmax=20.0,R1=0.1,RNmax=20.0)
np_reg = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp_reg,kmax_interpol=80,threshold=1e-10)

function solve2(interactions)
	pp = make_phys_params3B3D(;masses=[1.0,2.0,3.0], species=[:x,:y,:z], interactions=interactions)
	return ISGL_solve(pp, np_reg)[1:2]
end

# Baselines
e_single_ref = solve2([[vg_reg],[vg_reg],[vg_reg]])
e_double_ref = solve2([[vg_reg,vg_reg],[vg_reg],[vg_reg]])
e_double_fun_ref = solve2([[vg2_reg],[vg_reg],[vg_reg]])

# Single-potential representations are equivalent
@test all(isapprox.(solve2([[vgr_reg],[vgr_reg],[vgr_reg]]), e_single_ref; atol=1e-5))
@test all(isapprox.(solve2([[vg0_reg],[vg0_reg],[vg0_reg]]), e_single_ref; atol=1e-5))

# Multi-potential mixtures run and are equivalent to double-strength reference
multi_cases = [
	[[vgr_reg,vgr_reg],[vg_reg],[vg_reg]],
	[[vg0_reg,vg0_reg],[vg_reg],[vg_reg]],
	[[vg_reg,vgr_reg],[vg_reg],[vg_reg]],
	[[vgr_reg,vg_reg],[vg_reg],[vg_reg]],
	[[vg_reg,vg0_reg],[vg_reg],[vg_reg]],
	[[vg0_reg,vg_reg],[vg_reg],[vg_reg]],
	[[vgr_reg,vg0_reg],[vg_reg],[vg_reg]],
	[[vg0_reg,vgr_reg],[vg_reg],[vg_reg]],
]

for ints in multi_cases
	e = solve2(ints)
	@test all(isfinite.(e))
	@test all(isapprox.(e, e_double_ref; atol=1e-5))
	@test all(isapprox.(e, e_double_fun_ref; atol=1e-5))
end


# 11. PowerLawPotential: analytic treatment of V(r) = v0*r^p
# The analytic path must reproduce the numerical CentralPotential path, and for the
# harmonic oscillator (p=2) also the exact spectrum.

# 11a. constructor: p <= -3 must be rejected (radial integral diverges)
# 11a. validity of p is checked in sanity_checks (ISGL requires p > -3), not in the
# constructor, since the bound depends on dimension/angular momentum. ISGL_solve
# returns nothing instead of throwing (as for the other sanity-check failures).
@test PowerLawPotential(1.0, -3.0).p == -3.0   # construction itself is allowed
@test isnothing(ISGL_solve(make_phys_params3B3D(;masses=[m,m,m],species=[:x,:y,:z],interactions=[[PowerLawPotential(1.0,-3.0)],[],[]]),num_params))
@test isnothing(ISGL_solve(make_phys_params3B3D(;masses=[m,m,m],species=[:x,:y,:z],interactions=[[],[PowerLawPotential(1.0,-4.0)],[]]),num_params))
@test PowerLawPotential(-1.0, -1).p == -1.0   # integer exponents are accepted

# 11b. harmonic oscillator (p=2) against the exact spectrum
vho_pow = PowerLawPotential(1/a*(m)/2*omega^2, 2.0)
pp_pow = make_phys_params3B3D(;masses=[m,m,m], species=[:x,:y,:z], interactions=[[vho_pow],[vho_pow],[vho_pow]])
e_pow = ISGL_solve(pp_pow,num_params) /omega .- a;
@test all(isapprox.(e_pow[1:12], exact_results[1:12]; atol=1e-2))

# 11c. analytic vs numerical path for the very same potential (vcent_ho == vho_pow)
e_cent_ho = ISGL_solve(phys_params,num_params)
e_pow_ho  = ISGL_solve(pp_pow,num_params)
@test all(isapprox.(e_cent_ho[1:12], e_pow_ho[1:12]; rtol=1e-5)) # the interpolated path is only good to ~3e-6 here

# 11d. also with higher partial waves (exercises Lsum>0 in the shoulder solve)
e_pow_22 = ISGL_solve(pp_pow,num_params22) /omega .- a;
@test all(isapprox.(e_pow_22[1:12], exact_results[1:12]; atol=1e-3))

# 11e. Coulomb (p=-1): analytic vs numerical, on a Ps^- -like system
vatt_pow = PowerLawPotential(-1.0,-1.0); vrep_pow = PowerLawPotential(1.0,-1.0)
vatt_cent(r) = -1.0/r; vrep_cent(r) = 1.0/r
gp_coul = (;nmax=8,Nmax=8,r1=0.5,rnmax=25.0,R1=0.5,RNmax=25.0)
np_coul = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp_coul)
pp_coul_pow  = make_phys_params3B3D(;masses=[1.0,1.0,1.0], species=[:x,:y,:z], interactions=[[vrep_pow],[vatt_pow],[vatt_pow]])
pp_coul_cent = make_phys_params3B3D(;masses=[1.0,1.0,1.0], species=[:x,:y,:z], interactions=[[vrep_cent],[vatt_cent],[vatt_cent]])
e_coul_pow  = ISGL_solve(pp_coul_pow, np_coul)
e_coul_cent = ISGL_solve(pp_coul_cent, np_coul)
@test all(isapprox.(e_coul_pow[1:3], e_coul_cent[1:3]; atol=1e-6))
@test e_coul_pow[1] < -0.25   # Ps^- is bound below the Ps threshold

# 11f. non-integer exponent (exercises the general Gamma-function path)
vpl_pow  = PowerLawPotential(0.7, 1.5)
vpl_cent(r) = 0.7*r^1.5
pp_pl_pow  = make_phys_params3B3D(;masses=[m,m,m], species=[:x,:y,:z], interactions=[[vpl_pow],[vpl_pow],[vpl_pow]])
pp_pl_cent = make_phys_params3B3D(;masses=[m,m,m], species=[:x,:y,:z], interactions=[[vpl_cent],[vpl_cent],[vpl_cent]])
@test all(isapprox.(ISGL_solve(pp_pl_pow,num_params)[1:5], ISGL_solve(pp_pl_cent,num_params)[1:5]; rtol=1e-5))

# 11g. mixed interaction lists: checks index bookkeeping (powopt_arr padding, nint_arr)
ints_mix_pow  = [[vga,vho_pow],[vho_pow],[vho_pow]]
ints_mix_cent = [[vga,vcent_ho],[vcent_ho],[vcent_ho]]
pp_mix_pow  = make_phys_params3B3D(;masses=[m,m,m], species=[:x,:y,:z], interactions=ints_mix_pow)
pp_mix_cent = make_phys_params3B3D(;masses=[m,m,m], species=[:x,:y,:z], interactions=ints_mix_cent)
e_mix_pow  = ISGL_solve(pp_mix_pow,num_params)
e_mix_cent = ISGL_solve(pp_mix_cent,num_params)
@test all(isfinite.(e_mix_pow[1:5]))
@test all(isapprox.(e_mix_pow[1:5], e_mix_cent[1:5]; rtol=1e-5))

# 11h. complex scaling: theta=0 must reproduce the non-scaled result
np_pow_00csm = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp, theta_csm = 0.0)
@test all(isapprox.(ISGL_solve(pp_pow,np_pow_00csm,complex_scaling=false)[1:5],
                    ISGL_solve(pp_pow,np_pow_00csm,complex_scaling=true)[1:5]; atol=1e-3))

# 11i. complex scaling at finite angle: analytic path must match the numerical one.
# This is the actual check of the csmfac^(-p) factor used for power-law potentials.
np_pow_10csm = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp, theta_csm = 10.0)
e_csm_pow  = ISGL_solve(pp_pow,     np_pow_10csm, complex_scaling=true)
e_csm_cent = ISGL_solve(phys_params,np_pow_10csm, complex_scaling=true)
@test all(isapprox.(e_csm_pow[1:5], e_csm_cent[1:5]; rtol=1e-3))

# 11j. PowerLawPotential as a central observable (uses the analytic precompute_varr! method)
obs_pow  = (;stateindices=[1], centobs_arr=[PotentialFunction[PowerLawPotential(1.0,2.0)],PotentialFunction[],PotentialFunction[]], R2_arr=[1,0,0])
obs_cent = (;stateindices=[1], centobs_arr=[PotentialFunction[CentralPotential(r->r^2)],PotentialFunction[],PotentialFunction[]], R2_arr=[1,0,0])
_,_,co_pow,_  = ISGL_solve(phys_params, num_params; return_wavefunctions=true, observ_params=obs_pow)
_,_,co_cent,_ = ISGL_solve(phys_params, num_params; return_wavefunctions=true, observ_params=obs_cent)
@test isfinite(co_pow[1,1,1])
@test isapprox(co_pow[1,1,1], co_cent[1,1,1]; rtol=1e-6)


# 12. Complex-ranged basis functions: nu -> nu*(1 + i*omega), selectable per Jacobi coordinate.
# Supported only for the analytically treated potentials (GaussianPotential, PowerLawPotential),
# since the interpolated path for a generic central potential needs a real etaprc.

# 12a. option parsing
@test FewBodyToolkit.ISGL.parse_complex_ranged(:none) == (false,false)
@test FewBodyToolkit.ISGL.parse_complex_ranged(:r)    == (true,false)
@test FewBodyToolkit.ISGL.parse_complex_ranged(:R)    == (false,true)
@test FewBodyToolkit.ISGL.parse_complex_ranged(:both) == (true,true)
@test FewBodyToolkit.ISGL.parse_complex_ranged(true)  == (true,true)   # Bool alias, as in GEM2B_solve
@test FewBodyToolkit.ISGL.parse_complex_ranged(false) == (false,false)
@test_throws ErrorException FewBodyToolkit.ISGL.parse_complex_ranged(:rR)

# 12b. unsupported potentials and observables are rejected (sanity_checks returns nothing)
# Unit masses (not the m=1/40 of the tests above): with m=1/40 this Gaussian binds nothing, so the
# lowest eigenvalue is a discretised continuum state and cannot be compared across bases in 12f.
pp_cr_gauss = make_phys_params3B3D(;masses=[1.0,1.0,1.0], species=[:x,:y,:z], interactions=[[vga],[vga],[vga]])
pp_cr_cent  = make_phys_params3B3D(;masses=[1.0,1.0,1.0], species=[:x,:y,:z], interactions=[[vg],[vg],[vg]])
gp_cr = (;nmax=6,Nmax=6,r1=0.5,rnmax=8.0,R1=0.5,RNmax=7.0)
np_cr = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp_cr,omega_cr=0.9)

@test isnothing(ISGL_solve(pp_cr_cent,np_cr;complex_ranged=:both))   # plain function -> CentralPotential
@test isnothing(ISGL_solve(pp_cr_cent,np_cr;complex_ranged=:r))
obs_r2 = (;stateindices=[1], centobs_arr=[PotentialFunction[],PotentialFunction[],PotentialFunction[]], R2_arr=[1,0,0])
@test isnothing(ISGL_solve(pp_cr_gauss,np_cr;complex_ranged=:r,return_wavefunctions=true,observ_params=obs_r2))
# ... while the analytic potentials are accepted:
@test !isnothing(ISGL_solve(pp_cr_gauss,np_cr;complex_ranged=:r))

# 12c. basis-size doubling per coordinate
n_none = length(ISGL_solve(pp_cr_gauss,np_cr;complex_ranged=:none))
n_r    = length(ISGL_solve(pp_cr_gauss,np_cr;complex_ranged=:r))
n_R    = length(ISGL_solve(pp_cr_gauss,np_cr;complex_ranged=:R))
n_both = length(ISGL_solve(pp_cr_gauss,np_cr;complex_ranged=:both))
@test n_r == 2*n_none
@test n_R == 2*n_none
@test n_both == 4*n_none

# 12d. energies stay real for complex ranges alone (H and S are hermitian, not just symmetric)
e_cr_r = ISGL_solve(pp_cr_gauss,np_cr;complex_ranged=:r)
@test eltype(e_cr_r) <: Real
@test all(isfinite.(e_cr_r[1:3]))

# 12e. omega -> 0 regression. With omega=0 the conjugate half of the basis duplicates the original
# half exactly, so S is rank-deficient by construction and the redundant directions are removed by
# the threshold in cutSmallEV. The surviving spectrum must be the real-basis one. This is the
# sharpest check of the conjugation of the bra ranges: getting it wrong changes these numbers.
np_cr0 = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp_cr,omega_cr=0.0)
e_ref0 = ISGL_solve(pp_cr_gauss,np_cr0;complex_ranged=:none)
for mode in (:r,:R,:both)
    e0 = ISGL_solve(pp_cr_gauss,np_cr0;complex_ranged=mode)
    @test all(isapprox.(e0[1:3], e_ref0[1:3]; atol=1e-6))
end

# the same regression with l,L > 0. This is the check on the branch cuts: for l>0 the norm carries
# the fractional power (2*nu)^((2l+3)/4) and the prefactors carry (pi/zeta)^(3/2)*(pi/etapr)^(3/2),
# none of which is exercised at l=L=0. A smaller basis keeps the (4x doubled) :both case affordable.
gp_cr1 = (;nmax=4,Nmax=4,r1=0.5,rnmax=8.0,R1=0.5,RNmax=7.0)
np_cr1_0 = make_num_params3B3D(;lmax=1,Lmax=1,gem_params=gp_cr1,omega_cr=0.0)
e_ref1 = ISGL_solve(pp_cr_gauss,np_cr1_0;complex_ranged=:none)
for mode in (:r,:R,:both)
    e1 = ISGL_solve(pp_cr_gauss,np_cr1_0;complex_ranged=mode)
    @test all(isapprox.(e1[1:3], e_ref1[1:3]; atol=1e-6))
end

# 12f. finite omega: the complex-ranged basis is a different span of twice the size, not a superset
# of the real one (the real Gaussian is not in it), so no variational inequality is claimed here.
# Both must locate the same bound state. The tolerance is set by the basis, not by the method: with
# this gp_cr the two bases differ by ~2.5e-3, shrinking to ~8e-4 at nmax=Nmax=8 and ~2e-4 at 10.
e_ref = ISGL_solve(pp_cr_gauss,np_cr;complex_ranged=:none)
@test e_ref[1] < 0   # a genuine bound state, so that the comparison below is meaningful
for mode in (:r,:R,:both)
    ecr = ISGL_solve(pp_cr_gauss,np_cr;complex_ranged=mode)
    @test isapprox(ecr[1], e_ref[1]; atol=3e-3)   # same state, not a spurious one
end

# 12g. power-law interactions with complex ranges (Ps^- -like Coulomb system)
np_coul_cr = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp_coul,omega_cr=0.9)
e_coul_real = ISGL_solve(pp_coul_pow, np_coul_cr; complex_ranged=:none)
e_coul_cr   = ISGL_solve(pp_coul_pow, np_coul_cr; complex_ranged=:r)
@test isapprox(e_coul_cr[1], e_coul_real[1]; atol=1e-3)
@test e_coul_cr[1] < -0.25   # still bound below the Ps threshold

# 12h. complex ranges together with complex scaling. Neither symmetry survives, so the full matrix
# is filled; at theta=0 this must still reproduce the complex-ranged result.
np_cr_00csm = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp_cr,omega_cr=0.9,theta_csm=0.0)
e_cr_csm0 = ISGL_solve(pp_cr_gauss,np_cr_00csm;complex_ranged=:r,complex_scaling=true)
@test all(isapprox.(real.(e_cr_csm0[1:3]), e_cr_r[1:3]; atol=1e-6))
@test all(abs.(imag.(e_cr_csm0[1:3])) .< 1e-8)

# 12i. reduction to the two-body problem: switch off two of the three interactions and give the
# spectator coordinate a single, very broad Gaussian, so that <T_R> -> 0. The r-part of the basis
# then mirrors the GEM2B basis exactly, and complex_ranged=:r must reproduce GEM2B's complex-ranged
# result. Equal masses are used so that the reduced mass of the interacting pair is unambiguous.
# Complex ranges are applied to r only: with Nmax=1 the conjugate partner of a single very broad
# range is nearly linearly dependent on it, which is exactly the situation per-coordinate control
# is meant to avoid.
gp_spec = (;nmax=8,Nmax=1,r1=1.0,rnmax=10.0,R1=1.0e4,RNmax=1.0e4)
np_spec = make_num_params3B3D(;lmax=0,Lmax=0,gem_params=gp_spec,omega_cr=0.9)
pp_spec = make_phys_params3B3D(;masses=[1.0,1.0,1.0], species=[:x,:y,:z], interactions=[[vga],[],[]])

pp2B_spec = make_phys_params2B(;mur=0.5, interactions=[vga], dim=3)
np2B_spec = make_num_params2B(;gem_params=(;nmax=8,r1=1.0,rnmax=10.0),omega_cr=0.9)

e2_real = GEM2B_solve(pp2B_spec,np2B_spec)
e2_cr   = GEM2B_solve(pp2B_spec,np2B_spec;complex_ranged=true)
e3_real = ISGL_solve(pp_spec,np_spec;complex_ranged=:none)
e3_cr   = ISGL_solve(pp_spec,np_spec;complex_ranged=:r)

nb2 = count(<(0), real.(e2_real))   # compare the bound states only
@test nb2 >= 1
@test all(isapprox.(e3_real[1:nb2], e2_real[1:nb2]; atol=1e-3))
@test all(isapprox.(e3_cr[1:nb2],   e2_cr[1:nb2];   atol=1e-3))
