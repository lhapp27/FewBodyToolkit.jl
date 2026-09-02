# Tests for the module GEM3B1D (three-body, 1D)
using QuadGK
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

# same comparison at lmax=Lmax=1. This is the case that matters: the exponents nn = la+La+lb+Lb-2s
# reaching the potential's range-interpolation are {0} for lmax=Lmax=0 but {0,1,2,3,4} here, so an
# n-dependent error in the complex-scaling prefactor of precompute_varr! is invisible at l=0 and
# catastrophic at l=1. (It was: the prefactor read csmfac^(2n+1) instead of csmfac^(n+1), which put
# entries such as -896.2-4460.3im where the correct lowest state is -5.619-2.5e-7im.)
np4_10csm_l1 = make_num_params3B1D(;lmax=1,Lmax=1,gem_params=gp, theta_csm = 10.0)
e4_csm_basisfkt_l1 = GEM3B1D_solve(pp4, np4_10csm_l1, complex_scaling=true)
e4_csm_analytical_l1 = GEM3B1D_solve(pp4a, np4_10csm_l1, complex_scaling=true)
@test all(isapprox.(e4_csm_basisfkt_l1[1:5], e4_csm_analytical_l1[1:5]; atol=1e-6))

# and for a non-Gaussian potential the rotated-argument identity can be checked directly on the
# radial integrals that precompute_varr! interpolates: the value it must reproduce is
# int V(r e^{i th}) r^n exp(-alpha r^2) dr.
let th = 10.0, csmfac = exp(-im*th*pi/180), Vt = r -> -3.0*sech(r)^2
    for n in 0:4, alpha in (0.3, 2.0)
        want = quadgk(r -> Vt(r*exp(im*th*pi/180))*r^n*exp(-alpha*r^2), -Inf, 0, Inf)[1]
        got  = quadgk(r -> Vt(r)*r^n*exp(-alpha*csmfac^2*r^2), -Inf, 0, Inf)[1]*csmfac^(n+1)
        @test isapprox(got, want; atol=1e-10, rtol=1e-10)
    end
end

# 4. return_wavefunctions=true
energies_wf, wfs = GEM3B1D_solve(phys_params, num_params; return_wavefunctions=true)
@test length(energies_wf) > 0
@test all(isfinite.(energies_wf[1:4]))
@test size(wfs, 1) == size(wfs, 2)

# 5. deprecated keyword aliases
@test_logs (:warn, r"wf_bool is deprecated") GEM3B1D_solve(phys_params, num_params; wf_bool=false)
@test_logs (:warn, r"csm_bool is deprecated") GEM3B1D_solve(phys_params, num_params; csm_bool=false)
@test_logs (:warn, r"debug_bool is deprecated") GEM3B1D_solve(phys_params, num_params; debug_bool=false)

# 6. sanity check failure paths (sanity_checks3B throws ErrorException)
pp_badsize1d  = (masses=[1.0,1.0], species=[:b,:b,:b], interactions=[[],[],[]], parity=1)
@test_throws ErrorException GEM3B1D_solve(pp_badsize1d, num_params)  # wrong masses size

pp_bad1ident1d = (masses=[1.0,1.0,1.0], species=[:b,:x,:y], interactions=[[],[],[]], parity=1)
@test_throws ErrorException GEM3B1D_solve(pp_bad1ident1d, num_params) # 1 boson: invalid

pp_badmasses1d = (masses=[1.0,2.0,1.0], species=[:b,:b,:y], interactions=[[],[],[]], parity=1)
@test_throws ErrorException GEM3B1D_solve(pp_badmasses1d, num_params) # bosons with unequal masses


# 7. Complex-ranged basis functions: nu -> nu*(1 + i*omega), selectable per Jacobi coordinate.
# The effective Gaussian range etaeff of the interpolated interaction becomes complex, but stays in
# the sector |arg| <= atan(omega); the radial integrals are interpolated over that sector.
# kmax_interpol is reduced throughout: the setup cost grows with the number of angular nodes, and
# 150 radial nodes already resolve the mesh to well below the tolerances used here.
V0_pt = 6.0; a_pt = 1.0
v_pt(r) = -V0_pt/cosh(r/a_pt)^2   # Poeschl-Teller: a non-analytic central potential
pp_pt = make_phys_params3B1D(;interactions=[[v_pt],[v_pt],[v_pt]])
gp_cr = (;nmax=8,Nmax=8,r1=0.2,rnmax=8.0,R1=0.2,RNmax=7.0)

# 7a. omega -> 0 regression. With omega=0 the conjugate half of the basis duplicates the original
# half exactly, so S is rank-deficient by construction and the redundant directions are removed by
# the threshold in cutSmallEV. The surviving spectrum must be the real-basis one. This is the
# sharpest check of the conjugation of the bra ranges: getting it wrong changes these numbers.
np_cr0 = make_num_params3B1D(;gem_params=gp_cr,omega_cr=0.0,kmax_interpol=150)
e_pt_ref0 = GEM3B1D_solve(pp_pt,np_cr0)
for (mode,factor) in ((:r,2),(:R,2),(:both,4))
    e0 = GEM3B1D_solve(pp_pt,np_cr0;complex_ranged=mode)
    @test all(isapprox.(e0[1:4], e_pt_ref0[1:4]; atol=1e-8))
    @test length(e0) == factor*length(e_pt_ref0) # each complex-ranged coordinate doubles its ranges
end


# 8. Accuracy of the two-dimensional (sector) interpolation itself: compare the interpolant against
# a direct quadgk evaluation, at points that are not mesh nodes, for several exponents n and both
# signs of arg(alpha) (the negative half is reached by Schwarz reflection).
const G1D = FewBodyToolkit.GEM3B1D
obs_none = (;stateindices=[],centobs_arr=[[],[],[]],R2_arr=[0,0,0])
function cr1d_interpol_error(kmax_theta)
    np = make_num_params3B1D(;lmax=1,Lmax=1,gem_params=gp_cr,omega_cr=0.9,kmax_interpol=150,kmax_theta)
    sp = G1D.size_estimate(pp_pt,np,obs_none,true,false)
    precomp_arrs,interpol_arrs,_,_ = G1D.preallocate_data(pp_pt,np,obs_none,sp,false,true,false)
    G1D.precompute_3B1D(pp_pt,np,sp,precomp_arrs,true,false)
    G1D.interpolNshoulder(pp_pt,np,obs_none,sp,precomp_arrs,interpol_arrs,false,false,true)

    la, ha = log(interpol_arrs.alpha_arr[1]), log(interpol_arrs.alpha_arr[end])
    theta_max = atan(0.9)
    worst = 0.0
    for n in (0,2,4), fa in (0.05,0.3,0.7,0.95), ft in (-0.9,-0.3,0.15,0.5,1.0)
        alpha = exp(la+(ha-la)*fa)*cis(ft*theta_max)
        val = FewBodyToolkit.interpol_lookup(interpol_arrs.w_interpol_arr[1,1,n],alpha,true)
        ref = G1D.vcent_integration(v_pt,alpha,n,nothing)
        worst = max(worst, abs(val-ref)/abs(ref))
    end
    return worst
end
err_coarse = cr1d_interpol_error(9)
err_default = cr1d_interpol_error(25)
err_fine = cr1d_interpol_error(33)
@test err_default < 1e-5              # default angular resolution
@test err_fine < err_coarse/10        # and it converges in kmax_theta


# 9. Complex ranges against real ranges for a converged Poeschl-Teller well: the two bases differ
# (the real Gaussian is not in the complex-ranged span), but must locate the same bound states.
gp_ptc = (;nmax=12,Nmax=12,r1=0.15,rnmax=12.0,R1=0.15,RNmax=12.0)
np_ptc = make_num_params3B1D(;gem_params=gp_ptc,omega_cr=0.9,kmax_interpol=150)
e_pt_real = GEM3B1D_solve(pp_pt,np_ptc)
@test e_pt_real[1] < 0
for mode in (:r,:R)
    @test all(isapprox.(GEM3B1D_solve(pp_pt,np_ptc;complex_ranged=mode)[1:4], e_pt_real[1:4]; atol=1e-3))
end

# reduction to the two-body problem: switch off two of the three interactions and give the spectator
# coordinate a single, very broad Gaussian, so that <T_R> -> 0. The r-part of the basis then mirrors
# the GEM2B basis exactly, and complex_ranged=:r must reproduce GEM2B's complex-ranged result.
gp_spec1d = (;nmax=10,Nmax=1,r1=0.3,rnmax=12.0,R1=1.0e4,RNmax=1.0e4)
np_spec1d = make_num_params3B1D(;gem_params=gp_spec1d,omega_cr=0.9,kmax_interpol=150)
pp_spec1d = make_phys_params3B1D(;masses=[1.0,1.0,1.0],species=[:x,:y,:z],interactions=[[v_pt],[],[]])
pp2B_pt = make_phys_params2B(;mur=0.5,interactions=[v_pt],dim=1)
np2B_pt = make_num_params2B(;gem_params=(;nmax=10,r1=0.3,rnmax=12.0),omega_cr=0.9)

@test isapprox(GEM3B1D_solve(pp_spec1d,np_spec1d)[1],
               GEM2B_solve(pp2B_pt,np2B_pt)[1]; atol=1e-6)
@test isapprox(GEM3B1D_solve(pp_spec1d,np_spec1d;complex_ranged=:r)[1],
               GEM2B_solve(pp2B_pt,np2B_pt;complex_ranged=true)[1]; atol=1e-4)

# complex ranges together with complex scaling. Neither symmetry survives, so the full matrix is
# filled and the theta-mesh loses the Schwarz reflection; at theta=0 this must still reproduce the
# complex-ranged result.
np_crcsm = make_num_params3B1D(;gem_params=gp_cr,omega_cr=0.9,kmax_interpol=150,theta_csm=0.0)
e_cr_only = GEM3B1D_solve(pp_pt,np_crcsm;complex_ranged=:r)
e_cr_csm0 = GEM3B1D_solve(pp_pt,np_crcsm;complex_ranged=:r,complex_scaling=true)
@test all(isapprox.(real.(e_cr_csm0[1:3]), e_cr_only[1:3]; atol=1e-6))
@test all(abs.(imag.(e_cr_csm0[1:3])) .< 1e-8)

# At a finite scaling angle the two effects combine in the integrand: the radial integrals are taken
# at alpha*exp(-2i*theta_csm) with |arg(alpha)| <= atan(omega), so they only converge while
# atan(omega) + 2*theta_csm < 90 deg. Inside that range the bound state is barely moved; past it the
# combination is rejected with a message rather than failing inside quadgk.
np_csm20 = make_num_params3B1D(;gem_params=gp_cr,omega_cr=0.9,kmax_interpol=150,theta_csm=20.0)
e_csm20 = GEM3B1D_solve(pp_pt,np_csm20;complex_ranged=:r,complex_scaling=true)
@test isapprox(real(e_csm20[1]), e_cr_only[1]; atol=1e-3)
np_csm40 = make_num_params3B1D(;gem_params=gp_cr,omega_cr=0.9,kmax_interpol=150,theta_csm=40.0)
@test_throws ErrorException GEM3B1D_solve(pp_pt,np_csm40;complex_ranged=:r,complex_scaling=true)
# ... while an analytically treated interaction is unaffected, since it is never integrated
pp_ga = make_phys_params3B1D(;interactions=[[vga],[vga],[vga]])
@test !isempty(GEM3B1D_solve(pp_ga,np_csm40;complex_ranged=:r,complex_scaling=true))
