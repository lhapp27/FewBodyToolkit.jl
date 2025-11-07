# Tests for the module GEM3B1D (three-body, 1D)
# Tests via Gaussian potential

## Setup:
# physical parameters:
v0 = -25.0; mu_g = 1.0
v_gauss = GaussianPotential(v0,mu_g)
v_gauss2(r) = v0*exp(-mu_g*r^2)
phys_params = make_phys_params3B1D(;svals=["b","b","b"], vint_arr=[[v_gauss],[v_gauss],[v_gauss]])
phys_paramsN = make_phys_params3B1D(;svals=["b","b","b"], vint_arr=[[v_gauss2],[v_gauss2],[v_gauss2]])

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
energies_arr = GEM3B1D_solve(phys_params,num_params,csm_bool=1)
@test all(isapprox.(energies_arr[1:4], energies_arrA[1:4]; atol=1e-3)) # 5th state already belongs to the continuum

# 3. complex scaling
# check if complex scaling with zero angle yields the same as without complex scaling
pp4 = make_phys_params3B1D(;vint_arr=[[vg],[vg],[vg]])
np4_nocsm = make_num_params3B1D(;lmax=0,Lmax=0,gem_params=gp, theta_csm = 0.0)
np4_00csm = make_num_params3B1D(;lmax=0,Lmax=0,gem_params=gp, theta_csm = 0.0)
e4_noscm = GEM3B1D_solve(pp4, np4_nocsm, csm_bool=0)
e4_00csm = GEM3B1D_solve(pp4, np4_00csm, csm_bool=1)
@test all(isapprox.(e4_noscm[1:5], e4_00csm[1:5]; atol=1e-3))

# check if complex scaling in basis functions and potential (only for gaussian!) gives the same result
pp4a = make_phys_params3B1D(;vint_arr=[[vga],[vga],[vga]])
np4_10csm = make_num_params3B1D(;lmax=0,Lmax=0,gem_params=gp, theta_csm = 10.0)

e4_csm_basisfkt = GEM3B1D_solve(pp4, np4_10csm, csm_bool=1)
e4_csm_analytical = GEM3B1D_solve(pp4a, np4_10csm, csm_bool=1)
@test all(isapprox.(e4_csm_basisfkt[1:5], e4_csm_analytical[1:5]; atol=1e-3))