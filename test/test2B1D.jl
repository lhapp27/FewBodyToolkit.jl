# Tests for the module GEM2B (1D scenario)
# Tests via Poschl-Teller potential.

## Setup:
# physical parameters:
lambda=8.0
v_poschl(r) = -lambda*(lambda+1)/2/mur*1/cosh(r)^2
mur = 1.0
phys_params = make_phys_params2B(;interactions=[v_poschl],dim=1)

# numerical parameters:
gem_params = (;nmax=10,r1=0.3,rnmax=3.0) # gem_params
num_params = make_num_params2B(;gem_params)


## Tests:
# Exact results for the Pöschl-Teller potential
exact_results = [-(lambda-i)^2/2/mur for i=0:2:Int(floor(lambda-1))] #

# 1. Standard inputs: csm_bool = 0, cr_bool = 0
energies_arr = GEM2B.GEM2B_solve(phys_params,num_params)
@test all(isapprox.(energies_arr[1:4], exact_results; atol=1e-3))

# 2. Using complex-ranged basis functions: complex_scaling = false, complex_ranged = true
energies_arr = GEM2B.GEM2B_solve(phys_params,num_params;complex_ranged=true)
@test all(isapprox.(energies_arr[1:4], exact_results; atol=1e-3))

# 3. Complex scaling (angle 0°) should have no effect: complex_scaling = true
energies_arr = GEM2B.GEM2B_solve(phys_params,num_params;complex_scaling=true)
@test all(isapprox.(real.(energies_arr[1:4]), exact_results; atol=1e-3))

# 4. Finite complex scaling angle (5°) should have very little effect on the bound states
num_paramsC = make_num_params2B(;gem_params,complex_scaling_angle=5.0)
energies_arr = GEM2B.GEM2B_solve(phys_params,num_paramsC;complex_scaling=true)
@test all(isapprox.(real.(energies_arr[1:4]), exact_results; atol=1e-3))

# 5. Finite complex scaling angle (5°) , with simultaneous complex-ranged basis functions
num_paramsC = make_num_params2B(;gem_params,complex_scaling_angle=5.0)
energies_arr = GEM2B.GEM2B_solve(phys_params,num_paramsC;complex_scaling=true,complex_ranged=true)
@test all(isapprox.(real.(energies_arr[1:4]), exact_results; atol=1e-3))

# PowerLawPotential (1D): analytic treatment of V(x) = v0*|x|^p
# 1D harmonic oscillator; lmax=0 selects the even states E=(2n+1/2)*omega,
# lmax=1 the odd ones E=(2n+3/2)*omega
omega_1d = 0.6
v_ho_pow  = PowerLawPotential(0.5*mur*omega_1d^2, 2.0)
v_ho_cent(r) = 0.5*mur*omega_1d^2*r^2
gp_1d = (;nmax=14,r1=0.2,rnmax=12.0)
np_1d = make_num_params2B(;gem_params=gp_1d)
for (l,off) in [(0,0.5),(1,1.5)]
    pp_pow  = make_phys_params2B(;mur,interactions=[v_ho_pow], dim=1,lmax=l,lmin=l)
    pp_cent = make_phys_params2B(;mur,interactions=[v_ho_cent],dim=1,lmax=l,lmin=l)
    e_pow  = GEM2B.GEM2B_solve(pp_pow, np_1d)
    e_cent = GEM2B.GEM2B_solve(pp_cent,np_1d)
    exact_1d = [(2*n+off)*omega_1d for n=0:3]
    @test all(isapprox.(e_pow[1:4], exact_1d; atol=1e-4))
    @test all(isapprox.(e_pow[1:4], e_cent[1:4]; rtol=1e-8))
end

# |x|^p with non-integer p: analytic vs numerical (the numerical path needs abs() explicitly)
v_pl_pow = PowerLawPotential(0.8,1.5)
v_pl_cent(r) = 0.8*abs(r)^1.5
@test all(isapprox.(GEM2B.GEM2B_solve(make_phys_params2B(;mur,interactions=[v_pl_pow],dim=1),np_1d)[1:4],
                    GEM2B.GEM2B_solve(make_phys_params2B(;mur,interactions=[v_pl_cent],dim=1),np_1d)[1:4]; rtol=1e-8))

# validity check: in 1D with lmax=0 the bound is p > -1, so a pure 1/|x| already diverges
@test_throws ErrorException GEM2B.GEM2B_solve(make_phys_params2B(;mur,interactions=[PowerLawPotential(-1.0,-1.0)],dim=1),np_1d)
# for lmax=1 the bound is p > -3, so p=-1 is fine there
@test all(isfinite.(GEM2B.GEM2B_solve(make_phys_params2B(;mur,interactions=[PowerLawPotential(-1.0,-1.0)],dim=1,lmax=1,lmin=1),np_1d)[1:2]))
