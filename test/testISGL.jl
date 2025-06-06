# Tests for the module ISGL (three-body, 3D)
# Tests via harmonic oscillator potential

## Setup:
# physical parameters:
m=1/40.0;
omega = 15; a = 3.0;
vcent_ho(r) = 1/a*(m)/2*omega^2*r^2
phys_params = make_phys_params3B3D(;mass_arr=[m,m,m], svals=["x","y","z"], vint_arr=[[vcent_ho],[vcent_ho],[vcent_ho]])

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

# 3. svals = ["b","b","b"]
exact_resultsBBB = vcat([0],2*ones(1),4*ones(2),6*ones(3),8*ones(4),10*ones(5)) # without degeneracy
phys_paramsBBB = make_phys_params3B3D(;mass_arr=[m,m,m], svals=["b","b","b"], vint_arr=[[vcent_ho],[vcent_ho],[vcent_ho]])
energies_arrBBB = ISGL_solve(phys_paramsBBB,num_params22) /omega .- a;
@test all(isapprox.(energies_arrBBB[1:5], exact_resultsBBB[1:5]; atol=1e-3)) # improved accuracy

