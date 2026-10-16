# Benchmark cases. Only the public API is used, so that every release >= 0.5.4 can be benchmarked
# with the same definitions. Changing a case changes the benchmark: add a new name instead.
# All cases use the default overlap cutoff (threshold = 1e-8); at 1e-10 the larger complex-ranged
# bases give spurious deep states, which are not suitable for comparing energies between versions.

v_isgl(r) = -10.0*exp(-r^2)  # numerical potential: exercises the range interpolation
v_1d(r) = -4.0*exp(-r^2)

isgl_pp() = make_phys_params3B3D(;masses=[1.0,2.0,3.0],species=[:x,:y,:z],interactions=[[v_isgl],[v_isgl],[v_isgl]])
isgl_np(n) = make_num_params3B3D(;gem_params=(;nmax=n,Nmax=n,r1=0.1,rnmax=100.0,R1=0.1,RNmax=100.0),
                                  kmax_interpol=2000,theta_csm=10.0,omega_cr=0.9)

gem1d_pp() = make_phys_params3B1D(;masses=[1.0,5.0,5.0],species=[:x,:b,:b],interactions=[[],[v_1d],[v_1d]])
gem1d_np(n) = make_num_params3B1D(;gem_params=(;nmax=n,Nmax=n,r1=0.05,rnmax=15.0,R1=0.15,RNmax=17.0),
                                   lmax=1,Lmax=1,theta_csm=8.0,omega_cr=0.9)

# kind: :isgl / :gem1d (solve + matrix build), :inverse (inverse_solve on ISGL matrices)
const CASES = [
    (name="isgl_real",  kind=:isgl,    sizes=[20,30,40], kw=(;)),
    (name="isgl_csm",   kind=:isgl,    sizes=[16,24,32], kw=(;complex_scaling=true)),
    (name="isgl_cr",    kind=:isgl,    sizes=[12,18,24], kw=(;complex_ranged=:both)),
    (name="gem1d_real", kind=:gem1d,   sizes=[20,30,40], kw=(;)),
    (name="gem1d_cr",   kind=:gem1d,   sizes=[8,12,16],  kw=(;complex_ranged=:both)),
    (name="inverse",    kind=:inverse, sizes=[24,32],    kw=(;)),
]
const WARMUP_SIZE = 4
const INVERSE_ENERGIES = range(-2.5, 0.0, length=6)

params(kind, n) = kind == :gem1d ? (gem1d_pp(), gem1d_np(n)) : (isgl_pp(), isgl_np(n))
solvers(kind) = kind == :gem1d ? (GEM3B1D_solve, GEM3B1D_matrices) : (ISGL_solve, ISGL_matrices)
