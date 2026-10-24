# Benchmark cases. Only the public API is used, so that every registered release (>= 0.3.1) can be
# benchmarked with the same definitions. Changing a case changes the benchmark: add a new name instead.
# All cases use the default overlap cutoff (threshold = 1e-8); at 1e-10 the larger complex-ranged
# bases give spurious deep states, which are not suitable for comparing energies between versions.

v_isgl(r) = -10.0*exp(-r^2)  # numerical potential: exercises the range interpolation
v_1d(r) = -4.0*exp(-r^2)

# releases < 0.4.0 use mass_arr/svals/vint_arr and csm_bool instead of masses/species/interactions
# and complex_scaling; the definitions below translate, so that the cases are identical for all versions.
oldapi() = pkgversion(FewBodyToolkit) < v"0.4.0"
pp_kw(m, s, v) = oldapi() ? (;mass_arr=m, svals=string.(s), vint_arr=v) : (;masses=m, species=s, interactions=v)
solve_kw(kw) = oldapi() ? NamedTuple(k == :complex_scaling ? (:csm_bool => Int(x)) : (k => x) for (k, x) in pairs(kw)) : kw

isgl_pp() = make_phys_params3B3D(;pp_kw([1.0,2.0,3.0],[:x,:y,:z],[[v_isgl],[v_isgl],[v_isgl]])...)
isgl_np(n) = make_num_params3B3D(;gem_params=(;nmax=n,Nmax=n,r1=0.1,rnmax=100.0,R1=0.1,RNmax=100.0),
                                  kmax_interpol=2000,theta_csm=10.0,omega_cr=0.9)

gem1d_pp() = make_phys_params3B1D(;pp_kw([1.0,5.0,5.0],[:x,:b,:b],[[],[v_1d],[v_1d]])...)
gem1d_np(n) = make_num_params3B1D(;gem_params=(;nmax=n,Nmax=n,r1=0.05,rnmax=15.0,R1=0.15,RNmax=17.0),
                                   lmax=1,Lmax=1,theta_csm=8.0,omega_cr=0.9)

# kind: :isgl / :gem1d (solve + matrix build), :inverse (inverse_solve on ISGL matrices)
# since: first release with the needed features (complex ranges, *_matrices / inverse_solve)
const CASES = [
    (name="isgl_real",  kind=:isgl,    sizes=[20,30,40], since=v"0.3.1", kw=(;)),
    (name="isgl_csm",   kind=:isgl,    sizes=[16,24,32], since=v"0.3.1", kw=(;complex_scaling=true)),
    (name="isgl_cr",    kind=:isgl,    sizes=[12,18,24], since=v"0.5.1", kw=(;complex_ranged=:both)),
    (name="gem1d_real", kind=:gem1d,   sizes=[20,30,40], since=v"0.3.1", kw=(;)),
    (name="gem1d_cr",   kind=:gem1d,   sizes=[8,12,16],  since=v"0.5.2", kw=(;complex_ranged=:both)),
    (name="inverse",    kind=:inverse, sizes=[24,32],    since=v"0.5.4", kw=(;)),
    # large: only run when named with --cases= (N = 14700; up to ~50 GB and ~30 min per point for v0.5.4)
    (name="isgl_real_large", kind=:isgl, sizes=[70], since=v"0.3.1", kw=(;), large=true),
    (name="isgl_cr_large",   kind=:isgl, sizes=[35], since=v"0.5.1", kw=(;complex_ranged=:both), large=true),
]
const WARMUP_SIZE = 4
const INVERSE_ENERGIES = range(-2.5, 0.0, length=6)

params(kind, n) = kind == :gem1d ? (gem1d_pp(), gem1d_np(n)) : (isgl_pp(), isgl_np(n))
solvers(kind) = kind == :gem1d ? (:GEM3B1D_solve, :GEM3B1D_matrices) : (:ISGL_solve, :ISGL_matrices)
