# For testing the new feature CSM in the 3D case
#=
using Pkg; Pkg.activate(".")
using FewBodyToolkit.GEM2B =#

using BenchmarkTools,Printf

## general inputs:
# Physical parameters
mass_arr=[1.0,Inf]; mur = 1/(1/mass_arr[1]+1/mass_arr[2])
v0,mu_g = -10.0,1.0
lambda=8.0
# r-interaction:
v_coulomb(r) = -1/r
v_gauss(r) = v0*exp(-mu_g*r^2)
v_gauss2 = GEM2B.GaussianPotential(v0,mu_g) # define a GaussianPotential object
@show(typeof(v_gauss),typeof(v_gauss2))

phys_params = make_phys_params2B(;vint_arr=[v_coulomb])
phys_paramsg = make_phys_params2B(;vint_arr=[v_gauss])
phys_paramsg2 = make_phys_params2B(;vint_arr=[v_gauss2])
# numerical parameters:
gem_params = (;nmax=10,r1=0.1,rnmax=50.0) # gem_params
num_params = make_num_params2B(;gem_params)
wf_bool,cr_bool,csm_bool = 0,0,0

#comparison function
function comparison(energies_arr;gauss_bool=0,ex_arr = [])
    simax = min(lastindex(energies_arr),6); # max stateindex
    num_arr = real.(energies_arr[1:simax])
    isempty(ex_arr) && (ex_arr = [-1/(2*i^2) for i=1:simax])
    gauss_bool == 1 && (ex_arr .*= 0.0)
    diff_arr = ex_arr .- num_arr
    println("Comparison:")
    @printf("%-7s %-15s %-15s %-15s\n", "Index",  "Numeric Result", "Exact Result", "Difference")
    for i in 1:simax
        @printf("%-7d %-15.6f %-15.6f %-15.6f\n", i, num_arr[i], ex_arr[i], diff_arr[i])
    end
end

## Tests via Coulomb:
#1a. Is the output for csm_bool = 0 the same as before?
println("\n1a) csm_bool = 0, cr_bool = 0:")
energies_arr1a = GEM2B.GEM_solve(phys_params,num_params; cr_bool=0,csm_bool=0)
comparison(energies_arr1a)

#1b. Test of preallocation
println("\n1b) Test of preallocation:")
pa = GEM2B.PreallocStruct2B(num_params, cr_bool, csm_bool)
GEM2B.GEM_solve!(pa,phys_params,num_params,wf_bool,cr_bool,csm_bool)
energies_arr1b = pa.energies
comparison(energies_arr1b)

#1c. Is the output for csm_bool = 0 the same as before? with cr
println("\n1c) csm_bool = 0, cr_bool = 1:")
energies_arr1c = GEM2B.GEM_solve(phys_params,num_params; cr_bool=1,csm_bool=0)
comparison(energies_arr1c)

#2. Does  csm_bool = 1, theta_csm = 0.0 recover the old results?
println("\n2) csm_bool = 1, theta_csm = 0.0:")
energies_arr2 = GEM2B.GEM_solve(phys_params,num_params; cr_bool=0,csm_bool=1)
comparison(energies_arr2)

#3. Any bugs with theta_csm != 0.0?
num_params = (;gem_params,omega_cr=0.5,theta_csm=5.0, threshold = 10^-8)
println("\n3) csm_bool = 1, theta_csm = $(num_params.theta_csm):")
energies_arr3 = GEM2B.GEM_solve(phys_params,num_params; cr_bool=0,csm_bool=1)
comparison(real.(energies_arr3))

#4a. combination of cr and csm: csm_bool = 1, cr_bool = 1, but theta_csm = 0.0
num_params = (;gem_params,omega_cr=0.5,theta_csm=0.0, threshold = 10^-8)
println("\n4a) cr_bool = 1, csm_bool = 1, theta_csm = $(num_params.theta_csm):")
energies_arr4a = GEM2B.GEM_solve(phys_params,num_params; cr_bool=1,csm_bool=1)
comparison(real.(energies_arr4a))

#4b. combination of cr and csm: csm_bool = 1, cr_bool = 1, but theta_csm = 5.0
num_params = (;gem_params,omega_cr=0.5,theta_csm=5.0, threshold = 10^-8)
println("\n4b) cr_bool = 1, csm_bool = 1, theta_csm = $(num_params.theta_csm):")
energies_arr4b = GEM2B.GEM_solve(phys_params,num_params; cr_bool=1,csm_bool=1)
comparison(real.(energies_arr4b))

#5. check functionality of gaussopt
println("\nCheck of gaussopt:")
println("\n5a)gaussopt = 0: csm_bool = 0, theta_csm = 0.0")
energies_arr5a = GEM2B.GEM_solve(phys_paramsg,num_params; cr_bool=0,csm_bool=0)
comparison(energies_arr5a;gauss_bool=1)
println("\n5b)gaussopt = 1: csm_bool = 0, theta_csm = 0.0")
energies_arr5b = GEM2B.GEM_solve(phys_paramsg2,num_params; cr_bool=0,csm_bool=0)
comparison(energies_arr5b;gauss_bool=1)
num_params = (;gem_params,omega_cr=0.5,theta_csm=5.0, threshold = 10^-8)
println("\n5c)gaussopt = 0: csm_bool = 1, theta_csm = $(num_params.theta_csm)")
energies_arr5c = GEM2B.GEM_solve(phys_paramsg,num_params; cr_bool=0,csm_bool=1)
comparison(energies_arr5c;gauss_bool=1)
println("\n5d)gaussopt = 1: csm_bool = 1, theta_csm = $(num_params.theta_csm)")
energies_arr5d = GEM2B.GEM_solve(phys_paramsg2,num_params; cr_bool=0,csm_bool=1)
comparison(energies_arr5d;gauss_bool=1)

#6. check coupled channels
num_params = (;gem_params,omega_cr=0.5,theta_csm=10.0, threshold = 10^-8)
nullfun(r) = 0.0
phys_paramsn = make_phys_params2B(;vint_arr=[nullfun])
wfun(r) = v_coulomb(r);wfun2(r) = 0.0 # testcase: interaction only in wfun; should yield coulomb results with 2fold degeneracy
dfun(r) = 0.0; dfun2(r) = 0.0 # has massive effect (e.g. on ground state energy). any way to check/estimate?
WCC = [wfun wfun2; wfun2 wfun]
dor = 1; #derivative-order
DCC = reshape([ [dor, dfun], [dor, dfun2], [dor, dfun2], [dor, dfun] ], 2, 2)

println("\n6) Check of coupled channels:")
println("\ncsm_bool = 0, theta_csm = $(num_params.theta_csm)")
energies_arr6a = GEM2B.GEM_solveCC(phys_paramsn, num_params, WCC, DCC;cr_bool=0, csm_bool=0, diff_bool=1)
comparison(energies_arr6a;ex_arr=[-1/2,-1/2,-1/8,-1/8,-1/18,-1/18])
println("\ncsm_bool = 1, theta_csm = $(num_params.theta_csm)")
energies_arr6b = GEM2B.GEM_solveCC(phys_paramsn, num_params, WCC, DCC; cr_bool=0, csm_bool=1, diff_bool=1)
comparison(energies_arr6b;ex_arr=[-1/2,-1/2,-1/8,-1/8,-1/18,-1/18])


#0. Benchmark to check code performance
bench_bool = 0
if bench_bool == 1
println("\nBenchmarking:")
println("nmax = $(num_params.nmax), vgauss")
print("gaussopt = 0, csm_bool = 0: ")
@btime GEM2B.GEM_solve(phys_paramsg,num_params; cr_bool=0,csm_bool=0);
print("\ngaussopt = 1, csm_bool = 0: ")
@btime GEM2B.GEM_solve(phys_paramsg2,num_params; cr_bool=0,csm_bool=0);
print("\ngaussopt = 0, csm_bool = 1: ")
@btime GEM2B.GEM_solve(phys_paramsg,num_params; cr_bool=0,csm_bool=1); 
print("\ngaussopt = 0, csm_bool = 1, cr_bool = 1: ")
@btime GEM2B.GEM_solve(phys_paramsg,num_params; cr_bool=1,csm_bool=1); 
print("\ngaussopt = 1, csm_bool = 1: ")
@btime GEM2B.GEM_solve(phys_paramsg2,num_params; cr_bool=0,csm_bool=1); 
print("\ngaussopt = 1, csm_bool = 1, cr_bool = 1: ")
@btime GEM2B.GEM_solve(phys_paramsg2,num_params; cr_bool=1,csm_bool=1); 
print("\nPreallocation: gaussopt = 0, csm_bool = 0: ")
@btime GEM2B.GEM_solve!(pa,phys_paramsg,num_params,wf_bool,0,0)
print("\nPreallocation: gaussopt = 1, csm_bool = 0: ")
@btime GEM2B.GEM_solve!(pa,phys_paramsg2,num_params,wf_bool,0,0)
print("\ncoupled channels(2), csm_bool = 0:")
@btime GEM2B.GEM_solveCC(phys_paramsg, num_params, WCC, DCC; cr_bool=0, csm_bool=0, diff_bool=1)
print("\ncoupled channels(2), csm_bool = 1:")
@btime GEM2B.GEM_solveCC(phys_paramsg, num_params, WCC, DCC; cr_bool=0, csm_bool=1, diff_bool=1)
end



return nothing