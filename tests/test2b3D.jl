# For testing the new feature CSM in the 3D case

using BenchmarkTools
using Printf

include("../src/GEM.jl")
using .GEM

## general inputs:
# Physical parameters
mass_arr=[1.0,Inf]; mu = 1/(1/mass_arr[1]+1/mass_arr[2])
v0,mu_g = -10.0,1.0
lambda=8.0
# r-interaction:
v_coulomb(r) = -1/r
v_gauss(r) = v0*exp(-mu_g*r^2)
gaussopt=[[1,v0,mu_g]]
phys_params = (;hbar=1.0,mass_arr,vint_arr=[v_coulomb],lmax=0)
# numerical parameters:
gem_params = (;nmax=10,r1=0.1,rnmax=50.0) # gem_params
num_params = (;gem_params,omega_cr=0.5,theta_csm=0.0, threshold = 10^-8)
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
energies_arr1a = GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=0)
comparison(energies_arr1a)

#1b. Test of preallocation
println("\n1b) Test of preallocation:")
pa = GEM.PreallocStruct(num_params, cr_bool, csm_bool)
GEM.GEM_solve!(pa,phys_params,num_params,wf_bool,cr_bool,csm_bool,[0,1.0,1.0])
energies_arr1b = pa.energies
comparison(energies_arr1b)

#1c. Is the output for csm_bool = 0 the same as before? with cr
println("\n1c) csm_bool = 0, cr_bool = 1:")
energies_arr1c = GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=1,csm_bool=0)
comparison(energies_arr1c)

#2. Does  csm_bool = 1, theta_csm = 0.0 recover the old results?
println("\n2) csm_bool = 1, theta_csm = 0.0:")
energies_arr2 = GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=1)
comparison(energies_arr2)

#3. Any bugs with theta_csm != 0.0?
num_params = (;gem_params,omega_cr=0.5,theta_csm=5.0, threshold = 10^-8)
println("\n3) csm_bool = 1, theta_csm = $(num_params.theta_csm):")
energies_arr3 = GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=1)
comparison(real.(energies_arr3))

#4a. combination of cr and csm: csm_bool = 1, cr_bool = 1, but theta_csm = 0.0
num_params = (;gem_params,omega_cr=0.5,theta_csm=0.0, threshold = 10^-8)
println("\n4a) cr_bool = 1, csm_bool = 1, theta_csm = $(num_params.theta_csm):")
energies_arr4a = GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=1,csm_bool=1)
comparison(real.(energies_arr4a))

#4b. combination of cr and csm: csm_bool = 1, cr_bool = 1, but theta_csm = 5.0
num_params = (;gem_params,omega_cr=0.5,theta_csm=5.0, threshold = 10^-8)
println("\n4b) cr_bool = 1, csm_bool = 1, theta_csm = $(num_params.theta_csm):")
energies_arr4b = GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=1,csm_bool=1)
comparison(real.(energies_arr4b))

#5. check functionality of gaussopt
phys_params = (;hbar=1.0,mass_arr,vint_arr=[v_gauss],lmax=0)
println("\nCheck of gaussopt:")
println("\n5a)gaussopt = 0: csm_bool = 0, theta_csm = 0.0")
energies_arr5a = GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=0)
comparison(energies_arr5a;gauss_bool=1)
println("\n5b)gaussopt = 1: csm_bool = 0, theta_csm = 0.0")
energies_arr5b = GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=0, gaussopt=gaussopt)
comparison(energies_arr5b;gauss_bool=1)
num_params = (;gem_params,omega_cr=0.5,theta_csm=5.0, threshold = 10^-8)
println("\n5c)gaussopt = 0: csm_bool = 1, theta_csm = $(num_params.theta_csm)")
energies_arr5c = GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=1)
comparison(energies_arr5c;gauss_bool=1)
println("\n5d)gaussopt = 1: csm_bool = 1, theta_csm = $(num_params.theta_csm)")
energies_arr5d = GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=1, gaussopt=gaussopt)
comparison(energies_arr5d;gauss_bool=1)

#6. check coupled channels
num_params = (;gem_params,omega_cr=0.5,theta_csm=0.0, threshold = 10^-12)
wf_bool,cr_bool,csm_bool = 0,0,0
phys_params = (;hbar=1.0,mass_arr,vint_arr=[v_coulomb],lmax=0)
wfun(r) = exp(-r);wfun2(r) = exp(-r) # for wfun = wfun2, the ground state energy recovers the result of wfun = 0?!
dfun(r) = exp(-0.5*r); dfun2(r) = 1.0*exp(-0.5*r) # has massive effect (e.g. on ground state energy). any way to check/estimate?
WCC = [wfun wfun2; wfun2 wfun]
dor = 1; #derivative-order
DCC = reshape([ [dor, dfun], [dor, dfun2], [dor, dfun2], [dor, dfun] ], 2, 2)

println("\n6) Check of coupled channels:")
println("\ncsm_bool = 0, theta_csm = $(num_params.theta_csm)")
energies_arr6a = GEM.GEM_solveCC(phys_params, num_params, wf_bool, WCC, DCC; cr_bool=0, csm_bool=0, diff_bool=1)
comparison(energies_arr6a)
println("\ncsm_bool = 1, theta_csm = $(num_params.theta_csm)")
energies_arr6b = GEM.GEM_solveCC(phys_params, num_params, wf_bool, WCC, DCC; cr_bool=0, csm_bool=1, diff_bool=1)
comparison(energies_arr6b)


#0. Benchmark to check code performance
phys_params = (;hbar=1.0,mass_arr,vint_arr=[v_gauss],lmax=0)

bench_bool = 0
if bench_bool == 1
println("\nBenchmarking:")
println("nmax = num_params.nmax; vgauss")
print("gaussopt = 0, csm_bool = 0: ")
@btime GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=0); #2.409 ms (112758 allocations: 1.80 MiB)
print("\ngaussopt = 1, csm_bool = 0: ")
@btime GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=0,gaussopt=gaussopt); #43.900 μs (73 allocations: 22.38 KiB)
print("\ngaussopt = 0, csm_bool = 1: ")
@btime GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=1); #2.507 ms (112760 allocations: 2.48 MiB)
print("\ngaussopt = 0, csm_bool = 1, cr_bool = 1: ")
@btime GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=1,csm_bool=1); #21.309 ms (918828 allocations: 20.07 MiB)
print("\ngaussopt = 1, csm_bool = 1: ")
@btime GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=1,gaussopt=gaussopt); #62.100 μs (75 allocations: 33.27 KiB)
print("\ngaussopt = 1, csm_bool = 1, cr_bool = 1: ")
@btime GEM.GEM_solve(phys_params,num_params,wf_bool;cr_bool=1,csm_bool=1,gaussopt=gaussopt); #392.500 μs (88 allocations: 161.36 KiB)
print("\nPreallocation: gaussopt = 0, csm_bool = 0: ")
@btime GEM.GEM_solve!(pa,phys_params,num_params,wf_bool,0,0,[0,1.0,1.0]) #2.414 ms (112744 allocations: 1.80 MiB)
print("\nPreallocation: gaussopt = 1, csm_bool = 0: ")
@btime GEM.GEM_solve!(pa,phys_params,num_params,wf_bool,0,0,gaussopt) #41.400 μs (58 allocations: 17.94 KiB)
print("\ncoupled channels(2), csm_bool = 0:")
@btime GEM.GEM_solveCC(phys_params, num_params, wf_bool, WCC, DCC; cr_bool=0, csm_bool=0, diff_bool=1) #11.068 ms (545830 allocations: 8.73 MiB) ## big increase from 2.77 ms to 11?!
print("\ncoupled channels(2), csm_bool = 1:")
@btime GEM.GEM_solveCC(phys_params, num_params, wf_bool, WCC, DCC; cr_bool=0, csm_bool=1, diff_bool=1) #12.053 ms (545835 allocations: 14.44 MiB) ## similar big increase!
end



return nothing