# For testing the new feature CSM (only 1D code currently)

using BenchmarkTools
using Printf

include("../src/GEM1D.jl")
using .GEM1D

## general inputs:
# Physical parameters
mass_arr=[1.0,Inf]; mu = 1/(1/mass_arr[1]+1/mass_arr[2])
v0,mu_g = -10.0,1.0
lambda=8.0
# r-interaction:
v_poschl(r) = -lambda*(lambda+1)/2/mu*1/cosh(r)^2
v_gauss(r) = v0*exp(-mu_g*r^2)
gaussopt=[1,v0,mu_g]
phys_params = (;hbar=1.0,mass_arr,vint_arr=[v_poschl],lmax=0)
# numerical parameters:
gem_params = (;nmax=10,r1=0.1,rnmax=10.0) # gem_params
num_params = (;gem_params,omega_cr=0.5,theta_csm=0.0)
wf_bool,cr_bool,csm_bool = 0,0,0

#comparison function
function comparison(energies_arr,lambda;gauss_bool=0)
    neven = lastindex(0:2:Int(floor(lambda-1)));
    nodd = lastindex(1:2:Int(floor(lambda-1)));
    simax = min(lastindex(energies_arr),neven); # max stateindex
    num_arr = real.(energies_arr[1:simax]) ## only real part!
    ex_arr = [-(lambda-i)^2/2/mu for i=0:2:Int(floor(lambda-1))] # only for even states!
    gauss_bool == 1 && (ex_arr .*= 0.0)
    diff_arr = ex_arr .- num_arr
    println("Comparison:")
    @printf("%-7s %-15s %-15s %-15s\n", "Index",  "Numeric Result", "Exact Result", "Difference")
    for i in 1:simax
        @printf("%-7d %-15.6f %-15.6f %-15.6f\n", i, num_arr[i], ex_arr[i], diff_arr[i])
    end
end

## Tests via Poschl-Teller:
#1. Is the output for csm_bool = 0 the same as before?
println("\ncsm_bool = 0:")
energies_arr = GEM1D.GEM_solve1D(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=0)
comparison(energies_arr,lambda)

#1b. Is the output for csm_bool = 0 the same as before? with cr
println("\ncsm_bool = 0, cr_bool = 1:")
energies_arr = GEM1D.GEM_solve1D(phys_params,num_params,wf_bool;cr_bool=1,csm_bool=0)
comparison(energies_arr,lambda)

#2. Does  csm_bool = 1, theta_csm = 0.0 recover the old results?
println("\ncsm_bool = 1, theta_csm = 0.0:")
energies_arr = GEM1D.GEM_solve1D(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=1)
comparison(energies_arr,lambda)

#3. Any bugs with theta_csm != 0.0?
num_params = (;gem_params,omega_cr=0.5,theta_csm=5.0)
println("\ncsm_bool = 1, theta_csm = $(num_params.theta_csm):")
energies_arr = GEM1D.GEM_solve1D(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=1)
comparison(energies_arr,lambda)

#4. check functionality of gaussopt
phys_params = (;hbar=1.0,mass_arr,vint_arr=[v_gauss],lmax=0)
println("\n gaussopt = 0: csm_bool = 1, theta_csm = $(num_params.theta_csm)")
energies_arr = GEM1D.GEM_solve1D(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=1)
comparison(energies_arr,lambda;gauss_bool=1)
println("\n gaussopt = 1: csm_bool = 1, theta_csm = $(num_params.theta_csm)")
energies_arr = GEM1D.GEM_solve1D(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=1, gaussopt=gaussopt)
comparison(energies_arr,lambda;gauss_bool=1)

#4. Benchmark to check if code became slower

#= println("\n Benchmarking:")
print("gaussopt = 0, csm_bool = 0: ")
@btime GEM1D.GEM_solve1D(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=0); #3.372 ms (174747 allocations: 2.78 MiB)
print("\n gaussopt = 1, csm_bool = 0: ")
@btime GEM1D.GEM_solve1D(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=0,gaussopt=gaussopt); #124.600 μs (107 allocations: 25.70 KiB)
print("\n gaussopt = 0, csm_bool = 1: ")
@btime GEM1D.GEM_solve1D(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=1); #3.779 ms (174859 allocations: 4.45 MiB)
print("\n gaussopt = 1, csm_bool = 1: ")
@btime GEM1D.GEM_solve1D(phys_params,num_params,wf_bool;cr_bool=0,csm_bool=1,gaussopt=gaussopt); #71.600 μs (109 allocations: 42.14 KiB)
=#


return nothing