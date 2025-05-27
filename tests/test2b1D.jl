# For testing the new feature CSM (only 1D code currently)
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
v_poschl(r) = -lambda*(lambda+1)/2/mur*1/cosh(r)^2
v_gauss(r) = v0*exp(-mu_g*r^2)
v_gauss2 = GEM2B.GaussianPotential(v0,mu_g)
phys_params = make_phys_params2B(;vint_arr=[v_poschl],dim=1)
phys_paramsg = make_phys_params2B(;vint_arr=[v_gauss],dim=1)
phys_paramsg2 = make_phys_params2B(;vint_arr=[v_gauss2],dim=1)
# numerical parameters:
gem_params = (;nmax=10,r1=0.1,rnmax=10.0) # gem_params
num_params = make_num_params2B(;gem_params,theta_csm=5.0)
wf_bool,cr_bool,csm_bool = 0,0,0

#comparison function
function comparison(energies_arr,lambda;gauss_bool=0)
    neven = lastindex(0:2:Int(floor(lambda-1)));
    nodd = lastindex(1:2:Int(floor(lambda-1)));
    simax = min(lastindex(energies_arr),neven); # max stateindex
    num_arr = real.(energies_arr[1:simax]) ## only real part!
    ex_arr = [-(lambda-i)^2/2/mur for i=0:2:Int(floor(lambda-1))] # only for even states!
    gauss_bool == 1 && (ex_arr .*= 0.0)
    #diff_arr = ex_arr .- num_arr
    println("Comparison:")
    @printf("%-7s %-15s %-15s %-15s\n", "Index",  "Numeric Result", "Exact Result", "Difference")
    for i in 1:simax
        @printf("%-7d %-15.6f %-15.6f %-15.6f\n", i, num_arr[i], ex_arr[i], ex_arr[i] - num_arr[i])
    end
end

## Tests via Poschl-Teller:
#1. Is the output for csm_bool = 0 the same as before?
println("\ncsm_bool = 0:")
energies_arr = GEM2B.GEM_solve(phys_params,num_params;cr_bool=0,csm_bool=0)
comparison(energies_arr,lambda)


#1b. Is the output for csm_bool = 0 the same as before? with cr
println("\ncsm_bool = 0, cr_bool = 1:")
energies_arr = GEM2B.GEM_solve(phys_params,num_params;cr_bool=1,csm_bool=0)
comparison(energies_arr,lambda)

#2. Does  csm_bool = 1, theta_csm = 0.0 recover the old results?
println("\ncsm_bool = 1, theta_csm = 0.0:")
energies_arr = GEM2B.GEM_solve(phys_params,num_params;cr_bool=0,csm_bool=1)
comparison(energies_arr,lambda)

#3. Any bugs with theta_csm != 0.0?
println("\ncsm_bool = 1, theta_csm = $(num_params.theta_csm):")
energies_arr = GEM2B.GEM_solve(phys_params,num_params;cr_bool=0,csm_bool=1)
comparison(energies_arr,lambda)

#4. check functionality of gaussopt
println("\n gaussopt = 0: csm_bool = 1, theta_csm = $(num_params.theta_csm)")
energies_arr = GEM2B.GEM_solve(phys_paramsg,num_params;cr_bool=0,csm_bool=1)
comparison(energies_arr,lambda;gauss_bool=1)
println("\n gaussopt = 1: csm_bool = 1, theta_csm = $(num_params.theta_csm)")
energies_arr = GEM2B.GEM_solve(phys_paramsg2,num_params;cr_bool=0,csm_bool=1)
comparison(energies_arr,lambda;gauss_bool=1)

#4. Benchmark to check if code became slower
bench_bool = 0
if bench_bool == 1
    println("\n Benchmarking: nmax = $(num_params.nmax), vgauss") # value of nmax? probably 10
    print("gaussopt = 0, csm_bool = 0: ")
    @btime GEM2B.GEM_solve(phys_paramsg,num_params;cr_bool=0,csm_bool=0);
    print("\n gaussopt = 1, csm_bool = 0: ")
    @btime GEM2B.GEM_solve(phys_paramsg2,num_params;cr_bool=0,csm_bool=0);
    print("\n gaussopt = 0, csm_bool = 1: ")
    @btime GEM2B.GEM_solve(phys_paramsg,num_params;cr_bool=0,csm_bool=1);
    print("\n gaussopt = 1, csm_bool = 1: ")
    @btime GEM2B.GEM_solve(phys_paramsg2,num_params;cr_bool=0,csm_bool=1);
end 

return nothing