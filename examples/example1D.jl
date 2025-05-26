# Example file for using the 2-body GEM program
# these two lines should be called from the parent directory
#=
using Pkg; Pkg.activate(".")
using FewBodyToolkit.GEM2B =#
#println("module GEM2B loaded.")

using Printf, Interpolations, BenchmarkTools

# Example for two particles with Pöschl-Teller interparticle interaction

## Input parameters:
# Physical parameters
mass_arr=[1.0,Inf] # array of masses of particles [m1,m2]
mur = 1/(1/mass_arr[1]+1/mass_arr[2]) # reduced mass

# r-interaction:
lambda=8.0
function v_poschl(r)
    return -lambda*(lambda+1)/2/mur*1/cosh(r)^2
end

phys_params = make_phys_params2B(;mur,vint_arr=[v_poschl],dim=1)

#interpolation example:
r_arr = -10.0:0.5:10.0
v_arr = v_poschl.(r_arr)
v_interpol = cubic_spline_interpolation(r_arr,v_arr,extrapolation_bc=Line())
v_int(r) = v_interpol(r) # we have to make the interaction an object of type "function" again. interpolation objects dont work currently.
phys_params_interpol = make_phys_params2B(;mur,vint_arr=[v_int],dim=1)

# numerical parameters:
nmax=6 # number of Gaussian basis functions for r-variable
r1=0.1;rnmax=10.0; # r1 and rnmax for r-basis functions
gem_params = (;nmax,r1,rnmax) # gem_params
num_params = make_num_params2B(;gem_params)

# Print formatted output and compare with analytical results:
#comparison function
function comparison(num_arr,ex_arr,simax;s1="Numerical", s2="Exact")
    @printf("%-7s %-15s %-15s %-15s\n", "Index",  s1, s2, "Difference")
    for i in 1:simax
        @printf("%-7d %-15.6f %-15.6f %-15.6f\n", i, num_arr[i], ex_arr[i], ex_arr[i] - num_arr[i])
    end
end

# Calling the program:
energies_arr = GEM_solve(phys_params,num_params)

# for comparison with exact energies:
simax = min(findlast(energies_arr.<0),6); # max state index
ex_arr = [-(lambda-i)^2/2/mur for i=0:2:Int(floor(lambda-1))] # only for even states!

println("1. Numerical solution of the 1D problem:")
comparison(energies_arr, ex_arr,simax)

# For comparison with interpolation:
println("\n2. Numerical solution using an interpolated interaction:")
energies_interpol = GEM_solve(phys_params_interpol,num_params)
comparison(energies_interpol, energies_arr, simax;s1="Interpolated", s2="Numerical")


## optimization of GEM-parameters:
# Parameters are optimized for the state with index stateindex:
optim_bool = 1
if optim_bool == 1
    stateindex = 3 # which state to optimize for
    params_opt = GEM2B.GEM_Optim_2B(phys_params, num_params, stateindex)
    gem_params_opt = (;nmax, r1 = params_opt[1], rnmax = params_opt[2])
    num_params_opt = make_num_params2B(;gem_params=gem_params_opt)
    e2_opt = GEM2B.GEM_solve(phys_params,num_params_opt)
    
    println("\n3. Optimization of GEM parameters for E2[$stateindex]:")
    @printf("%-15s %-15s %-15s %-15s %-15s\n", "r1", "rnmax", "E2[$(stateindex-1)]", "E2[$stateindex]", "E2[$(stateindex+1)]")
    println("Before optimization:")
    @printf("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", gem_params.r1, gem_params.rnmax, energies_arr[stateindex-1], energies_arr[stateindex], energies_arr[stateindex+1])
    println("after optimization:")
    @printf("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", gem_params_opt.r1, gem_params_opt.rnmax, e2_opt[stateindex-1], e2_opt[stateindex], e2_opt[stateindex+1])
end
println("")
comparison(e2_opt, ex_arr, simax)

##### For finding an interaction globally scaled from vint, such that state with stateindex has energy target_e2: #####
# GEM-Parameters are also updated
vscale_bool = 1
if vscale_bool == 1
    stateindex = 3; target_e2 = -18.0;
    println("\n4. Scaling the potential such that E2[$stateindex] = $target_e2:")
    phys_params_scaled,num_params_scaled,vscale = GEM2B.v0GEMOptim(phys_params,num_params,stateindex,target_e2)
    e2_v0 = GEM2B.GEM_solve(phys_params_scaled,num_params_scaled)
    
    println("After scaling:")
    @printf("%-15s %-15s %-15s %-15s %-15s\n", "r1", "rnmax", "E2[$(stateindex-1)]", "E2[$stateindex]", "E2[$(stateindex+1)]")
    @printf("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", num_params_scaled.gem_params.r1, num_params_scaled.gem_params.rnmax, e2_v0[stateindex-1], e2_v0[stateindex], e2_v0[stateindex+1])
    println("vscale = $(round(vscale,digits=8)) should be approximately (λ+2)*(λ+2+1)/(λ*(λ+1)) = ", round((lambda+2)*(lambda+2+1)/(lambda*(lambda+1)),digits=8) )
end


## currently not working as intended:
#= #### Using complex-ranged Gaussian basis functions (Gbf):
println("\n5. Complex-ranged Gaussian basis functions:")
lambda = 18.0 # deeper potential with 18 states
function v_poschl_deep(r)
    return -lambda*(lambda+1)/2/mur*1/cosh(r)^2
end
phys_params = make_phys_params2B(;mur,vint_arr=[v_poschl_deep],dim=1)
ex_arr = [-(lambda-i)^2/2/mur for i=0:2:Int(floor(lambda-1))] # only for even states!

# using optimized normal (real-ranged) Gaussians:
nmax=10
num_params = make_num_params2B(;gem_params=(nmax,r1=0.1,rnmax=50.0))
params_opt = GEM2B.GEM_Optim_2B(phys_params, num_params, 7)
gem_params_opt = (;nmax, r1 = params_opt[1], rnmax = params_opt[2])
num_params_opt = make_num_params2B(;gem_params=gem_params_opt)
e2_opt = GEM2B.GEM_solve(phys_params,num_params_opt)
simax = min(findlast(e2_opt.<0),10); # max state index

#complex-ranged Gaussians:
gem_paramsCR = (;nmax=Int(nmax/2),r1=0.1,rnmax=50.0) # for complex-ranged Gbf, nmax_eff = nmax*2
num_paramsCR = make_num_params2B(;gem_params=gem_paramsCR)
params_optCR = GEM2B.GEM_Optim_2B(phys_params, num_paramsCR, 7; cr_bool=1)
gem_params_optCR = (;nmax=Int(nmax/2), r1 = params_optCR[1], rnmax = params_optCR[2])
num_params_optCR = make_num_params2B(;gem_params=gem_params_optCR)
e2_optCR = GEM2B.GEM_solve(phys_params,num_params_optCR;cr_bool=1)

comparison(e2_optCR, ex_arr, simax;s1="Complex-ranged")
comparison(e2_optCR, e2_opt, simax;s1="Complex-ranged", s2="Real-ranged") =#
