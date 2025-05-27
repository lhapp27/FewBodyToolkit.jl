# Example file for using the 2-body GEM program for a 3D system
#= these two lines should be called from the parent directory:
using Pkg; Pkg.activate(".")
using FewBodyToolkit.GEM2B
=#

using Printf, BenchmarkTools, Interpolations, Plots

## Input parameters:
# Physical parameters:
mass_arr = [1.0, Inf] # array of masses of particles [m1,m2]
mur = 1 / (1/mass_arr[1] + 1/mass_arr[2])
Z = 1.0

function v_coulomb(r)
    return -Z/r
end
phys_params = make_phys_params2B(;mur,vint_arr=[v_coulomb],dim=3)

#interpolation example:
r_arr = 0.001:0.01:20.001
v_arr = v_coulomb.(r_arr)
v_interpol = cubic_spline_interpolation(r_arr,v_arr,extrapolation_bc=Line())
v_int(r) = v_interpol(r) # we have to make the interaction an object of type "function" again. interpolation objects dont work currently.
phys_params_i = make_phys_params2B(;mur,vint_arr=[v_int],dim=3)

# numerical parameters:
nmax=10 # number of Gaussian basis functions for r-variable
r1=0.1;rnmax=10.0; # r1 and rnmax for r-basis functions
gem_params = (;nmax,r1,rnmax) # gem_params
num_params = make_num_params2B(;gem_params)

# Calling the program:
e2 = GEM2B.GEM_solve(phys_params,num_params)
e2_interpol = GEM2B.GEM_solve(phys_params_i,num_params)

simax = min(lastindex(e2),6); # max state index

#comparison function
function comparison(num_arr,ex_arr,simax;s1="Numerical", s2="Exact")
    @printf("%-7s %-15s %-15s %-15s\n", "Index",  s1, s2, "Difference")
    for i in 1:simax
        @printf("%-7d %-15.6f %-15.6f %-15.6f\n", i, num_arr[i], ex_arr[i], ex_arr[i] - num_arr[i])
    end
end

e2_exact = ex_arr = [-Z^2/(2*i^2) for i=1:simax]

println("1. Numerical solution of the 3D problem:")
comparison(e2,e2_exact,simax)

println("\n2. Numerical solution using an interpolated interaction:")
comparison(e2_interpol,e2,simax; s1="Interpolated", s2="Numerical")


##### For optimizing GEM-parameters: #####
# Parameters are optimized for state with index stateindex:
optim_bool = 1
if optim_bool == 1
    stateindex = 4
    params_opt = GEM2B.GEM_Optim_2B(phys_params, num_params, stateindex)
    gem_params_opt = (;nmax, r1 = params_opt[1], rnmax = params_opt[2])
    num_params_opt = make_num_params2B(;gem_params=gem_params_opt)
    e2_opt = GEM2B.GEM_solve(phys_params,num_params_opt)
    
    println("\n3. Optimization of GEM parameters for E2[$stateindex]:")
    @printf("%-15s %-15s %-15s %-15s %-15s\n", "r1", "rnmax", "E2[$(stateindex-1)]", "E2[$stateindex]", "E2[$(stateindex+1)]")
    println("Before optimization:")
    @printf("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", gem_params.r1, gem_params.rnmax, e2[stateindex-1], e2[stateindex], e2[stateindex+1])
    println("after optimization:")
    @printf("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", params_opt[1], params_opt[2], e2_opt[stateindex-1], e2_opt[stateindex], e2_opt[stateindex+1])
    comparison(e2_opt,e2_exact,simax; s1="Optimized")
end


##### For finding an interaction globally scaled from vint, such that state with stateindex has energy target_e2: #####
# GEM-Parameters are also updated
vscale_bool = 0
if vscale_bool == 1
    stateindex = 2; target_e2 = -0.5;
    println("\n4. Scaling the potential such that E2[$stateindex] = $target_e2:")
    phys_params_scaled,num_params_scaled,vscale = GEM2B.v0GEMOptim(phys_params,num_params,stateindex,target_e2)
    e2_v0 = GEM2B.GEM_solve(phys_params_scaled,num_params_scaled)

    @printf("%-15s %-15s %-15s %-15s %-15s\n", "r1", "rnmax", "E2[$(stateindex-1)]", "E2[$stateindex]", "E2[$(stateindex+1)]")
    @printf("%-15.6f %-15.6f %-15.6f %-15.6f %-15.6f\n", num_params.gem_params.r1, num_params.gem_params.rnmax, e2_v0[stateindex-1], e2_v0[stateindex], e2_v0[stateindex+1])
    println("vscale = $(round(vscale,digits=5)) should be approximately 2.0")
end

## Calculation of the wave function
nmax=2 # number of Gaussian basis functions for r-variable
r1=1.0;rnmax=20.0; # r1 and rnmax for r-basis functions
gem_params = (;nmax,r1,rnmax) # gem_params
num_params = make_num_params2B(;gem_params)
e2,wfs = GEM2B.GEM_solve(phys_params,num_params;wf_bool=1)

nu_arr = zeros(Float64,nmax)
GEM2B.buildnu(nu_arr,r1,rnmax,nmax)
wffun(r) = GEM2B.wavefun_point(r,nu_arr,wfs[:,1],phys_params.lmax,phys_params.dim)

dr = 0.1
r_arr = 0.0:dr:50.0
p = plot(xlabel="r", ylabel="\$ r^2\\,|\\psi(r)|^2 \$", title="2-body densities for 3D coulomb system")
density = zeros(length(r_arr),4)
for si = 1:4
    wf = wfs[:,si]
    psi_arr = GEM2B.wavefun_arr(r_arr,phys_params,num_params,wf;cr_bool=0)

    density[:,si] .= abs2.(psi_arr).*r_arr.^2

    plot!(r_arr,density[:,si],label="n=$(si-1), density", lw=2)
    #plot!(r_arr,psi_arr,label="n=$(si-1), wave function", lw=2)

    # check normalization:
    @show(sum(density[:,si])*dr)
end
display(p)
