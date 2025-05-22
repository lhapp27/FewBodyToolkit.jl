# not sure why this is an extra file.

using Printf
using BenchmarkTools
using HCubature
using Plots

# Loading the module:
include("../src/ISGL.jl")
using .ISGL
println("module ISGL loaded.")

## Define pair-interaction:
v0,mu_g = -10.0,1.0
v_gauss(r) = v0*exp(-mu_g*r^2)

# Observables:
rad(r) = r
rad2(r) = r^2

## Input parameters:
# Physical parameters
phys_params = (;mass_arr=[1.0,1.0,9.0],svals=["b","b","z"],vint_arr=[[v_gauss],[v_gauss],[v_gauss]],J_tot=0,parity=+1)

# numerical parameters:
gp = (;nmax=7,Nmax=7,r1=0.5,rnmax=10.0,R1=0.5,RNmax=10.0)
num_params = (;lmax=1,Lmax=1,gem_params=gp,theta_csm=0,omega_cr=0,mu0=0.08,c_shoulder=1.6,kmax_interpol=5000,threshold=10^-10)
gb=1
gaussopt3=[[[gb,v0,mu_g]],[[gb,v0,mu_g]],[[gb,v0,mu_g]]]

# Booleans:
wf_bool=0 # 0 -> dont calculate wave functions
csm_bool = 0

# Observables:
obs_arr=[[rad,rad2],[rad,rad2],[rad,rad2]] # old input info -> vmtl rad und rad 2 für alle Jacobi-Sets berechnen
stateindices = 1:3#[1,2]
observ_params = (stateindices,centobs_arr = [[rad,rad2],[rad,rad2],[rad,rad2]],R2_arr = [0,0,0])

## Calling the program:
function comparison(reference_arr,energies_arr)
    simax = min(lastindex(energies_arr),6); # max stateindex
    diff_arr = real.(reference_arr[1:simax] .- energies_arr[1:simax])
    #println("Comparison:")
    @printf("%-7s %-15s %-15s %-15s\n", "Index",  "Reference", "Test", "Difference")
    for i in 1:simax
        @printf("%-7d %-15.6f %-15.6f %-15.6f\n", i, reference_arr[i], real(energies_arr[i]), diff_arr[i])
    end
end

# testing gaussopt:
#println("\ngaussopt = 0, csm_bool = 0:")
#reference_arr = ISGL.ISGL_solve(phys_params,num_params,observ_params,wf_bool);

println("\ngaussopt = 1, csm_bool = 0:")
energies_arr = ISGL.ISGL_solve(phys_params,num_params,observ_params,wf_bool;gaussopt=gaussopt3);
comparison(reference_arr,energies_arr)


return nothing

