using Printf, BenchmarkTools

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
phys_params = (;mass_arr=[1.0,1.0,1.0],svals=["x","y","z"],vint_arr=[[v_gauss],[v_gauss],[v_gauss]],J_tot=0,parity=+1,spin_arr=[0,0,0])
phys_paramsg = (;mass_arr=[1.0,1.0,1.0],svals=["x","y","z"],vint_arr=[[],[],[]],J_tot=0,parity=+1,spin_arr=[0,0,0]) # for test of gaussopt

# numerical parameters:
gp = (;nmax=10,Nmax=10,r1=0.1,rnmax=25.0,R1=0.1,RNmax=25.0)
num_params = (;lmax=0,Lmax=0,gem_params=gp,theta_csm=0,omega_cr=0,mu0=0.08,c_shoulder=1.6,kmax_interpol=5000,threshold=10^-10,lmin=0,Lmin=0)
gb=1
gaussopt=[[[gb,v0,mu_g]],[[gb,v0,mu_g]],[[gb,v0,mu_g]]]

# Booleans:
wf_bool=0 # 0 -> dont calculate wave functions
csm_bool = 0

# Observables:
stateindices = 1:3#[1,2]
observ_params = (stateindices,centobs_arr = [[rad,rad2],[rad,rad2],[rad,rad2]],R2_arr = [0,0,0])

## Calling the program:
function comparison(reference_arr,energies_arr)
    simax = min(lastindex(energies_arr),6); # max stateindex
    diff_arr = real.(reference_arr .- energies_arr)
    #println("Comparison:")
    @printf("%-7s %-15s %-15s %-15s\n", "Index",  "Reference", "Test", "Difference")
    for i in 1:simax
        @printf("%-7d %-15.6f %-15.6f %-15.6f\n", i, reference_arr[i], real(energies_arr[i]), diff_arr[i])
    end
end

# testing gaussopt:
println("\ngaussopt = 0, csm_bool = 0:")
reference_arr = ISGL.ISGL_solve(phys_params,num_params,observ_params,wf_bool);

println("\ngaussopt = 1, csm_bool = 0:")
energies_arr = ISGL.ISGL_solve(phys_paramsg,num_params,observ_params,wf_bool;gaussopt=gaussopt);
comparison(reference_arr,energies_arr)

println("\ncsm_bool = 1, theta_csm = $(num_params.theta_csm):")
energies_arr = ISGL.ISGL_solve(phys_params,num_params,observ_params,wf_bool;csm_bool=1);
comparison(reference_arr,energies_arr)

println("\ngaussopt=1, csm_bool = 1, theta_csm = $(num_params.theta_csm):")
energies_arr = ISGL.ISGL_solve(phys_paramsg,num_params,observ_params,wf_bool;csm_bool=1,gaussopt=gaussopt);
comparison(reference_arr,energies_arr)

num_params = (;lmax=0,Lmax=0,gem_params=gp,theta_csm=10,omega_cr=0,mu0=0.08,c_shoulder=1.6,kmax_interpol=5000,threshold=10^-10,lmin=0,Lmin=0)
println("\ncsm_bool = 1, theta_csm = $(num_params.theta_csm): only real part: differences expected")
energies_arr = ISGL.ISGL_solve(phys_params,num_params,observ_params,wf_bool;csm_bool=1);
comparison(reference_arr,energies_arr)

println("\ngaussopt=1, csm_bool = 1, theta_csm = $(num_params.theta_csm): only real part: differences expected")
energies_arr = ISGL.ISGL_solve(phys_paramsg,num_params,observ_params,wf_bool;csm_bool=1,gaussopt=gaussopt);
comparison(reference_arr,energies_arr)


#4. Benchmark to check if code became slower
bench_bool = 1
if bench_bool == 1
    println("\n Benchmarking:")
    println("nmax = $(gp.nmax), Nmax = $(gp.Nmax)") # nmax = 10, Nmax = 10
    print("gaussopt = 0, csm_bool = 0: ")
    @btime ISGL.ISGL_solve(phys_params,num_params,observ_params,wf_bool); #893.461 ms (36241001 allocations: 590.97 MiB)
    print("\ngaussopt = 1, csm_bool = 0: ")
    @btime energies_arr = ISGL.ISGL_solve(phys_paramsg,num_params,observ_params,wf_bool;gaussopt=gaussopt); #108.708 ms (419 allocations: 17.22 MiB)
    print("\ngaussopt = 1, csm_bool = 1: ")
    @btime energies_arr = ISGL.ISGL_solve(phys_paramsg,num_params,observ_params,wf_bool;gaussopt=gaussopt, csm_bool = 1); #139.399 ms (423 allocations: 23.48 MiB)
    print("\ngaussopt = 0, csm_bool = 1: ")
    @btime ISGL.ISGL_solve(phys_params,num_params,observ_params,wf_bool;csm_bool = 1); #864.556 ms (36240985 allocations: 590.97 MiB)
end

# after first spin-orbit inclusion (13.05.25): code has become much slower, almost factor 2!

return nothing

