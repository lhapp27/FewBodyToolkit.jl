# # Three-body Coulomb problem in 3D: HD+ (proton, deuteron, electron)

# In this example we solve the three-body Coulomb problem of a proton, deuteron and electron in 3D using the ISGL module. The code calculates the bound states and observables like the mean radii between two of the three particles. The results are compared with the high-precision state-of-the-art literature values of https://doi.org/10.1063/1.1850905.


# ## Setup
using Printf, FewBodyToolkit.ISGL

# ## Input parameters:
# Define pair-interactions:
vde(r) = -1/r #deuteron-electron = V23 = V1
vep(r) = -1/r #electron-proton = V31 = V2
vpd(r) = +1/r #proton-deuteron = V12 = V3

# Physical parameters
mass_arr=[1836.15267343,3670.48296788,1.00] # array of masses of particles (m1,m2,m3), here: proton, deuteron, electron
phys_params = make_phys_params3B3D(;mass_arr,vint_arr=[[vde],[vep],[vpd]])

# numerical parameters:
gp = (;nmax=25,Nmax=25,r1=0.1,rnmax=25.0,R1=0.1,RNmax=25.0)
num_params = make_num_params3B3D(;gem_params=gp)

# In this example we also calculate central observables like the mean radii (and squared radii) between two of the three particles. With "central", we mean observables of any form that only depend on r, the direct distance between two particles. For this we define the observables as functions:
rad(r) = r # radius
rad2(r) = r^2 # squared radius

# The code also allows for the observable \\( \langle R^2\rangle \\) in the other Jacobi coordinate, via the input `R2_arr`.
stateindices = 1:3 # for which states to calculate observables
observ_params = (stateindices,centobs_arr = [[rad,rad2],[rad,rad2],[rad,rad2]],R2_arr = [1,1,1] ) # R2_arr=[0,0,0] means no R^2 calculation

# ## Obtaining energies and observables

# using the optional keyword arguement observ_params (together with wf_bool=1) we can obtain the results for the observables on-the-fly.
@time energies_arr,wf_arrs,co_out,R2_out = ISGL.ISGL_solve(phys_params,num_params;observ_params,wf_bool=1)

# ## Comparison with literature values

# Helper function:
function print_comparison_with_diff(label1, arr1, label2, arr2, indices)
    diff = arr1 .- arr2
    @printf("%-10s %-12s %-12s %-12s\n", "Index", label1, label2, "Diff")
    for (i, si) in enumerate(indices)
        @printf("%-10d %-12.3f %-12.3f %-12.3f\n", si, arr1[i], arr2[i], diff[i])
    end
    println("")
end

compmax = min(6,lastindex(energies_arr));
num_arr = energies_arr[1:compmax];
ex_arr = -[0.5978979685,0.5891818291,0.5809037001,0.5730505464,0.5656110418,0.5585755200][1:compmax]; # literature values for the energies

# #### Energy spectrum
println("HD+:  Energies:")
print_comparison_with_diff("Numeric", num_arr, "Literature", ex_arr, 1:compmax)
println("---------------------------------------------------")

# We find good agreement with the vaklues of the literature, with deviations between 10^-3 and a few 10^2.

# #### Mean radii \\( \langle r \rangle \\), and squared radii \\( \langle r^2 \rangle) \)

# radii:
rdp_lit = [2.055,2.171,2.292,2.417,2.547,2.683][stateindices];
rde_lit = [1.688,1.750,1.813,1.880,1.948,2.020][stateindices];
rpe_lit = [1.688,1.750,1.814,1.881,1.950,2.022][stateindices];
r_lit = [rde_lit,rpe_lit,rdp_lit]
# square radii:
rdp2_lit = [4.268,4.855,5.492,6.185,6.942,7.771][stateindices];
rde2_lit = [3.534,3.839,4.169,4.526,4.915,5.339][stateindices];
rpe2_lit = [3.537,3.843,4.173,4.531,4.921,5.346][stateindices];
r2_lit = [rde2_lit,rpe2_lit,rdp2_lit]


println("\nHD+:  radii  ⟨r⟩ with differences")
strings = ["r_de", "r_pe", "r_dp"]
for ii in [3,1,2]
    print_comparison_with_diff(strings[ii], co_out[ii,1,stateindices], string(strings[ii],"(lit)"), r_lit[ii], stateindices)
end
println("---------------------------------------------------")

println("\nHD+:  squared radii  ⟨r²⟩ with differences")
strings2 = ["r²_de", "r²_pe", "r²_dp"]
for ii in [3,1,2]
    print_comparison_with_diff(strings2[ii], co_out[ii,2,stateindices], string(strings2[ii],"(lit)"), r2_lit[ii], stateindices)
end
println("---------------------------------------------------")

# For the mean radii, the results for the ground state are still quite good, however accuracy decreases quickly for higher excited states. The repulsive Coulomb interaction between proton and deuteron poses a difficult problem for the centered Gaussian basis functions. The largest deviations are found for the (squared) radii of the proton-deuteron pair. Other basis functions might prove more suitable to capture the details of the system.

# #### Mean squared radii \\( \langle R^2 \rangle \\)

# For this observable the reference does not provide any values, so we just print the numerical results:
println("\nHD+:  Mean squared radii ⟨R²⟩")
row_labels = ["Proton rel. to (D+,e-) pair:  ", "Deuteron rel. to (e-,p+) pair:", "Electron rel. to (p+,D+) pair:"]
println(rpad(" ⟨R²⟩ ", 34), join(state_labels, "   "))
for i in 1:length(row_labels)
    print(rpad(row_labels[i], 30))
    for j in 1:length(state_labels)
        print(@sprintf("%10.4f", R2_out[i, j]))
    end
    println()
end

# Due to the almost negligible mass of the electron, the results for ⟨R²⟩ in the first two rows are almost identical and also almost match to the first column of the mean radii ⟨r²⟩ between proton and deuteron. The third row shows the results for the electron with respect to the center-of-ass of the proton-deuteron pair and are much smaller, since it feels an attractive Coulomb force of twice the strength (the effective pair of proton and deuteron has charge +2).