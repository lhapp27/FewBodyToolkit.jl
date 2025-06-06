# # Three-body Coulomb problem in 3D: HD+ (proton, deuteron, electron)

# In this example we solve the three-body Coulomb problem of a proton, deuteron and electron in 3D using the ISGL module. The code calculates the bound states and observables like the mean radii between two of the three particles. The results are compared with the high-precision state-of-the-art literature values of [bubin2005](@cite).


# ## Setup
using Printf, FewBodyToolkit#.ISGL

# ## Input parameters:
# Define pair-interactions:
vde(r) = -1/r #deuteron-electron = V23 = V1
vep(r) = -1/r #electron-proton = V31 = V2
vpd(r) = +1/r; #proton-deuteron = V12 = V3

# Physical parameters
mass_arr=[1836.15267343,3670.48296788,1.00] # array of masses of particles (m1,m2,m3), here: proton, deuteron, electron
phys_params = make_phys_params3B3D(;mass_arr,vint_arr=[[vde],[vep],[vpd]]);

# numerical parameters:
gp = (;nmax=25,Nmax=25,r1=0.1,rnmax=25.0,R1=0.1,RNmax=25.0)
num_params = make_num_params3B3D(;gem_params=gp);

# In this example we also calculate central observables like the mean radii (and squared radii) between two of the three particles. With "central", we mean observables of any form that only depend on r, the direct distance between two particles. For this we define the observables as functions:
rad(r) = r # radius
rad2(r) = r^2; # squared radius

# The code also allows for the observable \\( \langle R^2\rangle \\) in the other Jacobi coordinate, via the input `R2_arr`.
stateindices = 1:3 # for which states to calculate observables
observ_params = (stateindices,centobs_arr = [[rad,rad2],[rad,rad2],[rad,rad2]],R2_arr = [1,1,1] ) # R2_arr=[0,0,0] means no R^2 calculation

# ## Obtaining energies and observables

# Using the optional keyword argument `observ_params` (together with `wf_bool=1`) we can obtain the results for the observables on-the-fly. The energies are stored in the `energies` array, the eigenvectors in `wfs`, and the central observables in `co_out`. The mean squared radii for the R-coordinate are stored in `R2_out`.
energies,wfs,co_out,R2_out = ISGL.ISGL_solve(phys_params,num_params;observ_params,wf_bool=1);

# ## Comparison with literature values

compmax = min(6,lastindex(energies));
num_arr = energies[1:compmax];
ex_arr = -[0.5978979685,0.5891818291,0.5809037001,0.5730505464,0.5656110418,0.5585755200][1:compmax]; # literature values for the energies

# #### Energy spectrum
println("HD+:  Energies:")
comparison(num_arr, ex_arr, compmax; s1="Numeric", s2="Literature")
println("---------------------------------------------------")

# We find good agreement with the values of the literature, with deviations between 10^-3 and a few percent.

# #### Mean radii, and squared radii for the r-coordinate

# Literature values for radii (`dp` indicates deuteron-proton distance, `de` for deuteron-electron, and `pe` for proton-electron):
rdp_lit = [2.055,2.171,2.292,2.417,2.547,2.683][stateindices]
rde_lit = [1.688,1.750,1.813,1.880,1.948,2.020][stateindices]
rpe_lit = [1.688,1.750,1.814,1.881,1.950,2.022][stateindices]
r_lit = [rde_lit,rpe_lit,rdp_lit];

# ... and square radii:
rdp2_lit = [4.268,4.855,5.492,6.185,6.942,7.771][stateindices]
rde2_lit = [3.534,3.839,4.169,4.526,4.915,5.339][stateindices]
rpe2_lit = [3.537,3.843,4.173,4.531,4.921,5.346][stateindices]
r2_lit = [rde2_lit,rpe2_lit,rdp2_lit];

# Comparison of the numerical results with the literature:
println("\nHD+:  radii ⟨r⟩ with differences")
strings = ["r_de", "r_pe", "r_dp"]
for ii in [3,1,2]
    comparison(co_out[ii,1,stateindices], r_lit[ii], 3; s1=strings[ii], s2=string(strings[ii],"(lit)"))
end
println("---------------------------------------------------")

println("\nHD+:  squared radii ⟨r²⟩ with differences")
strings2 = ["r²_de", "r²_pe", "r²_dp"]
for ii in [3,1,2]
    comparison(co_out[ii,2,stateindices], r2_lit[ii], 3; s1=strings2[ii], s2=string(strings2[ii],"(lit)"))
end
println("---------------------------------------------------")

# For the mean radii, the results for the ground state are still quite good, however accuracy decreases quickly for higher excited states. The repulsive Coulomb interaction between proton and deuteron poses a difficult problem for the centered Gaussian basis functions. Accordingly, the largest deviations are found for the (squared) radii of the proton-deuteron pair. Other types of basis functions might prove more suitable to capture the details of the system.

# #### Mean squared radii for the R-coordinate

# For this observable the reference does not provide any values, so we just print the numerical results:
println("\nHD+:  Mean squared radii ⟨R²⟩")
row_labels = ["Proton rel. to (D+,e-) pair:  ", "Deuteron rel. to (e-,p+) pair:", "Electron rel. to (p+,D+) pair:"]
state_labels = ["State 1", "State 2", "State 3"]
println(rpad(" ⟨R²⟩ ", 34), join(state_labels, "   "))
for i in 1:length(row_labels)
    print(rpad(row_labels[i], 30))
    for j in 1:length(state_labels)
        print(@sprintf("%10.4f", R2_out[i, j]))
    end
    println()
end

# Due to the almost negligible mass of the electron, the results for ⟨R²⟩ in the first two rows are almost identical and also almost match to the first column of the mean radii ⟨r²⟩ between proton and deuteron. The third row shows the results for the electron with respect to the center-of-ass of the proton-deuteron pair and are clearly smaller, since it feels an attractive Coulomb force of twice the strength (the effective pair of proton and deuteron has charge +2).