```@meta
EditURL = "../../examples/ISGL_ps-.jl"
```

# Three-body Coulomb problem in 3D: ps- (electron, electron, positron)

In this example we solve the three-body Coulomb problem of two electrons and a positron in 3D using the ISGL module. The code calculates the bound states and observables like the mean radii between two of the three particles. The results are compared with the high-precision state-of-the-art literature values of [frolov1999](@cite).

## Setup

````@example ISGL_ps-
using Printf, FewBodyToolkit
````

## Input parameters:
Define pair-interactions:

````@example ISGL_ps-
vee(r) = +1/r #electron-electron: V12
vep(r) = -1/r #positron-electron: V31, V23
````

Physical parameters

````@example ISGL_ps-
mass_arr=[1.0,1.0,1.0] # array of masses of particles (m1,m2,m3), here: electron, electron, positron
svals=["b","b","z"] # since the two electrons are in asymmetric spin states, we can treat them as bosons for the spatial part
phys_params = make_phys_params3B3D(;mass_arr, svals, vint_arr=[[vep],[vep],[vee]]);
nothing #hide
````

numerical parameters:

````@example ISGL_ps-
gp = (;nmax=10,Nmax=10,r1=0.1,rnmax=25.0,R1=0.1,RNmax=25.0)
num_params = make_num_params3B3D(;gem_params=gp);
nothing #hide
````

In this example we also calculate central observables like the mean radii (and squared radii) between two of the three particles. With "central", we mean observables of any form that only depend on ``r``, the direct distance between two particles. For this we define the observables as functions:

````@example ISGL_ps-
rad(r) = r # radius
invrad(r) = 1/r # inverse radius
rad2(r) = r^2; # squared radius
nothing #hide
````

The code also allows for the observable `` \langle R^2\rangle `` in the other Jacobi coordinate, via the input `R2_arr`.

````@example ISGL_ps-
stateindices = [1] # for which states to calculate observables
observ_params = (stateindices,centobs_arr = [[rad,invrad,rad2],[rad,invrad,rad2],[rad,invrad,rad2]],R2_arr = [1,1,1] ) # R2_arr=[0,0,0] means no R^2 calculation
````

## Obtaining energies and observables

Using the optional keyword argument `observ_params` (together with `wf_bool=1`) we can obtain the results for the observables on-the-fly. The energies are stored in the `energies` array, the eigenvectors in `wfs`, and the central observables in `co_out`. The mean squared radii for the R-coordinate are stored in `R2_out`.

````@example ISGL_ps-
energies,wfs,co_out,R2_out = ISGL.ISGL_solve(phys_params,num_params;observ_params,wf_bool=1);
nothing #hide
````

## Comparison with literature values

````@example ISGL_ps-
simax = 1 # only 1 bound state
num_arr = energies[1:simax];
ex_arr = -[0.262005070232978][1:simax]; # literature value for the energy
nothing #hide
````

#### Energy spectrum

````@example ISGL_ps-
println("ps-:  binding energy:")
comparison(num_arr, ex_arr, simax; s1="Numeric", s2="Literature")
println("---------------------------------------------------")
````

We find good agreement with the literature value

#### Mean values for the radii, inverse radii, and squared radii for the r-coordinate

Literature values for radii (`ee` indicates the electron-electron distance, and `pe` the positron-electron distance):

````@example ISGL_ps-
ree_lit = [8.54858] #r21
rpe_lit = [5.48963] #r31
r_lit = [rpe_lit,rpe_lit,ree_lit];
nothing #hide
````

... inverse radii:

````@example ISGL_ps-
iree_lit = [0.15563] #1/r21
irpe_lit = [0.33982] #1/r31
ir_lit = [irpe_lit,irpe_lit,iree_lit];
nothing #hide
````

... and square radii:

````@example ISGL_ps-
ree2_lit = [93.1786] #r21^2
rpe2_lit = [48.4189] #r31^2
r2_lit = [rpe2_lit,rpe2_lit,ree2_lit];
nothing #hide
````

Comparison of the numerical results with the literature:

````@example ISGL_ps-
println("\nps-:  radii ⟨r⟩ with differences")
strings = ["r_pe","r_pe","r_ee"]
for ii in [1,2,3]
    comparison(co_out[ii,1,stateindices], r_lit[ii], simax; s1=strings[ii], s2=string(strings[ii],"(lit)"))
end
println("---------------------------------------------------")

println("\nps-:  inverse radii ⟨1/r⟩ with differences")
strings = ["1/r_pe","1/r_pe","1/r_ee"]
for ii in [1,2,3]
    comparison(co_out[ii,2,stateindices], ir_lit[ii], simax; s1=strings[ii], s2=string(strings[ii],"(lit)"))
end
println("---------------------------------------------------")

println("\nps-:  squared radii ⟨r²⟩ with differences")
strings2 = ["r²_pe","r²_pe","r²_ee"]
for ii in [1,2,3]
    comparison(co_out[ii,3,stateindices], r2_lit[ii], simax; s1=strings2[ii], s2=string(strings2[ii],"(lit)"))
end
println("---------------------------------------------------")
````

Since there is only one bound state, we can reproduce both energies and geometric properties with good accuracy. In contrast, the HD+ system supports several excited states, whose geometric properties are difficult to describe by the simple centered Gaussian basis functions.

#### Mean squared radii for the R-coordinate

For this observable the reference does not provide any direct comparison value, so we just print the numerical results:

````@example ISGL_ps-
println("\nps+:  Mean squared radii ⟨R²⟩")
row_labels = ["Electron 1 rel. to (e-_2,p+) pair:  ", "Electron 2 rel. to (p+,e-_1) pair:  ", "Positron rel. to (e-_1,e-_2) pair:  "]
state_labels = ["Ground state"]
println(rpad(" ⟨R²⟩ ", 37), join(state_labels, "   "))
for i in 1:length(row_labels)
    print(rpad(row_labels[i], 30))
    for j in 1:length(state_labels)
        print(@sprintf("%10.4f", R2_out[i, j]))
    end
    println()
end
````

Since particles 1 and 2 are identical electrons, the first two rows are identical. The third row shows the results for the positron with respect to the center-of-mass of the electron-electron pair and are clearly smaller, probably since it feels an attractive Coulomb force of twice the strength (the effective pair of the two electrons has charge -2).

## page References

```@bibliography
Pages = ["ISGL_ps-.md"]
Canonical = false
```

See also the [full bibliography](@ref References) for further references cited throughout this documentation.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

