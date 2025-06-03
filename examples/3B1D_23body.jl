# # Consistency check with two-body module GEM2B

# In this example we check if we can reproduce two-body results using the three-body code. This is a good test for new features of the three-body code and whenever there are no other known results to compare to. For that we use only one interaction between particles 1 and 2, i.e. particle 3 remains a spectator. Moreover, we employ only a single basis function to describe the relative motion of particle 3 relative to the center of mass of particles 1 and 2. This ensures minimal impact on the two-body subsystem and hence reproduces the two-body results.

# ## Setup
using Printf, Plots, FewBodyToolkit.GEM3B1D, FewBodyToolkit.GEM2B


# ## Input parameters:
# #### Physical parameters

# We use only a single interaction between particles 1 and 2. Particle 3 therefore does not interact and acts as a spectator.

v2(r) = -40 * 1/(1+r^4)
vint_arr=[[],[],[v2]] #[[v23],[v31],[v12]]

mass_arr = [1.0,20.0,20.0]# array of masses of particles (m1,m2,m3)
mur = 1/(1/mass_arr[1]+1/mass_arr[2]) # reduced mass

phys_params2B = make_phys_params2B(;mur,vint_arr=vint_arr[3],dim=1)
phys_params3B = make_phys_params3B1D(;mass_arr,vint_arr)

# #### Numerical parameters

# For the \\( R \\)- Jacobi coordinate we use only a single basis function with a very large range. This ensures the contribution from the kinetic energy operator corresponding to this Jacobi-coordinate is negligible. 

nmax = 8; r1=1.0;rnmax=10.0;
num_params2B = make_num_params2B(;gem_params=(;nmax, r1, rnmax))
num_params3B = make_num_params3B1D(;gem_params=(;nmax, r1, rnmax, Nmax= 1, R1=10000.0, RNmax=10000.0))

# ## Helper: comparison function

# We define a utility to compare two sets of eigenvalues:
function comparison(num_arr,ex_arr,simax;s1="Numerical", s2="Exact")
    @printf("%-7s %-15s %-15s %-15s\n", "Index",  s1, s2, "Difference")
    for i in 1:simax
        @printf("%-7d %-15.6f %-15.6f %-15.6f\n", i, num_arr[i], ex_arr[i], ex_arr[i] - num_arr[i])
    end
end;

# ## Numerical solution
e2 = GEM2B.GEM2B_solve(phys_params2B,num_params2B);
e3 = GEM3B1D.GEM3B1D_solve(phys_params3B,num_params3B);

# Determine the number of bound states
simax = findlast(e2.<0)

# Since particle 3 only acts as a spectator, the three-body system is effectively a two-body system. We can compare the results of the two codes. The three-body code can reproduce the two-body results, as it should:
println("Check if we can recover 2-body results with the 3-body code:")
comparison(e2,e3,simax;s1="2-body",s2="3-body")
