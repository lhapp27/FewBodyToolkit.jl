# Script for testing the GEM-3body-1D code:

using Printf, Plots, FewBodyToolkit.GEM3B1D

## Define pair-interaction:
v0 = -3.0
mu_g = 1.0
vgauss(r) = v0*exp(-mu_g*r^2)
vgc = CentralPotential(vgauss) # why do we need this?
vgauss2 = GaussianPotential(v0,mu_g) # we need to be using FewBodyToolkit for this. would be better if FBTK.GEM3B1D would also work.

vtest(r) = v0 * 1/(1 + r^6)
vtest2 = CentralPotential(vtest)

## Input parameters:
# Physical parameters
mass_arr = [1.0,20.0,20.0]# array of masses of particles (m1,m2,m3)
vint_arr=[[],[vtest2],[vtest2]] #[[v23],[v31],[v12]]
phys_params = make_phys_params3B1D(;mass_arr,vint_arr)

# numerical parameters:
num_params = make_num_params3B1D(;gem_params=(nmax=5, Nmax= 5, r1=1.0, rnmax=10.0, R1=1.0, RNmax=10.0))

# Calling the program:
energies_arr = GEM3B1D.GEM3B1D_solve(phys_params,num_params);

## Print output:
energies_arr[1:min(5,lastindex(energies_arr))] # print first 10 energies;
