# file for auxiliary functions

# functions for creating inputs
function make_phys_params3B3D(;hbar = 1.0, mass_arr=[1.0, 1.0, 1.0], svals=["x","y","z"], vint_arr=[[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)]], J_tot=0, parity=1, spin_arr=[0,0,0])
    return (;hbar, mass_arr, svals, vint_arr, J_tot, parity, spin_arr)
end

function make_num_params3B3D(; lmax=0, Lmax=0, gem_params=(nmax=3, Nmax=3, r1=1.0, rnmax=20.0, R1=1.0, RNmax=20.0),theta_csm=0.0, omega_cr=0.5, mu0=0.08, c_shoulder=1.6, kmax_interpol=1000, threshold=10^-8, lmin=0, Lmin=0)
    return (;lmax, Lmax, gem_params, theta_csm, omega_cr, mu0, c_shoulder, kmax_interpol, threshold, lmin, Lmin)
end