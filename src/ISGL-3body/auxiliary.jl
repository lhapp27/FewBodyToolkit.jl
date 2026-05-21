# file for auxiliary functions
# functions for creating inputs

"""
    make_phys_params3B3D(; hbar=1.0, masses=[1.0,1.0,1.0], species=[:x,:y,:z], interactions=[[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)]], J_tot=0, parity=+1, spins=[0,0,0])

Create and return a named tuple containing the physical parameters for a three-body system in 3D.

# Keyword arguments
- `hbar::Float64 = 1.0`: Reduced Planck constant used in calculations.
- `masses::Vector{Float64} = [1.0,1.0,1.0]`: 3-element vector containing the masses of the three particles.
- `species::Vector = [:x,:y,:z]`: 3-element vector containing particle labels used for automatic (anti-)symmetrization. Use Symbols (e.g. `:b`, `:f`).
- `interactions::Vector{Vector{Any}} = [[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)]]`: Vector of Vector of interaction potentials for each pair of particles: [[v23_1,v23_2,...],[v31_1,v31_2,...],[v12_1,v12_2,...]].
- `J_tot::Int = 0`: Total angular momentum of the system. Can take on half-integer values when spins are involved.
- `parity::Int = 1`: Parity `parity=(-1)^(l+L)` of the wave function. Possible values: +1,-1 for positive/negative parity, 0 for parity-violating potentials.
- `spins::Vector{Float64} = [0,0,0]`: 3-element vector containing the spins [z_1,z_2,z_3] of the three particles. Used for automatic (anti-)symmetrization.

# Returns
- `NamedTuple`: Named tuple with the specified physical parameters.

# Example
```julia
make_phys_params3B3D() # default: system of three different particles with the same mass and Gaussian interactions
make_phys_params3B3D(masses=[10.0,10.0,20.0], species=[:i,:j,:k], interactions=[[v23],[v31],[v12_1,v12_2]], J_tot=3/2, spins=[1/2,1/2,1/2],parity=-1, lmax=1, Lmax=1) # system of three different particles with half-integer spin, and interactions v23 (between particles 2 and 3), v31 (between particles 3 and 1), and v12_1, v12_2 (between particles 1 and 2). The interactions need to be defined before this call.
```
"""
function normalize_species(species)
    if any(s -> s isa String, species)
        @warn "String entries in `species` are deprecated, use Symbols instead (e.g. :b instead of \"b\")"
        return map(s -> s isa String ? Symbol(s) : s, species)
    end
    return species
end

function make_phys_params3B3D(;hbar = 1.0, masses=[1.0, 1.0, 1.0], species=[:x,:y,:z], interactions=[[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)]], J_tot=0, parity=1, spins=[0,0,0])
    species = normalize_species(species)
    return (;hbar, masses, species, interactions, J_tot, parity, spins)
end

"""
    make_num_params3B3D(;lmin=0, Lmin=0, lmax=0, Lmax=0, gem_params=(nmax=10, r1=0.2, rnmax=10.0, Nmax=10, R1=0.2, RNmax=20.0), complex_scaling_angle=0.0, complex_range_freq=0.9, mu0=0.08, c_shoulder=1.6, kmax_interpol=1000, threshold=10^-8)

Create and return a named tuple containing the numerical parameters for a three-body GEM calculation in 3D.

# Keyword arguments
- `lmin::Int = 0`: Minimum power `r^l` used in the basis functions of the `` r `` Jacobi coordinate.	
- `Lmin::Int = 0`: Minimum power `r^L` used in the basis functions of the `` R `` Jacobi coordinate.
- `lmax::Int = 0`: Maximum power `r^l` used in the basis functions of the `` r `` Jacobi coordinate.
- `Lmax::Int = 0`: Maximum power `r^L` used in the basis functions of the `` R `` Jacobi coordinate.
- `gem_params::NamedTuple = (nmax=5, r1=1.0, rnmax=10.0, Nmax=5, R1=1.0, RNmax=10.0)`: Parameters for the Gaussian Expansion Method (number of basis functions, smallest and largest range parameters for both Jacobi coordinates).
- `complex_scaling_angle::Float64 = 0.0`: Complex scaling angle (in degrees) for the Complex Scaling Method.
- `complex_range_freq::Float64 = 0.9`: Parameter controlling the frequency for complex-ranged basis functions. Currently unsupported.
- `mu0::Float64 = 0.08`: Parameter (prefactor) for the ISGL method.
- `c_shoulder::Float64 = 1.6`: Parameter (base) for the ISGL method.
- `kmax_interpol::Int = 1000`: Number of numerical integration with effective Gaussian ranges used for interpolation.
- `threshold::Float64 = 1e-8`: Numerical threshold for the generalized eigenvalue solver.

# Returns
- `NamedTuple`: Named tuple with the specified numerical parameters.

# Example
```julia
make_num_params3B3D() # for the default basis set
make_num_params3B3D(gem_params=(nmax=10, r1=0.5, rnmax=20.0, Nmax=15, R1=1.0, RNmax=100.0), kmax_interpol=5000) for a larger basis set, and a finer interpolation grid for matrix-element calculations of interactions with many features.
```
"""
function make_num_params3B3D(;lmin=0, Lmin=0, lmax=0, Lmax=0, gem_params=(nmax=10, Nmax=10, r1=0.2, rnmax=20.0, R1=0.2, RNmax=20.0),complex_scaling_angle=0.0, complex_range_freq=0.9, mu0=0.08, c_shoulder=1.6, kmax_interpol=1000, threshold=10^-8, theta_csm=nothing, omega_cr=nothing)
    complex_scaling_angle = isnothing(theta_csm) ? complex_scaling_angle : theta_csm
    complex_range_freq = isnothing(omega_cr) ? complex_range_freq : omega_cr
    return (;lmax, Lmax, gem_params, complex_scaling_angle, complex_range_freq, mu0, c_shoulder, kmax_interpol, threshold, lmin, Lmin)
end