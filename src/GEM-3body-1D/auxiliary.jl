# file for auxiliary functions of module GEM3B1D
# functions for creating inputs

"""
    make_phys_params3B1D(; hbar=1.0, masses=[1.0,1.0,1.0], species=[:x,:y,:z], interactions=[[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)]], parity=+1)

Create and return a named tuple containing the physical parameters for a three-body system in 1D.

# Keyword arguments
- `hbar::Float64 = 1.0`: Reduced Planck constant used in calculations.
- `masses::Vector{Float64} = [1.0,1.0,1.0]`: 3-element vector containing the masses of the three particles.
- `species::Vector = [:x,:y,:z]`: 3-element vector containing particle labels used for automatic (anti-)symmetrization. Use Symbols `:b`, (`:f`) for identical bosons or fermions
- `interactions::Vector{Vector{Any}} = [[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)]]`: Vector of Vector of interaction potentials for each pair of particles: [[v23_1,v23_2,...],[v31_1,v31_2,...],[v12_1,v12_2,...]].
- `parity::Int = 1`: Parity `parity=(-1)^(l+L)` of the wave function. Possible values: +1,-1 for positive/negative parity, 0 for parity-violating potentials.

# Returns
- `NamedTuple`: Named tuple with the specified physical parameters.

# Example
```julia
make_phys_params3B1D() # default: system of three different particles with the same mass and Gaussian interactions
make_phys_params3B1D(masses=[1.0,10.0,20.0], species=[:i,:j,:k], interactions=[[v23],[v31],[v12_1,v12_2]]) # system of three different particles with interactions v23 (between particles 2 and 3), v31 (between particles 3 and 1), and v12_1, v12_2 (between particles 1 and 2). The interactions need to be defined before this call.
```
"""
function normalize_species(species)
    if any(s -> s isa String, species)
        @warn "String entries in `species` are deprecated, use Symbols instead (e.g. :b instead of \"b\")"
        return map(s -> s isa String ? Symbol(s) : s, species)
    end
    return species
end

function make_phys_params3B1D(;hbar = 1.0, masses=[1.0,1.0,1.0], species=[:x,:y,:z], interactions=[[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)],[GaussianPotential(-1.0, 1.0)]], parity=+1)
    species = normalize_species(species)
    return (;hbar, masses, species, interactions, parity)
end


"""
    make_num_params3B1D(;lmin=0, Lmin=0, lmax=0, Lmax=0, gem_params=(nmax=5, r1=1.0, rnmax=10.0, Nmax=5, R1=1.0, RNmax=10.0), complex_scaling_angle=0.0, complex_range_freq=0.9, kmax_interpol=1000, threshold=10^-8)

Create and return a named tuple containing the numerical parameters for a three-body GEM calculation in 1D.

# Keyword arguments
- `lmin::Int = 0`: Minimum power `r^l` used in the basis functions of the ``r`` Jacobi coordinate.	
- `Lmin::Int = 0`: Minimum power `r^L` used in the basis functions of the ``R`` Jacobi coordinate.
- `lmax::Int = 0`: Maximum power `r^l` used in the basis functions of the ``r`` Jacobi coordinate.
- `Lmax::Int = 0`: Maximum power `r^L` used in the basis functions of the ``R`` Jacobi coordinate.
- `gem_params::NamedTuple = (nmax=5, r1=1.0, rnmax=10.0, Nmax=5, R1=1.0, RNmax=10.0)`: Parameters for the Gaussian Expansion Method (number of basis functions, smallest and largest range parameters for both Jacobi coordinates).
- `complex_scaling_angle::Float64 = 0.0`: Complex scaling angle (in degrees) for the Complex Scaling Method.
- `complex_range_freq::Float64 = 0.9`: Parameter controlling the frequency for complex-ranged basis functions. Currently unsupported.
- `kmax_interpol::Int = 1000`: Number of numerical integration with effective Gaussian ranges used for interpolation.
- `threshold::Float64 = 1e-8`: Numerical threshold for the generalized eigenvalue solver.

# Returns
- `NamedTuple`: Named tuple with the specified numerical parameters.

# Example
```julia
make_num_params3B1D() # for the default basis set
make_num_params3B1D(gem_params=(nmax=10, r1=0.5, rnmax=20.0, Nmax=15, R1=1.0, RNmax=100.0), complex_scaling_angle=10.0) for a larger basis set and non-zero rotation angle for complex scaling method
```
"""
function make_num_params3B1D(;lmin=0,Lmin=0, lmax=0, Lmax=0, gem_params=(nmax=5, r1=1.0, rnmax=10.0, Nmax=5, R1=1.0, RNmax=10.0), complex_scaling_angle=0.0, complex_range_freq=0.9, kmax_interpol=1000, threshold=10^-8, theta_csm=nothing, omega_cr=nothing)
    complex_scaling_angle = isnothing(theta_csm) ? complex_scaling_angle : theta_csm
    complex_range_freq = isnothing(omega_cr) ? complex_range_freq : omega_cr
    return (;lmin, Lmin, lmax, Lmax, gem_params, complex_scaling_angle, complex_range_freq, kmax_interpol, threshold)
end
