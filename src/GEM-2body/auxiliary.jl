# collection of auxiliary functions for module GEM2B

# functions for creating inputs

"""
    make_phys_params2B(; hbar=1.0, mur=1.0, vint_arr=[[GaussianPotential(-1.0, 1.0)]], lmin=0, lmax=0, dim=3)

Create and return a named tuple containing the physical parameters for a two-body system.

# Keyword arguments
- `hbar::Float64 = 1.0`: Reduced Planck constant used in calculations.
- `mur::Float64 = 1.0`: Reduced mass of the two-body system.
- `vint_arr::Vector{Any} = [[GaussianPotential(-1.0, 1.0)]]`: Array of interaction potentials or related parameters.
- `lmin::Int = 0`: Minimum orbital angular momentum quantum number.
- `lmax::Int = 0`: Maximum orbital angular momentum quantum number.
- `dim::Int = 3`: Spatial dimension of the system.

# Returns
- `NamedTuple`: Named tuple with the specified physical parameters.

# Example
```julia
make_phys_params2B(vint_arr=[GaussianPotential(-1.0, 0.5)], dim=1) # 1D system with Gaussian potential
make_phys_params2B(mur=0.5, vint_arr=[r -> -1/r], lmax=2)       # 3D Coulomb potential in d-wave (l=2) with reduced mass 0.5
```
"""
function make_phys_params2B(;hbar = 1.0, mur=1.0, vint_arr=[[GaussianPotential(-1.0, 1.0)]], lmin=0, lmax=0, dim=3)
    return (;hbar, mur, vint_arr, lmin, lmax, dim)
end

"""
    make_num_params2B(; gem_params=(nmax=5, r1=1.0, rnmax=10.0), theta_csm=0.0, omega_cr=0.9, threshold=1e-8)

Create and return a named tuple containing the numerical parameters for a two-body GEM calculation.

# Keyword arguments
- `gem_params::NamedTuple = (nmax=5, r1=1.0, rnmax=10.0)`: Parameters for the Gaussian Expansion Method (number of basis functions, smallest and largest range parameters).
- `theta_csm::Float64 = 0.0`: Complex scaling angle (in radians) for the Complex Scaling Method.
- `omega_cr::Float64 = 0.9`: Parameter controlling the frequency for complex-ranged basis functions.
- `threshold::Float64 = 1e-8`: Numerical threshold generalized eigenvalue solver.

# Returns
- `NamedTuple`: Named tuple with the specified numerical parameters.

# Example
```julia
make_num_params2B(gem_params=(nmax=20, r1=0.1, rnmax=50.0)) # for a larger basis set
make_num_params2B(gem_params=(nmax=20, r1=0.1, rnmax=50.0), theta_csm = 10.0) # non-zero rotation angle for complex scaling method (CSM)
```
"""
function make_num_params2B(; gem_params=(nmax=5, r1=1.0, rnmax=10.0),theta_csm=0.0, omega_cr=0.9, threshold=10^-8)
    return (;gem_params, theta_csm, omega_cr, threshold)
end

"""
    PreallocStruct2B{TTV, TS, TE}

A structure that holds preallocated arrays for the GEM-2B solver, used to store basis parameters, matrices, and results for two-body calculations.

# Fields
- `nu_arr::Vector{TTV}`: Array of nonlinear variational parameters (basis exponents).
- `S::Matrix{TS}`: Overlap matrix between basis functions.
- `T::Matrix{TTV}`: Kinetic energy matrix.
- `V::Matrix{TTV}`: Potential energy matrix.
- `energies::Vector{TE}`: Array to store computed eigenvalues (energies).
- `wavefunctions::Matrix{TTV}`: Matrix to store eigenvectors (wavefunctions).

# Type Parameters
- `TTV`: Element type for kinetic, potential, and basis parameter arrays (e.g., `Float64` or `ComplexF64`).
- `TS`: Element type for the overlap matrix (e.g., `Float64` or `ComplexF64`).
- `TE`: Element type for the energies array (e.g., `Float64` or `ComplexF64`).

# Description
This struct is designed to minimize memory allocations and improve performance by reusing arrays during repeated GEM-2B calculations. The types and sizes of the arrays are determined by the numerical parameters and whether complex rotation or complex scaling is used.

# Keyword arguments
- `num_params`: A named tuple containing numerical parameters, including the maximum number of basis functions (`nmax`).
- `cr_bool`: Boolean indicating whether complex range basis functions are used (1 for true, 0 for false).
- `csm_bool`: Boolean indicating whether complex scaling is used (1 for true, 0 for false).

# Example
```julia
PreallocStruct2B(num_params, cr_bool=0, csm_bool=0) # for real basis functions and no complex scaling
PreallocStruct2B(num_params, cr_bool=1, csm_bool=0) # for complex range basis functions and no complex scaling
PreallocStruct2B(num_params, cr_bool=0, csm_bool=1) # for real basis functions with complex scaling
```
"""
struct PreallocStruct2B{TTV, TS, TE}
    nu_arr::Vector{TTV}
    S::Matrix{TS}
    T::Matrix{TTV}
    V::Matrix{TTV}
    energies::Vector{TE}
    wavefunctions::Matrix{TTV}

    function PreallocStruct2B(num_params, cr_bool, csm_bool)
        nbasis = num_params.gem_params.nmax

        # Determine types: TTV = type of T and V matrix; TS = type of S matrix; TE = type of energies array
        if cr_bool == 0 && csm_bool == 0
            TTV = Float64;      TS = Float64;       TE = Float64
        elseif cr_bool == 1 && csm_bool == 0
            TTV = ComplexF64;   TS = ComplexF64;    TE = Float64        
        elseif cr_bool == 0 && csm_bool == 1
            TTV = ComplexF64;   TS = Float64;       TE = ComplexF64
        elseif cr_bool == 1 && csm_bool == 1
            TTV = ComplexF64;   TS = ComplexF64;    TE = ComplexF64
        else
            error("Error with (cr_bool = $cr_bool, csm_bool = $csm_bool). Only (0,0), (1,0), (0,1), (1,1) allowed.")
        end

        # Double nbasis if cr_bool == 1
        if cr_bool == 1
            nbasis *= 2
        end

        nu_arr = zeros(TTV, nbasis)
        S  = zeros(TS, nbasis, nbasis)
        T  = zeros(TTV, nbasis, nbasis)
        V  = zeros(TTV, nbasis, nbasis)
        energies = zeros(TE, nbasis)
        wavefunctions = zeros(TTV, nbasis, nbasis)

        new{TTV, TS, TE}(nu_arr, S, T, V, energies, wavefunctions)
    end
end
