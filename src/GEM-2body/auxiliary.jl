# collection of auxiliary functions for module GEM2B

# functions for creating inputs
function make_phys_params2B(;hbar = 1.0, mur=1.0, vint_arr=[[]], lmin=0, lmax=0, dim=3)
    return (;hbar, mur, vint_arr, lmin, lmax, dim)
end

function make_num_params2B(; gem_params=(nmax=5, r1=1.0, rnmax=10.0),theta_csm=0.0, omega_cr=0.5, threshold=10^-10)
    return (;gem_params, theta_csm, omega_cr, threshold)
end

# Structure which contains all preallocated arrays for the GEM-2B solver.
struct PreallocStruct2B
    nu_arr::Vector{T} where T
    S::Matrix{TS} where TS
    T::Matrix{T} where T
    V::Matrix{T} where T
    energies::Vector{TE} where TE
    wavefunctions::Matrix{T} where T
    
    function PreallocStruct2B(num_params, cr_bool, csm_bool)        
        nbasis = num_params.gem_params.nmax
        
        # Determine types: TTV = type of T and V matrix; TS = type of S matrix; TE = type of energies array
        if cr_bool == 0 && csm_bool == 0
            TTV = Float64; TS = Float64; TE = Float64
        elseif cr_bool == 1 && csm_bool == 0
            TTV = ComplexF64; TS = ComplexF64; TE = Float64; nbasis *= 2;        
        elseif cr_bool == 0 && csm_bool == 1
            TTV = ComplexF64; TS = Float64; TE = ComplexF64
        elseif cr_bool == 1 && csm_bool == 1
            TTV = ComplexF64; TS = ComplexF64; TE = ComplexF64; nbasis *= 2;
        else
            error("Error with (cr_bool = $cr_bool, csm_bool = $csm_bool). Only (0,0), (1,0), (0,1), (1,1) allowed.")
        end
        
        # Initialize arrays
        nu_arr = zeros(TTV, nbasis)
        S  = zeros(TS, nbasis, nbasis)
        T  = zeros(TTV, nbasis, nbasis)
        V  = zeros(TTV, nbasis, nbasis)
        energies = zeros(TE, nbasis)
        wavefunctions = zeros(TTV, nbasis, nbasis)
        
        new(nu_arr, S, T, V, energies, wavefunctions) # creates one instance of the struct. only works within an inner "constructor" function
    end
end
