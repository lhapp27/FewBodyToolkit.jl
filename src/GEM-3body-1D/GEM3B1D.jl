## Function/Module for using the Gaussian expansion method (GEM) for solving three-body problems in 1D(!)
## the code might seem unnecessary complicated at places. this is due to reusing the structure from the 3D ISGL code.


module GEM3B1D

using .. FewBodyToolkit
using ..FewBodyToolkit: parse_complex_ranged, hermitian_fill!, theta_mesh, interpol_lookup, check_cr_csm_sector # shared with ISGL, see common/auxiliary.jl
using SpecialFunctions, QuadGK, LinearAlgebra, StaticArrays, Roots, Interpolations, OffsetArrays
using Printf: @printf

include("auxiliary.jl")
include("sanitycheck.jl")
include("size_estimate.jl")
include("preallocate.jl")
include("precomputation.jl")
include("interpolationNshoulder.jl")
include("fillTVS.jl")
include("solveHS.jl")
include("../common/eigen2step.jl")

export GEM3B1D_solve
export make_phys_params3B1D,make_num_params3B1D

"""
    GEM3B1D_solve(phys_params, num_params; return_wavefunctions=false, complex_scaling=false, complex_ranged=:none, observ_params=(;stateindices=[],centobs_arr=[[],[],[]],R2_arr=[0,0,0]))

Solves the 1D three-body problem using the Gaussian Expansion Method (GEM).

# Arguments
- `phys_params`: Physical parameters for the three-body system (e.g., masses, interaction potentials, etc.).
- `num_params`: Numerical parameters for the GEM calculation (e.g., basis size, grid parameters, etc.).
- `return_wavefunctions`: (optional) If `true`, also returns wavefunction-related observables. Default is `false`.
- `complex_scaling`: (optional) If `true`, uses complex scaling method. Default is `false`.
- `complex_ranged`: (optional) Selects complex-ranged basis functions `nu -> nu*(1 + i*omega)` per Jacobi
  coordinate, with `omega = complex_range_freq` from `num_params`. One of `:none` (default), `:r`
  (only the internal coordinate `` r ``), `:R` (only `` R ``), or `:both`. A `Bool` is also accepted,
  with `true` meaning `:both`, for consistency with `GEM2B_solve`.
  Each complex-ranged coordinate doubles its number of basis functions, since the ranges enter as
  conjugate pairs; choosing `:both` therefore multiplies the total basis size by four.
  All interaction types are supported, including generic central potentials: their radial integrals
  are interpolated over the sector of the complex plane in which the effective Gaussian range lives,
  with `kmax_theta` angular nodes (see `make_num_params3B1D`).
- `observ_params`: (optional) Parameters for observable calculations. Currently not supported for 1D.

# Returns
- If `return_wavefunctions=false`: Returns an array of computed energies.
- If `return_wavefunctions=true`: Returns a tuple `(energies, wavefunctions)`.

# Example
```julia
phys_params = make_phys_params3B1D()
num_params = make_num_params3B1D()
energies = GEM3B3D_solve(phys_params, num_params) #solving with default parameters: three particles with the same mass and gaussian interaction
```
"""
function GEM3B1D_solve(phys_params, num_params;
    return_wavefunctions=false, complex_scaling=false, complex_ranged=:none, observ_params=(;stateindices=[],centobs_arr=[[],[],[]],R2_arr=[0,0,0]), debug=false,
    wf_bool=nothing, csm_bool=nothing, debug_bool=nothing)

    complex_ranged_r, complex_ranged_R = parse_complex_ranged(complex_ranged)
    cr_any = complex_ranged_r || complex_ranged_R

    if !isnothing(wf_bool)
        @warn "wf_bool is deprecated, use return_wavefunctions instead"
        return_wavefunctions = wf_bool
    end
    if !isnothing(csm_bool)
        @warn "csm_bool is deprecated, use complex_scaling instead"
        complex_scaling = csm_bool
    end
    if !isnothing(debug_bool)
        @warn "debug_bool is deprecated, use debug instead"
        debug = debug_bool
    end
    
    ## 1. interpretation of inputs
    show_details_bool = 0
    if show_details_bool == 1
        println("Inputs:")
        @show(phys_params)
        @show(num_params)
        println("")
    end
    
    ## 2. sanity checks:
    error_code = sanity_checks3B(phys_params);
    if error_code != 0
        error("Program stopped due to erroneous inputs. Error code: $error_code")
    end

    if observ_params != (;stateindices=[],centobs_arr=[[],[],[]],R2_arr=[0,0,0])
        error("Observables are currently not supported for 1D calculations.")
    end
    
    ## 3. computations to determine sizes of arrays for allocation:   
    size_params = size_estimate(phys_params,num_params,observ_params,complex_ranged_r,complex_ranged_R)
    
    ## 4. preallocation:
    precomp_arrs,interpol_arrs,fill_arrs,result_arrs = preallocate_data(phys_params,num_params,observ_params,size_params,complex_scaling,complex_ranged_r,complex_ranged_R)

    ## 5. precomputation:
    precompute_3B1D(phys_params,num_params,size_params,precomp_arrs,complex_ranged_r,complex_ranged_R)

    ## 6. preparation of interpolation & shoulder:
    interpolNshoulder(phys_params,num_params,observ_params,size_params,precomp_arrs,interpol_arrs,return_wavefunctions,complex_scaling,cr_any)
    
    ## 7. Calculation of matrix elements
    fill_TVS(num_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,complex_scaling,phys_params.hbar,debug,complex_ranged_r,complex_ranged_R)
    
    ## 8. Solving the generalized eigenproblem:
    solveHS(num_params,fill_arrs,result_arrs,return_wavefunctions)
    
    ## 9. Calculate observables: (currently not supported for 1D)
    if !return_wavefunctions
        return result_arrs.energies_arr
    elseif return_wavefunctions
        #calc_observables(num_params,observ_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,result_arrs)
        return result_arrs.energies_arr,result_arrs.wavefun_arr
    end
end


end ## end of module