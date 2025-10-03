## Function/Module for using the Gaussian expansion method (GEM) for solving three-body problems in 1D(!)
## the code might seem unnecessary complicated at places. this is due to reusing the structure from the 3D ISGL code.


module GEM3B1D

using .. FewBodyToolkit
using SpecialFunctions, QuadGK, LinearAlgebra, Optim, StaticArrays, Roots, Interpolations, OffsetArrays
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
    GEM3B1D_solve(phys_params, num_params; wf_bool=0, csm_bool=0, observ_params=(;stateindices=[],centobs_arr=[[],[],[]],R2_arr=[0,0,0]))

Solves the 1D three-body problem using the Gaussian Expansion Method (GEM).

# Arguments
- `phys_params`: Physical parameters for the three-body system (e.g., masses, interaction potentials, etc.).
- `num_params`: Numerical parameters for the GEM calculation (e.g., basis size, grid parameters, etc.).
- `wf_bool`: (optional) If `1`, also returns wavefunction-related observables. Default is `0`.
- `csm_bool`: (optional) If `1`, uses complex scaling method. Default is `0`.
- `observ_params`: (optional) Parameters for observable calculations. Currently not supported for 1D.

# Returns
- If `wf_bool == 0`: Returns an array of computed energies.
- If `wf_bool == 1`: Returns a tuple `(energies, wavefunctions)`.

# Example
```julia
phys_params = make_phys_params3B1D()
num_params = make_num_params3B1D()
energies = GEM3B3D_solve(phys_params, num_params) #solving with default parameters: three particles with the same mass and gaussian interaction
```
"""
function GEM3B1D_solve(phys_params, num_params; wf_bool=0, csm_bool=0, observ_params=(;stateindices=[],centobs_arr=[[],[],[]],R2_arr=[0,0,0]), debug_bool=0)
    
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
    size_params = size_estimate(phys_params,num_params,observ_params)
    
    ## 4. preallocation:
    precomp_arrs,interpol_arrs,fill_arrs,result_arrs = preallocate_data(phys_params,num_params,observ_params,size_params,csm_bool)

    ## 5. precomputation:
    precompute_3B1D(phys_params,num_params,size_params,precomp_arrs)

    ## 6. preparation of interpolation & shoulder:
    interpolNshoulder(phys_params,num_params,observ_params,size_params,precomp_arrs,interpol_arrs,wf_bool,csm_bool)
    
    ## 7. Calculation of matrix elements
    fill_TVS(num_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,csm_bool,phys_params.hbar,debug_bool)
    
    ## 8. Solving the generalized eigenproblem:
    solveHS(num_params,fill_arrs,result_arrs,wf_bool)
    
    ## 9. Calculate observables: (currently not supported for 1D)
    if wf_bool == 0
        return result_arrs.energies_arr
    elseif wf_bool == 1
        #calc_observables(num_params,observ_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,result_arrs)
        return result_arrs.energies_arr,result_arrs.wavefun_arr
    end
end


end ## end of module