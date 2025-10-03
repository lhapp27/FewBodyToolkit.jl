## Function/Module for using the Infinitesimally shifted Gaussian lobe functions (ISGL) within the Gaussian expansion method (GEM) for solving three-body problems

## see bottom for a documentation ("concept of the program"); might be a bit outdated

## 22.01.24: Change to structs
## 21.02.24: Change to module
## 14.03.25: Spin-implementation


module ISGL

using .. FewBodyToolkit
using LinearAlgebra,StaticArrays,OffsetArrays,Interpolations, SpecialFunctions,QuadGK,PartialWaveFunctions, WignerSymbols
using Printf: @printf

include("sanitycheck.jl")
include("size_estimate.jl")
include("auxiliary.jl")
include("preallocate.jl")
include("precomputation.jl")
include("interpolationNshoulder.jl")
include("fillTVS.jl")
include("solveHS.jl")
include("../common/eigen2step.jl")
include("observables.jl")

export ISGL_solve, make_phys_params3B3D, make_num_params3B3D

# default parameters for the observables
const DEFAULT_OBS = (
    stateindices = Int[],                           # Vector{Int}
    centobs_arr  = [Vector{PotentialFunction}() for _ in 1:3],  # Vector{Vector{PotentialFunction}}
    R2_arr       = Int[0, 0, 0],                    # Vector{Int}
)

"""
    ISGL_solve(phys_params, num_params; wf_bool=0, csm_bool=0, observ_params=(;stateindices=[],centobs_arr=[[],[],[]],R2_arr=[0,0,0]))

Solves the 3D three-body problem using the Gaussian Expansion Method (GEM).

# Arguments
- `phys_params`: Physical parameters for the three-body system (e.g., masses, interaction potentials, etc.).
- `num_params`: Numerical parameters for the GEM calculation (e.g., basis size, grid parameters, etc.).
- `wf_bool`: (optional) If `1`, also returns wavefunction-related observables. Default is `0`.
- `csm_bool`: (optional) If `1`, uses complex scaling method. Default is `0`.
- `observ_params`: (optional) Parameters for observable calculations.
    + `stateindices`: Indices of states for which observables are calculated.
    + `centobs_arr`: Array of central (only dependent on `` r ``; must be defined as functions) observables, for each Jacobi set (similar to `vint_arr` in `phys_params`).
    + `R2_arr`: Array which indicates whether the observable `` \\langle R^2 \\rangle `` should be calculated (1) for any of the three Jacobi sets, or not (0).

# Returns
- If `wf_bool == 0`: Returns an array of computed energies.
- If `wf_bool == 1`: Returns a tuple `(energies, wavefunctions, centobs_output, R2_output)`.
    + `energies`: Vector of computed energies.
    + `wavefunctions`: Matrix of eigenvectors (column-wise) which contain the coefficients of the basis functions.
    + `centobs_output`: Mean values of central observables for the specified states. The first dimension corresponds to the Jacobi sets, the second to the observables, and the third to the states.
    + `R2_output`: Mean squared radii for the R-coordinate. Rows corresond to the Jacobi sets, columns to the states.

# Example
```julia
phys_params = make_phys_params3B3D()
num_params = make_num_params3B3D()
energies = ISGL_solve(phys_params, num_params) #solving with default parameters: three particles with the same mass and gaussian interaction
```
"""
function ISGL_solve(phys_params, num_params; wf_bool = 0, csm_bool = 0, observ_params=DEFAULT_OBS, debug_bool = 0)
    
    ## 1. interpretation of inputs
    show_details_bool = 0
    if show_details_bool == 1
        println("Inputs:")
        @show(phys_params)
        @show(num_params)
        println("")
    end
    
    ## 2. sanity checks:
    error_code = sanity_checks(phys_params);
    if error_code != 0
        println("Erroneous inputs. Program stopped")
        return
    end
    
    ## 3. computations to determine sizes of arrays for allocation:   
    size_params = size_estimate(phys_params,num_params,observ_params,csm_bool)
    
    ## 4. preallocation: #is it really necessary? and/or can it not simply be done within the precomputation? better like this for performance analysis    
    precomp_arrs,temp_arrs,interpol_arrs,fill_arrs,result_arrs = preallocate_data(phys_params,num_params,observ_params,size_params,csm_bool)

    ## 5. precomputation:
    precompute_ISGL(phys_params,num_params,size_params,precomp_arrs,temp_arrs)    
    
    ## 6. preparation of interpolation & shoulder:
    interpolNshoulder(phys_params,num_params,observ_params,size_params,precomp_arrs,interpol_arrs,wf_bool,csm_bool)
    
    ## 7. Calculation of matrix elements
    fill_TVS(num_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,csm_bool,phys_params.hbar,debug_bool)
    
    ## 8. Solving the generalized eigenproblem:
    solveHS(num_params,fill_arrs,result_arrs,wf_bool)
    
    ## 9. Calculate observables: (limited to "central" observables at the moment)
    if wf_bool == 0
        return result_arrs.energies_arr
    elseif wf_bool == 1
        calc_observables(num_params,observ_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,result_arrs,csm_bool)
        return result_arrs.energies_arr,result_arrs.wavefun_arr,result_arrs.centobs_output,result_arrs.R2_output
    end
end



end ## end of module





### concept of the program:

## 1. interpretation of inputs
#   - phys_params = mass_arr,svals,v0_arr,muG_arr,J_tot,parity
#   - abcvals = avals,bvals,cvals
#   - num_params = lmax,Lmax,gem_params,theta_csm,omega_cr,mu0,c
#   - gem_params = nmax,Nmax,r1,rnmax,R1,RNmax
#   - ...

## 2. sanity check of inputs

## 3. computations to determine sizes of arrays for allocation:
#   - avals,bvals   : aus cvals,svals
#   - lL_arr        : J,parity,svals,bvals,lmax,Lmax
#   - imax_arr      : lL_arr
#   - kmax_dict     : lmax

# lLmax'    = lastindex(lL_arr)
# lmax'     = max(lastindex(L_list),lastindex(l_list))


## 4. allocations: (sizes)
#   - gamma_dict    : nmax?
#   - cleb_arr      : (nl+1) x (2*nl+1) x (nl+1) x (2*nl+1)

#   - jmat          : Matrix{SMatrix{2, 2, Float64, 4}}(undef, 3, 3);

#   - ranges        : nmax; Nmax
#   - norms         : lmax' x nmax ; lmax' x Nmax ?

#   - imax_arr      : lLmax' x lLmax'; (dict?)
#   - kmax_dict     : lmax' x (2*lmax' + 1)
#   - mij_arr       : lastindex(imax_arr)
#   - Clmk_arr      : lmax' x (2*lmax' +1) x max(k_max_arr)
#   - Dlmk_arr      : lmax' x (2*lmax' +1) x max(k_max_arr) x 3
#   - S_arr         : lLmax x lLmax x max(imax_arr)

#   - alpha_arr     : kmax_interpol
#   - v_arr         : kmax_interpol x (2*maxlmax+1)
#   - A_mat         : (2*maxlmax+1) x (2*maxlmax+1)
#   - w_arr         : kmax_interpol x (2*maxlmax+1) x [number of different interactions?!]

#   - T,V,S         : (lLmax'*nmax*Nmax) x (lLmax'*nmax*Nmax)

#   - energies_arr  : (lLmax'*nmax*Nmax)
#   - wavefun_arr   : (lLmax'*nmax*Nmax) x (lLmax'*nmax*Nmax)


## 5. Precomputations before loops over basis functions: (needs)
#   - gamma_dict    : nmax?
#   - cleb_arr      : J_tot,nl
#   - jmat          : mass_arr
#   - ranges        : gem_params
#   - norms         : lL_arr,ranges,gem_params
#   - imax_arr      : lL_arr
#   - kmax_dict     : l_list,L_list (use the one with more elements)
#   - mij_arr       : lL_arr
#   - Clmk_arr      : l_list,L_list (the one with more elements), kmax_dict
#   - Dlmk_arr      : l_list,L_list (the one with more elements), kmax_dict
#   - S_arr         : lL_arr,J_tot,imax_arr,mij_arr,cleb_arr

## 6. Preparations for range-interpolation and shoulder method:
#   - alpha_arr     : gem_params
#   - v_arr         : 
#   - A_mat         : 
#   - w_arr         : 

## 7. TVS matrices & diagonalization
#   - T,V,S         : jmat,ranges,norms,S_arr,phys_params,abcvals,lL_arr,gem_params
#   - energies,wf   : T,V,S

## 8. output