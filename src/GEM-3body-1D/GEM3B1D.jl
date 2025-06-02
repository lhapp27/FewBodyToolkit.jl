## Function/Module for using the Gaussian expansion method (GEM) for solving three-body problems in 1D(!)
## the code might seem unnecessary complicated at places. this is due to reusing the structure from the 3D ISGL code.

#TD:
#- 2Bdocs: Antique
#- 3B1D: finish cleanup, fix missing arguments, potential type dispatch and handling
#- 3B1D: example (and tests)
#- 3B1D: docs based on example



module GEM3B1D

using .. FewBodyToolkit
using SpecialFunctions, QuadGK, GSL, LinearAlgebra, Optim, StaticArrays, Roots, Interpolations, OffsetArrays
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
#include("../common/potentialtypes.jl")

export GEM3B1D_solve
export make_phys_params3B1D,make_num_params3B1D
#export wavefun

# gaussopt=[bool,v0,mu_g] for central gauss interaction: V(r) = v0*exp(-mu_g*r^2)
function GEM3B1D_solve(phys_params, num_params; wf_bool=0, csm_bool=0, observ_params=(;stateindices=[],centobs_arr=[[],[],[]],R2_arr=[0,0,0]))
    
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
    
    ## 3. computations to determine sizes of arrays for allocation:   
    size_params = size_estimate(phys_params,num_params,observ_params)
    
    ## 4. preallocation: #is it really necessary? and/or can it not simply be done within the precomputation? better like this for performance analysis    
    precomp_arrs,interpol_arrs,fill_arrs,result_arrs = preallocate_data(phys_params,num_params,observ_params,size_params,csm_bool)

    ## 5. precomputation:
    precompute_3B1D(phys_params,num_params,size_params,precomp_arrs)

    ## 6. preparation of interpolation & shoulder:
    interpolNshoulder(phys_params,num_params,observ_params,size_params,precomp_arrs,interpol_arrs,wf_bool,csm_bool)
    
    ## 7. Calculation of matrix elements
    fill_TVS(num_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,csm_bool)
    
    ## 8. Solving the generalized eigenproblem:
    solveHS(num_params,fill_arrs,result_arrs,wf_bool)
    
    ## 9. Calculate observables: (currently not supported for 1D)
    if wf_bool == 0
        return result_arrs.energies_arr
    elseif wf_bool == 1
        #calc_observables(num_params,observ_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,result_arrs)
        return result_arrs.energies_arr,result_arrs.centobs_output,result_arrs.R2_output
    end
end


end ## end of module