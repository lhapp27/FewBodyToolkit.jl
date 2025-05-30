## ISGL_solve is some kind of "main-method"
## Function/Module for using the Infinitesimally shifted Gaussian lobe functions (ISGL) within the Gaussian expansion method (GEM) for solving three-body problems

## see bottom for a documentation ("concept of the program"); might be a bit outdated

## 22.01.24: Change to structs
## 21.02.24: Change to module
## 14.03.25: Spin-implementation


module ISGL
using .. FewBodyToolkit
using LinearAlgebra,StaticArrays,OffsetArrays,Interpolations, SpecialFunctions,QuadGK,PartialWaveFunctions, WignerSymbols
using Printf: @printf

# try to make use of multiple dispatch for different potential types
## where should this part be placed ideally?
# types for the different potential functions
#= abstract type PotentialFunction end ##moved to common/potentialtypes.jl
struct CentralPotential <: PotentialFunction
    f::Function
end
struct SpinOrbitPotential <: PotentialFunction
    f::Function
end
# for evaluating these functions, e.g. in quadgk; used also for observables
function (v::PotentialFunction)(r)
    return v.f(r)
end =#


# functions for creating inputs
function make_phys_params(;hbar = 1.0, mass_arr=[1.0, 1.0, 1.0], svals=["b","b","b"], vint_arr=[[],[],[]], J_tot=0, parity=1, spin_arr=[0,0,0])
    return (;hbar, mass_arr, svals, vint_arr, J_tot, parity, spin_arr)
end

function make_num_params(; lmax=0, Lmax=0, gem_params=(nmax=3, Nmax=3, r1=1.0, rnmax=20.0, R1=1.0, RNmax=20.0),theta_csm=0.0, omega_cr=0.5, mu0=0.08, c_shoulder=1.6, kmax_interpol=1000, threshold=10^-7, lmin=0, Lmin=0)
    return (;lmax, Lmax, gem_params, theta_csm, omega_cr, mu0, c_shoulder, kmax_interpol, threshold, lmin, Lmin)
end

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


export ISGL_solve, make_phys_params, make_num_params
#export wavefun

# gaussopt=[bool,v0,mu_g] for central gauss interaction: V(r) = v0*exp(-mu_g*r^2)
function ISGL_solve(phys_params, num_params; wf_bool = 0, csm_bool = 0, observ_params=(;stateindices=[],centobs_arr=[[],[],[]],R2_arr=[0,0,0]), gaussopt=[[[0,-1.0,1.0]],[[0,-1.0,1.0]],[[0,-1.0,1.0]]]) # default inputs here look a bit cumbersome
    
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
    
    ## 2.b change gaussopt to complex for csm-Bool and "active gaussopt"
    gaussopt = csmgaussopt(gaussopt,csm_bool,num_params)
    
    ## 3. computations to determine sizes of arrays for allocation:   
    size_params = size_estimate(phys_params,num_params,observ_params,gaussopt)
    
    ## 4. preallocation: #is it really necessary? and/or can it not simply be done within the precomputation? better like this for performance analysis    
    precomp_arrs,temp_arrs,interpol_arrs,fill_arrs,result_arrs = preallocate_data(phys_params,num_params,observ_params,size_params,csm_bool)

    ## 5. precomputation:
    precompute_ISGL(phys_params,num_params,size_params,precomp_arrs,temp_arrs)    
    
    ## 6. preparation of interpolation & shoulder:
    interpolNshoulder(phys_params,num_params,observ_params,size_params,precomp_arrs,interpol_arrs,gaussopt,wf_bool,csm_bool)
    
    ## 7. Calculation of matrix elements
    fill_TVS(num_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,gaussopt,csm_bool,phys_params.hbar)
    
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