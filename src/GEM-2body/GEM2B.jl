## Function/Module for using the Gaussian expansion method (GEM) for solving two-body problems in 1D,2D, or 3D

## by Lucas Happ

module GEM2B

using .. FewBodyToolkit
using SpecialFunctions, QuadGK, LinearAlgebra, Optim, StaticArrays, Roots
using Printf: @printf

include("auxiliary.jl")
include("../common/eigen2step.jl")
#include("../common/potentialtypes.jl")
include("optimV0.jl")
include("MatrixElements123D.jl")
include("wavefunctions.jl")

export GEM2B_solve
export make_phys_params2B,make_num_params2B,PreallocStruct2B
export v0GEMOptim,GEM_Optim_2B
export wavefun_arr,wavefun_point

"""
    GEM2B_solve(phys_params, num_params; return_wavefunctions=0, complex_ranged=0, complex_scaling=0)

Solves two-body quantum mechanical problems using the Gaussian Expansion Method (GEM).

# Arguments
- `phys_params`: Physical parameters describing the two-body system:
    + `hbar::Float64`: reduced Planck constant 
    + `mur::Float64`: reduced mass 
    + `interactions=Vector{Any}`: a vector of interactions
    + `lmax::Int`: power of `r^lmax` in the basis functions; indicator for the angular momentum in 3D
    + `dim::Int`: the spatial dimension
- `num_params`: Numerical parameters struct containing information on the set of basis functions:
    + `gem_params::NamedTuple`: (number of basis functions, smallest and largest range parameters).
    + `complex_scaling_angle::Float64`: Complex scaling angle (in radians) for the Complex Scaling Method.
    + `complex_range_freq::Float64`: Parameter controlling the frequency for complex-ranged basis functions.
    + `threshold::Float64`: Numerical threshold generalized eigenvalue solver.

# Keywords
- `return_wavefunctions=0`: Whether to return wavefunctions (0: energies only, 1: energies and wavefunctions)
- `complex_ranged=0`: Whether to use complex rotation method (0: no, 1: yes)
- `complex_scaling=0`: Whether to use complex scaling method (0: no, 1: yes)

# Returns
- `energies`: Array of energy eigenvalues
- `wavefunctions`: (Optional) Array of eigenvectors if return_wavefunctions=1

# Example
```julia
phys_params = make_phys_params2B(hbar=1.0, mur=1.0, interactions=[GaussianPotential(-1.0, 0.5)], lmax=0, dim=3)
num_params = make_num_params2B(gem_params=(nmax=5, r1=1.0, rnmax=10.0), complex_range_freq=0.5, complex_scaling_angle=0.0, threshold=1e-10)
energies = GEM2B_solve(phys_params, num_params; return_wavefunctions=0, complex_ranged=0, complex_scaling=0)
# or with wavefunctions:
energies, wavefunctions = GEM2B_solve(phys_params, num_params; return_wavefunctions=1)
# Note: The function can handle 1D, 2D, or 3D problems based on the `dim` parameter in `phys_params`.
```
"""
function GEM2B_solve(phys_params, num_params;
    return_wavefunctions=false, complex_ranged=false, complex_scaling=false, debug=false, inverse=false, target_energy=0.0,
    wf_bool=nothing, cr_bool=nothing, csm_bool=nothing, debug_bool=nothing, inverse_bool=nothing)

    if !isnothing(wf_bool)
        @warn "wf_bool is deprecated, use return_wavefunctions instead"
        return_wavefunctions = wf_bool
    end
    if !isnothing(cr_bool)
        @warn "cr_bool is deprecated, use complex_ranged instead"
        complex_ranged = cr_bool
    end
    if !isnothing(csm_bool)
        @warn "csm_bool is deprecated, use complex_scaling instead"
        complex_scaling = csm_bool
    end
    if !isnothing(debug_bool)
        @warn "debug_bool is deprecated, use debug instead"
        debug = debug_bool
    end
    if !isnothing(inverse_bool)
        @warn "inverse_bool is deprecated, use inverse instead"
        inverse = inverse_bool
    end
        
    # preallocations:
    prealloc_arrs = PreallocStruct2B(num_params, complex_ranged, complex_scaling)
    
    # function call with preallocations:
    GEM2B_solve!(prealloc_arrs,phys_params,num_params,return_wavefunctions,complex_ranged,complex_scaling,debug,inverse,target_energy)
    
    if !return_wavefunctions
        if debug
            return prealloc_arrs
        end
        return prealloc_arrs.energies

    elseif return_wavefunctions
        if debug
            return prealloc_arrs
        end
        return prealloc_arrs.energies, prealloc_arrs.wavefunctions
    else
        error("error in return_wavefunctions: only boolean values allowed.")
    end
end

# function with preallocated arrays:
function GEM2B_solve!(prealloc_arrs,phys_params,num_params,return_wavefunctions::Bool,complex_ranged::Bool,complex_scaling::Bool,debug::Bool,inverse::Bool,target_energy)

    (;nu_arr,S,T,V,energies,wavefunctions) = prealloc_arrs
    
    # Destructuring struct:
    (;hbar,mur,interactions,lmax,dim) = phys_params # dimension is moved to a physical parameter
    (;gem_params,complex_range_freq,complex_scaling_angle,threshold) = num_params
    (;nmax,r1,rnmax) = gem_params    
    
    ## 1. Preliminaries:
    
    # gamma function: Dict due to half-integer arguments. Change to array possible by multiplication of argument by 2
    gamma_dict = Dict{Float64, Float64}()
    for i = 0.5:0.5:max(nmax,2*lmax+1)+0.5
        gamma_dict[i] = gamma(i)
    end
    
    # gaussian ranges:
    buildnu(nu_arr,r1,rnmax,nmax)
    if complex_ranged
        nu_arr .*= (1+complex_range_freq*im)
        nu_arr[nmax+1:2*nmax] .= conj.(nu_arr[1:nmax])
    end
    
    ## 2. Matrix elements
    # for numerical integration:
    if complex_ranged == 1 || complex_scaling == 1
        buf = alloc_segbuf(Float64,ComplexF64,Float64)
    else
        buf = alloc_segbuf(Float64,Float64,Float64)
    end

    MatrixS(S,lmax,nu_arr,dim)
    MatrixT(T,lmax,nu_arr,hbar,mur,complex_scaling,complex_scaling_angle,complex_ranged,dim)
    MatrixV(V,lmax,nu_arr,interactions,gamma_dict,buf,complex_scaling,complex_scaling_angle,complex_ranged,dim)
    
    if debug
        stp = min(9, size(T, 1))  # Adjust size_to_print as needed
        println("T:")
        display(T[1:stp,1:stp])
        println("V:")
        display(V[1:stp,1:stp])
        println("S:")
        display(S[1:stp,1:stp])
    end

    # symmetric fill:
    S .= Hermitian(S,:L) # if Hermitian or Symmetric: type-unstable? hermitian overall is ok, even if real-symmetric
    if !complex_scaling && !complex_ranged
        #T .+= V # T becomes H=T+V
        T .= Symmetric(T,:L)
        V .= Symmetric(V,:L)
    elseif complex_scaling && !complex_ranged
        #T .+= V # can be added before symmetry
        T .= Symmetric(T,:L) # T will be complex-symmetric
        V .= Symmetric(V,:L) # V will be complex-symmetric
    elseif !complex_scaling && complex_ranged
        #T .+= V 
        T .= Hermitian(T,:L) 
        V .= Hermitian(V,:L)
    elseif complex_scaling && complex_ranged # no symmetric filling possible for T and V !!!
        #T .+= V 
    end

    
    ## 3. Eigensolver
    if !return_wavefunctions
        if inverse
            # solve the inverse problem: find the values of v0 such that the first energy is close to target_energy
            eigen2step(energies,T .- target_energy.*S,-V;threshold=threshold) # The usual variable "energies" is used here for the eigenvalues which are the critical values of v0 leading to the target energy; a smaller threshold might be required!
        else
            eigen2step(energies,T.+V,S;threshold=threshold) # only energies
        end
        
        return energies
    else
        if inverse
            # solve the inverse problem: find basis parameters such that the first energy is close to target_energy
            eigen2step_valvec(energies,wavefunctions,T .- target_energy.*S,-V;threshold=threshold)
        else
            eigen2step_valvec(energies,wavefunctions,T.+V,S;threshold=threshold)
        end
        return energies,wavefunctions
    end    
end



## Gaussian ranges:
@views @inbounds function buildnu(nu_arr,r1,rnmax,nmax)
    nu_arr[1] = 1 /abs(r1)^2;
    nmax >1 && @. nu_arr[2:nmax] = 1/abs(r1)^2 * abs(r1/rnmax)^(2*((2:nmax)-1)/(nmax-1))
end







## Coupled-channel version:
function GEM2B_solveCC(phys_params, num_params, WCC, DCC;
    return_wavefunctions=false, complex_ranged=false, complex_scaling=false, return_diff=false,
    wf_bool=nothing, cr_bool=nothing, csm_bool=nothing, diff_bool=nothing)

    if !isnothing(wf_bool)
        @warn "wf_bool is deprecated, use return_wavefunctions instead"
        return_wavefunctions = wf_bool
    end
    if !isnothing(cr_bool)
        @warn "cr_bool is deprecated, use complex_ranged instead"
        complex_ranged = cr_bool
    end
    if !isnothing(csm_bool)
        @warn "csm_bool is deprecated, use complex_scaling instead"
        complex_scaling = csm_bool
    end
    if !isnothing(diff_bool)
        @warn "diff_bool is deprecated, use return_diff instead"
        return_diff = diff_bool
    end

    if complex_ranged && complex_scaling
        error("Currently only either CSM or CR are supported, but not both simultaneously.")
    end
    
    # return_diff denotes whether derivative terms in DCC should be included. it might be better to change this to diff_order = 0
    pa = GEM2B.PreallocStruct2B(num_params, complex_ranged, complex_scaling) # preallocate only once
    
    (;nu_arr,S,T,V,energies,wavefunctions) = pa # de-struct
    (;lmax,mur,dim) = phys_params
    (;gem_params,complex_range_freq,complex_scaling_angle,threshold) = num_params
    (;nmax,r1,rnmax) = gem_params 
    
    TT = typeof(T[1])
    TT2 = typeof(energies[1])
    
    # for convenience
    cmax = size(WCC,1) # number of channels
    nmax = lastindex(nu_arr) # incorporates already factor 2 in case of complex_ranged=1
    
    ## 1. Allocation of big matrices:
    WD = similar(V) # dummy matrix for WCC matrices
    Hbig = zeros(TT, cmax*nmax, cmax*nmax)
    Sbig = zeros(TT, cmax*nmax, cmax*nmax)
    energiesbig = zeros(TT2, cmax*nmax)
    wfbig = zeros(TT, cmax*nmax, cmax*nmax)
    
    ## 2. Preliminaries: 
    
    # gamma function:
    gamma_dict = Dict{Float64, Float64}()
    for i = 0.5:0.5:max(nmax,2*lmax+1)+0.5
        gamma_dict[i] = gamma(i)
    end
    
    # gaussian ranges:
    GEM2B.buildnu(nu_arr,r1,rnmax,nmax)
    if complex_ranged == 1
        nu_arr .*= (1+complex_range_freq*im)
        nu_arr[nmax+1:2*nmax] .= conj.(nu_arr[1:nmax])
    end
    
    ## 2. Matrix elements
    # for numerical integration:
    if complex_ranged || complex_scaling
        buf = alloc_segbuf(Float64,ComplexF64,Float64)
    else
        buf = alloc_segbuf(Float64,Float64,Float64)
    end
    
    # norm-overlap, kinetic energy, and channel-independent potential need to be calculated only once:
    GEM2B.MatrixS(S,lmax,nu_arr,dim)
    GEM2B.MatrixT(T,lmax,nu_arr,phys_params.hbar,mur,complex_scaling,complex_scaling_angle,complex_ranged,dim)
    GEM2B.MatrixV(V,lmax,nu_arr,phys_params.interactions,gamma_dict,buf,complex_scaling,complex_scaling_angle,complex_ranged,dim)
    
    # symmetric fill: seems correct even for coupled channels, since the same basis is used in each channel
    S .= Hermitian(S,:L) # if Hermitian or Symmetric: type-unstable?
    if !complex_scaling
        T .= Hermitian(T,:L)
        V .= Hermitian(V,:L)
        #T .+= V
    elseif complex_scaling
        T .= Symmetric(T,:L) # T will be complex-symmetric
        V .= Symmetric(V,:L) # V will be complex-symmetric
        #T .+= V
    end
    
    
    # loop over channels:
    for c = 1:cmax
        for cp = 1:c # sym_bool is always fulfilled, no?
            #sym_bool == 1 && cp > c && continue # only fill lower triangle of big matrix
            
            # WCC and DCC matrices:
            GEM2B.MatrixWD(WD,lmax,nu_arr,WCC[c,cp],DCC[c,cp],gamma_dict,buf,complex_scaling,complex_scaling_angle,return_diff,dim)
            
            if !complex_scaling
                WD .= Hermitian(WD,:L)
            elseif complex_scaling
                WD .= Symmetric(WD,:L) # matrix will be complex-symmetric
            end
            
            # fill Hbig:
            Hbig[(c-1)*nmax+1:c*nmax,(cp-1)*nmax+1:cp*nmax] .= WD
            c == cp && (Hbig[(c-1)*nmax+1:c*nmax,(cp-1)*nmax+1:cp*nmax] .+= T .+ V)
            
            # fill Sbig:
            c == cp && (Sbig[(c-1)*nmax+1:c*nmax,(cp-1)*nmax+1:cp*nmax] .= S)
        end
    end
    
    Sbig .= Symmetric(Sbig,:L)
    
    # not sure if correct. Can we always assume that WCC and DCC are symmetric or hermitian, and hence Hbig is so too?
    if complex_scaling == 0
        Hbig .= Hermitian(Hbig,:L)
    elseif complex_scaling == 1
        Hbig .= Symmetric(Hbig,:L) # matrix will be complex-symmetric
    end
    
    ## 3. Eigensolver
    if return_wavefunctions == 0
        eigen2step(energiesbig,Hbig,Sbig;threshold=num_params.threshold) # only energies
        return energiesbig
    else
        eigen2step_valvec(energiesbig,wfbig,Hbig,Sbig;threshold=num_params.threshold)
        return energiesbig,wfbig
    end
    
end



end ## end of module
