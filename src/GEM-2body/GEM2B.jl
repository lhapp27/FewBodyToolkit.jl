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
    GEM2B_solve(phys_params, num_params; wf_bool=0, cr_bool=0, csm_bool=0)

Solves two-body quantum mechanical problems using the Gaussian Expansion Method (GEM).

# Arguments
- `phys_params`: Physical parameters describing the two-body system:
    + `hbar::Float64`: reduced Planck constant 
    + `mur::Float64`: reduced mass 
    + `vint_arr=Vector{Any}`: a vector of interactions
    + `lmax::Int`: power of `r^lmax` in the basis functions; indicator for the angular momentum in 3D
    + `dim::Int`: the spatial dimension
- `num_params`: Numerical parameters struct containing information on the set of basis functions:
    + `gem_params::NamedTuple`: (number of basis functions, smallest and largest range parameters).
    + `theta_csm::Float64`: Complex scaling angle (in radians) for the Complex Scaling Method.
    + `omega_cr::Float64`: Parameter controlling the frequency for complex-ranged basis functions.
    + `threshold::Float64`: Numerical threshold generalized eigenvalue solver.

# Keywords
- `wf_bool=0`: Whether to return wavefunctions (0: energies only, 1: energies and wavefunctions)
- `cr_bool=0`: Whether to use complex rotation method (0: no, 1: yes)
- `csm_bool=0`: Whether to use complex scaling method (0: no, 1: yes)

# Returns
- `energies`: Array of energy eigenvalues
- `wavefunctions`: (Optional) Array of eigenvectors if wf_bool=1

# Example
```julia
phys_params = make_phys_params2B(hbar=1.0, mur=1.0, vint_arr=[GaussianPotential(-1.0, 0.5)], lmax=0, dim=3)
num_params = make_num_params2B(gem_params=(nmax=5, r1=1.0, rnmax=10.0), omega_cr=0.5, theta_csm=0.0, threshold=1e-10)
energies = GEM2B_solve(phys_params, num_params; wf_bool=0, cr_bool=0, csm_bool=0)
# or with wavefunctions:
energies, wavefunctions = GEM2B_solve(phys_params, num_params; wf_bool=1)
# Note: The function can handle 1D, 2D, or 3D problems based on the `dim` parameter in `phys_params`.
```
"""
function GEM2B_solve(phys_params, num_params; wf_bool=0, cr_bool=0, csm_bool=0, debug_bool=0, inverse_bool=0, target_energy=0.0)
        
    # preallocations:
    prealloc_arrs = PreallocStruct2B(num_params, cr_bool, csm_bool)
    
    # function call with preallocations:
    GEM2B_solve!(prealloc_arrs,phys_params,num_params,wf_bool,cr_bool,csm_bool,debug_bool,inverse_bool,target_energy)
    
    if wf_bool == 0
        if debug_bool == 1
            return prealloc_arrs
        end
        return prealloc_arrs.energies

    elseif wf_bool == 1
        if debug_bool == 1
            return prealloc_arrs
        end
        return prealloc_arrs.energies, prealloc_arrs.wavefunctions
    else
        error("error in wf_bool: only values of 0 or 1 allowed.")
    end
end

# function with preallocated arrays:
function GEM2B_solve!(prealloc_arrs,phys_params,num_params,wf_bool,cr_bool,csm_bool,debug_bool,inverse_bool,target_energy)

    (;nu_arr,S,T,V,energies,wavefunctions) = prealloc_arrs
    
    # Destructuring struct:
    (;hbar,mur,vint_arr,lmax,dim) = phys_params # dimension is moved to a physical parameter
    (;gem_params,omega_cr,theta_csm,threshold) = num_params
    (;nmax,r1,rnmax) = gem_params    
    
    ## 1. Preliminaries:
    
    # gamma function: Dict due to half-integer arguments. Change to array possible by multiplication of argument by 2
    gamma_dict = Dict{Float64, Float64}()
    for i = 0.5:0.5:max(nmax,2*lmax+1)+0.5
        gamma_dict[i] = gamma(i)
    end
    
    # gaussian ranges:
    buildnu(nu_arr,r1,rnmax,nmax)
    if cr_bool == 1
        nu_arr .*= (1+omega_cr*im)
        nu_arr[nmax+1:2*nmax] .= conj.(nu_arr[1:nmax])
    end
    
    ## 2. Matrix elements
    # for numerical integration:
    if cr_bool == 1 || csm_bool == 1
        buf = alloc_segbuf(Float64,ComplexF64,Float64)
    else
        buf = alloc_segbuf(Float64,Float64,Float64)
    end

    MatrixS(S,lmax,nu_arr,dim)
    MatrixT(T,lmax,nu_arr,hbar,mur,csm_bool,theta_csm,cr_bool,dim)
    MatrixV(V,lmax,nu_arr,vint_arr,gamma_dict,buf,csm_bool,theta_csm,cr_bool,dim)
    
    if debug_bool == 1
        minsize = min(10, size(T, 1))
        print_matrix("T", T, minsize)
        print_matrix("V", V, minsize)
        print_matrix("S", S, minsize)
    end

    # symmetric fill:
    S .= Hermitian(S,:L) # if Hermitian or Symmetric: type-unstable? hermitian overall is ok, even if real-symmetric
    if csm_bool == 0 && cr_bool == 0
        #T .+= V # T becomes H=T+V
        T .= Symmetric(T,:L)
        V .= Symmetric(V,:L)
    elseif csm_bool == 1 && cr_bool == 0
        #T .+= V # can be added before symmetry
        T .= Symmetric(T,:L) # T will be complex-symmetric
        V .= Symmetric(V,:L) # V will be complex-symmetric
    elseif csm_bool == 0 && cr_bool == 1
        #T .+= V 
        T .= Hermitian(T,:L) 
        V .= Hermitian(V,:L)
    elseif csm_bool == 1 && cr_bool == 1 # no symmetric filling possible for T and V !!!
        #T .+= V 
    end

    
    ## 3. Eigensolver
    if wf_bool == 0
        if inverse_bool == 1
            # solve the inverse problem: find the values of v0 such that the first energy is close to target_energy
            eigen2step(energies,T .- target_energy.*S,-V;threshold=threshold) # The usual variable "energies" is used here for the eigenvalues which are the critical values of v0 leading to the target energy; a smaller threshold might be required!
        else
            eigen2step(energies,T.+V,S;threshold=threshold) # only energies
        end
        
        return energies
    else
        if inverse_bool == 1
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

## for debugging:
function print_matrix(name, M, minsize)
    r, c = size(M)
    rmax = min(minsize, r)
    cmax = min(minsize, c)
    println("Matrix ", name, " (", r, "x", c, "):")
    for i in 1:rmax
        for j in 1:cmax
            if M[i, j] isa Complex
                @printf("(%8.4f,%8.4f) ", real(M[i,j]), imag(M[i,j]))
            else
                @printf("%8.4f ", M[i,j])
            end
        end
        println()
    end
    println()
end








## Coupled-channel version:
function GEM2B_solveCC(phys_params, num_params, WCC, DCC; wf_bool=0, cr_bool=0, csm_bool=0, diff_bool=0)
    if cr_bool == 1 && csm_bool == 1
        error("Currently only either CSM or CR are supported, but not both simultaneously.")
    end
    
    # diff_bool denotes whether derivative terms in DCC should be included. it might be better to change this to diff_order = 0
    pa = GEM2B.PreallocStruct2B(num_params, cr_bool, csm_bool) # preallocate only once
    
    (;nu_arr,S,T,V,energies,wavefunctions) = pa # de-struct
    (;lmax,mur,dim) = phys_params
    (;gem_params,omega_cr,theta_csm,threshold) = num_params
    (;nmax,r1,rnmax) = gem_params 
    
    TT = typeof(T[1])
    TT2 = typeof(energies[1])
    
    # for convenience
    cmax = size(WCC,1) # number of channels
    nmax = lastindex(nu_arr) # incorporates already factor 2 in case of cr_bool=1
    
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
    if cr_bool == 1
        nu_arr .*= (1+omega_cr*im)
        nu_arr[nmax+1:2*nmax] .= conj.(nu_arr[1:nmax])
    end
    
    ## 2. Matrix elements
    # for numerical integration:
    if cr_bool == 1 || csm_bool == 1
        buf = alloc_segbuf(Float64,ComplexF64,Float64)
    else
        buf = alloc_segbuf(Float64,Float64,Float64)
    end
    
    # norm-overlap, kinetic energy, and channel-independent potential need to be calculated only once:
    GEM2B.MatrixS(S,lmax,nu_arr,dim)
    GEM2B.MatrixT(T,lmax,nu_arr,phys_params.hbar,mur,csm_bool,theta_csm,cr_bool,dim)
    GEM2B.MatrixV(V,lmax,nu_arr,phys_params.vint_arr,gamma_dict,buf,csm_bool,theta_csm,cr_bool,dim)
    
    # symmetric fill: seems correct even for coupled channels, since the same basis is used in each channel
    S .= Hermitian(S,:L) # if Hermitian or Symmetric: type-unstable?
    if csm_bool == 0
        T .= Hermitian(T,:L)
        V .= Hermitian(V,:L)
        #T .+= V
    elseif csm_bool == 1
        T .= Symmetric(T,:L) # T will be complex-symmetric
        V .= Symmetric(V,:L) # V will be complex-symmetric
        #T .+= V
    end
    
    
    # loop over channels:
    for c = 1:cmax
        for cp = 1:c # sym_bool is always fulfilled, no?
            #sym_bool == 1 && cp > c && continue # only fill lower triangle of big matrix
            
            # WCC and DCC matrices:
            GEM2B.MatrixWD(WD,lmax,nu_arr,WCC[c,cp],DCC[c,cp],gamma_dict,buf,csm_bool,theta_csm,diff_bool,dim)
            
            if csm_bool == 0
                WD .= Hermitian(WD,:L)
            elseif csm_bool == 1
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
    if csm_bool == 0
        Hbig .= Hermitian(Hbig,:L)
    elseif csm_bool == 1
        Hbig .= Symmetric(Hbig,:L) # matrix will be complex-symmetric
    end
    
    ## 3. Eigensolver
    if wf_bool == 0
        eigen2step(energiesbig,Hbig,Sbig;threshold=num_params.threshold) # only energies
        return energiesbig
    else
        eigen2step_valvec(energiesbig,wfbig,Hbig,Sbig;threshold=num_params.threshold)
        return energiesbig,wfbig
    end
    
end



end ## end of module
