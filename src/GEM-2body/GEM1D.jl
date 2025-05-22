## 1D!
## Function/Module for using the Gaussian expansion method (GEM) for solving two-body problems
## 1-file version:

## by Lucas Happ,  27.04.2024

module GEM1D

using LinearAlgebra, SpecialFunctions, QuadGK, GSL, Optim

include("eigen2step.jl")
include("MatrixElements1D.jl")
include("optim_v0_scatt.jl")
include("wavefunctions.jl")

export GEM_solve
export PreallocStruct
export MatrixS1D,MatrixT1D,MatrixV1D
export v0GEMOptim, GEM_Optim_2B
export wavefun_arr, wavefun_point1D

#= struct PreallocStruct{T,T2}
    nu_arr::Vector{T}
    S::Matrix{T}
    T::Matrix{T}
    V::Matrix{T}
    energies::Vector{T2}
    wavefunctions::Matrix{T}
end =#

struct PreallocStruct
    nu_arr::Vector{T} where T
    S::Matrix{T} where T
    T::Matrix{T} where T
    V::Matrix{T} where T
    energies::Vector{T2} where T2
    wavefunctions::Matrix{T} where T
    
    function PreallocStruct(l_arr, num_params, cr_bool, csm_bool)        
        nbasis = num_params.gem_params.nmax *lastindex(l_arr)
        
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
        
        # maybe easier: but also less precise
        #=         if cr_bool == 0 && csm_bool == 0
            TTV = Float64; TS = Float64; TE = Float64
        elseif cr_bool == 1 || csm_bool == 1
            TTV = ComplexF64; TS = ComplexF64; TE = ComplexF64
        end
        if cr_bool == 1
            nmax *= 2 #needs change to nbasis
        end =#
        
        
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

#gaussopt: [boolean,v0,mu_g] #gaussian potential
#shiftgopt: [boolean,[v0_i],[mu_g_i],[z0_i]] #gaussian potential shifted by z0_i

function GEM_solve(phys_params, num_params, wf_bool; cr_bool=0, csm_bool=0, gaussopt=[0,-1.0,1.0], shiftgopt =[0,[-1.0],[1.0],[1.0]])
    if cr_bool == 1 && csm_bool == 1
        error("Currently only either CSM or CR are supported, but not both simultaneously.") # in prealloc it is already done, but not in the matrix elements.
    end
    
    (;lmin,lmax) = phys_params
    lmin != lmax && error("lmin != lmax. we allow for the input, but were too lazy to implement it yet")

    l_arr = phys_params.lmin:1:phys_params.lmax # angular momentum of the two-body system
    # preallocations: now similar to 3D case
    prealloc_arrs = PreallocStruct(l_arr, num_params, cr_bool, csm_bool)
    
    # function call with preallocations:
    GEM_solve1D!(prealloc_arrs,phys_params,num_params,wf_bool,cr_bool,csm_bool,gaussopt,shiftgopt)
    
    if wf_bool == 0
        return energies
    elseif wf_bool == 1
        return energies, wavefunctions
    else
        error("error in wf_bool: only values of 0 or 1 allowed.")
    end
end

# function with preallocated arrays:
function GEM_solve1D!(prealloc_arrs,phys_params,num_params,wf_bool,cr_bool,csm_bool,gaussopt,shiftgopt)
    
    (;nu_arr,S,T,V,energies,wavefunctions) = prealloc_arrs
    
    # Destructuring struct:
    (;hbar,mass_arr,vint_arr,lmax,lmin) = phys_params
    vint = vint_arr[1];
    (;gem_params,omega_cr,theta_csm,threshold) = num_params
    (;nmax,r1,rnmax) = gem_params    
    
    ## 1. Preliminaries:
    # angular momenta:
    l_arr = phys_params.lmin:1:phys_params.lmax
    nbasis = num_params.gem_params.nmax *lastindex(l_arr) # not adjusted for CR! this is only done in PreallocStruct
    
    # reduced mass:
    mur = 1 / (1/mass_arr[1] + 1/mass_arr[2]);
    
    # gamma function:
    gamma_dict = Dict{Float64, Float64}()
    for i = 0.5:0.5:max(nmax,2*lmax+1)+0.5
        gamma_dict[i] = gamma(i)
    end
    
    # gaussian ranges:
    buildnu(nu_arr,r1,rnmax,nbasis)
    if cr_bool == 1
        nu_arr .*= (1+omega_cr*im)
        nu_arr[nbasis+1:2*nbasis] .= conj.(nu_arr[1:nbasis])
    end
    
    ## 2. Matrix elements
    # for numerical integration:
    if cr_bool == 1 || csm_bool == 1
        buf = alloc_segbuf(Float64,ComplexF64,Float64)
    else
        buf = alloc_segbuf(Float64,Float64,Float64)
    end
    
    MatrixS1D(S,lmax,nu_arr)
    MatrixT1D(T,lmax,nu_arr,hbar,mur,csm_bool,theta_csm)
    MatrixV1D(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt,shiftgopt,csm_bool,theta_csm)
    
    # symmetric fill:
    S .= Hermitian(S,:L) # if Hermitian or Symmetric: type-unstable?
    if csm_bool == 0
        T .+= V
        T .= Hermitian(T,:L)
    elseif csm_bool == 1
        T .= Symmetric(T,:L) # T will be complex-symmetric
        V .= Symmetric(V,:L) # V will be complex-symmetric
        T .+= V
    end
    
    ## 3. Eigensolver    
    if wf_bool == 0
        eigen2step(energies,T,S;threshold=threshold) # only energies
        return energies
    else
        eigen2step_valvec(energies,wavefunctions,T,S;threshold=threshold)
        return energies,wavefunctions        
    end    
end




## Gaussian ranges:
@views @inbounds function buildnu(nu_arr,r1,rnmax,nmax)
    nu_arr[1] = 1 /abs(r1)^2;
    nmax >1 && @. nu_arr[2:nmax] = 1/abs(r1)^2 * abs(r1/rnmax)^(2*((2:nmax)-1)/(nmax-1))
end


end ## end of module
