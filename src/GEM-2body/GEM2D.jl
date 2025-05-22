## Function/Module for using the Gaussian expansion method (GEM) for solving two-body problems in 2D!
## 1-file version:

## by Lucas Happ,  26.01.2025

module GEM2D

using LinearAlgebra, SpecialFunctions, QuadGK, Optim, Roots, StaticArrays

include("eigen2step.jl")
include("optim_v0_scatt.jl")
include("MatrixElements2D.jl")
include("wavefunctions.jl")

export GEM_solve
export PreallocStruct
export v0GEMOptim
export GEM_Optim_2B
export buildnu
export MatrixS2D,MatrixT2D,MatrixV2D
export wavefun_arr,wavefun_point

struct PreallocStruct
    nu_arr::Vector{T} where T
    S::Matrix{T} where T
    T::Matrix{T} where T
    V::Matrix{T} where T
    energies::Vector{T2} where T2
    wavefunctions::Matrix{T} where T
    
    function PreallocStruct(num_params, cr_bool, csm_bool)
        if cr_bool == 1 && csm_bool == 1
            error("Currently only either CSM OR CR are supported: cr_bool & csm_bool cannot be 1 both")
        end
        
        nmax = num_params.gem_params.nmax
        
        # Determine types
        if cr_bool == 0 && csm_bool == 0
            TT = Float64; TT2 = Float64
        elseif cr_bool == 1 && csm_bool == 0
            TT = ComplexF64; TT2 = Float64; nmax *= 2;        
        elseif cr_bool == 0 && csm_bool == 1
            TT = ComplexF64; TT2 = ComplexF64
        else
            error("Error with (cr_bool = $cr_bool, csm_bool = $csm_bool). Only (0,0), (1,0), (0,1) allowed.")
        end
        
        # Initialize arrays
        nu_arr = zeros(TT, nmax)
        S  = zeros(TT, nmax, nmax)
        T  = zeros(TT, nmax, nmax)
        V  = zeros(TT, nmax, nmax)
        energies = zeros(TT2, nmax)
        wavefunctions = zeros(TT, nmax, nmax)
        
        new(nu_arr, S, T, V, energies, wavefunctions) # creates one instance of the struct. only works within an inner "constructor" function
    end
end


#gaussopt: [boolean,v0,mu_g]
function GEM_solve(phys_params, num_params, wf_bool; cr_bool=0, csm_bool=0, gaussopt=[[0,-1.0,1.0]])
    if cr_bool == 1 && csm_bool == 1
        error("Currently only either CSM OR CR are supported: cr_bool & csm_bool cannot be 1 both")
    end

    # if gaussopt ist used, all booleans must be the same. would make more sense to have only one gaussbool...
    # i.e. gaussopt=[[0,v0,mug],[1,v0,mug]] should not be allowed
    gblist = [gaussopt[i][1] for i = 1:lastindex(gaussopt)]
    if lastindex(unique(gblist)) > 1
        error("All gaussbooleans must be the same.")
    end
    
    # preallocations:
    prealloc_arrs = PreallocStruct(num_params, cr_bool, csm_bool)
    
    # function call with preallocations:
    GEM_solve!(prealloc_arrs,phys_params,num_params,wf_bool,cr_bool,csm_bool,gaussopt)
    
    if wf_bool == 0
        return prealloc_arrs.energies
    elseif wf_bool == 1
        return prealloc_arrs.energies, prealloc_arrs.wavefunctions
    else
        error("error in wf_bool=$(wf_bool): only values of 0 or 1 allowed.")
    end
end

# function with preallocated arrays:
function GEM_solve!(prealloc_arrs,phys_params,num_params,wf_bool,cr_bool,csm_bool,gaussopt)
    
    (;nu_arr,S,T,V,energies,wavefunctions) = prealloc_arrs
    
    # Destructuring struct:
    (;hbar,mur,vint_arr,lmax) = phys_params # lmax here plays the role of m as the angular contribution m^2/r^2
    vint = vint_arr[1];
    (;gem_params,omega_cr,theta_csm,threshold) = num_params
    (;nmax,r1,rnmax) = gem_params    
    
    ## 1. Preliminaries: 
    # reduced mass:
    #mur = 1 / (1/mass_arr[1] + 1/mass_arr[2]);
    # changed, such that mur is an input
    
    # gamma function:
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
    
    MatrixS2D(S,lmax,nu_arr)
    MatrixT2D(T,lmax,nu_arr,hbar,mur,csm_bool,theta_csm)
    MatrixV2D(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt,csm_bool,theta_csm)
    
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


## Coupled-channel version:
function GEM_solveCC2D(phys_params, num_params, wf_bool, WCC, DCC; cr_bool=0, csm_bool=0, diff_bool=0)
    
    error("Coupled-channel version not yet implemented for 2D")
    # some function names were changed to contain 2D, but MatrixElements need to be adapted, especially the derivatives
    
    # diff_bool denotes whether derivative terms in DCC should be included. changing this to diff_order = 0 would be better...
    pa = GEM2D.PreallocStruct(num_params, cr_bool, csm_bool) # preallocate only once
    
    (;nu_arr,S,T,V,energies,wavefunctions) = pa # de-struct
    (;gem_params,omega_cr,theta_csm,threshold) = num_params
    (;nmax,r1,rnmax) = gem_params 
    
    TT = typeof(T[1])
    TT2 = typeof(energies[1])
    
    gaussopt = SA[0,1.0,1.0] # i.e. no gaussopt
    
    # for convenience
    cmax = size(WCC,1) # number of channels
    nmax = lastindex(nu_arr) # incorporates already factor 2 in case of cr_bool=1
    lmax = phys_params.lmax
    
    ## 1. Allocation of big matrices:
    WD = similar(V) # dummy matrix for WCC matrices
    #D = similar(V) # dummy matrix for DCC matrices
    Hbig = zeros(TT, cmax*nmax, cmax*nmax)
    Sbig = zeros(TT, cmax*nmax, cmax*nmax)
    energiesbig = zeros(TT2, cmax*nmax)
    wfbig = zeros(TT, cmax*nmax, cmax*nmax)
    
    ## 2. Preliminaries: 
    # reduced mass:
    mur = 1 / (1/phys_params.mass_arr[1] + 1/phys_params.mass_arr[2]); # this should be the actual argument instead of mass_arr i think
    
    # gamma function:
    gamma_dict = Dict{Float64, Float64}()
    for i = 0.5:0.5:max(nmax,2*lmax+1)+0.5
        gamma_dict[i] = gamma(i)
    end
    
    # gaussian ranges:
    GEM2D.buildnu(nu_arr,r1,rnmax,nmax)
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
    GEM2D.MatrixS2D(S,lmax,nu_arr)
    GEM2D.MatrixT2D(T,lmax,nu_arr,phys_params.hbar,mur,csm_bool,theta_csm)
    GEM2D.MatrixV2D(V,lmax,nu_arr,phys_params.vint_arr[1],gamma_dict,buf,gaussopt,csm_bool,theta_csm)
    
    # symmetric fill: is this correct even for coupled channels? probably ok, since the same basis is used in each channel
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
            #GEM.MatrixV(W,lmax,nu_arr,WCC[c,cp],gamma_dict,buf,gaussopt,csm_bool,theta_csm) # can be combined into the function below
            GEM2D.MatrixWD2D(WD,lmax,nu_arr,WCC[c,cp],DCC[c,cp],gamma_dict,buf,csm_bool,theta_csm,diff_bool)
            
            if csm_bool == 0
                WD .= Hermitian(WD,:L)
            elseif csm_bool == 1
                WD .= Symmetric(WD,:L) # matrix will be complex-symmetric
            end
            
            # fill Hbig:
            Hbig[(c-1)*nmax+1:c*nmax,(cp-1)*nmax+1:cp*nmax] .= WD
            c == cp && (Hbig[(c-1)*nmax+1:c*nmax,(cp-1)*nmax+1:cp*nmax] .+= T .+ V)
            #@show([c,cp,WD,T,V,Hbig])

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
    
    #@show([Hbig,Sbig])
    
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
