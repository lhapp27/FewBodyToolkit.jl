## Function/Module for using the Gaussian expansion method (GEM) for solving two-body problems in 3D

## by Lucas Happ,  first version from 06.03.2024

module GEM3D

#using LinearAlgebra, SpecialFunctions, QuadGK, Optim, Roots, StaticArrays, Printf
using Printf: @printf

include("../common/eigen2step.jl")
include("optimV0.jl")
include("MatrixElements.jl")
include("wavefunctions.jl")

export GEM_solve
export PreallocStruct
export v0GEMOptim
export GEM_Optim_2B
export buildnu
export MatrixS,MatrixT,MatrixV
export wavefun_arr,wavefun_point

struct PreallocStruct
    nu_arr::Vector{T} where T
    S::Matrix{T} where T
    T::Matrix{T} where T
    V::Matrix{T} where T
    energies::Vector{T2} where T2
    wavefunctions::Matrix{T} where T
    
    function PreallocStruct(num_params, cr_bool, csm_bool)        
        nmax = num_params.gem_params.nmax
        
        # Determine types: TTV = type of T and V matrix; TS = type of S matrix; TE = type of energies array
        if cr_bool == 0 && csm_bool == 0
            TTV = Float64; TS = Float64; TE = Float64
        elseif cr_bool == 1 && csm_bool == 0
            TTV = ComplexF64; TS = ComplexF64; TE = Float64; nmax *= 2;        
        elseif cr_bool == 0 && csm_bool == 1
            TTV = ComplexF64; TS = Float64; TE = ComplexF64
        elseif cr_bool == 1 && csm_bool == 1
            TTV = ComplexF64; TS = ComplexF64; TE = ComplexF64; nmax *= 2;
        else
            error("Error with (cr_bool = $cr_bool, csm_bool = $csm_bool). Only (0,0), (1,0), (0,1), (1,1) allowed.")
        end                
        
        # Initialize arrays
        nu_arr = zeros(TTV, nmax)
        S  = zeros(TS, nmax, nmax)
        T  = zeros(TTV, nmax, nmax)
        V  = zeros(TTV, nmax, nmax)
        energies = zeros(TE, nmax)
        wavefunctions = zeros(TTV, nmax, nmax)
        
        new(nu_arr, S, T, V, energies, wavefunctions) # creates one instance of the struct. only works within an inner "constructor" function
    end
end


#gaussopt: [boolean,v0,mu_g]
function GEM_solve(phys_params, num_params, wf_bool; cr_bool=0, csm_bool=0, gaussopt=[[0,-1.0,1.0]])
    
    # if gaussopt is used, all booleans must be the same. would make more sense to have only one gaussbool...
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
        error("error in wf_bool: only values of 0 or 1 allowed.")
    end
end

# function with preallocated arrays:
function GEM_solve!(prealloc_arrs,phys_params,num_params,wf_bool,cr_bool,csm_bool,gaussopt)
    
    (;nu_arr,S,T,V,energies,wavefunctions) = prealloc_arrs
    
    # Destructuring struct:
    (;hbar,mass_arr,vint_arr,lmax) = phys_params
    vint = vint_arr[1];
    (;gem_params,omega_cr,theta_csm,threshold) = num_params
    (;nmax,r1,rnmax) = gem_params    
    
    ## 1. Preliminaries: 
    # reduced mass:
    mur = 1 / (1/mass_arr[1] + 1/mass_arr[2]);
    
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
        
    MatrixS(S,lmax,nu_arr)
    MatrixT(T,lmax,nu_arr,hbar,mur,csm_bool,theta_csm,cr_bool)
    MatrixV(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt,csm_bool,theta_csm,cr_bool)
    
    debug_bool = 0
    if debug_bool == 1
        #@show(nu_arr)
        #print_matrix("T", T, 3)
        #print_matrix("V", V, 3)
    end

    # symmetric fill:
    S .= Hermitian(S,:L) # if Hermitian or Symmetric: type-unstable? hermitian overall is ok, even if real-symmetric
    if csm_bool == 0 && cr_bool == 0
        T .+= V # t becomes H=T+V
        T .= Symmetric(T,:L)
    elseif csm_bool == 1 && cr_bool == 0
        #T .= Symmetric(T,:L) # T will be complex-symmetric
        #V .= Symmetric(V,:L) # V will be complex-symmetric
        T .+= V # can be added before symmetry
        T .= Symmetric(T,:L) # T will be complex-symmetric
    elseif csm_bool == 0 && cr_bool == 1
        T .+= V 
        T .= Hermitian(T,:L) 
    elseif csm_bool == 1 && cr_bool == 1 # no symmetric filling possible for T and V !!!
        T .+= V 
    end
    
    if debug_bool == 1
        #print_matrix("H", T, 2)
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

## for debugging:
function print_matrix(name, M, minsize)
    r, c = size(M)
    rmax = min(minsize, r)
    cmax = min(minsize, c)
    println("Matrix ", name, " (", r, "x", c, "):")
    for i in 1:rmax
        for j in 1:cmax
            if M[i, j] isa Complex
                @printf("(%8.3f,%8.3f) ", real(M[i,j]), imag(M[i,j]))
            else
                @printf("%8.3f ", M[i,j])
            end
        end
        println()
    end
    println()
end








## Coupled-channel version:
function GEM_solveCC(phys_params, num_params, wf_bool, WCC, DCC; cr_bool=0, csm_bool=0, diff_bool=0)
    if cr_bool == 1 && csm_bool == 1
        error("Currently only either CSM or CR are supported, but not both simultaneously.")
    end
    
    # diff_bool denotes whether derivative terms in DCC should be included. changing this to diff_order = 0 would be better...
    pa = GEM.PreallocStruct(num_params, cr_bool, csm_bool) # preallocate only once
    
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
    GEM.buildnu(nu_arr,r1,rnmax,nmax)
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
    GEM.MatrixS(S,lmax,nu_arr)
    GEM.MatrixT(T,lmax,nu_arr,phys_params.hbar,mur,csm_bool,theta_csm,cr_bool)
    GEM.MatrixV(V,lmax,nu_arr,phys_params.vint_arr[1],gamma_dict,buf,gaussopt,csm_bool,theta_csm,cr_bool)
    
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
            GEM.MatrixWD(WD,lmax,nu_arr,WCC[c,cp],DCC[c,cp],gamma_dict,buf,csm_bool,theta_csm,diff_bool)
            
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
