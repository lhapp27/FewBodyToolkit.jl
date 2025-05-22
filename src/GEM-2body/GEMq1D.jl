## quasi-1D!
## Function/Module for using the Gaussian expansion method (GEM) for solving two-body problems
## 1-file version:

## by Lucas Happ,  29.04.2024
## in principle this could probably be implemented as coupled-channel version of the 1D code

module GEMq1D

#= using LinearAlgebra
using SpecialFunctions
using QuadGK
using HypergeometricFunctions =#

include("../common/eigen2step.jl")
include("MatrixElementsq1D.jl")

export GEM_solveq1D
export PreallocStruct
export MatrixSq1D,MatrixTq1D,MatrixVq1D
export acoeff

struct PreallocStruct{T}
    nu_arr::Vector{T}
    S::Matrix{T}
    T::Matrix{T}
    V::Matrix{T}
    energies::Vector{Float64}
    wavefunctions::Matrix{T}
end

#gaussopt: [boolean,v0,mu_g]
function GEM_solveq1D(phys_params, num_params, ho_params, wf_bool; cr_bool = 0, gaussopt = [0,-1.0,1.0])

    if gaussopt[1] == 0
        error("q1D currently only works with gaussopt!") # but actually, the limitation only comes from the calculation of acoeff? -> should be possible with numerical integration
    end
    
    nmax = num_params.gem_params.nmax
    (nhomax,mho,omega_ho) = ho_params
    
    ntot = nmax*(nhomax+1);

    # preallocations:
    if cr_bool == 0
        nu_arr = zeros(nmax); #unaffected by ho!
        S  = zeros(ntot, ntot);
        T  = zeros(ntot, ntot);
        V  = zeros(ntot, ntot);
        energies = zeros(ntot);
        wavefunctions = zeros(ntot, ntot);
    elseif cr_bool == 1
        nu_arr = zeros(ComplexF64,2*nmax)
        S  = zeros(ComplexF64, 2*ntot, 2*ntot)
        T  = zeros(ComplexF64, 2*ntot, 2*ntot)
        V  = zeros(ComplexF64, 2*ntot, 2*ntot)
        energies = zeros(2*ntot);
        wavefunctions = zeros(ComplexF64, 2*ntot, 2*ntot);
    else
        error("Error with cr_bool = $cr_bool, only 0 or 1 allowed.")
    end
    
    prealloc_arrs = PreallocStruct(nu_arr,S,T,V,energies,wavefunctions)
    
    # function call with preallocations:
    GEM_solveq1D!(prealloc_arrs,phys_params,num_params,ho_params,wf_bool,cr_bool,gaussopt)
    
    if wf_bool == 0
        return energies
    elseif wf_bool == 1
        return energies, wavefunctions
    else
        error("error in wf_bool: only values of 0 or 1 allowed.")
    end
end

# function with preallocated arrays:
function GEM_solveq1D!(prealloc_arrs,phys_params,num_params,ho_params,wf_bool,cr_bool,gaussopt)
    
    (;nu_arr,S,T,V,energies,wavefunctions) = prealloc_arrs
    
    # Destructuring struct:
    (;hbar,mass_arr,vint_arr,lmax) = phys_params
    vint = vint_arr[1];
    (;gem_params,omega_cr,threshold) = num_params
    (;nmax,r1,rnmax) = gem_params

    (;nhomax,mho,omega_ho) = ho_params

    
    ## 1. Preliminaries: 
    # reduced mass:
    mur = 1 / (1/mass_arr[1] + 1/mass_arr[2]);
    # ho parameters:
    lambda = omega_ho*mur/hbar
    aho = 1/sqrt(2*lambda)
    
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
    buf = alloc_segbuf(Float64,typeof(nu_arr[1]),Float64) # besser mit if?
    
    # i think due to length information of nu_arr this fills only the matrix in the first ho-channel
    MatrixSq1D(S,lmax,nu_arr)
    MatrixTq1D(T,lmax,nu_arr,hbar,mur)
    MatrixVq1D(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt)
    
    # symmetric fill:
    S .= Hermitian(S,:L) # if Hermitian or Symmetric: type-unstable?
    T .= Hermitian(T,:L)
    V .= Hermitian(V,:L)

    # filling into HO-blocks; symmetry not used yet
    ngem = lastindex(nu_arr)
    for nhoC = 0:nhomax # for the columns
        for nhoR = 0:nhomax # for the rows
            ci = nhoC*ngem+1 # column-index in full matrix
            ri = nhoR*ngem+1 # row-index in full matrix

            if nhoC == nhoR
                T[ri:ri+ngem-1,ci:ci+ngem-1] .= T[1:ngem,1:ngem] .+ eho(nhoC,mho,hbar,omega_ho).* S[1:ngem,1:ngem] # adding ho energy to hamiltonian
                S[ri:ri+ngem-1,ci:ci+ngem-1] .= S[1:ngem,1:ngem]
            end

            V[ri:ri+ngem-1,ci:ci+ngem-1] .= V[1:ngem,1:ngem] .* acoeff(nhoC,nhoR,mho,gaussopt[3],lambda,gamma_dict)
        end
    end

    T .+= V # T becomes H

    
    ## 3. Eigensolver
    if wf_bool == 0
        eigen2step(energies,T,S,threshold=threshold) # only energies
        return energies
    else
        eigen2step_valvec(energies,wavefunctions,T,S,threshold=threshold)
        return energies,wavefunctions
    end    
end

## Gaussian ranges:
@views @inbounds function buildnu(nu_arr,r1,rnmax,nmax)
    nu_arr[1] = 1 /abs(r1)^2;
    nmax >1 && @. nu_arr[2:nmax] = 1/abs(r1)^2 * abs(r1/rnmax)^(2*((2:nmax)-1)/(nmax-1))
end

## ACoeff: copling coefficient for the different HO channels
#some tests indicate that its correct.
function acoeff(n,nprime,m,mu_g,lambda,gamma_dict)
    if lambda/mu_g == Inf
        if n==nprime
            return 1
        else
            return 0
        end
    else
        return (gamma_dict[n+abs(m)+1]*gamma_dict[nprime+abs(m)+1]/gamma_dict[n+1]/gamma_dict[nprime+1])^(1/2) * (mu_g/lambda)^(n+nprime)/(1+mu_g/lambda)^(n+nprime+abs(m)+1) * HypergeometricFunctions._₂F₁general(-nprime, -n, 1+abs(m), (lambda/mu_g)^2)/gamma_dict[1+abs(m)]
    end
end

# harmonic oscillator energy
function eho(nho,mho,hbar,omega_ho)
    return hbar*omega_ho*(2*nho+abs(mho)+1)
end


end ## end of module

