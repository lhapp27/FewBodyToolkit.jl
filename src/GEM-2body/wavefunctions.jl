# functions to output the 2-body wavefunction (not the coefficient arrays), based on the coefficients and num_params input

"""
    wavefun_arr(r_arr, phys_params, num_params, wf_arr; cr_bool=0)

Calculates the 2-body wavefunction at specified positions based on input parameters and coefficients.

# Arguments
- `r_arr::Vector{Float64}`: Array of positions where the wavefunction is evaluated.
- `phys_params::NamedTuple`: Physical parameters, here we only need `lmax`and `dim`
- `num_params::NamedTuple`: Numerical parameters containing the information about the Gaussian ranges
- `wf_arr::Vector`: Eigenvector from the diagonalization routine. Can be complex-valued.

# Keyword Arguments
- `cr_bool::Int=0`: Determines whether to use complex-ranged Gaussians (`1`) or not (`0`). Defaults to `0`.

# Returns
- `Vector{Float64}`: The wavefunction values at the specified positions r_arr.

# Example
```julia
wavefun_arr(0.0:0.1:10.0, phys_params, num_params, wf_arr; cr_bool=0)
```
"""
function wavefun_arr(r_arr, phys_params, num_params, wf_arr; cr_bool=0)
    # r_arr: array of positions
    # num_params: input to determine gaussian widths
    # cr_bool: 1 or 0 determining whether to use complex-ranged gaussians or not
    # wf_arr: eigenvector from diagonalizing routine
    
    (;gem_params,omega_cr) = num_params
    (;nmax,r1,rnmax) = gem_params
    (;lmax,dim) = phys_params
    
    # gaussian ranges:
    if cr_bool == 0
        nu_arr = zeros(nmax);
    elseif cr_bool == 1
        nu_arr = zeros(ComplexF64,2*nmax);
    end
    GEM2B.buildnu(nu_arr,r1,rnmax,nmax)
    
    # complex-ranged gaussians?
    if cr_bool == 1
        nu_arr .*= (1+omega_cr*im)
        nu_arr[nmax+1:2*nmax] .= conj.(nu_arr[1:nmax])
    end
    
    psi_arr = zeros(typeof(nu_arr[1]),lastindex(r_arr))
    for (ii,r) in enumerate(r_arr)
        psi_arr[ii] = wavefun_point(r,nu_arr,wf_arr,lmax,dim)
    end
    
    return psi_arr
end


"""
    wavefun_point(r, nu_arr, wf_arr, ll, dim)

Compute the value of the 2-body wavefunction at a given position `r`.

# Arguments
- `r::Float64`: The radial position where the wavefunction is evaluated.
- `nu_arr::Vector`: Array of Gaussian widths.
- `wf_arr::Vector`: Coefficients of the wavefunction for the given state.
- `ll::Int`: Orbital angular momentum quantum number.
- `dim::Int`: Dimensionality of the system.

# Returns
- `Float64`: The computed value of the wavefunction at the specified position `r`.

# Example
```julia
nu_arr = [1.0, 0.0025]
wf_arr = [0.8,0.6]
psi_value = wavefun_point(1.0, nu_arr, wf_arr, 0, 3)
```
"""
function wavefun_point(r,nu_arr,wf_arr,ll,dim) # wf_arr is already selected for a given stateindex
    psi = 0.0
    for nn = 1:lastindex(nu_arr)
        GaussNorm = (2*(nu_arr[nn] +nu_arr[nn]')^(ll+dim/2)/gamma(ll+dim/2))^(1/2) #it would be better to call gamma only once
        dim == 1 && (GaussNorm *= 1/sqrt(2)) # additional factor in 1D)
        z = (r .^ ll) .* exp.(-nu_arr[nn] * r .^ 2)
        #@show(GaussNorm,nu_arr[nn],z)
        psi += GaussNorm * wf_arr[nn] * (r .^ ll) .* exp.(-nu_arr[nn] * r .^ 2);
    end
    return psi    
end