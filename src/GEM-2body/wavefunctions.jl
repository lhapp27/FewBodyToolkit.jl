# functions to output the 2-body wavefunction (not the coefficient arrays), based on the coefficients and num_params input
# 3D only!!
# 25.11.2024: Updated compatibility for 1D (selection via moduleX)

# wavefunction at positions r_arr:
function wavefun_arr(r_arr, moduleX, phys_params, num_params, cr_bool, wf_arr)
    # r_arr: array of positions
    # num_params: input to determine gaussian widths
    # cr_bool: 1 or 0 determining whether to use complex-ranged gaussians or not
    # wf_arr: eigenvector from diagonalizing routine
    
    (;gem_params,omega_cr,theta_csm,threshold) = num_params
    (;nmax,r1,rnmax) = gem_params
    (;lmax) = phys_params #lmax is actually just l. and in 1D it represents the parity.
    
    # gaussian ranges:
    if cr_bool == 0
        nu_arr = zeros(nmax);
    elseif cr_bool == 1
        nu_arr = zeros(ComplexF64,2*nmax);
    end
    moduleX.buildnu(nu_arr,r1,rnmax,nmax)
    
    # complex-ranged gaussians?
    if cr_bool == 1
        nu_arr .*= (1+omega_cr*im)
        nu_arr[nmax+1:2*nmax] .= conj.(nu_arr[1:nmax])
    end
    
    #@show(nu_arr)

    psi_arr = zeros(typeof(nu_arr[1]),lastindex(r_arr))
    if nameof(moduleX) == :GEM
        for (ii,r) in enumerate(r_arr)
            psi_arr[ii] = wavefun_point(r,nu_arr,wf_arr,lmax)
        end        
    elseif nameof(moduleX) == :GEM1D
        for (ii,r) in enumerate(r_arr)
            psi_arr[ii] = wavefun_point1D(r,nu_arr,wf_arr,lmax) # need different function due to different normalization!
            #@show(ii,psi_arr[ii])
        end
    end
    
    return psi_arr
end

# For 3D:
function wavefun_point(r,nu_arr,wf_arr,ll) # wf_arr is already selected for a given stateindex
    psi = 0.0
    for nn = 1:lastindex(nu_arr)
        GaussNorm = (2*(nu_arr[nn] + conj(nu_arr[nn]))^(ll + 3 / 2) / gamma(ll + 3 / 2))^(1/2)
        psi += GaussNorm * wf_arr[nn] * (r .^ ll) .* exp.(-nu_arr[nn] * r .^ 2);
    end
    return psi    
end

# For 1D:
function wavefun_point1D(r,nu_arr,wf_arr,ll) # wf_arr is already selected for a given stateindex
    psi = 0.0
    for nn = 1:lastindex(nu_arr)
        GaussNorm = ((nu_arr[nn] + conj(nu_arr[nn]))^(ll + 1 / 2) / gamma(ll + 1 / 2))^(1/2)
        psi += GaussNorm * wf_arr[nn] * (r .^ ll) .* exp.(-nu_arr[nn] * r .^ 2);
    end
    return psi    
end