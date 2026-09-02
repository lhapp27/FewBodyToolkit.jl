# functions to precompute the w-arrays for the interaction via to be used in range-interpolation

function interpolNshoulder(phys_params,num_params,observ_params,size_params,precomp_arrs,interpol_arrs,return_wavefunctions::Bool,complex_scaling::Bool,complex_ranged::Bool=false)
    
    # Destruct Structs:
    (;interactions) = phys_params
    (;lmax,Lmax,gem_params,complex_scaling_angle,complex_range_freq,kmax_interpol) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;stateindices,centobs_arr,R2_arr) = observ_params
    (;cvals,maxlmax,lL_complete,central_indices) = size_params # add nint?
    (;gamma_dict,jmat,nu_arr,NU_arr) = precomp_arrs
    (;alpha_arr,theta_arr,alpha_grid,v_arr,w_interpol_arr) = interpol_arrs # ,v_obs_arr,w_obs_arr,w_obs_interpol_arr
    
    # effective Gaussian ranges for range-interpolation:
    # alpha_arr holds |alpha|; with complex ranges theta_arr adds the angular direction, so that the
    # mesh covers the sector in which etaeff lives (see theta_mesh). alpha_grid is the resulting mesh
    # of effective ranges at which the radial integrals are evaluated; for real ranges it is the
    # single column alpha_arr and everything below reduces to the previous, one-dimensional scheme.
    # only interpolated (central) potentials are integrated at complex alpha, so only they are subject
    # to the combined complex-range / complex-scaling sector limit
    any(!isempty(central_indices[cc]) for cc in cvals) &&
        check_cr_csm_sector(complex_ranged,complex_scaling,complex_range_freq,complex_scaling_angle)

    omega_cr = complex_ranged ? complex_range_freq : 0.0
    precompute_alpha_arr(alpha_arr,r1,rnmax,R1,RNmax,nu_arr,NU_arr,jmat,omega_cr)   
    precompute_alpha_grid!(alpha_grid,alpha_arr,theta_arr)
    
    # numerical integration and range-interpolation
    precompute_w(v_arr,alpha_arr,theta_arr,alpha_grid,w_interpol_arr,gamma_dict,maxlmax,cvals,interactions,centobs_arr,return_wavefunctions,complex_scaling,complex_scaling_angle,lmax,Lmax,central_indices,lL_complete) # ,w_obs_arr,v_obs_arr,w_obs_interpol_arr
    
end


# alpha_arr for range-interpolation method
function precompute_alpha_arr(alpha_arr,r1,rnmax,R1,RNmax,nu_arr,NU_arr,jmat,omega_cr=0.0)
    # returns alpha_arr [1:kmax_interpol] via buildnu function for geometric sequence
    
    # for calculating the minimum and maximum values of r and then effectively for etac which is necessary in the interpolation, in order to avoid extrapolation. maybe a bit heavy numerically?
    i=0;
    tempmin=tempmax=0.0
    j=0;
    tempmin2=tempmax2=0.0
    for a in 1:3
        for b in 1:3
            for c in 1:3
                j+=1
                tempA2 = transpose(jmat[a,c])*jmat[a,c] + transpose(jmat[b,c])*jmat[b,c]
                temp2 = det(tempA2)/tempA2[2,2]
                j==1 && (tempmin2 = temp2; tempmax2=temp2;)
                temp2 > tempmax2 && (tempmax2 = temp2)
                temp2 < tempmin2 && (tempmin2 = temp2)
                # testing showed that the min/max values are independent of the ranges
            end
        end
    end
    
    tempmin=tempmin2/max(rnmax,RNmax)^2;tempmax=tempmax2/min(r1,R1)^2;
    
    # complex ranges nu*(1 +/- i*omega) scale every range in modulus by sqrt(1+omega^2), and with it
    # the bounds on |etaeff| that this (range-independent) scan produces. Widening the window by that
    # factor at both ends covers it: checked numerically over the full basis for several mass sets,
    # where the upper bound grows by exactly sqrt(1+omega^2) and the lower one does not move at all.
    crfac = sqrt(1+omega_cr^2)
    tempmin /= crfac; tempmax *= crfac;
    
    # to avoid falling out of interpolating interval by numerical inaccuracy. alternatively we could use extrapolation.
    tempmin -= 100*eps();
    tempmax += 100*eps();
    
    min_r_eff = 1/sqrt(tempmax)
    max_r_eff = 1/sqrt(tempmin)
    
    buildnu!(max_r_eff,min_r_eff,lastindex(alpha_arr),alpha_arr) # alpha_min and _max interchanged, such that alpha_arr[1] is the smallest value!
end


# mesh of effective Gaussian ranges: |alpha| from alpha_arr, direction from theta_arr.
# For real ranges theta_arr == [0.0] and the mesh is exactly alpha_arr, kept real so that the
# real-range path is unchanged.
function precompute_alpha_grid!(alpha_grid::Matrix{Float64},alpha_arr,theta_arr)
    alpha_grid[:,1] .= alpha_arr
    return alpha_grid
end
function precompute_alpha_grid!(alpha_grid::Matrix{ComplexF64},alpha_arr,theta_arr)
    for kt = 1:lastindex(theta_arr)
        for k = 1:lastindex(alpha_arr)
            alpha_grid[k,kt] = alpha_arr[k]*cis(theta_arr[kt])
        end
    end
    return alpha_grid
end


# range-interpolation method:
# central interaction
#
# Under complex scaling the basis stays real (see MatrixT and element_VGauss) and only the
# Hamiltonian is rotated, so the required element is the potential at the rotated argument,
#
#     I(n,alpha) = int V(r*e^{i*theta}) r^n exp(-alpha*r^2) dr .
#
# Substituting r = u*e^{i*theta} and deforming the ray back onto the real axis (V analytic and
# decaying in the sector) turns this into the integral actually evaluated below,
#
#     I(n,alpha) = csmfac^(n+1) * int V(u) u^n exp(-alpha*csmfac^2*u^2) du ,   csmfac = e^{-i*theta},
#
# hence the prefactor csmfac^(n+1). Note that the integrand here carries r^n, not the r^(2n+2)
# of the 3D case in ISGL-3body, where the same derivation gives csmfac^(2n+3).
function precompute_varr!(v_arr,alpha_grid,nnlist,gamma_dict,vcent_fun::Union{Function,CentralPotential},buf,csmfac)
    for n in nnlist
        for kt = 1:size(alpha_grid,2)
            for k = 1:size(alpha_grid,1)
                v_arr[k,kt,n+1] = vcent_integration(vcent_fun,alpha_grid[k,kt]*csmfac^2,n,buf)*csmfac^(n+1) # v_arr is no offset-arr, hence n+1. v_arr[:,:,1] is for n=0, etc. 
            end
        end
    end
end
# for a symmetric interaction one could integrate only from 0 to Inf and multiply by 2. Not sure if we gain much from it.

function vcent_integration(vcent_fun,alpha,n,buf) #where {V}
    val = quadgk(r -> integrand(r,alpha,n,vcent_fun),-Inf,0,Inf;segbuf=buf)[1] #-Inf for 1D; intermediate 0
end

function integrand(r,alpha,n,vcent_fun) # this can be done similar as for 2-body for the different dimensions?
    return vcent_fun(r)*r^(n)*exp(-alpha*r^2) # removed +2 and factor 2 in front of n in exponent due to 1D!
end





### w_interpol_arr: range-interpolation
@views @inbounds function precompute_w(v_arr,alpha_arr,theta_arr,alpha_grid,w_interpol_arr,gamma_dict,maxlmax,cvals,interactions,centobs_arr,return_wavefunctions::Bool,complex_scaling::Bool,complex_scaling_angle,lmax,Lmax,central_indices,lL_complete) #,w_obs_arr,v_obs_arr,w_obs_interpol_arr
    # returns the Array w_interpol_arr[c in cvals,ivc in central_indices[c],nn in nnlist]
    
    log_alpha_range=range(log(alpha_arr[1]),log(alpha_arr[end]),lastindex(alpha_arr)) # log(alpha) to interpolate over uniform range of effective ranges alpha
    theta_range=range(theta_arr[1],theta_arr[end],lastindex(theta_arr))
    # complex ranges put etaeff into a sector of the complex plane, so the interpolant runs over
    # (log|alpha|, arg alpha) instead of over log(alpha) alone.
    complex_ranged = eltype(alpha_grid) <: Complex
    
    #Lsum_max = 2*maxlmax
    kmax = lastindex(alpha_arr)
    
    bufr = alloc_segbuf(Float64,Float64,Float64)
    bufc = alloc_segbuf(Float64,ComplexF64,Float64)
    
    csmfac = 1.0
    buf = bufr
    
    complex_ranged && (buf = bufc) # complex mesh -> complex integrand, even without complex scaling
    
    if complex_scaling
        csmfac = exp(-im*complex_scaling_angle*pi/180)
        buf = bufc
    end
    
    # determining nnlist: all necessary exponents nn = LL-2s; LL = la+La+lb+Lb; s = 0 : LL/2
    #LLlist = Int[] #not required
    nnlist = Int[]
    for (la,La) in lL_complete
        for (lb,Lb) in lL_complete
            LL = la+La+lb+Lb # LL = la+La+lb+Lb
            #push!(LLlist,LL)
            for nn in LL:-2:0 # nn = LL-2s
                push!(nnlist,nn)
            end
        end
    end
    #LLlist = unique(LLlist)
    nnlist = unique(nnlist)

    for cc in cvals
        
        for ivc in central_indices[cc]
            precompute_varr!(v_arr,alpha_grid,nnlist,gamma_dict,interactions[cc][ivc],buf,csmfac) #possible exponents via nnlist
            for nn in nnlist
                w_interpol_arr[cc,ivc,nn] = complex_ranged ?
                    cubic_spline_interpolation((log_alpha_range,theta_range), v_arr[:,:,nn+1]) :
                    cubic_spline_interpolation(log_alpha_range, v_arr[:,1,nn+1])
            end
        end
        
        # more loops could be added here for other PotentialTypes, in case they need numerical integration
        
        # no loop for Gaussian potentials, since they dont need numerical integration.        
        
        # observables are currently not supported in 1D        
    end
end

