# functions to precompute the w-arrays for the interaction via interpolation and the upon-the shoulder method

# complex_ranged is not read here: nu_arr already carries the complex ranges, so precompute_alpha_arr
# derives the |alpha| bounds of the sector from them, and theta_arr fixes its opening angle. The
# argument is kept for symmetry with GEM3B1D.interpolNshoulder and with the other stages.
function interpolNshoulder(phys_params,num_params,observ_params,size_params,precomp_arrs,interpol_arrs,return_wavefunctions::Bool,complex_scaling::Bool,complex_ranged::Bool=false)

    # Destruct Structs:
    (;interactions) = phys_params
    (;lmax,Lmax,gem_params,complex_scaling_angle,complex_range_freq,mu0,c_shoulder,kmax_interpol,complex_scaling_angle) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;stateindices,centobs_arr,R2_arr) = observ_params
    (;cvals,central_indices,so_indices,pow_indices,powopt_arr,maxlmax,nint_arr) = size_params
    (;gamma_dict,jmat,nu_arr,NU_arr) = precomp_arrs
    (;alpha_arr,theta_arr,alpha_grid,v_arr,A_mat,w_arr,w_interpol_arr,Ainv_arr_kine,v_pow,w_pow_arr,v_obs_arr,w_obs_arr,w_obs_interpol_arr) = interpol_arrs
    
    # range interpolation:
    # alpha_arr holds |alpha|; with complex ranges theta_arr adds the angular direction, so that the
    # mesh covers the sector in which etaprc lives (see theta_mesh). alpha_grid is the resulting mesh
    # of effective ranges at which the radial integrals are evaluated; for real ranges it is the
    # single column alpha_arr and everything below reduces to the previous, one-dimensional scheme.
    # only interpolated (central/spin-orbit) potentials are integrated at complex alpha, so only they
    # are subject to the combined complex-range / complex-scaling sector limit
    any(!isempty(central_indices[cc]) || !isempty(so_indices[cc]) for cc in cvals) &&
        check_cr_csm_sector(complex_ranged,complex_scaling,complex_range_freq,complex_scaling_angle)

    precompute_alpha_arr(alpha_arr,r1,rnmax,R1,RNmax,nu_arr,NU_arr,jmat)
    precompute_alpha_grid!(alpha_grid,alpha_arr,theta_arr)

    # wrap function types:
    # Auto-wrap plain functions as central, to ensure compatibility
    wrap_potential(f::Function) = CentralPotential(f)
    # Identity if already wrapped
    wrap_potential(p::PotentialFunction) = p
    vint_arr_wrapped = [wrap_potential.(potlist) for potlist in phys_params.interactions]
    centobs_arr_wrapped = [wrap_potential.(obslist) for obslist in centobs_arr]
    
    
    # upon-the-shoulder
    precompute_w(w_arr,v_arr,alpha_arr,theta_arr,alpha_grid,A_mat,w_interpol_arr,Ainv_arr_kine,v_pow,w_pow_arr,gamma_dict,maxlmax,mu0,c_shoulder,cvals,vint_arr_wrapped,centobs_arr_wrapped,w_obs_arr,v_obs_arr,w_obs_interpol_arr,return_wavefunctions,complex_scaling,complex_scaling_angle,central_indices,so_indices,pow_indices,powopt_arr,nint_arr)
    
end


# alpha_arr for range-interpolation method
function precompute_alpha_arr(alpha_arr,r1,rnmax,R1,RNmax,nu_arr,NU_arr,jmat)
    # returns alpha_arr [1:kmax_interpol] via buildnu function for geometric series
    max_r = max(rnmax,RNmax) # should be roughly ok, can be improved if the mass-coefficients are included!
    min_r = min(r1,R1) # maybe we need indeed mass-dependent max,min values in order to avoid extrapolation.
    
    # for calculating the minimum and maximum values of r and then effectively for etac which is necessary in the interpolation, in order to avoid extrapolation:
    # maybe a bit heavy numerically?
    i=0;
    tempmin=tempmax=0.0
    for a in 1:3
        for b in 1:3
            for c in 1:3
                for nua in nu_arr
                    for NUa in NU_arr
                        for nub in nu_arr
                            for NUb in NU_arr
                                i+=1
                                Aa = SA[nua 0.0 ; 0.0 NUa]
                                Ab = SA[nub 0.0 ; 0.0 NUb]
                                tempA = transpose(jmat[a,c])*Aa*jmat[a,c] + transpose(jmat[b,c])*Ab*jmat[b,c]
                                temp = abs(det(tempA)/tempA[2,2]) # abs is a no-op for real ranges (etaprc>0); for complex ranges it is the modulus of etaprc, i.e. the radial coordinate of the sector
                                i==1 && (tempmin = temp; tempmax=temp;)
                                temp > tempmax && (tempmax = temp)
                                temp < tempmin && (tempmin = temp)
                            end
                        end
                    end
                end
            end
        end
    end
    # to avoid falling out of interpolating interval by numerical inaccuracy
    tempmin -= 100*eps();
    tempmax += 100*eps();   
    
    min_r_eff = 1/sqrt(tempmax)-100*eps()#added eps for the same reason as above
    max_r_eff = 1/sqrt(tempmin)+100*eps()
    
    alpha_arr .= buildnu(max_r_eff,min_r_eff,lastindex(alpha_arr),alpha_arr) # alpha_min and _max interchanged, such that alpha_arr[1] is the smallest value!
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

# range-interpolation method

# central interaction
function precompute_varr!(v_arr,alpha_grid,Lsum,gamma_dict,vcent_fun::CentralPotential,buf,csmfac)
    for n = 0:Lsum
        for kt = 1:size(alpha_grid,2)
            for k = 1:size(alpha_grid,1)
                alpha = alpha_grid[k,kt]
                # alpha^(n+3/2) stays on the principal branch: |arg(alpha)| <= atan(omega) < pi/2
                norm_interpol = 1/2 * gamma_dict[n+1.5]/alpha^(n+3/2)
                v_arr[k,kt,n+1] = vcent_integration(vcent_fun,alpha*csmfac^2,n,buf)/gamma_dict[n+1.0]/norm_interpol * csmfac^(2*n+3)#584
            end
        end
    end
end

#= # spin-orbit interaction; postponed to future version
function precompute_varr!(v_arr,alpha_arr,Lsum,gamma_dict,vcent_fun::SpinOrbitPotential,buf,csmfac)
    for n = 0:Lsum # n=j
        nLS = n + 1 # for highlighting the difference due to LS
        for k=1:lastindex(alpha_arr)
            so_extrafac = 4*(1+n)/(2*n+3)/(2*n+2) # added missing factor 4!
            norm_interpol = 1/2 * gamma_dict[n+1.5]/alpha_arr[k]^(n+3/2)
            v_arr[k,n+1] = vcent_integration(vcent_fun,alpha_arr[k]*csmfac^2,nLS,buf)/gamma_dict[n+1.0]/norm_interpol * csmfac^(2*n+3) * so_extrafac
        end
    end
end =#

# power-law interaction: the integral is analytic, so no numerical integration is needed.
# The main interaction path does not use this (see the pow_indices block in precompute_w), but
# observables (centobs_arr) go through precompute_varr! as well, so a method is provided here.
function precompute_varr!(v_arr,alpha_grid,Lsum,gamma_dict,vcent_fun::PowerLawPotential,buf,csmfac)
    (;v0,p) = vcent_fun
    for n = 0:Lsum
        for kt = 1:size(alpha_grid,2)
            for k = 1:size(alpha_grid,1)
                alpha = alpha_grid[k,kt]
                norm_interpol = 1/2 * gamma_dict[n+1.5]/alpha^(n+3/2)
                vcent_analytic = v0*csmfac^(-p) * 1/2*gamma(n+(3+p)/2)/alpha^(n+(3+p)/2)
                v_arr[k,kt,n+1] = vcent_analytic/gamma_dict[n+1.0]/norm_interpol
            end
        end
    end
end

function vcent_integration(vcent_fun,alpha,n,buf) #where {V}
    val = quadgk(r -> integrand(r,alpha,n,vcent_fun),0,Inf;segbuf=buf)[1]
end
function integrand(r,alpha,n,vcent_fun)
    return vcent_fun(r)*r^(2*n+2)*exp(-alpha*r^2)
end


### w_arr: upon-the-shoulder method
@views @inbounds function precompute_w(w_arr,v_arr,alpha_arr,theta_arr,alpha_grid,A_mat,w_interpol_arr,Ainv_arr_kine,v_pow,w_pow_arr,gamma_dict,maxlmax,mu0,c_shoulder,cvals,interactions,centobs_arr,w_obs_arr,v_obs_arr,w_obs_interpol_arr,return_wavefunctions::Bool,complex_scaling::Bool,complex_scaling_angle,central_indices,so_indices,pow_indices,powopt_arr,nint_arr)
    # returns the Array w_arr[c in cvals,alpha=1:alphamax,theta=1:thetamax,Lsum = 1:2*maxlmax,n=1:Lsum+1]
    # note the +1 in the last argument: w_arr, v_arr are NOT offset-arrays due to problems with linear algebra package. 
    
    log_alpha_range=range(log(alpha_arr[1]),log(alpha_arr[end]),lastindex(alpha_arr))
    theta_range=range(theta_arr[1],theta_arr[end],lastindex(theta_arr))
    # complex ranges put etaprc into a sector of the complex plane, so the interpolant runs over
    # (log|alpha|, arg alpha) instead of over log(alpha) alone. Everything up to the construction of
    # the interpolants is common to both cases: the mesh alpha_grid simply has one column for real
    # ranges and ntheta of them otherwise.
    complex_ranged = eltype(alpha_grid) <: Complex
    ntheta = lastindex(theta_arr)
    alpha_grid_obs = reshape(alpha_arr,:,1) # observables stay on the real mesh (see sanity_checks)
    
    kmax = lastindex(alpha_arr)    # number of alpha values for interpolation
    
    bufr = alloc_segbuf(Float64,Float64,Float64)
    bufc = alloc_segbuf(Float64,ComplexF64,Float64)

    csmfac = 1.0; buf = bufr;
    complex_ranged && (buf = bufc) # complex mesh -> complex integrand, even without complex scaling
    if complex_scaling
        csmfac = exp(-im*complex_scaling_angle*pi/180)
        buf = bufc
    end
    
    for Lsum = 0:2*maxlmax
        
        # independent of cvals:
        for n = 0:Lsum
            for j = 0:Lsum
                A_mat[j+1,n+1] = (mu0*c_shoulder^n)^j/gamma_dict[j+1.0]
            end
        end
        
        A_mat_curr = A_mat[1:Lsum+1,1:Lsum+1] # currently (here, locally) used/necessary submatrix of A_mat
        Ainv = inv(A_mat_curr)
        
        # necessary for kine:
        if Lsum == 0
            Ainv_arr_kine[Lsum,0] = Ainv[0+1,0+1]
        else
            for n=0:Lsum
                Ainv_arr_kine[Lsum,n] = Ainv[n+1,0+1] + Ainv[n+1,1+1]
            end
        end
        
        # necessary for power-law interactions:
        # For V(r)=v0*r^p the radial integral of the range-interpolation method is closed-form,
        # Vtilde_n(alpha) = Gamma(n+(3+p)/2)/Gamma(n+3/2) * alpha^(-p/2), i.e. the alpha-dependence is a
        # single power alpha^(-p/2), the same for every n. It therefore factors out of the shoulder solve
        # and is applied later in element_VPow (via etaprc^(-p/2)); no interpolation over alpha is needed.
        # v0 is applied in element_VPow as well (as for the Gaussian), so only p enters here.
        for cc in cvals
            for iv in pow_indices[cc]
                p_pow = powopt_arr[cc][iv][2]
                for j = 0:Lsum
                    v_pow[j+1] = csmfac^(-p_pow) * gamma(j+(3+p_pow)/2)/gamma_dict[j+1.5]/gamma_dict[j+1.0]
                end
                for n = 0:Lsum # w = Ainv*v_pow, using the Ainv already formed above (as for Ainv_arr_kine)
                    temp_pow = zero(eltype(v_pow))
                    for j = 0:Lsum
                        temp_pow += Ainv[n+1,j+1]*v_pow[j+1]
                    end
                    w_pow_arr[cc,iv,Lsum,n] = temp_pow
                end
            end
        end
        
        for cc in cvals
            #performance: precompute_varr needs more time than the interpolation procedure for w_arr below. this is mostly due to the use of quadgk
            
            for iv = 1:nint_arr[cc]

                # numerical integration only necessary for central and spin-orbit interactions, (and more if added)
                if !(iv in central_indices[cc] || iv in so_indices[cc])
                    continue
                end

                # step 1: calculate v_arr via integration. v_arr is only temporary and can be reused for each interaction
                precompute_varr!(v_arr,alpha_grid,Lsum,gamma_dict,interactions[cc][iv],buf,csmfac)
                
                # step 2: the shoulder solve is a linear map of v at fixed alpha, so it is done
                # node by node and the angular direction just comes along
                for kt = 1:ntheta
                    for k = 1:kmax
                        w_arr[cc,k,kt,Lsum+1,1:Lsum+1] .= A_mat_curr\v_arr[k,kt,1:Lsum+1] # Lsum+1, because w_arr is no OffsetArray!; also consider only subarray 0:Lsum = 1:Lsum+1
                    end
                end
                
                for n = 0:Lsum
                    w_interpol_arr[cc,iv,Lsum,n] = complex_ranged ?
                        cubic_spline_interpolation((log_alpha_range,theta_range), w_arr[cc,:,:,Lsum+1,n+1]) :
                        cubic_spline_interpolation(log_alpha_range, w_arr[cc,:,1,Lsum+1,n+1])
                end
            end            
        end
        
        
        ## for the observables: (separate loop, as there could be observables for Jacobi sets without interaction!)
        for cco = 1:3
            !return_wavefunctions && continue # no need for observables if return_wavefunctions = false
            for (jj,obs) in enumerate(centobs_arr[cco]) # for each c there might be several observables -> loop over them
                
                # cal v_arr for observables:
                precompute_varr!(v_obs_arr,alpha_grid_obs,Lsum,gamma_dict,obs,bufr,1.0)
                
                # w_arr for observables:
                for k = 1:kmax
                    w_obs_arr[cco,jj,k,Lsum+1,1:Lsum+1] .= A_mat_curr\v_obs_arr[k,1,1:Lsum+1]
                end
                
                # w_interpol_arr for observables:
                for n = 0:Lsum
                    w_obs_interpol_arr[cco,jj,Lsum,n] = cubic_spline_interpolation(log_alpha_range, w_obs_arr[cco,jj,:,Lsum+1,n+1])
                end
                
            end
        end
        
    end
    
end

