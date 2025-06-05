# functions to precompute the w-arrays for the interaction via interpolation and the upon-the shoulder method

function interpolNshoulder(phys_params,num_params,observ_params,size_params,precomp_arrs,interpol_arrs,wf_bool,csm_bool)
    
    # Destruct Structs:
    (;vint_arr) = phys_params
    (;lmax,Lmax,gem_params,theta_csm,omega_cr,mu0,c_shoulder,kmax_interpol,theta_csm) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;stateindices,centobs_arr,R2_arr) = observ_params
    (;cvals,central_indices,so_indices,maxlmax,nint_arr) = size_params
    (;gamma_dict,jmat,nu_arr,NU_arr) = precomp_arrs
    (;alpha_arr,v_arr,A_mat,w_arr,w_interpol_arr,Ainv_arr_kine,v_obs_arr,w_obs_arr,w_obs_interpol_arr) = interpol_arrs
    
    # range interpolation:
    precompute_alpha_arr(alpha_arr,r1,rnmax,R1,RNmax,nu_arr,NU_arr,jmat)
    
    # wrap function types:
    # Auto-wrap plain functions as central, for keeping backwards compatibility
    wrap_potential(f::Function) = CentralPotential(f)
    # Identity if already wrapped
    wrap_potential(p::PotentialFunction) = p
    vint_arr_wrapped = [wrap_potential.(potlist) for potlist in phys_params.vint_arr]
    centobs_arr_wrapped = [wrap_potential.(obslist) for obslist in centobs_arr]
    
    
    # upon-the-shoulder
    precompute_w(w_arr,v_arr,alpha_arr,A_mat,w_interpol_arr,Ainv_arr_kine,gamma_dict,maxlmax,mu0,c_shoulder,cvals,vint_arr_wrapped,centobs_arr_wrapped,w_obs_arr,v_obs_arr,w_obs_interpol_arr,wf_bool,csm_bool,theta_csm,central_indices,so_indices,nint_arr)
    
end


# alpha_arr for range-interpolation method
function precompute_alpha_arr(alpha_arr,r1,rnmax,R1,RNmax,nu_arr,NU_arr,jmat)
    # returns alpha_arr [1:kmax_interpol] via buildnu function for geometric series
    max_r = max(rnmax,RNmax) # should be roughly ok, can be improved if the mass-coefficients are included!
    min_r = min(r1,R1) # maybe we need indeed mass-depenendet max,min values in order to avoid extrapolation.
    
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
                                temp = det(tempA)/tempA[2,2]
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

# range-interpolation method

# central interaction
function precompute_varr!(v_arr,alpha_arr,Lsum,gamma_dict,vcent_fun::CentralPotential,buf,csmfac)
    for n = 0:Lsum
        for k=1:lastindex(alpha_arr)
            norm_interpol = 1/2 * gamma_dict[n+1.5]/alpha_arr[k]^(n+3/2)
            v_arr[k,n+1] = vcent_integration(vcent_fun,alpha_arr[k]*csmfac^2,n,buf)/gamma_dict[n+1.0]/norm_interpol * csmfac^(2*n+3)#584
        end
    end
end

# spin-orbit interaction
function precompute_varr!(v_arr,alpha_arr,Lsum,gamma_dict,vcent_fun::SpinOrbitPotential,buf,csmfac)
    for n = 0:Lsum # n=j
        nLS = n + 1 # for highlighting the difference due to LS
        for k=1:lastindex(alpha_arr)
            so_extrafac = 4*(1+n)/(2*n+3)/(2*n+2) # added missing factor 4!
            norm_interpol = 1/2 * gamma_dict[n+1.5]/alpha_arr[k]^(n+3/2)
            v_arr[k,n+1] = vcent_integration(vcent_fun,alpha_arr[k]*csmfac^2,nLS,buf)/gamma_dict[n+1.0]/norm_interpol * csmfac^(2*n+3) * so_extrafac
        end
    end
end

function vcent_integration(vcent_fun,alpha,n,buf) #where {V}
    val = quadgk(r -> integrand(r,alpha,n,vcent_fun),0,Inf;segbuf=buf, rtol=1e-6, atol=1e-10)[1]
end
function integrand(r,alpha,n,vcent_fun)
    return vcent_fun(r)*r^(2*n+2)*exp(-alpha*r^2)
end


### w_arr: upon-the-shoulder method
@views @inbounds function precompute_w(w_arr,v_arr,alpha_arr,A_mat,w_interpol_arr,Ainv_arr_kine,gamma_dict,maxlmax,mu0,c_shoulder,cvals,vint_arr,centobs_arr,w_obs_arr,v_obs_arr,w_obs_interpol_arr,wf_bool,csm_bool,theta_csm,central_indices,so_indices,nint_arr)
    # returns the Array w_arr[c in cvals,alpha=1:alphamax,Lsum = 1:2*maxlmax,n=1:Lsum+1]
    # note the +1 in the last argument: w_arr, v_arr ar NOT offset-arrays due to problems with linear algebra package. 
    
    log_alpha_range=range(log(alpha_arr[1]),log(alpha_arr[end]),lastindex(alpha_arr))
    
    kmax = lastindex(v_arr[:,1])    # number of alpha values for interpolation
    
    bufr = alloc_segbuf(Float64,Float64,Float64)
    bufc = alloc_segbuf(Float64,ComplexF64,Float64)

    csmfac = 1.0; buf = bufr;
    if csm_bool == 1
        csmfac = exp(im*theta_csm*pi/180)
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
        
        for cc in cvals
            #performance: precompute_varr needs more time than the interpolation procedure for w_arr below. this is mostly due to the use of quadgk
            
            for iv in nint_arr[cc]

                # numerical integration only necessary for central and spin-orbit interactions, (and more if added)
                if !(iv in central_indices[cc] || iv in so_indices[cc])
                    continue
                end

                # step 1: calculate v_arr via integration. v_arr is only temporary and can be reused for each interaction
                precompute_varr!(v_arr,alpha_arr,Lsum,gamma_dict,vint_arr[cc][iv],buf,csmfac)
                
                for k = 1:kmax
                    w_arr[cc,k,Lsum+1,1:Lsum+1] .= A_mat_curr\v_arr[k,1:Lsum+1] # Lsum+1, because w_arr is no OffsetArray!; also consider only subarray 0:Lsum = 1:Lsum+1
                end
                
                for n = 0:Lsum
                    w_interpol_arr[cc,iv,Lsum,n] = cubic_spline_interpolation(log_alpha_range, w_arr[cc,:,Lsum+1,n+1])
                end
            end            
        end
        
        
        ## for the observables: (separate loop, as there could be observables for Jacobi sets without interaction!)
        for cco = 1:3
            wf_bool == 0 && continue # no need for observables if wf_bool = 0
            for (jj,obs) in enumerate(centobs_arr[cco]) # for each c there might be several observables -> loop over them
                
                # cal v_arr for observables:
                precompute_varr!(v_obs_arr,alpha_arr,Lsum,gamma_dict,obs,bufr,1.0)
                
                # w_arr for observables:
                for k = 1:kmax
                    w_obs_arr[cco,jj,k,Lsum+1,1:Lsum+1] .= A_mat_curr\v_obs_arr[k,1:Lsum+1]
                end
                
                # w_interpol_arr for observables:
                for n = 0:Lsum
                    w_obs_interpol_arr[cco,jj,Lsum,n] = cubic_spline_interpolation(log_alpha_range, w_obs_arr[cco,jj,:,Lsum+1,n+1])
                end
                
            end
        end
        
    end
    
end

