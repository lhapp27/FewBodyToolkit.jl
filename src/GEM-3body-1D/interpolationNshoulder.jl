# functions to precompute the w-arrays for the interaction via interpolation and the upon-the shoulder method

function interpolNshoulder(phys_params,num_params,observ_params,size_params,precomp_arrs,interpol_arrs,gaussopt,wf_bool,csm_bool,coulopt)
    #interpolNshoulder(alpha_arr,v_arr,A_mat,w_arr,w_interpol_arr,Ainv_arr_kine,gamma_dict,maxlmax,cvals,num_params,aac,gac,abc,gbc,nu_arr,NU_arr,jmat,murR_arr,vint_arr)
    
    # Destruct Structs:
    (;vint_arr) = phys_params
    (;lmax,Lmax,gem_params,theta_csm,omega_cr,mu0,c_shoulder,kmax_interpol,theta_csm) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;stateindices,centobs_arr,R2_arr) = observ_params
    (;cvals,maxlmax) = size_params
    (;gamma_dict,jmat,nu_arr,NU_arr) = precomp_arrs
    (;alpha_arr,v_arr,A_mat,w_arr,w_interpol_arr,Ainv_arr_kine,v_obs_arr,w_obs_arr,w_obs_interpol_arr) = interpol_arrs
    
    # range interpolation:
    precompute_alpha_arr(alpha_arr,r1,rnmax,R1,RNmax,nu_arr,NU_arr,jmat)
    #@time precompute_alpha_arr(alpha_arr,r1,rnmax,R1,RNmax,nu_arr,NU_arr,jmat) # ca 1% of precompute_w
    
    
    # upon-the-shoulder
    precompute_w(w_arr,v_arr,alpha_arr,A_mat,w_interpol_arr,Ainv_arr_kine,gamma_dict,maxlmax,mu0,c_shoulder,cvals,vint_arr,centobs_arr,w_obs_arr,v_obs_arr,w_obs_interpol_arr,gaussopt,wf_bool,csm_bool,theta_csm,coulopt,lmax,Lmax)
    #@time precompute_w(w_arr,v_arr,alpha_arr,A_mat,w_interpol_arr,Ainv_arr_kine,gamma_dict,maxlmax,mu0,c_shoulder,cvals,vint_arr,centobs_arr,w_obs_arr,v_obs_arr,w_obs_interpol_arr,gaussopt,wf_bool,csm_bool,theta_csm)
    
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
                #@show([a,b,c,tempmin2,tempmax2,temp2])
                #it turned out that nu is completely irrelevant --> big speedup. maybe it's possible to determine the min and max value even completely analytically, but this should be fast enough.
#=                  for nua in nu_arr
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
                end =#
            end
        end
    end

    tempmin=tempmin2/max(rnmax,RNmax)^2;tempmax=tempmax2/min(r1,R1)^2;

    # to avoid falling out of interpolating interval by numerical inaccuracy
    tempmin -= 20*eps();
    tempmax += 20*eps();
    
    #@show(log.([tempmin,tempmax]))
    
    # attempt for analytical approximation: fails, as it is too inaccurate
    #tempmax2 = (284*maximum(m_arr)^6*max(nu_arr[1],NU_arr[1])^2 - 24*minimum(m_arr)^6*min(nu_arr[end],NU_arr[end])^2)/(160*minimum(m_arr)*min(nu_arr[end],NU_arr[end]))
    #tempmin2 = (284*minimum(m_arr)^6*min(nu_arr[end],NU_arr[end])^2 - 24*maximum(m_arr)^6*max(nu_arr[1],NU_arr[1])^2)/(160*maximum(m_arr)*max(nu_arr[1],NU_arr[1]))
    

    min_r_eff = 1/sqrt(tempmax)#0.005;#
    max_r_eff = 1/sqrt(tempmin)#50.0;#
    #println("Achtung: min/max_r_eff hardcoded")
    
    alpha_arr .= buildnu(max_r_eff,min_r_eff,lastindex(alpha_arr),alpha_arr) # alpha_min and _max interchanged, such that alpha_arr[1] is the smallest value!
    #@show([tempmin,tempmax,tempmin2,tempmax2])
end

# range-interpolation method
# maybe we can change this to not have vcent_fun as argument, but vcent_integration. in that case we can provide another vcent_integration function for the cases when we have an analytical formula!

## new functions with reduced allocations:
vint_csm(r,vint,theta_csm) = vint(r*exp(im*theta_csm*pi/180))
function precompute_varr!(v_arr,alpha_arr,Lsum,gamma_dict,vcent_fun::V,buf) where {V}
    for n = 0:Lsum
        for k=1:lastindex(alpha_arr)
            #norm_interpol = 1/2 * gamma_dict[n+1.5]/alpha_arr[k]^(n+3/2)
            v_arr[k,n+1] = vcent_integration(vcent_fun,alpha_arr[k],n,buf)#/gamma_dict[n+1.0]/norm_interpol #normlization unnecessary in 1D, as we dont do "upon-the-shoulder".
            #@show([n,alpha_arr[k]],v_arr[k,n+1])
        end
    end
end

function integrand(r,alpha,n,vcent_fun)
    return vcent_fun(r)*r^(n)*exp(-alpha*r^2) # removed +2 and factor 2 in front of n in exponent due to 1D!
end

function vcent_integration(vcent_fun,alpha,n,buf) #where {V}
    val = quadgk(r -> integrand(r,alpha,n,vcent_fun),-Inf,0,Inf;segbuf=buf, rtol=1e-6, atol=1e-10)[1] #-Inf for 1D; intermediate 0?
end
## end of new functions


#= @views @inbounds function precompute_varr(v_arr,alpha_arr,Lsum,gamma_dict,vcent_fun)
    # returns the array v_arr[1:alphamax,0:Lsum] which contains the result of radial integration of vcent together with gaussian of range alpha and r^2n divided by the norm_interpol factor and n!
    
    for n = 0:Lsum
        for k=1:lastindex(alpha_arr)
            norm_interpol = 1/2 * gamma_dict[n+1.5]/alpha_arr[k]^(n+3/2)
            v_arr[k,n+1] = vcent_integration(vcent_fun,alpha_arr[k],n)/gamma_dict[n+1.0]/norm_interpol
        end
    end
    #return v_arr
end

function vcent_integration(vcent_fun,alpha,n)
    # returns the scalar result (for a single alpha and a single n) of integration of vcent times r^2n and times a gaussian of range alpha
    
    #println("alpha=",alpha,", n=",n)
    
    #integrand(r) = vcent_fun(r)*r^(2*n+2)*exp(-alpha*r^2)
    #val =  quadgk(r -> integrad(r),0,Inf)[1]#
    
    #val,err =  quadgk(r -> vcent_fun(r)*r^(2*n+2)*exp(-alpha*r^2),0,Inf)
    val =  quadgk(r -> vcent_fun(r)*r^(2*n+2)*exp(-alpha*r^2),0,Inf)[1]
    
    return val#,err
end

# how to include potential parameters?! global constants? function should have only 1 argument due to how its called in vcent_integration...

=#


### w_arr: upon-the-shoulder method
@views @inbounds function precompute_w(w_arr,v_arr,alpha_arr,A_mat,w_interpol_arr,Ainv_arr_kine,gamma_dict,maxlmax,mu0,c_shoulder,cvals,vint_arr,centobs_arr,w_obs_arr,v_obs_arr,w_obs_interpol_arr,gaussopt,wf_bool,csm_bool,theta_csm,coulopt,lmax,Lmax)
    # returns the Array w_arr[c in cvals,alpha=1:alphamax,Lsum = 1:2*maxlmax,n=1:Lsum+1]
    # note the +1 in the last argument: w_arr, v_arr ar NOT offset-arrays due to problems with linear algebra package. 
    
    log_alpha_range=range(log(alpha_arr[1]),log(alpha_arr[end]),lastindex(alpha_arr))
    
    #Lsum_max = 2*maxlmax
    kmax = lastindex(v_arr[:,1])
    
    bufr = alloc_segbuf(Float64,Float64,Float64)
    bufc = alloc_segbuf(Float64,ComplexF64,Float64)
    
    for Lsum = 0:2*(lmax+Lmax) # this should rather be 2*(lmax+Lmax) instead of 2*maxlmax??
        
        # Ainv is unnecessary in 1D:
#=         # independent of cvals:
        # needs numerical parameters: mu0,c_shoulder
        for n = 0:Lsum
            #mun_arr[n] = mu0*c_shoulder^n
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
        end =#
        
        for cc in cvals
            #performance: precompute_varr needs more time than the interpolation procedure for w_arr below. this is mainly/maybe due to the use of quadgk. at the moment i dont know how to improve on that.

println("summe über iv, also über alle v in vint_arr fehlt!!!")

            # in principle: if gaussopt[cc][1][1] == 1, then the interaction is gaussian and we dont need to do the numerical integration, nor the interpolation procedure. #second index in case we have several gaussians for one cc
            (gaussopt[cc][1][1] == 1 || coulopt[cc][1][1] == 1) && continue
            
            # step 1: calculate v_arr via integration. v_arr is only temporary and can be reused for each interaction
            if csm_bool == 0
                precompute_varr!(v_arr,alpha_arr,Lsum,gamma_dict,vint_arr[cc][1],bufr)
            elseif csm_bool == 1
                precompute_varr!(v_arr,alpha_arr,Lsum,gamma_dict,r->vint_csm(r,vint_arr[cc][1],theta_csm),bufc)
            end

            for n = 0:Lsum
                w_interpol_arr[cc,n] = cubic_spline_interpolation(log_alpha_range, v_arr[:,n+1]) # achtung, struktur von w_interpol_arr geändert!
            end
            
            ## FALLS wir interpolation verwenden wollen brauchen wir glaube ich nur v_arr für alle n = 0...Lsum --> in neuem array speichern? w_arr unnötig
            
#=             for k = 1:kmax
                #alternativ:
                #w_arr[cc,k,Lsum+1,1:Lsum+1] .= Ainv*v_arr[k,1:Lsum+1], but its slower
                w_arr[cc,k,Lsum+1,1:Lsum+1] .= A_mat_curr\v_arr[k,1:Lsum+1] # Lsum+1, because w_arr is no OffsetArray!; also consider only subarray 0:Lsum = 1:Lsum+1
            end
            
            for n = 0:Lsum
                w_interpol_arr[cc,Lsum,n] = cubic_spline_interpolation(log_alpha_range, w_arr[cc,:,Lsum+1,n+1])
            end   =#          
        end
        
        # observables are currently not supported in 1D
#=         ## for the observables: (separate loop, as there could be observables for Jacobi sets without interaction!)
        for cco = 1:3
            wf_bool == 0 && continue # no need for observables if wf_bool = 0
            for (jj,obs) in enumerate(centobs_arr[cco]) # for each c there might be several observables -> loop over them
                
                # cal v_arr for observables:
                precompute_varr!(v_obs_arr,alpha_arr,Lsum,gamma_dict,obs,bufr)
                
                # w_arr for observables:
                for k = 1:kmax
                    w_obs_arr[cco,jj,k,Lsum+1,1:Lsum+1] .= A_mat_curr\v_obs_arr[k,1:Lsum+1]
                end
                
                # w_interpol_arr for observables:
                for n = 0:Lsum
                    w_obs_interpol_arr[cco,jj,Lsum,n] = cubic_spline_interpolation(log_alpha_range, w_obs_arr[cco,jj,:,Lsum+1,n+1])
                end
                
            end
        end =#
        
    end

end

