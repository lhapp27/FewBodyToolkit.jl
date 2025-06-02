# functions to precompute the w-arrays for the interaction via interpolation and the upon-the shoulder method

function interpolNshoulder(phys_params,num_params,observ_params,size_params,precomp_arrs,interpol_arrs,wf_bool,csm_bool)
    
    # Destruct Structs:
    (;vint_arr) = phys_params
    (;lmax,Lmax,gem_params,theta_csm,omega_cr,mu0,c_shoulder,kmax_interpol,theta_csm) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;stateindices,centobs_arr,R2_arr) = observ_params
    (;cvals,maxlmax,nint_arr) = size_params # add nint?
    (;gamma_dict,jmat,nu_arr,NU_arr) = precomp_arrs
    (;alpha_arr,v_arr,A_mat,w_arr,w_interpol_arr,Ainv_arr_kine) = interpol_arrs # ,v_obs_arr,w_obs_arr,w_obs_interpol_arr
    
    # range interpolation:
    precompute_alpha_arr(alpha_arr,r1,rnmax,R1,RNmax,nu_arr,NU_arr,jmat)   
    
    # upon-the-shoulder: dwnti 1D?
    precompute_w(w_arr,v_arr,alpha_arr,A_mat,w_interpol_arr,Ainv_arr_kine,gamma_dict,maxlmax,mu0,c_shoulder,cvals,vint_arr,centobs_arr,wf_bool,csm_bool,theta_csm,lmax,Lmax,nint_arr) # ,w_obs_arr,v_obs_arr,w_obs_interpol_arr
    
end


# alpha_arr for range-interpolation method
function precompute_alpha_arr(alpha_arr,r1,rnmax,R1,RNmax,nu_arr,NU_arr,jmat)
    # returns alpha_arr [1:kmax_interpol] via buildnu function for geometric series
    
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
    
    # to avoid falling out of interpolating interval by numerical inaccuracy. alternatively we could use extrapolation.
    tempmin -= 100*eps();
    tempmax += 100*eps();
    
    min_r_eff = 1/sqrt(tempmax)
    max_r_eff = 1/sqrt(tempmin)
    
    alpha_arr .= buildnu(max_r_eff,min_r_eff,lastindex(alpha_arr),alpha_arr) # alpha_min and _max interchanged, such that alpha_arr[1] is the smallest value!
    #@show([tempmin,tempmax,tempmin2,tempmax2])
end

# range-interpolation method
# maybe we can change this to not have vcent_fun as argument, but vcent_integration. in that case we can provide another vcent_integration function for the cases when we have an analytical formula!

## new functions with reduced allocations:

# central interaction # 1D change needs to be confirmed
function precompute_varr!(v_arr,alpha_arr,nnlist,gamma_dict,vcent_fun::PotentialFunction,buf,csmfac)
    for n in nnlist
        for k=1:lastindex(alpha_arr)
            v_arr[k,n+1] = vcent_integration(vcent_fun,alpha_arr[k]*csmfac^2,n,buf)*csmfac^(2*n+1) # v_arr is no offset-arr, hence n+1. v_arr[:,1] is for n=0, etc. 
        end
    end
end
# for central interaction one could integrate only from 0 to Inf?

function vcent_integration(vcent_fun,alpha,n,buf) #where {V}
    val = quadgk(r -> integrand(r,alpha,n,vcent_fun),-Inf,0,Inf;segbuf=buf, rtol=1e-6, atol=1e-10)[1] #-Inf for 1D; intermediate 0?
end

function integrand(r,alpha,n,vcent_fun) # this can be done similar as for 2-body for the different dimensions?
    return vcent_fun(r)*r^(n)*exp(-alpha*r^2) # removed +2 and factor 2 in front of n in exponent due to 1D!
end

## old 2:
#= vint_csm(r,vint,theta_csm) = vint(r*exp(im*theta_csm*pi/180)) # should not be used anymore. did we update csm in the basis functions for 1D?
# this is not yet fully dispatched as was done for the 3D case. needs update.
function precompute_varr!(v_arr,alpha_arr,Lsum,gamma_dict,vcent_fun::V,buf) where {V}
    for n = 0:Lsum
        for k=1:lastindex(alpha_arr)
            #norm_interpol = 1/2 * gamma_dict[n+1.5]/alpha_arr[k]^(n+3/2)
            v_arr[k,n+1] = vcent_integration(vcent_fun,alpha_arr[k],n,buf)#/gamma_dict[n+1.0]/norm_interpol #normlization unnecessary in 1D, as we dont do "upon-the-shoulder".
            #@show([n,alpha_arr[k]],v_arr[k,n+1])
        end
    end
end =#

## end of new functions





### w_arr: upon-the-shoulder method
@views @inbounds function precompute_w(w_arr,v_arr,alpha_arr,A_mat,w_interpol_arr,Ainv_arr_kine,gamma_dict,maxlmax,mu0,c_shoulder,cvals,vint_arr,centobs_arr,wf_bool,csm_bool,theta_csm,lmax,Lmax,nint_arr) #,w_obs_arr,v_obs_arr,w_obs_interpol_arr
    # returns the Array w_arr[c in cvals,alpha=1:alphamax,Lsum = 1:2*maxlmax,n=1:Lsum+1]
    # note the +1 in the last argument: w_arr, v_arr ar NOT offset-arrays due to problems with linear algebra package. 
    
    log_alpha_range=range(log(alpha_arr[1]),log(alpha_arr[end]),lastindex(alpha_arr))
    
    #Lsum_max = 2*maxlmax
    kmax = lastindex(v_arr[:,1])
    
    bufr = alloc_segbuf(Float64,Float64,Float64)
    bufc = alloc_segbuf(Float64,ComplexF64,Float64)
    
    for cc in cvals
        
        for iv in 1:nint_arr[cc]
            
            csmfac = 1.0
            buf = bufr

            if csm_bool == 1
                csmfac = exp(im*theta_csm*pi/180)
                buf = bufc
            end
            
            nnlist = 2*(lmax+Lmax):-2:0 # all necessary exponents LL-2s: 2*(lmax+Lmax), 2*(lmax+Lmax)-2, ..., 0

            precompute_varr!(v_arr,alpha_arr,nnlist,gamma_dict,vint_arr[cc][iv],bufr,csmfac) # via Lsum we get all possible powers of LL-2s?
            
            for nn in nnlist # possible exponents LL-2s: 2*(lmax+Lmax), 2*(lmax+Lmax)-2, ..., 0
                # since we dont need shoulder for 1D, we directly store the interpolation into w_interpol_arr
                w_interpol_arr[cc,iv,nn] = cubic_spline_interpolation(log_alpha_range, v_arr[:,nn+1]) # wild guess... careful with +1? reduce one dimension of w_interopl_arr
            end            
        end
        
        # observables are currently not supported in 1D
        
    end    
end

