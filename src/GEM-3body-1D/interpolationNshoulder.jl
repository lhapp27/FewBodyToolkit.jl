# functions to precompute the w-arrays for the interaction via to be used in range-interpolation

function interpolNshoulder(phys_params,num_params,observ_params,size_params,precomp_arrs,interpol_arrs,wf_bool,csm_bool)
    
    # Destruct Structs:
    (;vint_arr) = phys_params
    (;lmax,Lmax,gem_params,theta_csm,kmax_interpol) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;stateindices,centobs_arr,R2_arr) = observ_params
    (;cvals,maxlmax,lL_complete,central_indices) = size_params # add nint?
    (;gamma_dict,jmat,nu_arr,NU_arr) = precomp_arrs
    (;alpha_arr,v_arr,w_interpol_arr) = interpol_arrs # ,v_obs_arr,w_obs_arr,w_obs_interpol_arr
    
    # effective Gaussian ranges alpha_arr for range-interpolation:
    precompute_alpha_arr(alpha_arr,r1,rnmax,R1,RNmax,nu_arr,NU_arr,jmat)   
    
    # numerical integration and range-interpolation
    precompute_w(v_arr,alpha_arr,w_interpol_arr,gamma_dict,maxlmax,cvals,vint_arr,centobs_arr,wf_bool,csm_bool,theta_csm,lmax,Lmax,central_indices,lL_complete) # ,w_obs_arr,v_obs_arr,w_obs_interpol_arr
    
end


# alpha_arr for range-interpolation method
function precompute_alpha_arr(alpha_arr,r1,rnmax,R1,RNmax,nu_arr,NU_arr,jmat)
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
    
    # to avoid falling out of interpolating interval by numerical inaccuracy. alternatively we could use extrapolation.
    tempmin -= 100*eps();
    tempmax += 100*eps();
    
    min_r_eff = 1/sqrt(tempmax)
    max_r_eff = 1/sqrt(tempmin)
    
    buildnu!(max_r_eff,min_r_eff,lastindex(alpha_arr),alpha_arr) # alpha_min and _max interchanged, such that alpha_arr[1] is the smallest value!
end


# range-interpolation method:
# central interaction # 1D change needs to be confirmed
function precompute_varr!(v_arr,alpha_arr,nnlist,gamma_dict,vcent_fun::Union{Function,CentralPotential},buf,csmfac)
    for n in nnlist
        for k=1:lastindex(alpha_arr)
            v_arr[k,n+1] = vcent_integration(vcent_fun,alpha_arr[k]*csmfac^2,n,buf)*csmfac^(2*n+1) # v_arr is no offset-arr, hence n+1. v_arr[:,1] is for n=0, etc. 
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
@views @inbounds function precompute_w(v_arr,alpha_arr,w_interpol_arr,gamma_dict,maxlmax,cvals,vint_arr,centobs_arr,wf_bool,csm_bool,theta_csm,lmax,Lmax,central_indices,lL_complete) #,w_obs_arr,v_obs_arr,w_obs_interpol_arr
    # returns the Array w_interpol_arr[c in cvals,ivc in central_indices[c],nn in nnlist]
    
    log_alpha_range=range(log(alpha_arr[1]),log(alpha_arr[end]),lastindex(alpha_arr)) # log(alpha) to interpolate over uniform range of effective ranges alpha
    
    #Lsum_max = 2*maxlmax
    kmax = lastindex(v_arr[:,1])
    
    bufr = alloc_segbuf(Float64,Float64,Float64)
    bufc = alloc_segbuf(Float64,ComplexF64,Float64)
    
    csmfac = 1.0
    buf = bufr
    
    if csm_bool == 1
        csmfac = exp(im*theta_csm*pi/180)
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
            precompute_varr!(v_arr,alpha_arr,nnlist,gamma_dict,vint_arr[cc][ivc],bufr,csmfac) #possible exponents via nnlist
            for nn in nnlist
                w_interpol_arr[cc,ivc,nn] = cubic_spline_interpolation(log_alpha_range, v_arr[:,nn+1])
            end
        end
        
        # more loops could be added here for other PotentialTypes, in case they need numerical integration
        
        # no loop for Gaussian potentials, since they dont need numerical integration.        
        
        # observables are currently not supported in 1D        
    end
end

