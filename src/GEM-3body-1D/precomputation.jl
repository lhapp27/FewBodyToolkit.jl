# functions to precompute repetedly used arrays within the ISGL program

# contains:
# - 

#using PartialWaveFunctions, SpecialFunctions, QuadGK

function precompute_ISGL(phys_params,num_params,size_params,precomp_arrs)
    
    #Destructing Structs:    
    (;mass_arr) = phys_params
    (;lmax,Lmax,gem_params,kmax_interpol) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;lL_complete,l_complete,nl) = size_params
    (;gamma_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr) = precomp_arrs
    
    # how to update the values?!
    #precompute_gamma(precomp_arrs.gamma_dict,max(nmax,2*max(lmax,Lmax)+1)) #?!
    #gamma_dict = precompute_gamma(gamma_dict,max(nmax,2*max(lmax,Lmax)+1))
    precompute_gamma(gamma_dict,max(nmax,3*max(lmax,Lmax)+3))
    
    #jmat = precompute_jmat(jmat,mass_arr)
    precompute_jmat(jmat,mass_arr)
    
    #murR_arr = precompute_murR(murR_arr,mass_arr)
    precompute_murR(murR_arr,mass_arr)
    
    #nu_arr,NU_arr = precompute_ranges(nu_arr,NU_arr,r1,rnmax,nmax,R1,RNmax,Nmax)
    precompute_ranges(nu_arr,NU_arr,r1,rnmax,nmax,R1,RNmax,Nmax)
    
    #norm_arr,NORM_arr = precompute_norms(norm_arr,NORM_arr,nu_arr,NU_arr,nl,l_complete,gamma_dict)
    precompute_norms1D(norm_arr,NORM_arr,nu_arr,NU_arr,nl,l_complete,gamma_dict)
    
end

### gamma functions
function precompute_gamma(gamma_dict,ende)
    # returns the dict gamma_dict with keys "n" for each value of n between 0.5 and ende in steps of 0.5
    for i=0.5:0.5:ende+0.5
        gamma_dict[i] = gamma(i)
    end
    #return gamma_dict
end


### jacobi transformation matrices:
function precompute_jmat(jmat,m_arr)
    for i = 1:3 # just calculate all combinations... better: i,j in cvals?
        for j=1:3
            jmat[i,j] = jcbtr(i,j,m_arr)
        end
    end
    #return jmat
end

function jcbtr(i::Int64,f::Int64,m_arr::Array{Float64}) # Für fixes m1,m2,m3. Alternativ: ma,mb,mc als arugment, dann muss es aber jedesmal entsprechend aufgerufen werden!    
    ma = m_arr[i];mb = m_arr[mod(i,3)+1];mc = m_arr[mod(i+1,3)+1]
    if f==i
        return SA[1.0 0.0;0.0 1.0]
    elseif mod(f-i,3)==1    # trafo a->b (+1), also if f = i+1
        return SA[-ma/(ma+mc) 1.0; -mc*(ma+mb+mc)/(ma+mc)/(mb+mc) -mb/(mb+mc)]
    elseif mod(f-i,3)==2    # trafo a->c (+2), also if f = i+2
        return SA[-ma/(ma+mb) -1.0; +mb*(ma+mb+mc)/(ma+mb)/(mb+mc) -mc/(mb+mc)]
    end
    error("something went wrong in jtrafo")    
end

### reduced masses:
function precompute_murR(murR_arr,m_arr)
    for b = 1:3# just calculate all combinations...
        ma = m_arr[mod(b+1,3)+1]#circshift(m_arr,1)[b];#m_arr[mod(b+1,3)+1]; # Ergibt b-1 in Möglichkeiten 1,2,3; mappt 1,2,3 auf 3,1,2.
        mb = m_arr[b];
        mc = m_arr[mod(b,3)+1]#circshift(m_arr,-1)[b];#m_arr[mod(b,3)+1]; # Ergibt b+1 in Möglichkeiten 1,2,3; mappt 1,2,3 auf 2,3,1
        murR_arr[1,b] = mc*ma/(mc+ma); # mur
        murR_arr[2,b] = mb*(ma+mc)/(ma+mb+mc); # muR
    end
    #return murR_arr
end


### ranges ###
function precompute_ranges(nu_arr,NU_arr,r1,rnmax,nmax,R1,RNmax,Nmax)
    # returns the array nu_arr[n=1:nmax] for each value of n. same for NU_arr[N=1:Nmax]
    nu_arr .= buildnu(r1,rnmax,nmax,nu_arr)
    NU_arr .= buildnu(R1,RNmax,Nmax,NU_arr)
    
    #return nu_arr,NU_arr
end

@views function buildnu(r1,rnmax,nmax,nu_arr)
    nu_arr[1] = 1 /r1^2;
    nmax >1 && @. nu_arr[2:nmax] = 1/r1^2 * (r1/rnmax)^(2*((2:nmax)-1)/(nmax-1))
    return nu_arr
end


### norms ###
function precompute_norms1D(norm_arr,NORM_arr,nu_arr,NU_arr,nl,l_complete,gamma_dict)
    # returns the array norm_arr[l,n] for each combination of n,l. same for NORM_arr
    
    #norm(nu,l) = (2*(2*nu)^(l+3/2)/gamma_dict[l+3/2])^(1/2);
    norm(nu,l) = ((2*nu)^(l+1/2)/gamma_dict[l+1/2])^(1/2); #adopted to 1D
    
    for n = 1:lastindex(nu_arr)
        for l in l_complete #lindex = 1:nl
            #l = l_complete[lindex]
            norm_arr[l,n] = norm(nu_arr[n],l)
        end
    end
    
    for N = 1:lastindex(NU_arr)
        for l in l_complete #lindex = 1:nl
            #l = l_complete[lindex]
            NORM_arr[l,N] = norm(NU_arr[N],l)
        end
    end
end
