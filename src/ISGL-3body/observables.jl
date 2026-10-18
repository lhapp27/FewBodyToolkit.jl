## Function for calculating the observables within the ISGL program (similar to fillTVS)

@views @inbounds function calc_observables(num_params,observ_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,result_arrs,complex_scaling::Bool)

    if complex_scaling
        error("Observables are not implemented for CSM: complex_scaling = $complex_scaling != false")
    end
    
    (;gem_params,mu0,c_shoulder) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;stateindices,centobs_arr,R2_arr) = observ_params
    (;abvals_arr,cvals,groupindex_arr,abI,factor_bf,box_size_arr,nbasis_total,starts,ends,bvalsdiag,maxlmax,mij_arr_dict) = size_params
    (;spintrafo_dict,facsymm_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr,Clmk_arr,Dlmk_arr,S_arr) = precomp_arrs
    (;alpha_arr,v_arr,A_mat,Ainv_arr_kine,w_obs_interpol_arr) = interpol_arrs
    (;T,V,S,basis) = fill_arrs
    (;energies_arr,wavefun_arr,centobs_output,R2_output) = result_arrs
    
    # error check:
    error_bool = 0
    if error_bool == 1
        for i = 1:lastindex(energies_arr)
            emcHc = energies_arr[i] .- wavefun_arr[:,i]'*T*wavefun_arr[:,i]
            if abs(emcHc) > 10^(-10)
                error("Error for eigenvector $i: E - c*H*c = $emcHc > 10^-10")
            end
        end
    end
    
    ## Norm:
    norms = zeros(lastindex(stateindices)); #now has length 2 for stateindex =[1,4]! # old: zeros(stateindices[end]) # stateindex 1-4: length 4!
    for (sik,si) in enumerate(stateindices) # e.g. stateindices = [1,4]: sik = 1,2; si = 1,4
        norms[sik] = real.(wavefun_arr[:,si]' * S * wavefun_arr[:,si])
    end
    
    ## Calculation of matrix elements for CENTRAL observables, one loop over the lower triangle for all states:
    temp_obs = zeros(lastindex(stateindices))
    for cco = 1:3 # cvals for centobs_arr
        isempty(centobs_arr[cco]) && continue
        for (obsind,obs) in enumerate(centobs_arr[cco]) # obs might not be needed here, obsind should be sufficient.
            temp_obs .= 0.0
            for coli in 1:nbasis_total, rowi in coli:nbasis_total # sum over contributions of all basis functions
                args = element_args(rowi,coli,basis,abvals_arr,bvalsdiag,nu_arr,NU_arr,norm_arr,NORM_arr,mij_arr_dict)
                
                obs_elem = 2*oab(jmat,mu0,c_shoulder,w_obs_interpol_arr,args,cco,obsind,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict)
                
                rowi == coli && (obs_elem *= 1/2) # only lower triangular matrix, therefore global factor 2, but not on the diagonal!
                
                for (sik,si) in enumerate(stateindices) # for which states should the observables be calculated?
                    temp_obs[sik] += real(wavefun_arr[rowi,si]' * obs_elem * wavefun_arr[coli,si])
                end
            end
            centobs_output[cco,obsind,:] .= temp_obs ./ norms
        end
    end
    
    
    ## For (non-central) R^2 observables:
    temp_R2 = zeros(lastindex(stateindices))
    for cco = 1:3
        R2_arr[cco] == 0 && continue # calculate <R^2> only if indicated by R2_arr[cco] = 1
        
        temp_R2 .= 0.0
        for coli in 1:nbasis_total, rowi in coli:nbasis_total
            args = element_args(rowi,coli,basis,abvals_arr,bvalsdiag,nu_arr,NU_arr,norm_arr,NORM_arr,mij_arr_dict)
            
            R2_elem = 2*R2ab(jmat,Ainv_arr_kine,mu0,c_shoulder,args,cco,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict)
            
            rowi == coli && (R2_elem *= 1/2) # only lower triangular matrix, therefore global factor 2, but not on the diagonal!
            
            for (sik,si) in enumerate(stateindices)
                temp_R2[sik] += real(wavefun_arr[rowi,si]' * R2_elem * wavefun_arr[coli,si])
            end
        end
        R2_output[cco,:] .= temp_R2 ./ norms
        
    end    
end


## functions to calculate one observables matrix element, summing over the necessary a,b values:
function oab(jmat,mu0,c_shoulder,w_obs_interpol_arr,temp_args_i,cco,obsind,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict)
    (;rowi,coli,avals,bvals,factor_ab,ranges,norm4,mij_arr,sa,JsSa,sb,JsSb,la,La,lb,Lb,Lsum,JlLa,JlLb) = temp_args_i

    factor_ab = 1; # for observables it is more difficult to make use of the symmetry, hence: factor_ab -> 1, avals_new -> avals; factor_symm is still necessary for the correct symmetrization of the WF

    tempO = 0.0

    (JsSa != JsSb || JlLa != JlLb) && return tempO #immediately skip if either one is violated

    for a in avals
        for b in bvals
            #factor_symm = facsymm(a,b,abI,la,lb,factor_bf)
            factor_symm = facsymm_dict[a,b,la,lb,sa,sb]
            # for spin-independent interactions,observables: overlap of spin functions can be calculated in Jacobi-set b:
            
            uabc = spintrafo_dict[a,b,JsSa,sa,sb]
            
            tempO += factor_ab*factor_symm*uabc*element_O(cco,ranges,norm4,jmat[a,cco],jmat[b,cco],mij_arr,S_arr,la,La,lb,Lb,mu0,c_shoulder,w_obs_interpol_arr,Lsum,obsind,JlLa);
        end
    end
    return tempO
end

## same as oab, but for R^2 (noncentral) observable
function R2ab(jmat,Ainv_arr_kine,mu0,c_shoulder,temp_args_i,cco,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict)
    (;rowi,coli,avals,bvals,factor_ab,ranges,norm4,mij_arr,sa,JsSa,sb,JsSb,la,La,lb,Lb,Lsum,JlLa,JlLb) = temp_args_i

    factor_ab = 1;
    
    tempR2 = 0.0    
    for a in avals
        for b in bvals
            #factor_symm = facsymm(a,b,abI,la,lb,factor_bf)
            factor_symm = facsymm_dict[a,b,la,lb,sa,sb]
            # for spin-independent interactions,observables: overlap of spin functions can be calculated in Jacobi-set b:
            if JsSa == JsSb
                uabc = spintrafo_dict[a,b,JsSa,sa,sb]
            else
                uabc = 0.0; continue # we dont need to calculate the matrix element any further since the spin parts are orthogonal
            end
        
            tempR2 += factor_ab*factor_symm*uabc*element_R2(ranges,norm4,jmat[a,cco],jmat[b,cco],mij_arr,S_arr,la,La,lb,Lb,Ainv_arr_kine,mu0,c_shoulder,Lsum,JlLa);
        end
    end
    return tempR2
end


# calculation of a single matrix element: central observable O(r_c)
function element_O(cco,ranges,norm4,jac,jbc,mij_arr_i,S_arr,la,La,lb,Lb,mu0,c_shoulder,w_obs_interpol_arr,Lsum,obsind,JlL)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    
    etac  = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2
    zetac = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2 
    xic = nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC
    etaprc = etac - xic^2/zetac
    
    # p, pprime, f, q:
    p = SA[nua*betaAC,NUa*deltaAC,nub*betaBC,NUb*deltaBC]
    ppr = SA[nua*alphaAC,NUa*gammaAC,nub*alphaBC,NUb*gammaBC]
    f = xic/zetac
    #fpr = xic/etac
    q = ppr .- f*p
    
    prefac = norm4 * (pi^2/zetac/etaprc)^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);
    
    # per n: interpolated weight wn (once per n, outside of the ii-loop) and the pair values gn (see IJ_PAIRS)
    summe = 0.0
    for n = 0:Lsum
        mun = mu0*c_shoulder^n
        wn = w_obs_interpol_arr[cco,obsind,Lsum,n](log(etaprc))
        gn = map(((i,j),) -> 2*p[i]*p[j]/zetac + 2*mun*q[i]*q[j]/etaprc, IJ_PAIRS)
        for ii = 1:lastindex(mij_arr_i)
            (m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]
            summe += S_arr[la,La,lb,Lb,JlL,ii]*wn*gn[1]^m12*gn[2]^m13*gn[3]^m14*gn[4]^m23*gn[5]^m24*gn[6]^m34
        end
    end
    
    return prefac*summe
end



# calculation of a single matrix element for R_c^2
function element_R2(ranges,norm4,jac,jbc,mij_arr_i,S_arr,la,La,lb,Lb,Ainv_arr_kine,mu0,c_shoulder,Lsum,JlL)
    ## copy from kinetic energy!
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    
    etac  = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2
    zetac = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2 
    xic = nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC
    etaprc = etac - xic^2/zetac
    zetaprc = zetac - xic^2/etac
    
    # p, pprime, f, q:
    # unclear how to relate to AetaAB
    p = SA[nua*betaAC,NUa*deltaAC,nub*betaBC,NUb*deltaBC]
    ppr = SA[nua*alphaAC,NUa*gammaAC,nub*alphaBC,NUb*gammaBC]
    f = xic/zetac
    fpr = xic/etac
    q = ppr .- f*p
    qpr = p .- fpr*ppr    # fixed typo!
    
    #we can reuse here some stuff from the kinetic energy
    K0 = 3/2/zetaprc    
    Kij_arr = 2/zetaprc^2*qpr.*transpose(qpr);
    
    prefac = norm4 * (pi^2/zetac/etaprc)^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);
    
    # per n: shoulder weight wn and the pair values kn (see IJ_PAIRS)
    sum = 0.0
    for n = 0:Lsum
        mun = mu0*c_shoulder^n
        wn = K0*Ainv_arr_kine[Lsum,n]
        kn = map(((i,j),) -> 2*p[i]*p[j]/zetac + 2*q[i]*q[j]/etaprc + mun*Kij_arr[i,j]/K0, IJ_PAIRS)
        for ii = 1:lastindex(mij_arr_i)
            (m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]
            sum += S_arr[la,La,lb,Lb,JlL,ii]*wn*kn[1]^m12*kn[2]^m13*kn[3]^m14*kn[4]^m23*kn[5]^m24*kn[6]^m34
        end
    end
    
    return prefac*sum
end