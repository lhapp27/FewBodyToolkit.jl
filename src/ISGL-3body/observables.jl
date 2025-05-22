## Function for calculating the observables within the ISGL program (similar to fillTVS)

@views @inbounds function calc_observables(num_params,observ_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,result_arrs,csm_bool)

    if csm_bool == 1
        error("Observables are not implemented for CSM: csm_bool = $csm_bool != 0")
    end
    
    (;gem_params,mu0,c_shoulder) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;stateindices,centobs_arr,R2_arr) = observ_params
    (;abvals_arr,cvals,groupindex_arr,abI,factor_bf,box_size_arr,starts,ends,bvalsdiag,maxlmax,mij_arr_dict) = size_params
    (;spintrafo_dict,facsymm_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr,Clmk_arr,Dlmk_arr,S_arr) = precomp_arrs
    (;alpha_arr,v_arr,A_mat,Ainv_arr_kine,w_obs_interpol_arr) = interpol_arrs
    (;w_arr_kine,kij_arr,gij_arr,T,V,S,temp_args_arr,temp_fill_mat,wn_obs_interpol_arr) = fill_arrs
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
    
    flati = lastindex(temp_args_arr)
    
    ## Norm:
    norms = zeros(lastindex(stateindices)); #now has length 2 for stateindex =[1,4]! # old: zeros(stateindices[end]) # stateindex 1-4: length 4!
    for (sik,si) in enumerate(stateindices) # e.g. stateindices = [1,4]: sik = 1,2; si = 1,4
        norms[sik] = real.(wavefun_arr[:,si]' * S * wavefun_arr[:,si])
    end
    
    ## Calculation of matrix elements and matrix filling via 1d loop for CENTRAL observables:
    # no parallelization:
    for cco = 1:3 # cvals for centobs_arr
        isempty(centobs_arr[cco]) && continue
        for (obsind,obs) in enumerate(centobs_arr[cco]) # ich glaube obsind ist ausreichend, oder? obs steckt ja in interpolation drin.
            for (sik,si) in enumerate(stateindices) # for which states should the observables be calculated?
                
                temp_obs = 0.0
                for index in 1:flati # sum over contributions of all basis functions
                    rowi,coli=temp_args_arr[index]
                    
                    obs_elem = 2*oab(jmat,gij_arr,mu0,c_shoulder,w_obs_interpol_arr,wn_obs_interpol_arr,temp_args_arr[index],cco,obsind,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict)
                    #@show([cco,rowi,coli,obs_elem])
                    
                    rowi == coli && (obs_elem *= 1/2) # only lower triangular matrix, therefore global factor 2, but not on the diagonal!
                    
                    temp_obs += real(wavefun_arr[rowi,si]' * obs_elem * wavefun_arr[coli,si])
                    
                end
                centobs_output[cco,obsind,sik] = temp_obs/norms[sik] 
            end
            
        end
    end
    
    
    ## For (non-central) R^2 observables:
    for cco = 1:3
        R2_arr[cco] == 0 && continue # calculate <R^2> only if indicated by R2_arr[cco] = 1
        
        for (sik,si) in enumerate(stateindices)
            temp_R2 = 0.0
            for index in 1:flati
                rowi,coli=temp_args_arr[index]
                
                R2_elem = 2*R2ab(jmat,w_arr_kine,Ainv_arr_kine,kij_arr,mu0,c_shoulder,temp_args_arr[index],cco,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict)                            
                
                rowi == coli && (R2_elem *= 1/2) # only lower triangular matrix, therefore global factor 2, but not on the diagonal!
                
                temp_R2 += real(wavefun_arr[rowi,si]' * R2_elem * wavefun_arr[coli,si])
            end
            R2_output[cco,sik] = temp_R2/norms[sik]
        end
        
    end    
end


## functions to calculate one observables matrix element, summing over the necessary a,b values:
function oab(jmat,gij_arr,mu0,c_shoulder,w_obs_interpol_arr,wn_obs_interpol_arr,temp_args_i,cco,obsind,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict)
    #rowi,coli,ranges,norm4,mij_arr,S_arr,la,La,lb,Lb,Lsum,avals_new,bvals_new,cvals,factor_ab,avals,bvals = temp_args_i
    (;rowi,coli,avals,bvals,factor_ab,ranges,norm4,mij_arr,sa,JsSa,sb,JsSb,la,La,lb,Lb,Lsum,JlLa,JlLb) = temp_args_i
    
    factor_ab = 1; # für observablen, Symmetrieausnutzung schwieriger, deshalb: factor_ab -> 1, avals_new -> avals; factor_symm ist trotzdem wichtig für die korrekte symmetrisierung der WF
    
    tempO = 0.0

    (JsSa != JsSb || JlLa != JlLb) && return tempO #immediately skip if either one is violated

    for a in avals
        for b in bvals
            #factor_symm = facsymm(a,b,abI,la,lb,factor_bf)
            factor_symm = facsymm_dict[a,b,la,lb,sa,sb]
            # for spin-independent interactions,observables: overlap of spin functions can be calculated in Jacobi-set b:
            
            uabc = spintrafo_dict[a,b,JsSa,sa,sb]
            
            tempO += factor_ab*factor_symm*uabc*element_O(cco,ranges,norm4,jmat[a,cco],jmat[b,cco],mij_arr,S_arr,la,La,lb,Lb,gij_arr,mu0,c_shoulder,w_obs_interpol_arr,Lsum,wn_obs_interpol_arr,obsind,JlLa);
            #@show([cco,avals,bvals,a,b,la,lb,factor_symm,tempO])
        end
    end
    return tempO
end

## same as oab, but for R^2 (noncentral) observable
function R2ab(jmat,w_arr_kine,Ainv_arr_kine,kij_arr,mu0,c_shoulder,temp_args_i,cco,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict)
    #rowi,coli,ranges,norm4,mij_arr,S_arr,la,La,lb,Lb,Lsum,avals_new,bvals_new,cvals,factor_ab,avals,bvals = temp_args_i
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
            
            # is this correct even for the potential, which is evaluated in c-set?
            if JlLa == JlLb
                JlL = JlLa
            else
                #println("JlLa!=JlLb. this can actually happen!")
                continue
            end
            tempR2 += factor_ab*factor_symm*uabc*element_R2(ranges,norm4,jmat[a,cco],jmat[b,cco],mij_arr,S_arr,la,La,lb,Lb,w_arr_kine,Ainv_arr_kine,kij_arr,mu0,c_shoulder,Lsum,JlLa);
        end
    end
    return tempR2
end


# calculation of a single matrix element: central observable O(r_c)
function element_O(cco,ranges,norm4,jac,jbc,mij_arr_i,S_arr,la,La,lb,Lb,gij_arr,mu0,c_shoulder,w_obs_interpol_arr,Lsum,wn_obs_interpol_arr,obsind,JlL)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    # alpha <-> gamma; beta <-> delta    
    
    etac  = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2
    zetac = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2 
    xic = nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC
    etaprc = etac - xic^2/zetac
    
    # p, pprime, f, q:
    # unclear how to relate to AetaAB
    p = SA[nua*betaAC,NUa*deltaAC,nub*betaBC,NUb*deltaBC]
    ppr = SA[nua*alphaAC,NUa*gammaAC,nub*alphaBC,NUb*gammaBC]
    f = xic/zetac
    #fpr = xic/etac
    q = ppr .- f*p
    
    #gij_arr:
    for n = 0:Lsum
        mun = mu0*c_shoulder^n
        for i=1:3
            for j=(i+1):4
                gij_arr[i,j,n] = 2*p[i]*p[j]/zetac + 2*mun*q[i]*q[j]/etaprc
            end
        end
    end
    
    #gij_arrS = @SVector [@SVector [2*p[i]*p[j]/zetac + 2*mu0*c_shoulder^n*q[i]*q[j]/etaprc for (i,j) in ij_vals] for n=0:Lsum] # same problem as in T: size (i.e. Lsum) not known at compile-time!
    
    prefac = norm4 * (pi^2/zetac/etaprc)^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);
    
    # here, interpolation is called (outside of ii-loop):
    for n = 0:Lsum
        wn_obs_interpol_arr[n] = w_obs_interpol_arr[cco,obsind,Lsum,n](log(etaprc))
    end
    
    summe = 0.0
    for ii = 1:lastindex(mij_arr_i)
        (m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]
        
        sum2 = 0.0
        for n=0:Lsum
            sum2 += wn_obs_interpol_arr[n]*gij_arr[1,2,n]^m12*gij_arr[1,3,n]^m13*gij_arr[1,4,n]^m14*gij_arr[2,3,n]^m23*gij_arr[2,4,n]^m24*gij_arr[3,4,n]^m34 # abh von ii nur in S, welches unabh von n ist.... nein, mij hängen von ii ab!
        end        
        summe += S_arr[la,La,lb,Lb,JlL,ii]*sum2
    end
    
    return prefac*summe
end



# calculation of a single matrix element for R_c^2
function element_R2(ranges,norm4,jac,jbc,mij_arr_i,S_arr,la,La,lb,Lb,w_arr_kine,Ainv_arr_kine,kij_arr,mu0,c_shoulder,Lsum,JlL)
    ## copy from kinetic energy!
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    #alpha <-> gamma; beta <-> delta    
    
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
    
    for n = 0:Lsum
        mun = mu0*c_shoulder^n
        w_arr_kine[n] = K0*Ainv_arr_kine[Lsum,n]
        
        for i=1:3
            for j=(i+1):4
                kij_arr[i,j,n] = 2*p[i]*p[j]/zetac + 2*q[i]*q[j]/etaprc + mun*Kij_arr[i,j]/K0
            end
        end
    end
    
    #w_arr_kineS = @SVector [K0*Ainv_arr_kine[Lsum,n] for n=0:Lsum]    
    #kij_arrSa = @SVector [2*p[i]*p[j]/zeta + 2*q[i]*q[j]/etapr for (i,j) in ij_vals]
    #kij_arrS = @SVector[kij_arrSa .+ mu0*c_shoulder^n*Kij_arrS/K0 for n=0:Lsum] # problem: unable to know size (Lsum+1) at compile-time. Try with maxlmax dosnt help.
    
    prefac = norm4 * (pi^2/zetac/etaprc)^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);
    
    sum = 0.0
    for ii = 1:lastindex(mij_arr_i)
        (m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]
        
        sum2 = 0.0
        for n=0:Lsum
            sum2 +=w_arr_kine[n]*kij_arr[1,2,n]^m12*kij_arr[1,3,n]^m13*kij_arr[1,4,n]^m14*kij_arr[2,3,n]^m23*kij_arr[2,4,n]^m24*kij_arr[3,4,n]^m34
        end
        sum += S_arr[la,La,lb,Lb,JlL,ii]*sum2        
    end
    
    return prefac*sum
end