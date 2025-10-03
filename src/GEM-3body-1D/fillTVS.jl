## Function for calculating the matrix elements and filling the matrices T,V,S within the GEM3B1D program

@views @inbounds function fill_TVS(num_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,csm_bool,hbar,debug_bool)
    
    (;gem_params,theta_csm) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;abvals_arr,cvals,groupindex_arr,abI,factor_bf,box_size_arr,starts,ends,bvalsdiag,lL_arr,maxlmax,gauss_indices,central_indices,contact1D_indices,gaussopt_arr,contact1Dopt_arr) = size_params
    (;gamma_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr) = precomp_arrs
    (;alpha_arr,v_arr,w_interpol_arr) = interpol_arrs
    (;wn_interpol_arr,T,V,S,temp_args_arr,temp_fill_mat) = fill_arrs
    

    # create a 1D array of NamedTuples to hold the arguments for the matrix element calculation. This also carries the matrix structure via row- and column indices.
    flati = flattento1Dloop(temp_args_arr,groupindex_arr,factor_bf,bvalsdiag,abvals_arr,lL_arr,nmax,Nmax,nu_arr,NU_arr,norm_arr,NORM_arr,starts,cvals)
    
    ## Calculation of matrix elements and filling via 1D loop:
    # Norm-Overlap S
    for index in 1:flati
        rowi,coli=temp_args_arr[index]
        temp_fill_mat[rowi,coli] = sab(jmat,temp_args_arr[index],abI,factor_bf,gamma_dict)
    end
    S .= Symmetric(temp_fill_mat,:L); # transpose fill
    

    # Kinetic energy T
    for index in 1:flati
        rowi,coli=temp_args_arr[index]
        temp_fill_mat[rowi,coli] = tab(jmat,murR_arr,temp_args_arr[index],abI,factor_bf,gamma_dict,hbar) # we can reuse the same temp_fill_mat
    end
    T .= Symmetric(temp_fill_mat,:L); # transpose fill:
    
    if csm_bool == 1
        T .*= exp(-2*im*theta_csm*pi/180)
    end
    

    # Interaction V
    for index in 1:flati
        rowi,coli=temp_args_arr[index]
        temp_fill_mat[rowi,coli] = vab(jmat,w_interpol_arr,wn_interpol_arr,temp_args_arr[index],abI,factor_bf,gamma_dict,gauss_indices,central_indices,contact1D_indices,gaussopt_arr,contact1Dopt_arr,csm_bool,theta_csm)
    end
    V .= Symmetric(temp_fill_mat,:L); # transpose fill:
    
    # Example usage for debugging
    if debug_bool == 1
        stp = min(9, size(T, 1))  # Adjust size_to_print as needed
        println("T:")
        display(T[1:stp,1:stp])
        println("V:")
        display(V[1:stp,1:stp])
        println("S:")
        display(S[1:stp,1:stp])
    end
    
    T .+= V # T becomes the full Hamiltonian matrix H = T + V

end


# factor to ensure proper (anti-)symmetrization in case of identical particles. Only needed for 2+1 systems
# Maybe we should precompute this into a small array or dict instead of calling it repeatedly?
function facsymm(a,b,abI,la,lb,factor_bf)
    #abI == 0 && return 1
    factor_symm = 1
    a == abI && (factor_symm *= factor_bf*(-1)^la)
    b == abI && (factor_symm *= factor_bf*(-1)^lb)
    return factor_symm
end



# calculation of a single matrix element: norm-overlap S. Summing over the necessary a,b,c values:
function sab(jmat,temp_args_i,abI,factor_bf,gamma_dict)
    (;avals_new,bvals_new,factor_ab,ranges,norm4,la,La,lb,Lb) = temp_args_i
    tempS = 0.0
    for a in avals_new # only for identical particles we sum over several values of a and b here.
        for b in bvals_new
            factor_symm = facsymm(a,b,abI,la,lb,factor_bf)
            tempS += factor_ab*factor_symm*element_S(ranges,norm4,jmat[a,b],la,La,lb,Lb,gamma_dict)
        end
    end
    return tempS
end

# calculation of a single matrix element: kinetic energy T
function tab(jmat,murR_arr,temp_args_i,abI,factor_bf,gamma_dict,hbar)
    (;avals_new,bvals_new,factor_ab,ranges,norm4,la,La,lb,Lb) = temp_args_i
    tempT = 0.0                                
    for a in avals_new
        for b in bvals_new
            factor_symm = facsymm(a,b,abI,la,lb,factor_bf)
            tempT += factor_ab*factor_symm*element_T(ranges,norm4,jmat[a,b],murR_arr,la,La,lb,Lb,b,gamma_dict,hbar)
        end
    end
    return tempT
end

# calculation of a single matrix element: interaction V(r_c)
function vab(jmat,w_interpol_arr,wn_interpol_arr,temp_args_i,abI,factor_bf,gamma_dict,gauss_indices,central_indices,contact1D_indices,gaussopt_arr,contact1Dopt_arr,csm_bool,theta_csm)
    (;avals_new,bvals_new,factor_ab,ranges,norm4,la,La,lb,Lb,cvals) = temp_args_i
    tempV = 0.0                                    
    for a in avals_new
        for b in bvals_new
            factor_symm = facsymm(a,b,abI,la,lb,factor_bf)
            for c in cvals
                
                for ivg in gauss_indices[c] # we need to provide these indices to the function.
                    gaussopt = gaussopt_arr[c][ivg] # Gaussian parameters are stored in gaussopt_arr
                    tempV += factor_ab*factor_symm*element_VGauss(c,ranges,norm4,jmat[a,c],jmat[b,c],la,La,lb,Lb,gaussopt,gamma_dict,csm_bool,theta_csm)
                end
                
                for ivc in central_indices[c]
                    tempV += factor_ab*factor_symm*element_V(c,ivc,ranges,norm4,jmat[a,c],jmat[b,c],la,La,lb,Lb,w_interpol_arr,wn_interpol_arr,gamma_dict)
                end

                for ivc1D in contact1D_indices[c]
                    contactopt = contact1Dopt_arr[c][ivc1D] # Contact potential parameters are stored in contact1Dopt_arr
                    tempV += factor_ab*factor_symm*element_VContact(c,ranges,norm4,jmat[a,c],jmat[b,c],la,La,lb,Lb,contactopt,gamma_dict,csm_bool,theta_csm)
                end
                
            end
        end
    end
    return tempV
end

# function for a single matrix element: norm-overlap S
function element_S(ranges,norm4,jab,la,La,lb,Lb,gamma_dict)
    
    mod(la+La+lb+Lb,2) == 1 && return 0.0 #parity of <a| and |b> must be the same, otherwise return 0
    
    (;nua,nub,NUa,NUb) = ranges;# named tuple unnecessary slow?
    (alpha,gamma,beta,delta) = jab # careful with order (gamma before beta)
    
    eta1 = nua*alpha^2 + NUa*gamma^2 + nub
    eta2 = 2*(nua*alpha*beta + NUa*gamma*delta)
    eta3 = nua*beta^2 + NUa*delta^2 + NUb
    
    sumk = 0.0
    for k = 0:la
        sumK = 0.0
        for K = 0:La
            sums = 0.0
            for s = 0:floor((k+K+Lb)/2)
                sums += gamma_dict[(la+La+lb+Lb-2*s+1)/2]/gamma_dict[s+1]/gamma_dict[k+K+Lb-2*s+1] * eta3^s*eta2^(k+K+Lb-2*s)/((eta1-eta2^2/4/eta3)^((la+La+lb+Lb-2*s+1)/2))
            end
            sumK += sums*gamma_dict[la+1]/(gamma_dict[la-k+1]*gamma_dict[k+1]) * gamma_dict[La+1]/(gamma_dict[La-K+1]*gamma_dict[K+1]) * alpha^(la-k) * beta^k * gamma^(La-K) * delta^K * gamma_dict[k+K+Lb+1]*(-1/2/eta3)^(k+K+Lb)
        end
        sumk += sumK
    end
    sumk *= norm4 * (pi/eta3)^(1/2)
    
    return sumk
end


# function for a single matrix element: kinetic energy T
function element_T(ranges,norm4,jab,murR_arr,la,La,lb,Lb,b,gamma_dict,hbar)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alpha,gamma,beta,delta) = jab # careful with order (gamma before beta)
    
    mur= murR_arr[1,b] # reduced mass of particles a,c
    muR= murR_arr[2,b] # reduced mass of particle b wrt. the pair (a,c)
    
    tr = -hbar^2/2/mur
    tR = -hbar^2/2/muR
    
    temp = 4*nub^2*element_S(ranges,norm4,jab,la,La,lb+2,Lb,gamma_dict) - 2*(1+2*lb)*nub*element_S(ranges,norm4,jab,la,La,lb,Lb,gamma_dict)
    lb >= 2 && ( temp += lb*(lb-1)*element_S(ranges,norm4,jab,la,La,lb-2,Lb,gamma_dict) )
    temp *= tr
    
    Temp = 4*NUb^2*element_S(ranges,norm4,jab,la,La,lb,Lb+2,gamma_dict) - 2*(1+2*Lb)*NUb*element_S(ranges,norm4,jab,la,La,lb,Lb,gamma_dict)
    Lb >= 2 && ( Temp += Lb*(Lb-1)*element_S(ranges,norm4,jab,la,La,lb,Lb-2,gamma_dict) )
    Temp *= tR
    
    return temp+Temp
end


# function for a single matrix element: interaction V(r_c) based on numerical integration and interpolation over effective Gaussian range.
function element_V(c,iv,ranges,norm4,jac,jbc,la,La,lb,Lb,w_interpol_arr,wn_interpol_arr,gamma_dict)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    #GEM review: alpha <-> gamma; beta <-> delta
    
    eta5 = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2
    eta6 = 2*(nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC)
    eta7 = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2
    
    etaeff = eta5 - eta6^2/4/eta7
    Lsum = la+La+lb+Lb

    for n = Lsum:-2:0
        wn_interpol_arr[n] = w_interpol_arr[c,iv,n](log(etaeff)) # interpolated numerical integration value for the effective Gaussian range etaeff
    end
    
    sumk = 0.0
    for k = 0:la
        sumK = 0.0
        for K = 0:La
            sumkp = 0.0
            for kp = 0:lb
                sumKp = 0.0
                for Kp = 0:Lb                    
                    KK = k+K+kp+Kp; LL = la+La+lb+Lb;
                    
                    sums = 0.0
                    for s = 0:floor(KK/2)
                        nn = Int(LL -2*s)
                        sums += eta7^s*eta6^(KK-2*s)/(gamma_dict[s+1]*gamma_dict[KK-2*s+1])*wn_interpol_arr[nn]
                    end
                    
                    # possibly costly due to many dict calls?
                    sumKp += sums * gamma_dict[la+1]/(gamma_dict[la-k+1]*gamma_dict[k+1]) * gamma_dict[La+1]/(gamma_dict[La-K+1]*gamma_dict[K+1]) * gamma_dict[lb+1]/(gamma_dict[lb-kp+1]*gamma_dict[kp+1]) * gamma_dict[Lb+1]/(gamma_dict[Lb-Kp+1]*gamma_dict[Kp+1]) * alphaAC^(la-k)*betaAC^k*gammaAC^(La-K)*deltaAC^K * alphaBC^(lb-kp)*betaBC^kp*gammaBC^(Lb-Kp)*deltaBC^Kp * gamma_dict[KK+1] * (pi/eta7)^(1/2)*(-1/2/eta7)^KK
                    
                end
                sumkp += sumKp
            end
            sumK += sumkp
        end
        sumk += sumK
    end
    sumk *= norm4
    
    return sumk
end



# calculation of a single matrix element: interaction VGauss(r_c)
function element_VGauss(c,ranges,norm4,jac,jbc,la,La,lb,Lb,gaussopt,gamma_dict,csm_bool,theta_csm)

    mod(la+La+lb+Lb,2) == 1 && return 0.0 #Gaussian potential is symmetric, hence parity of <a| and |b> must be the same, otherwise return 0
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    #GEM review: alpha <-> gamma; beta <-> delta
    
    v0,mu_g = gaussopt
    csm_bool == 1 && (mu_g *= exp(2*im*theta_csm*pi/180)) # apply CSM factor to the range of the Gaussian potential
    
    LL = la+La+lb+Lb;
        
    eta5 = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2 + mu_g
    eta6 = 2*(nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC)
    eta7 = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2
    
    sumk = 0.0
    for k = 0:la
        sumK = 0.0
        for K = 0:La
            sumkp = 0.0
            for kp = 0:lb
                sumKp = 0.0
                for Kp = 0:Lb                    
                    KK = k+K+kp+Kp;
                    
                    sums = 0.0
                    for s = 0:floor(KK/2)
                        sums += eta7^s*eta6^(KK-2s)/(gamma_dict[s+1]*gamma_dict[KK-2*s+1]) * gamma_dict[(LL-2*s+1)/2]/((eta5-eta6^2/4/eta7)^((LL-2*s+1)/2))
                    end
                    
                    sumKp += sums * gamma_dict[la+1]/(gamma_dict[la-k+1]*gamma_dict[k+1]) * gamma_dict[La+1]/(gamma_dict[La-K+1]*gamma_dict[K+1]) * gamma_dict[lb+1]/(gamma_dict[lb-kp+1]*gamma_dict[kp+1]) * gamma_dict[Lb+1]/(gamma_dict[Lb-Kp+1]*gamma_dict[Kp+1]) * alphaAC^(la-k)*betaAC^k*gammaAC^(La-K)*deltaAC^K * alphaBC^(lb-kp)*betaBC^kp*gammaBC^(Lb-Kp)*deltaBC^Kp * gamma_dict[KK+1] * (pi/eta7)^(1/2)*(-1/2/eta7)^KK
                    
                end
                sumkp += sumKp
            end
            sumK += sumkp
        end
        sumk += sumK
    end
    sumk *= norm4 * v0
    
    return sumk    
end

# calculation of a single matrix element: interaction VContact(r_c)
function element_VContact(c,ranges,norm4,jac,jbc,la,La,lb,Lb,contactopt,gamma_dict,csm_bool,theta_csm)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    
    v0,z0 = contactopt
    csm_bool == 1 && (z0 *= exp(2*im*theta_csm*pi/180)) # apply CSM factor to the position of the contact interaction
    
    LL = la+La+lb+Lb;
        
    eta5 = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2
    eta6 = 2*(nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC)
    eta7 = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2
    
    sumk = 0.0
    for k = 0:la
        sumK = 0.0
        for K = 0:La
            sumkp = 0.0
            for kp = 0:lb
                sumKp = 0.0
                for Kp = 0:Lb                    
                    KK = k+K+kp+Kp;

                    if z0 == 0
                        mod(KK,2) == 1 && continue
                        KK-LL != 0 && continue
                    end

                    sums = 0.0
                    for s = 0:floor(KK/2) #expression via 1F1 exists too
                        sums += eta7^s*eta6^(KK-2s)/(gamma_dict[s+1]*gamma_dict[KK-2*s+1]) * z0^(LL-2*s) *exp(-(eta5-eta6^2/4/eta7)*z0^2) #eta7^(-KK) in external factor, like for Gaussian
                    end
                    
                    sumKp += sums * gamma_dict[la+1]/(gamma_dict[la-k+1]*gamma_dict[k+1]) * gamma_dict[La+1]/(gamma_dict[La-K+1]*gamma_dict[K+1]) * gamma_dict[lb+1]/(gamma_dict[lb-kp+1]*gamma_dict[kp+1]) * gamma_dict[Lb+1]/(gamma_dict[Lb-Kp+1]*gamma_dict[Kp+1]) * alphaAC^(la-k)*betaAC^k*gammaAC^(La-K)*deltaAC^K * alphaBC^(lb-kp)*betaBC^kp*gammaBC^(Lb-Kp)*deltaBC^Kp * gamma_dict[KK+1] * (pi/eta7)^(1/2)*(-1/2/eta7)^KK
                    
                end
                sumkp += sumKp
            end
            sumK += sumkp
        end
        sumk += sumK
    end
    sumk *= norm4 * v0
    
    return sumk    
end


# returns flati and fills temp_args_arr
function flattento1Dloop(temp_args_arr,groupindex_arr,factor_bf,bvalsdiag,abvals_arr,lL_arr,nmax,Nmax,nu_arr,NU_arr,norm_arr,NORM_arr,starts,cvals)
    # Keep loop-structure and write necessary functions arguments for matrix-element-calculation into 1-dim array temp_args_arr
    flati = 0
    # Iterate over boxes:
    for boxC in groupindex_arr
        for boxR in groupindex_arr
            boxR < boxC && continue # fill only lower-triangular (boxes, not elements!) only works for real- or complex-symmetric, or hermitian matrices; needs some work for CSM and CR!
            
            # if there are some identical particles: we can ignore the sum over a-values and simply multiply by a factor which is equal to the number of a-values. ONLY ON THE BOX-DIAGONAL! (boxC = boxR)
            if boxC == boxR
                bvals_new = bvalsdiag[boxC] # for changing the role of a,b to be in line with lower triangular!
                factor_ab = lastindex(abvals_arr[boxR]) # factor for amount of a-values normally
                diag_bool = 1 # for skipping lower-triangular calculation within each box!
            else
                bvals_new = abvals_arr[boxC]
                factor_ab = 1
                diag_bool = 0
            end
            avals_new = abvals_arr[boxR] # a bit unneccessary, but for more consistent naming. effectively only used here in fillTVS
            avals = abvals_arr[boxR] # necessary for observables
            bvals = abvals_arr[boxC] # necessary for observables
            
            alphab = 0;
            for (lb,Lb) in lL_arr[bvals_new[1]] # why bvals_new[1]? -> because either only 1 value in it, or the particles are identical and all lL_arr[bvals_new[1,2,...]] are the same!
                for nb = 1:nmax
                    nub = nu_arr[nb]
                    normb = norm_arr[lb,nb]
                    for Nb = 1:Nmax
                        NUb = NU_arr[Nb]
                        NORMb = NORM_arr[Lb,Nb]
                        
                        alphab += 1                        
                        alpha = 0
                        
                        for (la,La) in lL_arr[avals_new[1]]
                            Lsum=0 #Int64((la+La+lb+Lb)/2) #unnecessary in 1D
                            
                            for na = 1:nmax
                                nua = nu_arr[na]
                                norma = norm_arr[la,na]
                                for Na = 1:Nmax
                                    NUa = NU_arr[Na]
                                    NORMa = NORM_arr[La,Na]
                                    
                                    alpha += 1
                                    diag_bool == 1 && alpha < alphab && continue # skip upper-triangular only on diagonal boxes!
                                    
                                    norm4 = norma*normb*NORMa*NORMb;
                                    ranges = (;nua,nub,NUa,NUb);
                                    
                                    # alpha, alphab are row and column indices within a given box defined by boxR,boxC. Now transforming into the row- and column indices of the big matrix:
                                    rowi = starts[boxR] + alpha - 1
                                    coli = starts[boxC] + alphab - 1
                                    
                                    flati += 1
                                    temp_args_arr[flati] = (;rowi,coli,ranges,norm4,la,La,lb,Lb,Lsum,avals_new,bvals_new,cvals,factor_ab,avals,bvals) # write all loop-index-depending function-arguments into a temporary 1D array
                                end
                            end
                        end
                    end
                end
            end                
        end
    end
    return flati
end