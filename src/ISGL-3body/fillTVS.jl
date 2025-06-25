## Function for calculating the matrix elements and filling the matrices T,V,S within the ISGL program

@views @inbounds function fill_TVS(num_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,csm_bool,hbar)
    
    (;gem_params,mu0,c_shoulder,theta_csm) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;abvals_arr,cvals,gauss_indices,central_indices,so_indices,groupindex_arr,abI,factor_bf,box_size_arr,starts,ends,bvalsdiag,s_arr,JsS_arr,s_complete,JsS_complete,JlL_arr,lL_nested,maxlmax,mij_arr_dict,mijSO_arr_dict,gaussopt_arr) = size_params
    (;gamma_dict,spintrafo_dict,spinoverlap_dict,global6j_dict,facsymm_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr,Clmk_arr,Dlmk_arr,S_arr,SSO_arr) = precomp_arrs
    (;alpha_arr,v_arr,A_mat,w_interpol_arr,Ainv_arr_kine) = interpol_arrs
    (;w_arr_kine,wn_interpol_arr,kij_arr,gij_arr,T,V,S,temp_args_arr,temp_fill_mat) = fill_arrs
    
    # reducing complicated many-loop structure to a single 1D loop:
    flati = flattento1Dloop(temp_args_arr,groupindex_arr,bvalsdiag,abvals_arr,s_arr,JsS_arr,JlL_arr,lL_nested,nmax,Nmax,nu_arr,NU_arr,norm_arr,NORM_arr,mij_arr_dict,starts)
    
    ## Calculation of matrix elements and matrix filling via 1d loop (so we dont have redundant loops over all spin configs, lL combinations, and nmax)
    for index in 1:flati
        (;rowi,coli) = temp_args_arr[index]
        temp_fill_mat[rowi,coli] = sab(jmat,temp_args_arr[index],abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict)
    end
    # transpose fill:
    S .= Symmetric(temp_fill_mat,:L);
    
    #combined t and v for some speedup: removed combination to cope for CSM
    for index in 1:flati
        (;rowi,coli) = temp_args_arr[index]
        temp_fill_mat[rowi,coli] = tab(jmat,murR_arr,w_arr_kine,Ainv_arr_kine,kij_arr,mu0,c_shoulder,temp_args_arr[index],abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict,hbar)
    end
    # transpose fill:
    T .= Symmetric(temp_fill_mat,:L);
    
    if csm_bool == 1
        T .*= exp(-2*im*theta_csm*pi/180)
    end
    
    #v
    for index in 1:flati
        (;rowi,coli) = temp_args_arr[index]
        temp_fill_mat[rowi,coli] = vab(jmat,gij_arr,mu0,c_shoulder,w_interpol_arr,wn_interpol_arr,temp_args_arr[index],abI,factor_bf,S_arr,SSO_arr,cvals,spintrafo_dict,spinoverlap_dict,facsymm_dict,gauss_indices,central_indices,so_indices,s_arr,global6j_dict,mijSO_arr_dict,gaussopt_arr,csm_bool) # do we need hbar^2 for SO?
    end
    # transpose fill:
    V .= Symmetric(temp_fill_mat,:L);
    
    debug_bool = 0
    if debug_bool == 1
        stp = min(9, size(T, 1))  # Adjust size_to_print as needed
        println("T:")
        display(T[1:stp,1:stp])
        println("V:")
        display(V[1:stp,1:stp])
        #println("S:")
        #print_matrices(S, size_to_print)
    end
    
    T .+= V
    
end

# Debugging: Print matrices in a formatted way
function print_matrices(M, size_to_print)
    for i in 1:size_to_print, j in 1:size_to_print
        @printf("%10.3f ", real(M[i, j]))
        j == size_to_print && println()
    end
end

## functions to calculate one matrix element, summing over the necessary a,b,c values:
function sab(jmat,temp_args_i,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict)
    (;rowi,coli,avals_new,bvals_new,factor_ab,ranges,norm4,mij_arr,sa,JsSa,sb,JsSb,la,La,lb,Lb,JlLa,JlLb) = temp_args_i
    tempS = 0.0
    
    # maybe easier:
    (JsSa != JsSb || JlLa != JlLb) && return tempS #immediately skip if either one is violated
    
    for a in avals_new # this does not indicate the "boxes"!! just if in one box several a and b values are summed over (mostly due to identical particles)!
        for b in bvals_new
            #factor_symm = facsymm(a,b,abI,la,lb,factor_bf,spin_arr,sa,sb)
            factor_symm = facsymm_dict[a,b,la,lb,sa,sb]
            
            uab = spintrafo_dict[a,b,JsSa,sa,sb]
            JlL = JlLa
            
            tempS += factor_ab*factor_symm*uab*element_S(ranges,norm4,jmat[a,b],mij_arr,S_arr,la,La,lb,Lb,JlL)
            #@show(factor_ab,factor_symm,uab,sss)
            #println("sab: a=$a, b=$b, uab=$uab, factor_ab=$factor_ab, factor_symm=$factor_symm, sss=$sss")
        end
    end
    return tempS
end

function tab(jmat,murR_arr,w_arr_kine,Ainv_arr_kine,kij_arr,mu0,c_shoulder,temp_args_i,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict,hbar)
    #rowi,coli,ranges,norm4,mij_arr,S_arr,la,La,lb,Lb,Lsum,avals_new,bvals_new,cvals,factor_ab,factor_symm = temp_args_i # is this faster?
    (;avals_new,bvals_new,factor_ab,ranges,norm4,mij_arr,sa,JsSa,sb,JsSb,la,La,lb,Lb,Lsum,JlLa,JlLb) = temp_args_i
    tempT = 0.0
    
    # maybe easier:
    (JsSa != JsSb || JlLa != JlLb) && return tempT #immediately skip if either one is violated
    
    for a in avals_new
        for b in bvals_new
            #factor_symm = facsymm(a,b,abI,la,lb,factor_bf,spin_arr,sa,sb)
            factor_symm = facsymm_dict[a,b,la,lb,sa,sb]
            
            uab = spintrafo_dict[a,b,JsSa,sa,sb]
            JlL = JlLa
            
            tempT += hbar^2*factor_ab*factor_symm*uab*element_T(ranges,norm4,jmat[a,b],murR_arr,mij_arr,S_arr,la,La,lb,Lb,w_arr_kine,Ainv_arr_kine,kij_arr,b,mu0,c_shoulder,Lsum,JlL)
        end
    end
    return tempT
end

function vab(jmat,gij_arr,mu0,c_shoulder,w_interpol_arr,wn_interpol_arr,temp_args_i,abI,factor_bf,S_arr,SSO_arr,cvals,spintrafo_dict,spinoverlap_dict,facsymm_dict,gauss_indices,central_indices,so_indices,s_arr,global6j_dict,mijSO_arr_dict,gaussopt_arr,csm_bool)
    (;rowi,coli,avals_new,bvals_new,factor_ab,ranges,norm4,mij_arr,sa,JsSa,sb,JsSb,la,La,lb,Lb,Lsum,JlLa,JlLb) = temp_args_i
    #println("vab: $rowi, $coli, JsSa=$JsSa, JsSb=$JsSb")
    tempV = 0.0
    
    
    for a in avals_new
        for b in bvals_new
            #factor_symm = facsymm(a,b,abI,la,lb,factor_bf,spin_arr,sa,sb)
            factor_symm = facsymm_dict[a,b,la,lb,sa,sb]
            
            # for spin-independent interactions: overlap of spin functions can be calculated in advance:
            uab = spintrafo_dict[a,b,JsSa,sa,sb]
            
            for c in cvals
                
                spinoverlap = spinoverlap_dict[a,b,c,JsSa,JsSb,sa,sb]
                
                #@show([uab,spinoverlap,global6jfac])
                
                for ivg in gauss_indices[c] #loop over the gaussian interactions for this c.
                    (JsSa != JsSb || JlLa != JlLb) && return tempV #immediately skip if it is violated.
                    v0,mu_g = gaussopt_arr[c][ivg]

                    tempV += factor_ab*factor_symm*uab*element_VGauss(c,ranges,norm4,jmat[a,c],jmat[b,c],mij_arr,S_arr,la,La,lb,Lb,gij_arr,mu0,c_shoulder,Lsum,wn_interpol_arr,v0,mu_g,JlLa)
                end
                for ivc in central_indices[c] #loop over the central interactions for this c.
                    (JsSa != JsSb || JlLa != JlLb) && return tempV #immediately skip if it is violated.
                    tempV += factor_ab*factor_symm*uab*element_V(c,ranges,norm4,jmat[a,c],jmat[b,c],mij_arr,S_arr,la,La,lb,Lb,gij_arr,mu0,c_shoulder,w_interpol_arr,Lsum,wn_interpol_arr,JlLa,ivc)
                end
                for ivso in so_indices[c] # loop over the spin-orbit interactions for this c.
                    (abs(JlLa-JlLb) <= 1 <= JlLa+JlLb) == false && return tempV #immediately skip if it is violated.
                    (abs(JsSa-JsSb) <= 1 <= JsSa+JsSb) == false && return tempV #immediately skip if it is violated.
                    #only relevant for SO interactions
                    global6jfac = global6j_dict[JsSa,JsSb,JlLa,JlLb]
                    tempV += factor_ab*factor_symm*global6jfac*spinoverlap*element_VSO(c,ranges,norm4,jmat[a,c],jmat[b,c],mijSO_arr_dict,SSO_arr,la,La,lb,Lb,gij_arr,mu0,c_shoulder,w_interpol_arr,Lsum,wn_interpol_arr,JlLa,JlLb,ivso)
                end
                
            end
            
        end
    end
    #@show(tempV)
    return tempV
end


function element_S(ranges,norm4,jab,mij_arr_i,S_arr,la,La,lb,Lb,JlL)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alpha,gamma,beta,delta) = jab # careful with order (gamma before beta)
    
    eta  = nua*alpha^2 + NUa*gamma^2 + nub*1.0^2 + NUb*0.0
    zeta = nua*beta^2 + NUa*delta^2 + nub*0.0 + NUb*1.0^2 
    xi = nua*alpha*beta + NUa*gamma*delta + nub*1.0*0.0 + NUb*0.0*1.0        
    etapr = eta - xi^2/zeta
    
    # p, pprime, f, q:
    # unclear how to relate to AetaAB
    p = SA[nua*beta,NUa*delta,0.0,NUb]
    ppr = SA[nua*alpha,NUa*gamma,nub,0.0]
    f = xi/zeta
    q = ppr .- f*p
    
    fij_arr = @SMatrix[2*p[i]*p[j]/zeta + 2*q[i]*q[j]/etapr for i=1:4, j=1:4]    
    prefac = norm4 * (pi^2/zeta/etapr)^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);
    
    sum = 0.0
    for ii = 1:lastindex(mij_arr_i)
        (m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]      
        sum += S_arr[la,La,lb,Lb,JlL,ii]*fij_arr[1,2]^m12*fij_arr[1,3]^m13*fij_arr[1,4]^m14*fij_arr[2,3]^m23*fij_arr[2,4]^m24*fij_arr[3,4]^m34 * prefac
    end   
    
    return sum
end


# calculation of a single matrix element: kinetic energy T
function element_T(ranges,norm4,jab,murR_arr,mij_arr_i,S_arr,la,La,lb,Lb,w_arr_kine,Ainv_arr_kine,kij_arr,b,mu0,c_shoulder,Lsum,JlL)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alpha,gamma,beta,delta) = jab # careful with order (gamma before beta)
    
    mur= murR_arr[1,b]
    muR= murR_arr[2,b]
    
    # global factor hbar^2 outside of the function.
    tr = -1/2/mur
    tR = -1/2/muR
    
    eta  = nua*alpha^2 + NUa*gamma^2 + nub*1.0^2 + NUb*0.0^2
    zeta = nua*beta^2 + NUa*delta^2 + nub*0.0^2 + NUb*1.0^2 
    xi = nua*alpha*beta + NUa*gamma*delta + nub*1.0*0.0 + NUb*0.0*1.0        
    etapr = eta - xi^2/zeta
    zetapr = zeta - xi^2/eta    
    
    # p, pprime, f, q:
    # unclear how to relate to AetaAB
    p = SA[nua*beta,NUa*delta,0.0,NUb]
    ppr = SA[nua*alpha,NUa*gamma,nub,0.0]
    f = xi/zeta
    fpr = xi/eta
    q = ppr .- f*p
    qpr = p .- fpr*ppr    # fixed typo!
    
    K0 = 6*tr*nub*(nub/etapr - 1.0) + 6*tR*NUb*(NUb/zetapr - 1.0)    
    
    Kij_arr = @SMatrix[8*tr*nub^2/etapr*(q[i]*q[j]/etapr - q[1]*I[i,1]*I[j,3] - q[2]*I[i,2]*I[j,3] - q[4]*I[i,3]*I[j,4]) + 8*tR*NUb^2/zetapr*(qpr[i]*qpr[j]/zetapr - qpr[1]*I[i,1]*I[j,4] - qpr[2]*I[i,2]*I[j,4] - qpr[3]*I[i,3]*I[j,4]) for i=1:4, j=1:4]
    
    for n = 0:Lsum
        mun = mu0*c_shoulder^n
        w_arr_kine[n] = K0*Ainv_arr_kine[Lsum,n]
        
        for i=1:3
            for j=(i+1):4
                kij_arr[i,j,n] = 2*p[i]*p[j]/zeta + 2*q[i]*q[j]/etapr + mun*Kij_arr[i,j]/K0
            end
        end
    end
        
    prefac = norm4 * (pi^2/zeta/etapr)^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);
    
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

# calculation of a single matrix element: interaction V(r_c)
function element_V(c,ranges,norm4,jac,jbc,mij_arr_i,S_arr,la,La,lb,Lb,gij_arr,mu0,c_shoulder,w_interpol_arr,Lsum,wn_interpol_arr,JlL,iv)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    #there exist other definitions! alpha <-> gamma; beta <-> delta    
    
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
        
    prefac = norm4 * (pi^2/zetac/etaprc)^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);
    
    # here, interpolation is called (outside of ii-loop):
    for n = 0:Lsum
        wn_interpol_arr[n] = w_interpol_arr[c,iv,Lsum,n](log(etaprc))
    end
    
    summe = 0.0
    for ii = 1:lastindex(mij_arr_i) # the correct mij_arr is already selected in the flatento1Dloop function
        (m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]
        
        sum2 = 0.0
        for n=0:Lsum
            sum2 += wn_interpol_arr[n]*gij_arr[1,2,n]^m12*gij_arr[1,3,n]^m13*gij_arr[1,4,n]^m14*gij_arr[2,3,n]^m23*gij_arr[2,4,n]^m24*gij_arr[3,4,n]^m34 # abh von ii nur in S, welches unabh von n ist.... nein, mij hängen von ii ab!
        end
        
        summe += S_arr[la,La,lb,Lb,JlL,ii]*sum2
    end
    
    return prefac*summe
end

# I really think we need a second function. I dont see how we can combine it efficiently into a single one.
function element_VSO(c,ranges,norm4,jac,jbc,mijSO_arr_dict,SSO_arr,la,La,lb,Lb,gij_arr,mu0,c_shoulder,w_interpol_arr,Lsum,wn_interpol_arr,JlLa,JlLb,ivso)
    
    #prechecks: they should ideally never trigger, due to proper handling before. not sure if they cost performance
    Lsum < 1 && return 0.0
    #(abs(JlLa-JlLb) <= 1 <= JlLa+JlLb) == false && return 0.0 # done outside already?
    LsumSO = Lsum - 1
    
    (;nua,nub,NUa,NUb) = ranges; # is this usage of named tuple slow?
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    #careful of other definitions! 
    
    etac  = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2
    zetac = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2 
    xic = nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC
    etaprc = etac - xic^2/zetac
    
    # p, pprime, f, q:
    # unclear how to relate to AetaAB
    p = SA[nua*betaAC,NUa*deltaAC,nub*betaBC,NUb*deltaBC]
    ppr = SA[nua*alphaAC,NUa*gammaAC,nub*alphaBC,NUb*gammaBC]
    f = xic/zetac
    fso = (nub*betaBC*alphaBC + NUb*deltaBC*gammaBC)/zetac
    
    q = ppr .- f*p
    kron34 = SA[0,0,1,1]
    gg = ppr.*kron34 .- fso*p
    
    
    #gij_arr:
    for n = 0:LsumSO
        mun = mu0*c_shoulder^n
        for i=1:3
            for j=(i+1):4
                gij_arr[i,j,n] = 2*p[i]*p[j]/zetac + 2*mun*q[i]*q[j]/etaprc
            end
        end
    end
    
    #hij_arr:
    #hij_arr0 = @SMatrix[p[i]*ppr[j] for i=1:4, j=1:4]*fso # first version. incorrect.
    hij_arr = @SMatrix[q[i]*gg[j] - q[j]*gg[i] for i=1:4, j=1:4]*(-2) #corrected
        
    prefac = norm4 * (pi^2/zetac/etaprc)^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);
    
    # here, interpolation is called (outside of ii-loop):
    for n = 0:LsumSO
        wn_interpol_arr[n] = w_interpol_arr[c,ivso,LsumSO,n](log(etaprc))
    end
    
    ivjv_arr = @SVector[(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)] # can also be used for iuju
    
    summe = 0.0 # total sum, could also be named sumv, since sum over v (or iv,jv)
    for (v,(iv,jv)) in enumerate(ivjv_arr)
        mij_arr_i = mijSO_arr_dict[((la,La),(lb,Lb),v)] # this is a different approach. for Vcent, the mij_arr_i is already selected within the flattento1Dloop function.
        
        sumi = 0.0 # sum over i
        for ii = 1:lastindex(mij_arr_i) # imax
            #(m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]
            
            sumn = 0.0 # sum over n
            for n=0:LsumSO
                
                # a bit closer to the formulae:
                produ = 1.0 # product over u (or iu,ju)
                for (u,(iu,ju)) in enumerate(ivjv_arr)
                    produ *= gij_arr[iu,ju,n]^mij_arr_i[ii][u]
                end
                sumn += wn_interpol_arr[n]*produ
            end
            sumi += SSO_arr[la,La,lb,Lb,JlLa,JlLb,v,ii]*sumn
        end
        summe += sumi * hij_arr[iv,jv]
    end
    return prefac*summe
end





# calculation of a single matrix element: interaction VGauss(r_c)
function element_VGauss(c,ranges,norm4,jac,jbc,mij_arr_i,S_arr,la,La,lb,Lb,gij_arr,mu0,c_shoulder,Lsum,wn_interpol_arr,v0,mu_g,JlL)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    
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
    ggij_arr = @SMatrix[2*p[i]*p[j]/zetac + 2*q[i]*q[j]/(etaprc + mu_g) for i=1:4, j=1:4]
    
    prefac = norm4 * (pi^2/zetac/(etaprc + mu_g))^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);
    
    summe = 0.0
    for ii = 1:lastindex(mij_arr_i)
        (m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]
        
        summe += S_arr[la,La,lb,Lb,JlL,ii]*ggij_arr[1,2]^m12*ggij_arr[1,3]^m13*ggij_arr[1,4]^m14*ggij_arr[2,3]^m23*ggij_arr[2,4]^m24*ggij_arr[3,4]^m34
    end
    
    return v0*prefac*summe
end


# returns flati and fills temp_args_arr
function flattento1Dloop(temp_args_arr,groupindex_arr,bvalsdiag,abvals_arr,s_arr,JsS_arr,JlL_arr,lL_nested,nmax,Nmax,nu_arr,NU_arr,norm_arr,NORM_arr,mij_arr_dict,starts)
    # Keep loop-structure and write necessary functions arguments for matrix-element-calculation into 1-dim array temp_args_arr
    flati = 0
    # Iterate over boxes:
    for boxC in groupindex_arr
        for boxR in groupindex_arr
            boxR < boxC && continue # fill only lower-triangular (boxes, not elements!) only works for real-symmetric or hermitian matrices; NOT anymore for CSM? Also works for CSM -> complex symmetric matrices
            
            # if there are some identical particles: we can ignore the sum over a-values and simply multiply by a factor which is equal to the number of a-values. ONLY ON THE BOX-DIAGONAL! (boxC = boxR)
            if boxC == boxR
                bvals_new = bvalsdiag[boxC] # for changing the role of a,b to be in line with lower triangular!
                factor_ab = lastindex(abvals_arr[boxR]) # factor for amount of a-values normally
                diag_bool = 1 # for skipping lower-triangular calculation within each box!
            else
                #avals_new = abvals_arr[boxR]
                bvals_new = abvals_arr[boxC]
                factor_ab = 1
                diag_bool = 0
            end
            #bvals_new = abvals_arr[boxC] # for more consistent naming...
            avals_new = abvals_arr[boxR] # a bit unneccessary, but for more consistent naming. effectively only used here in fillTVS
            avals = abvals_arr[boxR] # necessary for observables
            bvals = abvals_arr[boxC] # necessary for observables
            
            
            alphab = 0; # this is actually beta
            for (isb,sb) in enumerate(s_arr[bvals_new[1]]) # why bvals_new[1]? -> cuz either only 1 value in it, or the particles are identical and all entries are the same!
                for (ijssb,JsSb) in enumerate(JsS_arr[bvals_new[1]][isb])
                    for (ijllb,JlLb) in enumerate(JlL_arr[bvals_new[1]][isb][ijssb])
                        for (lb,Lb) in lL_nested[bvals_new[1]][isb][ijssb][ijllb]
                            for nb = 1:nmax
                                nub = nu_arr[nb]
                                normb = norm_arr[lb,nb]
                                for Nb = 1:Nmax
                                    NUb = NU_arr[Nb]
                                    NORMb = NORM_arr[Lb,Nb]
                                    
                                    alphab += 1                        
                                    alpha = 0
                                    
                                    for (isa,sa) in enumerate(s_arr[avals_new[1]])
                                        for (ijssa,JsSa) in enumerate(JsS_arr[avals_new[1]][isa])
                                            for (ijlla,JlLa) in enumerate(JlL_arr[avals_new[1]][isa][ijssa])
                                                for (la,La) in lL_nested[avals_new[1]][isa][ijssa][ijlla]
                                                    Lsum=Int64((la+La+lb+Lb)/2)
                                                    mij_arr = mij_arr_dict[(la,La),(lb,Lb)]
                                                    for na = 1:nmax
                                                        nua = nu_arr[na]
                                                        norma = norm_arr[la,na]
                                                        for Na = 1:Nmax
                                                            NUa = NU_arr[Na]
                                                            NORMa = NORM_arr[La,Na]
                                                            
                                                            alpha += 1
                                                            diag_bool == 1 && alpha < alphab && continue # skip upper triangular only on diagonal boxes
                                                            
                                                            norm4 = norma*normb*NORMa*NORMb;
                                                            ranges = (;nua,nub,NUa,NUb);                                    
                                                            
                                                            # alpha, alphab are row and column within a given box defined by boxR,boxC. Now transforming into the row- and column indices of the big matrix:
                                                            rowi = starts[boxR] + alpha - 1
                                                            coli = starts[boxC] + alphab - 1
                                                            
                                                            flati += 1

                                                            temp_args_arr[flati] = TempArgs(rowi,coli,ranges,norm4,mij_arr,sa,JsSa,sb,JsSb,JlLa,JlLb,la,La,lb,Lb,Lsum,avals_new,bvals_new,factor_ab,avals,bvals) # write all loop-index-depending function-arguments into a temporary 1D array
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end                                    
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