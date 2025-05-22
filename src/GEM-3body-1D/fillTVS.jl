## Function for calculating the matrix elements and filling the matrices T,V,S within the ISGL program

@views @inbounds function fill_TVS(num_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,gaussopt,csm_bool)
    
    (;gem_params,mu0,c_shoulder,theta_csm) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;abvals_arr,cvals,groupindex_arr,abI,factor_bf,box_size_arr,starts,ends,bvalsdiag,lL_arr,maxlmax) = size_params
    (;gamma_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr) = precomp_arrs
    (;alpha_arr,v_arr,A_mat,w_interpol_arr,Ainv_arr_kine) = interpol_arrs
    (;w_arr_kine,wn_interpol_arr,kij_arr,gij_arr,T,V,S,temp_args_arr,temp_fill_mat) = fill_arrs
    
    flati = flattento1Dloop(temp_args_arr,groupindex_arr,factor_bf,bvalsdiag,abvals_arr,lL_arr,nmax,Nmax,nu_arr,NU_arr,norm_arr,NORM_arr,starts,cvals)
    
    #removed parallelization entirely, as mostly not useful. the parallelization is much better done outside, e.g. when calling the 3-body code in a loop
    
    ## Calculation of matrix elements and matrix filling via 1d loop:
    # the idea is that this should be relatively easy to parallelize. also it gets automatically adjusted to different loop structure from above
    for index in 1:flati
        rowi,coli=temp_args_arr[index]
        temp_fill_mat[rowi,coli] = sab(jmat,temp_args_arr[index],abI,factor_bf,gamma_dict)
    end
    # transpose fill:
    S .= Symmetric(temp_fill_mat,:L);
    
    #combined t and v for some speedup: removed combination to cope for CSM
    for index in 1:flati
        rowi,coli=temp_args_arr[index]
        temp_fill_mat[rowi,coli] = tab(jmat,murR_arr,w_arr_kine,Ainv_arr_kine,kij_arr,mu0,c_shoulder,temp_args_arr[index],abI,factor_bf,gamma_dict)
    end
    # transpose fill:
    T .= Symmetric(temp_fill_mat,:L);
    
    ## unclear type-stability
    #@show(gaussopt,typeof(gaussopt))
    if csm_bool == 1
        T .*= exp(-2*im*theta_csm*pi/180)
        #gaussoptc = gaussopt .*(1.0+0.0*im)
        #for cc in 1:lastindex(gaussopt)
        #    for vi in 1:lastindex(gaussopt[cc])
        #        gaussoptc[cc][vi][3] *= exp(2*im*theta_csm*pi/180)
        #    end
        #end
    else
        #gaussoptc = gaussopt
    end
    #gaussopt=gaussoptc
    #@show(gaussopt,typeof(gaussopt))
    
    #v
    for index in 1:flati
        rowi,coli=temp_args_arr[index]
        temp_fill_mat[rowi,coli] = vab(jmat,gij_arr,mu0,c_shoulder,w_interpol_arr,wn_interpol_arr,temp_args_arr[index],gaussopt,abI,factor_bf,gamma_dict)
    end
    # transpose fill:
    V .= Symmetric(temp_fill_mat,:L);
    

    # Example usage for debugging
    size_to_print = min(10, size(T, 1))  # Adjust size_to_print as needed
    println("T:")
    print_matrices(T, size_to_print)
    println("V:")
    print_matrices(V, size_to_print)
    #println("S:")
    #print_matrices(S, size_to_print)
    
    T .+= V
    
end

# Debugging: Print matrices in a formatted way
function print_matrices(M, size_to_print)
    for i in 1:size_to_print, j in 1:size_to_print
        @printf("%10.4f ", real(M[i, j]))
        j == size_to_print && println()
    end
end

# will it cost a lot of performance?
# this was implemented in the 3D code to fix some problem with symmetrization. but what was it?
function facsymm(a,b,abI,la,lb,factor_bf)
    #abI == 0 && return 1
    factor_symm = 1
    a == abI && (factor_symm *= factor_bf*(-1)^la)
    b == abI && (factor_symm *= factor_bf*(-1)^lb)
    return factor_symm
end



## functions to calculate one matrix element, summing over the necessary a,b,c values:
function sab(jmat,temp_args_i,abI,factor_bf,gamma_dict)
    (;avals_new,bvals_new,factor_ab,ranges,norm4,la,La,lb,Lb) = temp_args_i
    tempS = 0.0
    for a in avals_new
        for b in bvals_new
            factor_symm = facsymm(a,b,abI,la,lb,factor_bf)
            tempS += factor_ab*factor_symm*element_S(ranges,norm4,jmat[a,b],la,La,lb,Lb,gamma_dict)
        end
    end
    return tempS
end

function tab(jmat,murR_arr,w_arr_kine,Ainv_arr_kine,kij_arr,mu0,c_shoulder,temp_args_i,abI,factor_bf,gamma_dict)
    #rowi,coli,ranges,norm4,mij_arr,S_arr,la,La,lb,Lb,Lsum,avals_new,bvals_new,cvals,factor_ab,factor_symm = temp_args_i
    (;avals_new,bvals_new,factor_ab,ranges,norm4,la,La,lb,Lb,Lsum) = temp_args_i
    tempT = 0.0                                
    for a in avals_new
        for b in bvals_new
            factor_symm = facsymm(a,b,abI,la,lb,factor_bf)
            tempT += factor_ab*factor_symm*element_T(ranges,norm4,jmat[a,b],murR_arr,la,La,lb,Lb,b,gamma_dict)
        end
    end
    return tempT
end

function vab(jmat,gij_arr,mu0,c_shoulder,w_interpol_arr,wn_interpol_arr,temp_args_i,gaussopt,abI,factor_bf,gamma_dict)
    #rowi,coli,ranges,norm4,mij_arr,S_arr,la,La,lb,Lb,Lsum,avals_new,bvals_new,cvals,factor_ab,factor_symm = temp_args_i
    (;avals_new,bvals_new,factor_ab,ranges,norm4,la,La,lb,Lb,Lsum,cvals) = temp_args_i
    tempV = 0.0                                    
    for a in avals_new
        for b in bvals_new
            factor_symm = facsymm(a,b,abI,la,lb,factor_bf)
            for c in cvals
                # maybe not the fastest way to have the if-condition within this loop, but thats what it is currently.
                if gaussopt[c][1][1] == 1
                    tempV += factor_ab*factor_symm*element_VGauss(c,ranges,norm4,jmat[a,c],jmat[b,c],la,La,lb,Lb,gaussopt,gamma_dict)
                #elseif coulopt[c][1][1] == 1
                #    tempV += factor_ab*factor_symm*element_VCoulomb(c,ranges,norm4,jmat[a,c],jmat[b,c],la,La,lb,Lb,coulopt,gamma_dict)
                else
                    tempV += factor_ab*factor_symm*element_V(c,ranges,norm4,jmat[a,c],jmat[b,c],la,La,lb,Lb,w_interpol_arr,wn_interpol_arr,gamma_dict)
                    #error("Error in fillTVS: gaussopt[c][1][1] = $(gaussopt[c][1][1]) has invalid value. Only 0 or 1 allowed.")
                end                        
            end
        end
    end
    #@show([la,La,lb,Lb],tempV)
    return tempV
end

function element_S(ranges,norm4,jab,la,La,lb,Lb,gamma_dict)
    
    mod(la+La+lb+Lb,2) == 1 && return 0.0

    (;nua,nub,NUa,NUb) = ranges;
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
                #z = gamma_dict[(la+La+lb+Lb-2*s+1)/2]/gamma_dict[s+1]/gamma_dict[k+K+Lb-2*s+1] * eta3^s*eta2^(k+K+Lb-2*s)/((eta1-eta2^2/4/eta3)^((la+La+lb+Lb-2*s+1)/2))
                #@show([k,K,Lb,s,eta2,eta2^(k+K+Lb-2*s)])
            end
            #@show([k,K,Lb,sums])
            sumK += sums*gamma_dict[la+1]/(gamma_dict[la-k+1]*gamma_dict[k+1]) * gamma_dict[La+1]/(gamma_dict[La-K+1]*gamma_dict[K+1]) * alpha^(la-k) * beta^k * gamma^(La-K) * delta^K * gamma_dict[k+K+Lb+1]*(-1/2/eta3)^(k+K+Lb)
        end
        sumk += sumK
    end
    sumk *= norm4 * (pi/eta3)^(1/2)
    
    return sumk
end


# calculation of a single matrix element: kinetic energy T
function element_T(ranges,norm4,jab,murR_arr,la,La,lb,Lb,b,gamma_dict)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alpha,gamma,beta,delta) = jab # careful with order (gamma before beta)
    
    mur= murR_arr[1,b]
    muR= murR_arr[2,b]
    
    # -hbar^2/2/mu # hbar=1 here. trb,tRb:
    tr = -1/2/mur
    tR = -1/2/muR

    temp = 4*nub^2*element_S(ranges,norm4,jab,la,La,lb+2,Lb,gamma_dict) - 2*(1+2*lb)*nub*element_S(ranges,norm4,jab,la,La,lb,Lb,gamma_dict)
    lb >= 2 && ( temp += lb*(lb-1)*element_S(ranges,norm4,jab,la,La,lb-2,Lb,gamma_dict) )
    temp *= tr

    Temp = 4*NUb^2*element_S(ranges,norm4,jab,la,La,lb,Lb+2,gamma_dict) - 2*(1+2*Lb)*NUb*element_S(ranges,norm4,jab,la,La,lb,Lb,gamma_dict)
    Lb >= 2 && ( Temp += Lb*(Lb-1)*element_S(ranges,norm4,jab,la,La,lb,Lb-2,gamma_dict) )
    Temp *= tR

    return temp+Temp
end


# calculation of a single matrix element: interaction V(r_c)
function element_V(c,ranges,norm4,jac,jbc,la,La,lb,Lb,w_interpol_arr,wn_interpol_arr,gamma_dict)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    #EH: alpha <-> gamma; beta <-> delta    
    
    eta5 = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2
    eta6 = 2*(nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC)
    eta7 = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2
    
    etaeff = eta5 - eta6^2/4/eta7
    Lsum = la+La+lb+Lb

    for n = 0:Lsum
        wn_interpol_arr[n] = w_interpol_arr[c,n](log(etaeff))
    end
    #@show([la,La,lb,Lb],wn_interpol_arr)

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
                        #zz = eta7^s/(eta6^(2*s)*gamma_dict[s+1]*gamma_dict[KK-2*s+1])*wn_interpol_arr[nn]
                        #z1 = wn_interpol_arr[nn]
                        #z2 = mod(LL-2*s+1,2)*gamma_dict[(LL-2*s+1)/2]/((eta5+1.0-eta6^2/4/eta7)^((LL-2*s+1)/2))
#                        @show([eta7^s,eta6^(2*s),gamma_dict[s+1],gamma_dict[KK-2*s+1],wn_interpol_arr[nn]])
                        #@show([s,eta6,1/eta6^(2s)])
                    end

                    # this seems a bit inefficient to compute due to many dict calls?
                    sumKp += sums * gamma_dict[la+1]/(gamma_dict[la-k+1]*gamma_dict[k+1]) * gamma_dict[La+1]/(gamma_dict[La-K+1]*gamma_dict[K+1]) * gamma_dict[lb+1]/(gamma_dict[lb-kp+1]*gamma_dict[kp+1]) * gamma_dict[Lb+1]/(gamma_dict[Lb-Kp+1]*gamma_dict[Kp+1]) * alphaAC^(la-k)*betaAC^k*gammaAC^(La-K)*deltaAC^K * alphaBC^(lb-kp)*betaBC^kp*gammaBC^(Lb-Kp)*deltaBC^Kp * gamma_dict[KK+1] * (pi/eta7)^(1/2)*(-1/2/eta7)^KK

                    #@show([la,La,lb,Lb,k,K,kp,Kp],sums,sumKp)
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
function element_VGauss(c,ranges,norm4,jac,jbc,la,La,lb,Lb,gaussopt,gamma_dict)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    #EH: alpha <-> gamma; beta <-> delta

    v0 = gaussopt[c][1][2]
    mu_g = gaussopt[c][1][3]

    LL = la+La+lb+Lb;

    sumk = 0.0
    mod(LL+1,2) == 0 && return sumk

    etac  = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2
    zetac = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2 
    xic = nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC
    etaprc = etac - xic^2/zetac

    eta5 = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2 + mu_g
    eta6 = 2*(nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC)
    eta7 = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2
    
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
                        #zz    = eta7^s*eta6^(KK-2s)/(gamma_dict[s+1]*gamma_dict[KK-2*s+1]) * gamma_dict[(LL-2*s+1)/2]/((eta5-eta6^2/4/eta7)^((LL-2*s+1)/2));
                    end
                    #@show(sums)

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

# calculation of a single matrix element: interaction VCoulomb(r_c)
# BEWARE! pure Coulomb is singular in 1D! -> the "s"-wave component makes problems, and it cannot be excluded, as it is always implicitly included via the Jacobi-transformations!
function element_VCoulomb(c,ranges,norm4,jac,jbc,la,La,lb,Lb,coulopt,gamma_dict)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    #EH: alpha <-> gamma; beta <-> delta

    ZaZb = coulopt[c][1][2] # product of the two relevant charges
    LL = la+La+lb+Lb-1;

    sumk = 0.0
    mod(LL+1+1,2) == 0 && return sumk

    eta5 = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2
    eta6 = 2*(nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC)
    eta7 = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2
    

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
                        @show([LL,KK,s])
                        sums += eta7^s*eta6^(KK-2*s)*gamma_dict[(LL-2*s+1)/2]/(gamma_dict[s+1]*gamma_dict[KK-2*s+1]*(eta5-eta6^2/4/eta7)^((LL-2*s+1)/2))
                    end

                    sumKp += sums * gamma_dict[la+1]/(gamma_dict[la-k+1]*gamma_dict[k+1]) * gamma_dict[La+1]/(gamma_dict[La-K+1]*gamma_dict[K+1]) * gamma_dict[lb+1]/(gamma_dict[lb-kp+1]*gamma_dict[kp+1]) * gamma_dict[Lb+1]/(gamma_dict[Lb-Kp+1]*gamma_dict[Kp+1]) * alphaAC^(la-k)*betaAC^k*gammaAC^(La-K)*deltaAC^K * alphaBC^(lb-kp)*betaBC^kp*gammaBC^(Lb-Kp)*deltaBC^Kp * gamma_dict[KK+1] * (pi/eta7)^(1/2)*(-1/2/eta7)^KK
                end
                sumkp += sumKp
            end
            sumK += sumkp
        end
        sumk += sumK
    end
    sumk *= norm4 * ZaZb
    
    return sumk    
end


# returns flati and fills temp_args_arr
function flattento1Dloop(temp_args_arr,groupindex_arr,factor_bf,bvalsdiag,abvals_arr,lL_arr,nmax,Nmax,nu_arr,NU_arr,norm_arr,NORM_arr,starts,cvals)
    # Keep loop-structure and write necessary functions arguments for matrix-element-calculation into 1-dim array temp_args_arr
    flati = 0
    # Iterate over boxes:
    for boxC in groupindex_arr
        for boxR in groupindex_arr
            boxR < boxC && continue # fill only lower-triangular (boxes, not elements!) only works for real-symmetric or hermitian matrices; NOT anymore for CSM? Also works for CSM -> complex symmetric matrices
            
            # if there are some identical particles: we can ignore the sum over a-values and simply multiply by a factor which is equal to the number of a-values. ONLY ON THE BOX-DIAGONAL! (boxC = boxR)
            if boxC == boxR
                #avals_new = avalsdiag[boxR]#[abvals_arr[boxR][1]] # take only the first of possible values
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
                            Lsum=0 #Int64((la+La+lb+Lb)/2) #completely useless in 1D.
                            #mij_arr = mij_arr_dict[(la,La),(lb,Lb)]

                            for na = 1:nmax
                                nua = nu_arr[na]
                                norma = norm_arr[la,na]
                                for Na = 1:Nmax
                                    NUa = NU_arr[Na]
                                    NORMa = NORM_arr[La,Na]
                                    
                                    alpha += 1
                                    diag_bool == 1 && alpha < alphab && continue # skip lower-triangular only on diagonal boxes! -> actually skips upper triangular!
                                    
                                    norm4 = norma*normb*NORMa*NORMb;
                                    ranges = (;nua,nub,NUa,NUb);
                                    
                                    # alpha, alphab are row and column within a given box defined by boxR,boxC. Now transforming into the row- and column indices of the big matrix:
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