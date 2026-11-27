## Function for calculating the matrix elements and filling the matrices T,V,S within the ISGL program

@views @inbounds function fill_TVS(num_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,complex_scaling::Bool,hbar,complex_ranged_r::Bool=false,complex_ranged_R::Bool=false,add_V::Bool=true,parallel::Bool=false)

    (;gem_params,mu0,c_shoulder,complex_scaling_angle) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;abvals_arr,cvals,gauss_indices,central_indices,so_indices,pow_indices,groupindex_arr,abI,factor_bf,box_size_arr,nbasis_total,starts,ends,bvalsdiag,s_arr,JsS_arr,s_complete,JsS_complete,JlL_arr,lL_nested,maxlmax,mij_arr_dict,mijSO_arr_dict,gaussopt_arr,powopt_arr) = size_params
    (;gamma_dict,spintrafo_dict,spinoverlap_dict,global6j_dict,facsymm_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr,Clmk_arr,Dlmk_arr,S_arr,SSO_arr) = precomp_arrs
    (;alpha_arr,v_arr,A_mat,w_interpol_arr,Ainv_arr_kine,w_pow_arr) = interpol_arrs
    (;T,V,S,basis) = fill_arrs
    
    complex_ranged = complex_ranged_r || complex_ranged_R
    # With complex ranges the bra ranges are conjugated, so S, T and V become hermitian instead of
    # symmetric. Combined with complex scaling neither symmetry survives and the full matrix is filled.
    fill_full = complex_ranged && complex_scaling
    add_V || size(V) == size(T) || error("fill_TVS: add_V=false needs a separate V, see preallocate_data")

    ## Calculation of matrix elements and matrix filling: one loop over the elements of the lower triangle (full matrix
    # for fill_full), mirrored in place afterwards. The arguments of each element are generated from the basis list.
    # V is accumulated into T, which becomes the Hamiltonian H = T + V. Kept separately only for add_V=false (inverse problem, matrix export)
    csmfac = exp(-2*im*complex_scaling_angle*pi/180)
    function element(rowi,coli)
        args = element_args(rowi,coli,basis,abvals_arr,bvalsdiag,nu_arr,NU_arr,norm_arr,NORM_arr,mij_arr_dict)
        sel = sab(jmat,args,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict)
        tel = convert(eltype(T), tab(jmat,murR_arr,Ainv_arr_kine,mu0,c_shoulder,args,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict,hbar))
        complex_scaling && (tel *= csmfac)
        vel = vab(jmat,mu0,c_shoulder,w_interpol_arr,w_pow_arr,args,abI,factor_bf,S_arr,SSO_arr,cvals,spintrafo_dict,spinoverlap_dict,facsymm_dict,gauss_indices,central_indices,so_indices,pow_indices,s_arr,global6j_dict,mijSO_arr_dict,gaussopt_arr,powopt_arr,complex_scaling) # do we need hbar^2 for SO?
        return sel, tel, vel
    end
    # complex ranges without complex scaling: conjugate-paired basis functions (see fill_partners)
    partner = complex_scaling ? nothing : fill_partners(element,nbasis_total,gem_params,complex_ranged_r,complex_ranged_R)
    function fill_columns(cols)
        for coli in cols, rowi in (fill_full ? 1 : coli):nbasis_total
            if !isnothing(partner)
                pr, pc, cj = pair_position(rowi,coli,partner)
                (pc, pr) < (coli, rowi) && continue # already filled together with its partner element
            end
            sel, tel, vel = element(rowi,coli)
            store_element!(S,T,V,rowi,coli,sel,tel,vel,add_V)
            !isnothing(partner) && (pr, pc) != (rowi, coli) && store_element!(S,T,V,pr,pc,sel,tel,vel,add_V,cj)
        end
    end
    fill_parallel(fill_columns,nbasis_total,parallel)
    # transpose fill: the overlap is hermitian for complex ranges (also together with complex scaling,
    # since the ISGL complex-scaling rotates the potential and leaves the ranges alone)
    hermitian_fill!(S,complex_ranged)
    if !fill_full
        hermitian_fill!(T,complex_ranged)
        add_V || hermitian_fill!(V,complex_ranged)
    end
    
end


# pairs (i,j) of the shoulder sums, in the order of the exponents (m12,m13,m14,m23,m24,m34) in mij_arr
const IJ_PAIRS = ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4))

## functions to calculate one matrix element, summing over the necessary a,b,c values:
function sab(jmat,temp_args_i,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict)
    (;rowi,coli,avals_new,bvals_new,factor_ab,ranges,norm4,mij_arr,sa,JsSa,sb,JsSb,la,La,lb,Lb,JlLa,JlLb) = temp_args_i
    tempS = zero(norm4) # complex for complex-ranged basis functions
    
    # maybe easier:
    (JsSa != JsSb || JlLa != JlLb) && return tempS #immediately skip if either one is violated
    
    for a in avals_new # this does not indicate the "boxes"!! just if in one box several a and b values are summed over (mostly due to identical particles)!
        for b in bvals_new
            factor_symm = facsymm_dict[a,b,la,lb,sa,sb]
            
            uab = spintrafo_dict[a,b,JsSa,sa,sb]
            JlL = JlLa
            
            tempS += factor_ab*factor_symm*uab*element_S(ranges,norm4,jmat[a,b],mij_arr,S_arr,la,La,lb,Lb,JlL)
        end
    end
    return tempS
end

function tab(jmat,murR_arr,Ainv_arr_kine,mu0,c_shoulder,temp_args_i,abI,factor_bf,S_arr,spintrafo_dict,facsymm_dict,hbar)
    (;avals_new,bvals_new,factor_ab,ranges,norm4,mij_arr,sa,JsSa,sb,JsSb,la,La,lb,Lb,Lsum,JlLa,JlLb) = temp_args_i
    tempT = zero(norm4)
    
    # maybe easier:
    (JsSa != JsSb || JlLa != JlLb) && return tempT #immediately skip if either one is violated
    
    for a in avals_new
        for b in bvals_new
            factor_symm = facsymm_dict[a,b,la,lb,sa,sb]
            
            uab = spintrafo_dict[a,b,JsSa,sa,sb]
            JlL = JlLa
            
            tempT += hbar^2*factor_ab*factor_symm*uab*element_T(ranges,norm4,jmat[a,b],murR_arr,mij_arr,S_arr,la,La,lb,Lb,Ainv_arr_kine,b,mu0,c_shoulder,Lsum,JlL)
        end
    end
    return tempT
end

function vab(jmat,mu0,c_shoulder,w_interpol_arr,w_pow_arr,temp_args_i,abI,factor_bf,S_arr,SSO_arr,cvals,spintrafo_dict,spinoverlap_dict,facsymm_dict,gauss_indices,central_indices,so_indices,pow_indices,s_arr,global6j_dict,mijSO_arr_dict,gaussopt_arr,powopt_arr,complex_scaling)
    (;rowi,coli,avals_new,bvals_new,factor_ab,ranges,norm4,mij_arr,sa,JsSa,sb,JsSb,la,La,lb,Lb,Lsum,JlLa,JlLb) = temp_args_i
    # only the mesh over arg(etaprc) >= 0 is tabulated unless complex scaling made the potential
    # complex and broke the Schwarz reflection (see theta_mesh)
    reflect = !complex_scaling
    # complex for complex ranges (via norm4) and/or complex scaling (via the interpolated and power-law weights)
    tempV = zero(promote_type(typeof(norm4),eltype(w_pow_arr)))
    
    for a in avals_new
        for b in bvals_new
            factor_symm = facsymm_dict[a,b,la,lb,sa,sb]
            
            # for spin-independent interactions: overlap of spin functions can be calculated in advance:
            uab = spintrafo_dict[a,b,JsSa,sa,sb]
            
            for c in cvals
                
                spinoverlap = spinoverlap_dict[a,b,c,JsSa,JsSb,sa,sb]
                
                for ivg in gauss_indices[c] #loop over the gaussian interactions for this c.
                    (JsSa != JsSb || JlLa != JlLb) && return tempV #immediately skip if it is violated.
                    v0,mu_g = gaussopt_arr[c][ivg]

                    tempV += factor_ab*factor_symm*uab*element_VGauss(c,ranges,norm4,jmat[a,c],jmat[b,c],mij_arr,S_arr,la,La,lb,Lb,mu0,c_shoulder,Lsum,v0,mu_g,JlLa)
                end
                for ivp in pow_indices[c] #loop over the power-law interactions for this c.
                    (JsSa != JsSb || JlLa != JlLb) && return tempV #immediately skip if it is violated.
                    v0,p_pow = powopt_arr[c][ivp]

                    tempV += factor_ab*factor_symm*uab*element_VPow(c,ranges,norm4,jmat[a,c],jmat[b,c],mij_arr,S_arr,la,La,lb,Lb,mu0,c_shoulder,Lsum,w_pow_arr,v0,p_pow,JlLa,ivp)
                end
                for ivc in central_indices[c] #loop over the central interactions for this c.
                    (JsSa != JsSb || JlLa != JlLb) && return tempV #immediately skip if it is violated.
                    tempV += factor_ab*factor_symm*uab*element_V(c,ranges,norm4,jmat[a,c],jmat[b,c],mij_arr,S_arr,la,La,lb,Lb,mu0,c_shoulder,w_interpol_arr,Lsum,JlLa,ivc,reflect)
                end
                #= commented out; spin-orbit is moved to future development
                for ivso in so_indices[c] # loop over the spin-orbit interactions for this c.
                    (abs(JlLa-JlLb) <= 1 <= JlLa+JlLb) == false && return tempV #immediately skip if it is violated.
                    (abs(JsSa-JsSb) <= 1 <= JsSa+JsSb) == false && return tempV #immediately skip if it is violated.
                    #only relevant for SO interactions
                    global6jfac = global6j_dict[JsSa,JsSb,JlLa,JlLb]
                    tempV += factor_ab*factor_symm*global6jfac*spinoverlap*element_VSO(c,ranges,norm4,jmat[a,c],jmat[b,c],mijSO_arr_dict,SSO_arr,la,La,lb,Lb,gij_arr,mu0,c_shoulder,w_interpol_arr,Lsum,wn_interpol_arr,JlLa,JlLb,ivso,reflect)
                end =#
                
            end
            
        end
    end
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
    p = SA[nua*beta,NUa*delta,0.0,NUb]
    ppr = SA[nua*alpha,NUa*gamma,nub,0.0]
    f = xi/zeta
    q = ppr .- f*p
    
    fij_arr = @SMatrix[2*p[i]*p[j]/zeta + 2*q[i]*q[j]/etapr for i=1:4, j=1:4]
    # the two 3/2-powers are kept separate on purpose: for complex ranges each factor is then
    # continued from its own (small) argument, instead of from the argument of the merged product.
    prefac = norm4 * (pi/zeta)^(3/2)*(pi/etapr)^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);

    sum = zero(prefac)
    for ii = 1:lastindex(mij_arr_i)
        (m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]      
        sum += S_arr[la,La,lb,Lb,JlL,ii]*fij_arr[1,2]^m12*fij_arr[1,3]^m13*fij_arr[1,4]^m14*fij_arr[2,3]^m23*fij_arr[2,4]^m24*fij_arr[3,4]^m34 * prefac
    end   
    
    return sum
end


# calculation of a single matrix element: kinetic energy T
function element_T(ranges,norm4,jab,murR_arr,mij_arr_i,S_arr,la,La,lb,Lb,Ainv_arr_kine,b,mu0,c_shoulder,Lsum,JlL)
    
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
    p = SA[nua*beta,NUa*delta,0.0,NUb]
    ppr = SA[nua*alpha,NUa*gamma,nub,0.0]
    f = xi/zeta
    fpr = xi/eta
    q = ppr .- f*p
    qpr = p .- fpr*ppr    # fixed typo!
    
    K0 = 6*tr*nub*(nub/etapr - 1.0) + 6*tR*NUb*(NUb/zetapr - 1.0)    
    
    Kij_arr = @SMatrix[8*tr*nub^2/etapr*(q[i]*q[j]/etapr - q[1]*I[i,1]*I[j,3] - q[2]*I[i,2]*I[j,3] - q[4]*I[i,3]*I[j,4]) + 8*tR*NUb^2/zetapr*(qpr[i]*qpr[j]/zetapr - qpr[1]*I[i,1]*I[j,4] - qpr[2]*I[i,2]*I[j,4] - qpr[3]*I[i,3]*I[j,4]) for i=1:4, j=1:4]
    
    prefac = norm4 * (pi/zeta)^(3/2)*(pi/etapr)^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);

    # per n: shoulder weight wn and the pair values kn (see IJ_PAIRS)
    sum = zero(promote_type(typeof(prefac),typeof(K0)))
    for n = 0:Lsum
        mun = mu0*c_shoulder^n
        wn = K0*Ainv_arr_kine[Lsum,n]
        kn = map(((i,j),) -> 2*p[i]*p[j]/zeta + 2*q[i]*q[j]/etapr + mun*Kij_arr[i,j]/K0, IJ_PAIRS)
        for ii = 1:lastindex(mij_arr_i)
            (m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]
            sum += S_arr[la,La,lb,Lb,JlL,ii]*wn*kn[1]^m12*kn[2]^m13*kn[3]^m14*kn[4]^m23*kn[5]^m24*kn[6]^m34
        end
    end
    
    return prefac*sum
end

# calculation of a single matrix element: interaction V(r_c)
function element_V(c,ranges,norm4,jac,jbc,mij_arr_i,S_arr,la,La,lb,Lb,mu0,c_shoulder,w_interpol_arr,Lsum,JlL,iv,reflect::Bool)
    
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
    
    prefac = norm4 * (pi/zetac)^(3/2)*(pi/etaprc)^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);

    # per n: interpolated weight wn (once per n, outside of the ii-loop) and the pair values gn (see IJ_PAIRS).
    # For real ranges this is the usual lookup at log(etaprc). For complex ranges etaprc is complex
    # but confined to the sector |arg| <= atan(omega), over which the interpolant is tabulated in
    # two dimensions; interpol_lookup picks the right variant (see common/auxiliary.jl).
    summe = zero(promote_type(typeof(prefac),eltype(eltype(w_interpol_arr))))
    for n = 0:Lsum
        mun = mu0*c_shoulder^n
        wn = interpol_lookup(w_interpol_arr[c,iv,Lsum,n],etaprc,reflect)
        gn = map(((i,j),) -> 2*p[i]*p[j]/zetac + 2*mun*q[i]*q[j]/etaprc, IJ_PAIRS)
        for ii = 1:lastindex(mij_arr_i)
            (m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]
            summe += S_arr[la,La,lb,Lb,JlL,ii]*wn*gn[1]^m12*gn[2]^m13*gn[3]^m14*gn[4]^m23*gn[5]^m24*gn[6]^m34
        end
    end
    
    return prefac*summe
end


# calculation of a single matrix element: power-law interaction V(r) = v0*r^p (analytic)
# Same structure as element_V, but the interpolated weights w_interpol_arr[c,iv,Lsum,n](log(etaprc))
# are replaced by the exact w_n(etaprc) = etaprc^(-p/2) * w_pow_arr[c,iv,Lsum,n]. The common factor
# etaprc^(-p/2) is independent of n and is therefore absorbed into prefac.
function element_VPow(c,ranges,norm4,jac,jbc,mij_arr_i,S_arr,la,La,lb,Lb,mu0,c_shoulder,Lsum,w_pow_arr,v0,p_pow,JlL,iv)
    
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
    
    # etaprc^(-p/2) uses the principal branch. For complex ranges arg(etaprc) is bounded by
    # atan(omega) < pi/2, so the principal branch is the continuous continuation of the real-range
    # result; x^y is evaluated as exp(y*log(x)) with the principal log of the *base*, so no wrapping
    # occurs even for large |p|.
    prefac = norm4 * (pi/zetac)^(3/2)*(pi/etaprc)^(3/2) * etaprc^(-p_pow/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);

    # no interpolation necessary here: w_pow_arr is exact and alpha-independent. Per n: pair values gn (see IJ_PAIRS)
    summe = zero(promote_type(typeof(prefac),eltype(w_pow_arr)))
    for n = 0:Lsum
        mun = mu0*c_shoulder^n
        wn = w_pow_arr[c,iv,Lsum,n]
        gn = map(((i,j),) -> 2*p[i]*p[j]/zetac + 2*mun*q[i]*q[j]/etaprc, IJ_PAIRS)
        for ii = 1:lastindex(mij_arr_i)
            (m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]
            summe += S_arr[la,La,lb,Lb,JlL,ii]*wn*gn[1]^m12*gn[2]^m13*gn[3]^m14*gn[4]^m23*gn[5]^m24*gn[6]^m34
        end
    end
    
    return v0*prefac*summe
end

# matrix element for Spin-Orbit interaction; postponed to future version
#=
function element_VSO(c,ranges,norm4,jac,jbc,mijSO_arr_dict,SSO_arr,la,La,lb,Lb,gij_arr,mu0,c_shoulder,w_interpol_arr,Lsum,wn_interpol_arr,JlLa,JlLb,ivso,reflect::Bool)
    
    #prechecks: they should ideally never trigger, due to proper handling before.
    Lsum < 1 && return 0.0
    #(abs(JlLa-JlLb) <= 1 <= JlLa+JlLb) == false && return 0.0 # done outside already?
    LsumSO = Lsum - 1
    
    (;nua,nub,NUa,NUb) = ranges; # is this usage of named tuple slow?
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
        wn_interpol_arr[n] = interpol_lookup(w_interpol_arr[c,ivso,LsumSO,n],etaprc,reflect)
    end
    
    ivjv_arr = @SVector[(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)] # can also be used for iuju
    
    summe = 0.0 # total sum, could also be named sumv, since sum over v (or iv,jv)
    for (v,(iv,jv)) in enumerate(ivjv_arr)
        mij_arr_i = mijSO_arr_dict[((la,La),(lb,Lb),v)] # this is a different approach. for Vcent, the mij_arr_i is already selected in element_args.
        
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
end =#





# calculation of a single matrix element: interaction VGauss(r_c)
function element_VGauss(c,ranges,norm4,jac,jbc,mij_arr_i,S_arr,la,La,lb,Lb,mu0,c_shoulder,Lsum,v0,mu_g,JlL)
    
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
    
    #gij_arr:
    ggij_arr = @SMatrix[2*p[i]*p[j]/zetac + 2*q[i]*q[j]/(etaprc + mu_g) for i=1:4, j=1:4]
    
    prefac = norm4 * (pi/zetac)^(3/2)*(pi/(etaprc + mu_g))^(3/2)/(nua^la*NUa^La*nub^lb*NUb^Lb);
    
    summe = zero(eltype(ggij_arr))
    for ii = 1:lastindex(mij_arr_i)
        (m12,m13,m14,m23,m24,m34) = mij_arr_i[ii]

        summe += S_arr[la,La,lb,Lb,JlL,ii]*ggij_arr[1,2]^m12*ggij_arr[1,3]^m13*ggij_arr[1,4]^m14*ggij_arr[2,3]^m23*ggij_arr[2,4]^m24*ggij_arr[3,4]^m34
    end
    
    return v0*prefac*summe
end


# Basis functions in matrix order, one entry per row/column: box, spin and angular-momentum channel, and the
# range indices n,N (effective ones, i.e. doubled for a complex-ranged coordinate). Same loop order as the matrix.
function basis_list(size_params,nmax,Nmax)
    (;groupindex_arr,abvals_arr,s_arr,JsS_arr,JlL_arr,lL_nested,starts,nbasis_total) = size_params
    basis = Vector{BasisFun}(undef,nbasis_total)
    for box in groupindex_arr
        a = abvals_arr[box][1] # either only 1 value in it, or the particles are identical and all entries are the same
        alpha = 0
        for (is,s) in enumerate(s_arr[a])
            for (ijss,JsS) in enumerate(JsS_arr[a][is])
                for (ijll,JlL) in enumerate(JlL_arr[a][is][ijss])
                    for (l,L) in lL_nested[a][is][ijss][ijll]
                        for n = 1:nmax
                            for N = 1:Nmax
                                alpha += 1
                                basis[starts[box] + alpha - 1] = (;box,s,JsS,JlL,l,L,n,N)
                            end
                        end
                    end
                end
            end
        end
    end
    return basis
end

# arguments of the matrix element between the basis functions rowi (bra) and coli (ket), see TempArgs
function element_args(rowi,coli,basis,abvals_arr,bvalsdiag,nu_arr,NU_arr,norm_arr,NORM_arr,mij_arr_dict)
    a, b = basis[rowi], basis[coli]
    # if there are some identical particles: we can ignore the sum over a-values and simply multiply by a factor which is equal to the number of a-values. ONLY ON THE BOX-DIAGONAL!
    if a.box == b.box
        bvals_new = bvalsdiag[b.box]
        factor_ab = lastindex(abvals_arr[a.box])
    else
        bvals_new = abvals_arr[b.box]
        factor_ab = 1
    end
    # the bra ranges enter conjugated. This is a no-op for real ranges and is what makes S,T,V hermitian for complex ones.
    nua, norma = conj(nu_arr[a.n]), conj(norm_arr[a.l,a.n])
    NUa, NORMa = conj(NU_arr[a.N]), conj(NORM_arr[a.L,a.N])
    nub, normb = nu_arr[b.n], norm_arr[b.l,b.n]
    NUb, NORMb = NU_arr[b.N], NORM_arr[b.L,b.N]
    norm4 = norma*normb*NORMa*NORMb
    Lsum = Int64((a.l+a.L+b.l+b.L)/2)
    return TempArgs(rowi,coli,(;nua,nub,NUa,NUb),norm4,mij_arr_dict[a.l,a.L,b.l,b.L],a.s,a.JsS,b.s,b.JsS,a.JlL,b.JlL,a.l,a.L,b.l,b.L,Lsum,abvals_arr[a.box],bvals_new,factor_ab,abvals_arr[a.box],abvals_arr[b.box])
end