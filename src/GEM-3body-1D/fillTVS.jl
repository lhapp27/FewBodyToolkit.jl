## Function for calculating the matrix elements and filling the matrices T,V,S within the GEM3B1D program

@views @inbounds function fill_TVS(num_params,size_params,precomp_arrs,interpol_arrs,fill_arrs,complex_scaling::Bool,hbar,complex_ranged_r::Bool=false,complex_ranged_R::Bool=false,add_V::Bool=true)
    
    (;gem_params,complex_scaling_angle) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;abvals_arr,cvals,groupindex_arr,abI,factor_bf,box_size_arr,nbasis_total,starts,ends,bvalsdiag,lL_arr,maxlmax,gauss_indices,central_indices,contact1D_indices,gaussopt_arr,contact1Dopt_arr) = size_params
    (;gamma_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr) = precomp_arrs
    (;alpha_arr,v_arr,w_interpol_arr) = interpol_arrs
    (;wn_interpol_arr,T,V,S,basis) = fill_arrs

    complex_ranged = complex_ranged_r || complex_ranged_R
    # With complex ranges the bra ranges are conjugated, so S, T and V become hermitian instead of
    # symmetric. Combined with complex scaling neither symmetry survives and the full matrix is filled.
    fill_full = complex_ranged && complex_scaling
    add_V || size(V) == size(T) || error("fill_TVS: add_V=false needs a separate V, see preallocate_data")

    ## Calculation of matrix elements and filling: one loop over the elements of the lower triangle (full matrix
    # for fill_full), mirrored in place afterwards. The arguments of each element are generated from the basis list.
    # V is accumulated into T, which becomes the Hamiltonian H = T + V. Kept separately only for add_V=false (inverse problem, matrix export)
    csmfac = exp(-2*im*complex_scaling_angle*pi/180)
    function element(rowi,coli)
        args = element_args(rowi,coli,basis,abvals_arr,bvalsdiag,cvals,nu_arr,NU_arr,norm_arr,NORM_arr)
        sel = sab(jmat,args,abI,factor_bf,gamma_dict)
        tel = convert(eltype(T), tab(jmat,murR_arr,args,abI,factor_bf,gamma_dict,hbar))
        complex_scaling && (tel *= csmfac)
        vel = vab(jmat,w_interpol_arr,wn_interpol_arr,args,abI,factor_bf,gamma_dict,gauss_indices,central_indices,contact1D_indices,gaussopt_arr,contact1Dopt_arr,complex_scaling,complex_scaling_angle)
        return sel, tel, vel
    end
    # complex ranges without complex scaling: conjugate-paired basis functions (see fill_partners)
    partner = complex_scaling ? nothing : fill_partners(element,nbasis_total,gem_params,complex_ranged_r,complex_ranged_R)
    for coli in 1:nbasis_total, rowi in (fill_full ? 1 : coli):nbasis_total
        if !isnothing(partner)
            pr, pc, cj = pair_position(rowi,coli,partner)
            (pc, pr) < (coli, rowi) && continue # already filled together with its partner element
        end
        sel, tel, vel = element(rowi,coli)
        store_element!(S,T,V,rowi,coli,sel,tel,vel,add_V)
        !isnothing(partner) && (pr, pc) != (rowi, coli) && store_element!(S,T,V,pr,pc,sel,tel,vel,add_V,cj)
    end
    # transpose fill: the overlap is hermitian for complex ranges (also together with complex scaling,
    # since the complex-scaling rotates the potential and leaves the ranges alone)
    hermitian_fill!(S,complex_ranged)
    # transpose fill:
    if !fill_full
        hermitian_fill!(T,complex_ranged)
        add_V || hermitian_fill!(V,complex_ranged)
    end

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
    tempS = zero(norm4) # complex for complex-ranged basis functions
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
    tempT = zero(norm4)                                
    for a in avals_new
        for b in bvals_new
            factor_symm = facsymm(a,b,abI,la,lb,factor_bf)
            tempT += factor_ab*factor_symm*element_T(ranges,norm4,jmat[a,b],murR_arr,la,La,lb,Lb,b,gamma_dict,hbar)
        end
    end
    return tempT
end

# calculation of a single matrix element: interaction V(r_c)
function vab(jmat,w_interpol_arr,wn_interpol_arr,temp_args_i,abI,factor_bf,gamma_dict,gauss_indices,central_indices,contact1D_indices,gaussopt_arr,contact1Dopt_arr,complex_scaling,complex_scaling_angle)
    (;avals_new,bvals_new,factor_ab,ranges,norm4,la,La,lb,Lb,cvals) = temp_args_i
    tempV = zero(promote_type(typeof(norm4),eltype(wn_interpol_arr)))
    # only the mesh over arg(etaeff) >= 0 is tabulated unless complex scaling made the potential
    # complex and broke the Schwarz reflection (see theta_mesh)
    reflect = !complex_scaling
    for a in avals_new
        for b in bvals_new
            factor_symm = facsymm(a,b,abI,la,lb,factor_bf)
            for c in cvals
                
                for ivg in gauss_indices[c] # we need to provide these indices to the function.
                    gaussopt = gaussopt_arr[c][ivg] # Gaussian parameters are stored in gaussopt_arr
                    tempV += factor_ab*factor_symm*element_VGauss(c,ranges,norm4,jmat[a,c],jmat[b,c],la,La,lb,Lb,gaussopt,gamma_dict,complex_scaling,complex_scaling_angle)
                end
                
                for ivc in central_indices[c]
                    tempV += factor_ab*factor_symm*element_V(c,ivc,ranges,norm4,jmat[a,c],jmat[b,c],la,La,lb,Lb,w_interpol_arr,wn_interpol_arr,gamma_dict,reflect)
                end

                for ivc1D in contact1D_indices[c]
                    contactopt = contact1Dopt_arr[c][ivc1D] # Contact potential parameters are stored in contact1Dopt_arr
                    tempV += factor_ab*factor_symm*element_VContact(c,ranges,norm4,jmat[a,c],jmat[b,c],la,La,lb,Lb,contactopt,gamma_dict,complex_scaling,complex_scaling_angle)
                end
                
            end
        end
    end
    return tempV
end

# function for a single matrix element: norm-overlap S
function element_S(ranges,norm4,jab,la,La,lb,Lb,gamma_dict)
    
    mod(la+La+lb+Lb,2) == 1 && return zero(norm4) #parity of <a| and |b> must be the same, otherwise return 0
    
    (;nua,nub,NUa,NUb) = ranges;# named tuple unnecessary slow?
    (alpha,gamma,beta,delta) = jab # careful with order (gamma before beta)
    
    eta1 = nua*alpha^2 + NUa*gamma^2 + nub
    eta2 = 2*(nua*alpha*beta + NUa*gamma*delta)
    eta3 = nua*beta^2 + NUa*delta^2 + NUb
    sqeta = sqrt(eta1-eta2^2/4/eta3) # x^(n/2) == sqrt(x)^n, also for complex x (principal branch); integer powers are faster
    
    sumk = zero(norm4)
    for k = 0:la
        sumK = zero(norm4)
        for K = 0:La
            sums = zero(norm4)
            for s = 0:div(k+K+Lb,2)
                sums += gamma_dict[(la+La+lb+Lb-2*s+1)/2]/gamma_dict[s+1]/gamma_dict[k+K+Lb-2*s+1] * eta3^s*eta2^(k+K+Lb-2*s)/sqeta^(la+La+lb+Lb-2*s+1)
            end
            sumK += sums*gamma_dict[la+1]/(gamma_dict[la-k+1]*gamma_dict[k+1]) * gamma_dict[La+1]/(gamma_dict[La-K+1]*gamma_dict[K+1]) * alpha^(la-k) * beta^k * gamma^(La-K) * delta^K * gamma_dict[k+K+Lb+1]*(-1/2/eta3)^(k+K+Lb)
        end
        sumk += sumK
    end
    sumk *= norm4 * sqrt(pi/eta3)
    
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
    
    S0 = element_S(ranges,norm4,jab,la,La,lb,Lb,gamma_dict)
    temp = 4*nub^2*element_S(ranges,norm4,jab,la,La,lb+2,Lb,gamma_dict) - 2*(1+2*lb)*nub*S0
    lb >= 2 && ( temp += lb*(lb-1)*element_S(ranges,norm4,jab,la,La,lb-2,Lb,gamma_dict) )
    temp *= tr
    
    Temp = 4*NUb^2*element_S(ranges,norm4,jab,la,La,lb,Lb+2,gamma_dict) - 2*(1+2*Lb)*NUb*S0
    Lb >= 2 && ( Temp += Lb*(Lb-1)*element_S(ranges,norm4,jab,la,La,lb,Lb-2,gamma_dict) )
    Temp *= tR
    
    return temp+Temp
end


# function for a single matrix element: interaction V(r_c) based on numerical integration and interpolation over effective Gaussian range.
function element_V(c,iv,ranges,norm4,jac,jbc,la,La,lb,Lb,w_interpol_arr,wn_interpol_arr,gamma_dict,reflect::Bool)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    #GEM review: alpha <-> gamma; beta <-> delta
    
    eta5 = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2
    eta6 = 2*(nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC)
    eta7 = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2
    
    etaeff = eta5 - eta6^2/4/eta7
    Lsum = la+La+lb+Lb

    # For real ranges this is the usual lookup at log(etaeff). For complex ranges etaeff is complex
    # but confined to the sector |arg| <= atan(omega), over which the interpolant is tabulated in
    # two dimensions; interpol_lookup picks the right variant (see common/auxiliary.jl).
    for n = Lsum:-2:0
        wn_interpol_arr[n] = interpol_lookup(w_interpol_arr[c,iv,n],etaeff,reflect) # interpolated numerical integration value for the effective Gaussian range etaeff
    end
    
    zz = zero(promote_type(typeof(norm4),eltype(wn_interpol_arr)))
    sumk = zz
    for k = 0:la
        sumK = zz
        for K = 0:La
            sumkp = zz
            for kp = 0:lb
                sumKp = zz
                for Kp = 0:Lb                    
                    KK = k+K+kp+Kp; LL = la+La+lb+Lb;
                    
                    sums = zz
                    for s = 0:div(KK,2)
                        nn = Int(LL -2*s)
                        sums += eta7^s*eta6^(KK-2*s)/(gamma_dict[s+1]*gamma_dict[KK-2*s+1])*wn_interpol_arr[nn]
                    end
                    
                    # possibly costly due to many dict calls?
                    sumKp += sums * gamma_dict[la+1]/(gamma_dict[la-k+1]*gamma_dict[k+1]) * gamma_dict[La+1]/(gamma_dict[La-K+1]*gamma_dict[K+1]) * gamma_dict[lb+1]/(gamma_dict[lb-kp+1]*gamma_dict[kp+1]) * gamma_dict[Lb+1]/(gamma_dict[Lb-Kp+1]*gamma_dict[Kp+1]) * alphaAC^(la-k)*betaAC^k*gammaAC^(La-K)*deltaAC^K * alphaBC^(lb-kp)*betaBC^kp*gammaBC^(Lb-Kp)*deltaBC^Kp * gamma_dict[KK+1] * sqrt(pi/eta7)*(-1/2/eta7)^KK
                    
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
function element_VGauss(c,ranges,norm4,jac,jbc,la,La,lb,Lb,gaussopt,gamma_dict,complex_scaling,complex_scaling_angle)

    mod(la+La+lb+Lb,2) == 1 && return zero(norm4) #Gaussian potential is symmetric, hence parity of <a| and |b> must be the same, otherwise return 0
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    #GEM review: alpha <-> gamma; beta <-> delta
    
    v0,mu_g = gaussopt
    complex_scaling && (mu_g *= exp(2*im*complex_scaling_angle*pi/180)) # apply CSM factor to the range of the Gaussian potential
    
    LL = la+La+lb+Lb;
        
    eta5 = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2 + mu_g
    eta6 = 2*(nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC)
    eta7 = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2
    sqeta = sqrt(eta5-eta6^2/4/eta7) # see element_S
    
    zz = zero(promote_type(typeof(norm4),typeof(eta5)))
    sumk = zz
    for k = 0:la
        sumK = zz
        for K = 0:La
            sumkp = zz
            for kp = 0:lb
                sumKp = zz
                for Kp = 0:Lb                    
                    KK = k+K+kp+Kp;
                    
                    sums = zz
                    for s = 0:div(KK,2)
                        sums += eta7^s*eta6^(KK-2s)/(gamma_dict[s+1]*gamma_dict[KK-2*s+1]) * gamma_dict[(LL-2*s+1)/2]/sqeta^(LL-2*s+1)
                    end
                    
                    sumKp += sums * gamma_dict[la+1]/(gamma_dict[la-k+1]*gamma_dict[k+1]) * gamma_dict[La+1]/(gamma_dict[La-K+1]*gamma_dict[K+1]) * gamma_dict[lb+1]/(gamma_dict[lb-kp+1]*gamma_dict[kp+1]) * gamma_dict[Lb+1]/(gamma_dict[Lb-Kp+1]*gamma_dict[Kp+1]) * alphaAC^(la-k)*betaAC^k*gammaAC^(La-K)*deltaAC^K * alphaBC^(lb-kp)*betaBC^kp*gammaBC^(Lb-Kp)*deltaBC^Kp * gamma_dict[KK+1] * sqrt(pi/eta7)*(-1/2/eta7)^KK
                    
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
function element_VContact(c,ranges,norm4,jac,jbc,la,La,lb,Lb,contactopt,gamma_dict,complex_scaling,complex_scaling_angle)
    
    (;nua,nub,NUa,NUb) = ranges;
    (alphaAC,gammaAC,betaAC,deltaAC) = jac
    (alphaBC,gammaBC,betaBC,deltaBC) = jbc # careful with order (gamma before beta)
    
    v0,z0 = contactopt
    complex_scaling && (z0 *= exp(2*im*complex_scaling_angle*pi/180)) # apply CSM factor to the position of the contact interaction
    
    LL = la+La+lb+Lb;
        
    eta5 = nua*alphaAC^2 + NUa*gammaAC^2 + nub*alphaBC^2 + NUb*gammaBC^2
    eta6 = 2*(nua*alphaAC*betaAC + NUa*gammaAC*deltaAC + nub*alphaBC*betaBC + NUb*gammaBC*deltaBC)
    eta7 = nua*betaAC^2 + NUa*deltaAC^2 + nub*betaBC^2 + NUb*deltaBC^2
    
    zz = zero(promote_type(typeof(norm4),typeof(eta5),typeof(z0)))
    sumk = zz
    for k = 0:la
        sumK = zz
        for K = 0:La
            sumkp = zz
            for kp = 0:lb
                sumKp = zz
                for Kp = 0:Lb                    
                    KK = k+K+kp+Kp;

                    if z0 == 0
                        mod(KK,2) == 1 && continue
                        KK-LL != 0 && continue
                    end

                    sums = zz
                    for s = 0:div(KK,2) #expression via 1F1 exists too
                        sums += eta7^s*eta6^(KK-2s)/(gamma_dict[s+1]*gamma_dict[KK-2*s+1]) * z0^(LL-2*s) *exp(-(eta5-eta6^2/4/eta7)*z0^2) #eta7^(-KK) in external factor, like for Gaussian
                    end
                    
                    sumKp += sums * gamma_dict[la+1]/(gamma_dict[la-k+1]*gamma_dict[k+1]) * gamma_dict[La+1]/(gamma_dict[La-K+1]*gamma_dict[K+1]) * gamma_dict[lb+1]/(gamma_dict[lb-kp+1]*gamma_dict[kp+1]) * gamma_dict[Lb+1]/(gamma_dict[Lb-Kp+1]*gamma_dict[Kp+1]) * alphaAC^(la-k)*betaAC^k*gammaAC^(La-K)*deltaAC^K * alphaBC^(lb-kp)*betaBC^kp*gammaBC^(Lb-Kp)*deltaBC^Kp * gamma_dict[KK+1] * sqrt(pi/eta7)*(-1/2/eta7)^KK
                    
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


# Basis functions in matrix order, one entry per row/column: box, angular-momentum channel, and the range
# indices n,N (effective ones, i.e. doubled for a complex-ranged coordinate). Same loop order as the matrix.
function basis_list(size_params,nmax,Nmax)
    (;groupindex_arr,abvals_arr,lL_arr,starts,nbasis_total) = size_params
    basis = Vector{BasisFun}(undef,nbasis_total)
    for box in groupindex_arr
        a = abvals_arr[box][1] # either only 1 value in it, or the particles are identical and all lL_arr[a] are the same
        alpha = 0
        for (l,L) in lL_arr[a]
            for n = 1:nmax
                for N = 1:Nmax
                    alpha += 1
                    basis[starts[box] + alpha - 1] = (;box,l,L,n,N)
                end
            end
        end
    end
    return basis
end

# arguments of the matrix element between the basis functions rowi (bra) and coli (ket), see TempArgs1D
function element_args(rowi,coli,basis,abvals_arr,bvalsdiag,cvals,nu_arr,NU_arr,norm_arr,NORM_arr)
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
    Lsum = 0 # unnecessary in 1D
    return TempArgs1D(rowi,coli,(;nua,nub,NUa,NUb),norm4,a.l,a.L,b.l,b.L,Lsum,abvals_arr[a.box],bvals_new,cvals,factor_ab,abvals_arr[a.box],abvals_arr[b.box])
end