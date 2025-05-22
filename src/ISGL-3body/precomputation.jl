# functions to precompute repetedly used arrays within the ISGL program

# contains:
# - 

# for spins we need to add transformation coefficients from spins in jacobi-set a,b to c (and a to b)

function precompute_ISGL(phys_params,num_params,size_params,precomp_arrs,temp_arrs)    
    #Destructing Structs:
    #evtl einfacher alles zu destructen?
    #(;lmax,Lmax,gem_params,theta_csm,omega_cr,mu0,c_shoulder,kmax_interpol) = num_params
    
    (;mass_arr,J_tot,spin_arr) = phys_params
    (;lmax,Lmax,gem_params,kmax_interpol) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;cvals,s_arr,abI,factor_bf,s_complete,JsS_arr,JsS_complete,JlL_arr,JlL_complete,lL_complete,l_complete,nl,imax_dict,mij_arr_dict,kmax_dict,imaxSO_dict,mijSO_arr_dict,so_indices) = size_params
    (;gamma_dict,cleb_arr,spintrafo_dict,spinoverlap_dict,global6j_dict,facsymm_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr,Clmk_arr,Dlmk_arr,S_arr,SSO_arr) = precomp_arrs
    (;temp_clmk,temp_dlmk,temp_S,temp_D1,temp_D2) = temp_arrs
    
    # maximum JlL that needs to be considered:
    JlLmax = 0
    for c in cvals
        for is in keys(s_arr[c])
            for iss in keys(JsS_arr[c][is])
                JlLmax = Int64(max(JlLmax,maximum(JlL_arr[c][is][iss])))
            end
        end
    end
    
    #@show(J_tot,JlLmax)
    
    
    # how to update the values?!
    #precompute_gamma(precomp_arrs.gamma_dict,max(nmax,2*max(lmax,Lmax)+1)) #?!
    #gamma_dict = precompute_gamma(gamma_dict,max(nmax,2*max(lmax,Lmax)+1))
    precompute_gamma(gamma_dict,max(nmax,2*max(lmax,Lmax)+1+2)) #spin-orbit needs higher values, so let's give it +2...
    
    #cleb_arr = precompute_cleb(cleb_arr,J_tot,nl,l_complete)
    precompute_cleb(cleb_arr,JlL_complete,nl,l_complete)
    #@show(cleb_arr)
    
    #jmat = precompute_jmat(jmat,mass_arr)
    precompute_jmat(jmat,mass_arr)
    
    #murR_arr = precompute_murR(murR_arr,mass_arr)
    precompute_murR(murR_arr,mass_arr)
    
    #nu_arr,NU_arr = precompute_ranges(nu_arr,NU_arr,r1,rnmax,nmax,R1,RNmax,Nmax)
    precompute_ranges(nu_arr,NU_arr,r1,rnmax,nmax,R1,RNmax,Nmax)
    
    #norm_arr,NORM_arr = precompute_norms(norm_arr,NORM_arr,nu_arr,NU_arr,nl,l_complete,gamma_dict)
    precompute_norms(norm_arr,NORM_arr,nu_arr,NU_arr,nl,l_complete,gamma_dict)
    precompute_spintrafo(spin_arr,s_arr,JsS_arr,spintrafo_dict)
    precompute_facsymm(facsymm_dict,abI,factor_bf,spin_arr,l_complete,s_complete,s_arr)
    precompute_global6jfac(global6j_dict,J_tot,JsS_complete,JlL_complete)# also only for LS interaction, but necessary fur current fillTVS implementation
    precompute_spinoverlap(spinoverlap_dict,J_tot,spin_arr,s_arr,JsS_arr,spintrafo_dict) # same as above
    
    #Clmk_arr = precompute_Clmk(Clmk_arr,l_complete,nl,kmax_dict,gamma_dict,temp_clmk)
    precompute_Clmk(Clmk_arr,l_complete,nl,kmax_dict,gamma_dict,temp_clmk)
    
    #Dlmk_arr = precompute_Dlmk(Dlmk_arr,l_complete,nl,kmax_dict,temp_dlmk)
    precompute_Dlmk(Dlmk_arr,l_complete,nl,kmax_dict,temp_dlmk)
    
    #S_arr = precompute_S(S_arr,Clmk_arr,Dlmk_arr,J_tot,lL_complete,l_complete,imax_dict,mij_arr_dict,kmax_dict,gamma_dict,cleb_arr,temp_S,temp_D1,temp_D2)
    precompute_S(S_arr,Clmk_arr,Dlmk_arr,JlL_complete,lL_complete,l_complete,imax_dict,mij_arr_dict,kmax_dict,gamma_dict,cleb_arr,temp_S,temp_D1,temp_D2)
    
    # necessary only for spin-orbit interactions:
    if prod(isempty.(so_indices)) == false
        precompute_SSO(SSO_arr,Clmk_arr,Dlmk_arr,JlL_complete,lL_complete,l_complete,imaxSO_dict,mijSO_arr_dict,kmax_dict,gamma_dict,cleb_arr,temp_S,temp_D1,temp_D2)
    end
    
end

### gamma functions
function precompute_gamma(gamma_dict,ende)
    # returns the dict gamma_dict with keys "n" for each value of n between 0.5 and ende in steps of 0.5
    for i=0.5:0.5:ende+0.5
        gamma_dict[i] = gamma(i)
    end
    #return gamma_dict
end


### clebsch-gordan coefficients: for l-L coupling
@views @inbounds function precompute_cleb(cleb_arr,JlL_complete,nl,l_complete)
    # returns the (Offset-)array cleb_arr[la,ma,La,Ma,JlL] of Float64s for each combination of la,ma,La,Ma,JlL
    
    for JlL in JlL_complete
        for la in l_complete #lai = 1:nl
            #la = l_complete[lai]
            for La in l_complete #Lai = 1:nl
                #La = l_complete[Lai]
                for J = 0:JlL
                    if (JlL < abs(la-La)) || (JlL > la+La)
                        #cleb_arr[la,:,La,:,J,:] .= 0.0
                        continue
                    end
                    for ma = -la:la
                        for Ma = -La:La
                            M_tot = ma+Ma
                            if (abs(M_tot) > J)
                                #cleb_arr[la,ma,La,Ma,J,:] .= 0.0
                                continue
                            end
                            
                            #@show([la,ma,La,Ma,J,M_tot])
                            cleb_arr[la,ma,La,Ma,JlL] = PartialWaveFunctions.clebschgordan(la,ma,La,Ma,J,M_tot)
                            #println("la,ma,La,Ma,cleb=",la,", ",ma,", ",La,", ",Ma,", ",", ",cleb_arr[la,ma,La,Ma])
                        end
                    end
                end
            end
        end
    end
    
    #return cleb_arr    
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

function jcbtr(i::Int64,f::Int64,m_arr::Array{<:Real}) # Für fixes m1,m2,m3. Alternativ: ma,mb,mc als arugment, dann muss es aber jedesmal entsprechend aufgerufen werden!    
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
function precompute_norms(norm_arr,NORM_arr,nu_arr,NU_arr,nl,l_complete,gamma_dict)
    # returns the array norm_arr[l,n] for each combination of n,l. same for NORM_arr
    
    norm(nu,l) = (2*(2*nu)^(l+3/2)/gamma_dict[l+3/2])^(1/2);
    
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
    
    #return norm_arr,NORM_arr
end

#global 6j prefactor for spin-orbit interaction
function precompute_global6jfac(global6j_dict,J_tot,JsS_complete,JlL_complete)
    
    for JsSa in JsS_complete
        for JsSb in JsS_complete
            for JlLa in JlL_complete
                for JlLb in JlL_complete
                    global6j_dict[JsSa,JsSb,JlLa,JlLb] = (-1)^(J_tot+JsSa+JlLb)*WignerSymbols.wigner6j(J_tot,JsSa,JlLa,1,JlLb,JsSb) # 1 due to spin-orbit interaction. probably 2 for tensor interaction
                end
            end
        end
    end
    
end


#spintrafo is for the U-functions
function precompute_spintrafo(spin_arr,s_arr,JsS_arr,spintrafo_dict)
    # returns an array for all possible spin transformation coefficients
    # there should be some simplifications such that we dont need to calculate them all (440-444)
    for c = 1:3 # for simplicity all c-values, even if not used
        for cprime = 1:3
            for (is,sc) in enumerate(s_arr[c])
                for (isp,scp) in enumerate(s_arr[cprime])
                    for JsSc in JsS_arr[c][is]
                        # here we implicitly assume that JsSc is the same as JsScprime. Otherwise the spin-overlap is zero. This is however treated explicitly only in fill_TVS's sab,tab,vab
                        spintrafo_dict[c,cprime,JsSc,sc,scp] = spintrafo_fun(c,cprime,JsSc,sc,scp,spin_arr)
                        #@show(c,cprime,JsSc,JsS_arr[cprime][isp],sc,scp,spintrafo_dict[c,cprime,JsSc,sc,scp])
                    end
                end
            end
        end
    end
end

function spintrafo_fun(c,cp,JsS,sc,scp,spin_arr)
    a,b = mod(c,3)+1,mod(c+1,3)+1 # a=c+1,b=c-1, for given c^
    za,zb,zc = spin_arr[[a,b,c]]
    
    # cyclic permutation of jacobi sets needs to be treated separately:
    if cp == c # racah coefficients cannot cope for coupling in the same basis?
        #println("cp=$cp == c=$c")
        if sc == scp
            return 1.0
        else
            return 0.0
        end
    elseif cp == mod(c,3)+1 # cp = c+1
        #println("cp=$cp == c+1=$(mod((c+1),3)+1)")
        return ((2*sc+1)*(2*scp+1))^(1/2) * WignerSymbols.racahW(za,zb,JsS,zc,sc,scp) *(-1.0)^(za+scp-JsS)
    elseif cp == mod(c+1,3)+1 # cp = c-1
        #println("cp=$cp == c-1=$(mod(c,3)+1)")
        return ((2*sc+1)*(2*scp+1))^(1/2) * WignerSymbols.racahW(zb,za,JsS,zc,sc,scp) *(-1.0)^(2*za+2*zb+zc-sc-JsS)
    end    
end

# matrix element of the spin part for l-s interactions
# this corresponds only to the double-bar matrix element of the spin part (Eq. 60, if we dont change the numbering.) The "global" 6j-symbol carrying the JlL(a,b) and JsS(a,b) dependenceies should be calculated within TVSfill. Therefore this function is independent of JlL!
function precompute_spinoverlap(spinoverlap_dict,J_tot,spin_arr,s_arr,JsS_arr,spintrafo_dict)
    # Eq. (63) in new notes
    J = J_tot
    for a in 1:3
        for b in 1:3
            for c in 1:3
                for (isa,sa) in enumerate(s_arr[a])
                    for (isb,sb) in enumerate(s_arr[b])
                        for (issa,JsSa) in enumerate(JsS_arr[a][isa])
                            for (issb,JsSb) in enumerate(JsS_arr[b][isb]) # is there a better way to use only JsSa values that are also allowed for b?
                                #@show([a,b,c,JsS,sa,sb,JlL])
                                
                                zc = spin_arr[c]
                                uabc = 0.0
                                for sc in s_arr[c]
                                    uac = spintrafo_dict[a,c,JsSa,sa,sc]
                                    ubc = conj(spintrafo_dict[b,c,JsSb,sb,sc])
                                    uabc += uac*ubc * (-1)^(sc+zc+JsSb+1) * sqrt(sc*(sc+1)*(sc+2)) * WignerSymbols.wigner6j(sc,JsSa,zc,JsSb,sc,1)
                                    
                                    #check if this cures comparison to 2-body system:
                                    #zc == 0 && ( uabc *= sqrt((2*sc+1)/(sc+2)) ) ## THIS IS NOT SATISFACTORY!, and not completely correct either...

                                    #@show([a,b,c,sa,sb,sc,JsSa,JsSb,zc])
                                    #@show([uac,ubc,uabc])
                                    #this should be checked again:
                                    #= #check: indeed no contribution for sc = 0
                                    if sc == 0.0
                                        zzz = sqrt((sc+2)/(2sc+1)) * (sc*(sc+1) + JsS*(JsS+1) - zc*(zc+1))
                                        uac*ubc*zzz != 0 && @show([uac,ubc,zzz])
                                    end =#
                                    
                                end
                                uabc *= sqrt((2*JsSa + 1)*(2*JsSb + 1)) # *hbar? hbar set to 1.0 here...
                                
                                x = sqrt(sa*(sa+1)*(sa+2))
                                y = sqrt(sa*(sa+1)*(2*sa+1))


                                #@show([uabc,x,y])

                                spinoverlap_dict[a,b,c,JsSa,JsSb,sa,sb] = uabc
                            end
                        end
                    end
                end    
            end
        end
    end
end



# for a 3-element vector v this outputs the sum of the two entries, whose index does not match "a"
# in this case it is the sum of zb and zc
function ex_sum(v, a)
    sum(v[setdiff(1:3, a)])
end

function facsymm(a,b,abI,la,lb,factor_bf,spin_arr,sa,sb)
    #abI == 0 && return 1
    factor_symm = 1
    a == abI && (factor_symm *= factor_bf*(-1)^la*(-1)^(ex_sum(spin_arr, a)-sa))#;println("a = $a == abI = $abI"))
    b == abI && (factor_symm *= factor_bf*(-1)^lb*(-1)^(ex_sum(spin_arr, b)-sb))#;println("b = $b == abI = $abI"))
    #@show([a,b,abI,la,lb,factor_bf,sa,sb,factor_symm])
    return factor_symm
end

function precompute_facsymm(facsymm_dict,abI,factor_bf,spin_arr,l_complete,s_complete,s_arr)
    # solution only via dict?!
    for a = 1:3
        for b = 1:3
            for la in l_complete
                for lb in l_complete
                    for sa in s_arr[a]
                        for sb in s_arr[b]
                            facsymm_dict[a,b,la,lb,sa,sb] = facsymm(a,b,abI,la,lb,factor_bf,spin_arr,sa,sb)
                        end
                    end
                end
            end
        end
    end
    #return facsymm_arr
end


### Clmk
@views @inbounds function precompute_Clmk(Clmk_arr,l_complete,nl,kmax_dict,gamma_dict,temp_arr)
    # returns the array of C[l,m,k] for each necessary l,m,k combination
    #temp_arr = MVector{findmax(kmax_dict)[1], Float64}(zeros(findmax(kmax_dict)[1]))
    for l in l_complete #lindex = 1:nl
        #l = l_complete[lindex]
        for m = 0:l
            Clmk_arr[l,m,1:kmax_dict[l,m]] .= clmk_fun(l,m,kmax_dict[l,m],gamma_dict,temp_arr)
            m>0 && (@. Clmk_arr[l,-m,1:kmax_dict[l,-m]] = (-1)^m*Clmk_arr[l,m,1:kmax_dict[l,m]])
        end
    end
    #return Clmk_arr
end


@views @inbounds function clmk_fun(l,m,kmax,gamma_dict,temp_arr)
    if m<0
        minus_coeff = (-1)^m;
    else
        minus_coeff = 1.0;
    end
    m = abs(m);
    
    k = 0
    bino(n,m) = gamma_dict[n+1]/(gamma_dict[m+1]*gamma_dict[n-m+1])
    
    for j = 0:Int64(floor((l-m)/2))
        Almj = ((2*l+1)*(gamma_dict[l-m+1]/4/pi/gamma_dict[l+m+1]))^(1/2)*gamma_dict[l+m+1]/2.0^m*(-1.0)^j/4.0^j/gamma_dict[j+1]/gamma_dict[m+j+1]/gamma_dict[l-m-2*j+1];
        for s = 0:(l-m-2*j)
            for t = 0:(j+m)
                for u = 0:j
                    k=k+1
                    #println(k)
                    k>kmax && break
                    temp_arr[k] = minus_coeff*bino(l-m-2*j,s)*bino(j+m,t)*bino(j,u)*(-1.0)^(l-s-t-u)*Almj/4^l;
                end
            end
        end
    end
    
    return temp_arr[1:kmax];
end



### Dlmk
@views @inbounds function precompute_Dlmk(Dlmk_arr,l_complete,nl,kmax_dict,temp_arr)
    # returns a big array D[l,m,k,1:3] of [Dx,Dy,Dz] for each l,m,k combination necessary
    for l in l_complete #lindex = 1:nl
        #l = l_complete[lindex]
        for m = 0:l
            Dlmk_arr[l,m,1:kmax_dict[l,m],:] .= dlmk_fun(l,m,kmax_dict[l,m],temp_arr)
            m>0 && (Dlmk_arr[l,-m,1:kmax_dict[l,-m],:] .= conj.(Dlmk_arr[l,m,1:kmax_dict[l,m],:]))
        end
    end
    #return Dlmk_arr
end

@views @inbounds function dlmk_fun(l,m,kmax,temp_arr)
    #temp_arr = Array{ComplexF64}(undef,kmax,3)*0.0
    #temp_arr = ComplexF64(zeros(kmax,3))
    k = 0
    for j = 0:Int64(floor((l-m)/2))
        for s = 0:(l-m-2j)
            for t = 0:(j+m)
                for u = 0:j
                    k=k+1
                    k>kmax && break
                    temp_arr[k,1] = 2*t-(j+m)+2*u-j;
                    temp_arr[k,2] = im*(2*t-(j+m)-2*u+j); # needs always real() wrapping after usage. alternatively: code here as real and keep -1 in mind!
                    temp_arr[k,3] = 2*s - (l-m-2*j);                        
                end
            end
        end
    end
    return temp_arr[1:kmax,:]
end


## for central interaction:
### S_coeff
@views @inbounds function precompute_S(S_arr,Clmk_arr,Dlmk_arr,JlL_complete,lL_complete,l_complete,imax_dict,mij_arr_dict,kmax_dict,gamma_dict,cleb_arr,temp_arr,tempD1,tempD2)
    # returns the OffsetArray S_arr[la,La,lb,Lb,1:maximax]. only 1:imax(la,La,lb,Lb) is filled with nonzeros.
    
    for JlL in JlL_complete
        for i in keys(lL_complete)
            (la,La) = lL_complete[i] # reihenfolge wichtig?!
            for j in keys(lL_complete)
                (lb,Lb) = lL_complete[j] # reihenfolge wichtig?!
                #=             if lb < la && Lb < La
                    S_arr[la,La,lb,Lb,:] .= S_arr[lb,Lb,la,La,:]
                    println("hi: la,La,lb,Lb=",la,",",La,",",lb,",",Lb)
                    continue
                end
                println("la,La,lb,Lb=",la,",",La,",",lb,",",Lb) =#
                #(lb >la || Lb > La) && continue #!!?? Achtung: Vereinfachung zum debuggen!!            
                S_arr[la,La,lb,Lb,JlL,:] .= Si_fun(JlL,la,La,lb,Lb,imax_dict[((la,La),(lb,Lb))],mij_arr_dict[((la,La),(lb,Lb))],kmax_dict,gamma_dict,cleb_arr,Clmk_arr,Dlmk_arr,temp_arr,tempD1,tempD2)
            end
        end
    end
    #return S_arr
end



# mit Wigner Eckart Vereinfachung: -> bis zu Faktor 100 schneller!?
@views @inbounds function Si_fun(JlL,la,La,lb,Lb,imax,mij_arr_i,kmax_dict,gamma_dict,cleb_arr,Clmk_arr,Dlmk_arr,temp_arr,D1,D2)
    # returns a vector S(la,La,lb,Lb,J)[i=1,i=2,...,i=maximax] only i=1:imax(la,La,lb,LB) is filled with nonzeros
    
    M_tot = JlL # why?! should be independent of M -> choose as we like
    
    # reset to zero (the same array is used in each function call)
    for i = 1:lastindex(temp_arr)
        temp_arr[i] = 0.0
    end
    
    for ma = la:la # statt -la:la (Wigner-Eckart!)
        Ma = M_tot - ma
        abs(Ma) > La && continue # its a guess atm: not sure if I fully understand this... dependency on M_tot?
        
        for mb = -lb:lb
            Mb = M_tot - mb
            abs(Mb) > Lb && continue # its a guess atm: not sure if I fully understand this... dependency on M_tot?
            
            cleb_fac = 1/cleb_arr[la,ma,La,Ma,JlL]*cleb_arr[lb,mb,Lb,Mb,JlL]
            
            for ka=1:kmax_dict[la,ma]
                for Ka=1:kmax_dict[La,Ma]
                    for kb=1:kmax_dict[lb,mb]
                        for Kb=1:kmax_dict[Lb,Mb]
                            # ccleb_fac = RCFAC
                            ccleb_fac = Clmk_arr[la,ma,ka]*Clmk_arr[La,Ma,Ka]*Clmk_arr[lb,mb,kb]*Clmk_arr[Lb,Mb,Kb] * cleb_fac
                            
                            D1 .= conj.(Dlmk_arr[la,ma,ka,:])
                            D2 .= conj.(Dlmk_arr[La,Ma,Ka,:])
                            D3 = Dlmk_arr[lb,mb,kb,:]
                            D4 = Dlmk_arr[Lb,Mb,Kb,:]
                            
                            D12 = real(transpose(D1)*D2)
                            D13 = real(transpose(D1)*D3)
                            D14 = real(transpose(D1)*D4)
                            D23 = real(transpose(D2)*D3)
                            D24 = real(transpose(D2)*D4)
                            D34 = real(transpose(D3)*D4)
                            
                            # calculation of S for all values of i, so it can be accessed later in the sum over i in fillTVS.
                            for i = 1:imax # loop over i inside, as the rest can be done independently
                                m12 = mij_arr_i[i][1]
                                m13 = mij_arr_i[i][2]
                                m14 = mij_arr_i[i][3]
                                m23 = mij_arr_i[i][4]
                                m24 = mij_arr_i[i][5]
                                m34 = mij_arr_i[i][6]
                                
                                # mij_fac = FACTNO:
                                mij_fac = 1/(gamma_dict[m12+1.0]*gamma_dict[m13+1.0]*gamma_dict[m14+1.0]*gamma_dict[m23+1.0]*gamma_dict[m24+1.0]*gamma_dict[m34+1.0])
                                
                                temp_arr[i] += D12^m12*D13^m13*D14^m14*D23^m23*D24^m24*D34^m34*ccleb_fac*mij_fac
                            end
                            
                        end
                    end
                end
            end
        end
    end
    
    return temp_arr
end



## for spin-orbit interaction:
### S_coeff
@views @inbounds function precompute_SSO(SSO_arr,Clmk_arr,Dlmk_arr,JlL_complete,lL_complete,l_complete,imaxSO_dict,mijSO_arr_dict,kmax_dict,gamma_dict,cleb_arr,temp_arr,tempD1,tempD2)
    # returns the OffsetArray S_arr[la,La,lb,Lb,1:maximax]. only 1:imax(la,La,lb,Lb) is filled with nonzeros.
    
    
    for JlLa in JlL_complete
        for JlLb in JlL_complete
            (abs(JlLa - JlLb) <= 1 <= JlLa + JlLb) == false && continue # only for allowed values of JlLa,JlLb
            for (la,La) in lL_complete
                (abs(la - La) <= JlLa <= la + La) == false && continue # only for allowed values of la,La,JlLa
                for (lb,Lb) in lL_complete
                    (abs(lb - Lb) <= JlLb <= lb + Lb) == false && continue # only for allowed values of lb,Lb,JlLb
                    for v = 1:6
                        imaxSO_dict[((la,La),(lb,Lb),v)] == 0 && continue # skip if no combinations of mijSO exist
                        
                        SSO_arr[la,La,lb,Lb,JlLa,JlLb,v,:] .= Si_funSO(JlLa,JlLb,la,La,lb,Lb,imaxSO_dict[((la,La),(lb,Lb),v)],mijSO_arr_dict[((la,La),(lb,Lb),v)],kmax_dict,gamma_dict,cleb_arr,Clmk_arr,Dlmk_arr,temp_arr,tempD1,tempD2,v)
                    end
                end
            end
        end
    end
    
end


# to assign pairs of Di,Dj based on v
function get_Ds(v, D1, D2, D3, D4)
    if v == 1
        return D1, D2
    elseif v == 2
        return D1, D3
    elseif v == 3
        return D1, D4
    elseif v == 4
        return D2, D3
    elseif v == 5
        return D2, D4
    else  # v == 6
        return D3, D4
    end
end


@views @inbounds function Si_funSO(JlLa,JlLb,la,La,lb,Lb,imax,mij_arr_i,kmax_dict,gamma_dict,cleb_arr,Clmk_arr,Dlmk_arr,temp_arr,D1,D2,v)
    # returns a vector S(la,La,lb,Lb,J)[i=1,i=2,...,i=maximax] only i=1:imax(la,La,lb,LB) is filled with nonzeros
    
    # reset to zero (the same array is used in each function call)
    for i = 1:lastindex(temp_arr)
        temp_arr[i] = 0.0
    end
    
    #done outside already
    #=     if (abs(JlLa - JlLb) <= 1 <= JlLa + JlLb) == false
        return temp_arr
    end =#
    
    M_tot = 0 # not really Mtot, but the index? (not rank) of the tensor operator
    MlLa = min(JlLa, JlLb)
    MlLb = MlLa
    
    mamax = min(la,MlLa+La) # from legacy code, origin not clear
    mamin = max(-la,MlLa-La)
    mbmax = min(lb,MlLb+Lb)
    mbmin = max(-lb,MlLb-Lb)
    if mamax < mamin || mbmax < mbmin
        return temp_arr
    end
    
    # MlLa = MlLb+0 is always true

    # for the cases when cleb_arr was not precomputed (mostly, if lmin,Lmin > 0)
    if cleb_arr[1,M_tot,JlLb,MlLb,JlLa] == 0
        clebJba = PartialWaveFunctions.clebschgordan(1,M_tot,JlLb,MlLb,JlLa,MlLa)
    else
        clebJba = cleb_arr[1,M_tot,JlLb,MlLb,JlLa]
    end
    clebjllabJ = (-1)^(-1+JlLa+JlLb)*sqrt(2*JlLa+1) / clebJba # this is my formula. legacy code has some weird factor, which i cannot reproduce since CLEBSH does not output clebsch-gordan...
    
    for ma = mamin:mamax # no easy truncation from Wigner-Eckart
        Ma = MlLa - ma
        abs(Ma) > La && continue
        
        for mb = mbmin:mbmax
            Mb = MlLb - mb
            abs(Mb) > Lb && continue 
            
            cleb_fac = cleb_arr[la,ma,La,Ma,JlLa]*cleb_arr[lb,mb,Lb,Mb,JlLb]*clebjllabJ
            
            #=             if la == La == lb == Lb == v == JlLa == JlLb == 1
                ca = cleb_arr[la,ma,La,Ma,JlLa]
                cb = cleb_arr[lb,mb,Lb,Mb,JlLb]
                @show([la,ma,La,Ma,JlLa,MlLa])
                @show([lb,mb,Lb,Mb,JlLb,MlLb])
                @show([ca,cb])
                @show(ca*cb*clebjllabJ)
            end =#
            
            for ka=1:kmax_dict[la,ma]
                for Ka=1:kmax_dict[La,Ma]
                    for kb=1:kmax_dict[lb,mb]
                        for Kb=1:kmax_dict[Lb,Mb]
                            # ccleb_fac = RCFAC
                            ccleb_fac = Clmk_arr[la,ma,ka]*Clmk_arr[La,Ma,Ka]*Clmk_arr[lb,mb,kb]*Clmk_arr[Lb,Mb,Kb] * cleb_fac
                            
                            # D_i, the i is only determined by la,ma,ka?
                            D1 .= conj.(Dlmk_arr[la,ma,ka,:])
                            D2 .= conj.(Dlmk_arr[La,Ma,Ka,:])
                            D3 = Dlmk_arr[lb,mb,kb,:]
                            D4 = Dlmk_arr[Lb,Mb,Kb,:]
                            
                            D12 = real(transpose(D1)*D2)
                            D13 = real(transpose(D1)*D3)
                            D14 = real(transpose(D1)*D4)
                            D23 = real(transpose(D2)*D3)
                            D24 = real(transpose(D2)*D4)
                            D34 = real(transpose(D3)*D4)
                            
                            # cleaner call, moved if-tree to extra function get_Ds
                            Di, Dj = get_Ds(v, D1, D2, D3, D4)                            
                            Dcrossz = im*(Di[1]*Dj[2]-Di[2]*Dj[1]) # DixDjy - DiyDjx
                            
                            # calculation of S for all values of i, so it can be accessed later in the sum over i in fillTVS.
                            for i = 1:imax # loop over i inside, as the rest can be done independently
                                m12 = mij_arr_i[i][1] # automatically has the correct v
                                m13 = mij_arr_i[i][2]
                                m14 = mij_arr_i[i][3]
                                m23 = mij_arr_i[i][4]
                                m24 = mij_arr_i[i][5]
                                m34 = mij_arr_i[i][6]
                                
                                # mij_fac = FACTNO:
                                mij_fac = 1/(gamma_dict[m12+1.0]*gamma_dict[m13+1.0]*gamma_dict[m14+1.0]*gamma_dict[m23+1.0]*gamma_dict[m24+1.0]*gamma_dict[m34+1.0])
                                
                                sss = D12^m12*D13^m13*D14^m14*D23^m23*D24^m24*D34^m34*ccleb_fac*mij_fac*Dcrossz
                                temp_arr[i] += sss
                            end
                            
                        end
                    end
                end
            end
        end
    end
    
    return temp_arr
end