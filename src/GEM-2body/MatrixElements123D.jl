# Common file for all dimensions dim=1,2,3

element_S(lmax,nu1,nu2,dim) = (2*(nu1*nu2)^(1/2)/(nu1+nu2))^(lmax+dim/2) # norm-overlap

# kinetic energy: no common formula for all dimensions
function element_T(lmax,nu1,nu2,hbar,mur,dim) # maybe slow
    if dim == 1
        return element_T1D(lmax,nu1,nu2,hbar,mur)
    elseif dim == 2
        return element_T2D(lmax,nu1,nu2,hbar,mur)
    elseif dim == 3
        return element_T3D(lmax,nu1,nu2,hbar,mur)
    else
        error("Invalid dimension: $dim")
    end
end
element_T1D(lmax,nu1,nu2,hbar,mur) = -hbar^2/2/mur * element_S(lmax,nu1,nu2,1)*(-2*(1+2*lmax)*nu1+4*nu1^2*(lmax+1/2)/(nu1+nu2) + lmax*(lmax-1)*2*(nu1+nu2)/(2*lmax-1))# kinetic energy
element_T2D(lmax,nu1,nu2,hbar,mur) = -hbar^2/2/mur * (-2^(3+lmax)*(1+lmax)*(nu1*nu2)^((3+lmax)/2)*(nu1+nu2)^(-2-lmax)) # kinetic energy
element_T3D(lmax,nu1,nu2,hbar,mur) = hbar^2/2/mur * (2*lmax+3)*(nu1*nu2)^(1/2)*(2*(nu1*nu2)^(1/2)/(nu1+nu2))^(lmax+5/2) # kinetic energy


## old functions with if-conditions for dimensions
#= element_VGauss(lmax,nu1,nu2,v0,mu_g,dim) = v0*(2*(nu1*nu2)^(1/2)/(nu1+nu2+mu_g))^(lmax+dim/2) # gaussian interaction

function element_V(lmax,nu1,nu2,vint,gamma_dict,buf,dim) # maybe slow
    if dim == 1
        return element_V1D(lmax,nu1,nu2,vint,gamma_dict,buf,dim)
    elseif (dim == 2) || (dim == 3)
        return element_V23D(lmax,nu1,nu2,vint,gamma_dict,buf,dim)
    else
        error("Invalid dimension: $dim")
    end
end
function element_V1D(lmax,nu1,nu2,vint,gamma_dict,buf,dim) # 1D needs different integration limits
    return quadgk(r -> integrand(r,lmax,nu1,nu2,vint,dim),-Inf,0,Inf;segbuf=buf)[1]*norm(nu1,lmax,gamma_dict,dim)*norm(nu2,lmax,gamma_dict,dim)/2 # 1D needs additional factor 1/2
end

function element_V23D(lmax,nu1,nu2,vint,gamma_dict,buf,dim)
    return quadgk(r -> integrand(r,lmax,nu1,nu2,vint,dim),0,Inf;segbuf=buf)[1]*norm(nu1,lmax,gamma_dict,dim)*norm(nu2,lmax,gamma_dict,dim)
end =#


# for central interaction via numerical integration:
norm(nu,lmax,gamma_dict,dim) = (2*(2*nu)^(lmax+dim/2)/gamma_dict[lmax+dim/2])^(1/2) # result for integration over [0,\infty); in 1D we therefore need an additional factor 1/sqrt(2), but this is taken care of in element_V
integrand(r,lmax,nu1,nu2,vint,dim) = r^(2*lmax+(dim-1))*exp(-(nu1+nu2)*r^2)*vint(r)

function element_V(vint::GaussianPotential,lmax,nu1,nu2,gamma_dict,buf,dim,lower_limit,dimfac) #gaussian interaction
    v0 = vint.v0
    mu_g = vint.mu_g
    return v0*(2*(nu1*nu2)^(1/2)/(nu1+nu2+mu_g))^(lmax+dim/2)
end

function element_V(vint::CentralPotential,lmax,nu1,nu2,gamma_dict,buf,dim,lower_limit,dimfac) # 1D; needs different integration limits
    return quadgk(r -> integrand(r,lmax,nu1,nu2,vint,dim),lower_limit,Inf;segbuf=buf)[1]*norm(nu1,lmax,gamma_dict,dim)*norm(nu2,lmax,gamma_dict,dim) * dimfac
    # for 1D an additional segmentation at 0 can be useful. Not sure how to implement this in quadgk since providing the domain as a single argument [-Inf,0,Inf] yields an error
end

# if central potential is simply defined as a function: wrap into CentralPotential type
function element_V(vint::Function,lmax,nu1,nu2,gamma_dict,buf,dim,lower_limit,dimfac)
    return element_V(CentralPotential(vint),lmax,nu1,nu2,gamma_dict,buf,dim,lower_limit,dimfac)
end

## functions for the full (lower-triangular) matrices. requires already preallocated matrices S,T,V
# norm-overlap
function MatrixS(S,lmax,nu_arr,dim)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr) # only lower triangular
            S[nrow,ncol] = element_S(lmax,nu_arr[nrow]',nu_arr[ncol],dim)
        end
    end
end


## T: csm via global factor;
function MatrixT(T, lmax, nu_arr, hbar, mur, csm_bool, theta_csm, cr_bool, dim)
    fill_full = (csm_bool == 1 && cr_bool == 1) #for simultaneous csm and cr: need to fill full matrix
    
    for ncol in 1:lastindex(nu_arr)
        row_start = fill_full ? 1 : ncol # full fill if csm AND cr, otherwise only lower triangular
        for nrow in row_start:lastindex(nu_arr)
            T[nrow, ncol] = element_T(lmax, nu_arr[nrow]', nu_arr[ncol], hbar, mur, dim)
        end
    end
    
    csm_bool == 1 && (T .*= exp(-2*im*theta_csm*pi/180))
end

## V:
function MatrixV(V, lmax, nu_arr, vint_arr, gamma_dict, buf, csm_bool, theta_csm, cr_bool, dim)
    fill_full = (csm_bool == 1 && cr_bool == 1) #for simultaneous csm and cr: need to fill full matrix
    if csm_bool == 1 # to apply the csm in the basis functions:
        csmfacnu = exp(-2*im*theta_csm*pi/180)
    elseif csm_bool == 0
        csmfacnu = 1.0
    end
    
    if dim == 1
        lower_limit = -Inf # for 1D, the integration limits are different
        dimfac = 1/2 # additional factor for 1D to cope for the norm
    elseif (dim == 2) || (dim == 3)
        lower_limit = 0.0 # for 2D and 3D, the integration limits are the same
        dimfac = 1.0 # no additional factor for 2D and 3D
    else
        error("Invalid dimension: $dim")
    end

    # sum over all interactions
    for vint in vint_arr

        for ncol in 1:lastindex(nu_arr)
            row_start = fill_full ? 1 : ncol # full fill if csm and cr, otherwise only lower triangular
            nucol = nu_arr[ncol]*csmfacnu
            for nrow in row_start:lastindex(nu_arr)
                nurow = nu_arr[nrow]'*csmfacnu
                V[nrow, ncol] += element_V(vint, lmax, nurow, nucol, gamma_dict, buf, dim, lower_limit, dimfac)
            end
        end

    end
    
end

#= # central potential
function MatrixV_central(V, lmax, nu_arr, vint, gamma_dict, buf, fill_full, csmfacnu, dim)
    for ncol in 1:lastindex(nu_arr)
        row_start = fill_full ? 1 : ncol # full fill if csm and cr, otherwise only lower triangular
        nucol = nu_arr[ncol]*csmfacnu
        for nrow in row_start:lastindex(nu_arr)
            nurow = nu_arr[nrow]'*csmfacnu
            V[nrow, ncol] += element_V(lmax, nurow, nucol, vint, gamma_dict, buf,dim)
        end
    end
end

# Gaussian potential
function MatrixV_gauss(V, lmax, nu_arr, gaussopt, gamma_dict, buf, fill_full, csmfacnu, dim)
    gmax = lastindex(gaussopt)
    for gi in 1:gmax # there could be several gaussians
        v0, mu_g = gaussopt[gi][2:3]
        for ncol in 1:lastindex(nu_arr)
            row_start = fill_full ? 1 : ncol # full fill if csm and cr
            nucol = nu_arr[ncol]*csmfacnu
            for nrow in row_start:lastindex(nu_arr)
                nurow = nu_arr[nrow]'*csmfacnu
                V[nrow, ncol] += element_VGauss(lmax, nurow, nucol, v0, mu_g, dim)
            end
        end
    end
end =#




####  for coupled-channels:

## D: Terms which include derivatives
function MatrixWD(WD,lmax,nu_arr,WCC,DCC,gamma_dict,buf,csm_bool,theta_csm,diff_bool,dim)
    # DCC = [n,fun], where n denotes the derivative order and fun the functional form of the prefactor
    
    n = DCC[1]
    # Dictionary to map derivative order to corresponding function
    # this needs better treatment for vanishing terms due to lmax < n. maybe its not a problem...
    Dfun_dict = Dict(
    1 => (r, l, nu) -> l/r - 2*nu*r,
    2 => (r, l, nu) -> l*(l-1)/r^2 - 2*(1+2l)*nu + 4*nu^2*r^2,
    3 => (r, l, nu) -> l*(l-1)*(l-2)/r^3 - 6l^2*nu/r + 12*(l+1)*nu^2*r - 8*nu^3*r^3,
    4 => (r, l, nu) -> l*(l-1)*(l-2)*(l-3)/r^4 - 4*nu*l*(l-1)*(2l-1)/r^2 + 12*nu^2*(1+2l+2l^2) - 16*(3+2l)*nu^3*r^2 + 16*nu^4*r^4
    )    
    if !haskey(Dfun_dict, n)
        error("Derivative order n=$n > 4 not implemented")
    end
    
    ##maybe use n==0 instead of diff_bool?
    
    if diff_bool == 1
        Dnfun = Dfun_dict[n] # function that carries the effect of the n-th derivative
        wdfun = (r,l,nu) -> Dnfun(r,l,nu)*DCC[2](r)  + WCC(r) # D times prefactor function + W
    elseif diff_bool == 0
        wdfun = (r,l,nu) -> WCC(r)
    end

    if dim == 1
        lower_limit = -Inf # for 1D, the integration limits are different
        dimfac = 1/2 # additional factor for 1D
    elseif (dim == 2) || (dim == 3)
        lower_limit = 0.0 # for 2D and 3D, the integration limits are the same
        dimfac = 1.0 # no additional factor for 2D and 3D
    else
        error("Invalid dimension: $dim")
    end
    
    ## not adopted yet to support simultaneous CSM and CR!
    if csm_bool == 0
        MatrixWD_cr(WD,lmax,nu_arr,wdfun,gamma_dict,buf,dim,lower_limit,dimfac)
    elseif csm_bool == 1
        MatrixWD_csm(WD,lmax,nu_arr,wdfun,gamma_dict,buf,theta_csm,dim,lower_limit,dimfac)
    end
end


# no csm
function MatrixWD_cr(WD,lmax,nu_arr,wdfun,gamma_dict,buf,dim,lower_limit,dimfac)
    #roots,weights = rootsweights(1.0)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)
            wdfunr(r) = wdfun(r,lmax,nu_arr[ncol])# function definition in each loop iteration... slow?
            WD[nrow,ncol] = element_V(wdfunr,lmax,nu_arr[nrow]',nu_arr[ncol],gamma_dict,buf,dim,lower_limit,dimfac)
            #WD[nrow,ncol] = element_V_Num2(lmax,nu_arr[nrow]',nu_arr[ncol],wdfunr,gamma_dict,buf,roots,weights)
        end
    end
end

#csm
function MatrixWD_csm(WD,lmax,nu_arr,wdfun,gamma_dict,buf,theta_csm,dim,lower_limit,dimfac)
    csmfac=exp(-2*im*theta_csm*pi/180)
    #roots,weights = rootsweights(1.0)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)
            wdfunr(r) = wdfun(r,lmax,nu_arr[ncol]*csmfac)
            WD[nrow,ncol] = element_V(wdfunr,lmax,nu_arr[nrow]'*csmfac,nu_arr[ncol]*csmfac,gamma_dict,buf,dim,lower_limit,dimfac)
            #WD[nrow,ncol] = element_V_Num2(lmax,nu_arr[nrow]',nu_arr[ncol],r->vint_csm(r,wdfunr,theta_csm),gamma_dict,buf,roots,weights)
        end
    end
end