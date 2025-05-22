# Separate file for the matrix elements. Includes convenience functions for the full matrices T,V,S (only lower-triangular)

## functions to calculate single matrix elements (normalization prefactor is already included):
element_S2D(lmax,nu1,nu2) = (2*(nu1*nu2)^(1/2)/(nu1+nu2))^(lmax+1) # norm-overlap
element_T2D(lmax,nu1,nu2,hbar,mur) = -hbar^2/2/mur * (-2^(3+lmax)*(1+lmax)*(nu1*nu2)^((3+lmax)/2)*(nu1+nu2)^(-2-lmax)) # kinetic energy
element_VGauss2D(lmax,nu1,nu2,v0,mu_g) = v0*(2*(nu1*nu2)^(1/2)/(nu1+nu2+mu_g))^(lmax+1) # gaussian interaction;
# central interaction via numerical integration:
norm(nu,lmax,gamma_dict) = 2*(2^lmax*nu^(lmax+1)/gamma_dict[lmax+1])^(1/2);

vint_csm(r,vint,theta_csm) = vint(r*exp(im*theta_csm*pi/180))

integrand(r,lmax,nu1,nu2,vint) = r^(2*lmax+1)*exp(-(nu1+nu2)*r^2)*vint(r)
function element_V2D(lmax,nu1,nu2,vint,gamma_dict,buf)
    return quadgk(r -> integrand(r,lmax,nu1,nu2,vint),0,Inf;segbuf=buf)[1]*norm(nu1,lmax,gamma_dict)*norm(nu2,lmax,gamma_dict) # providing additional point of interval-segmentation can help!
end

# hard-coded Gauss-Quadrature for possibly higher speed: currently NOT WORKING!
#= using FastGaussQuadrature
function rootsweights(alpha)
    alpha = 1.0
    roots11, weights11 = gausslegendre(75); # using FastGaussQuadrature.jl
    
    # 24-point Gauss-Legendre quadrature for the interval [-1,1]
    #roots11 = SA[-0.9951872199970214,-0.9747285559713095,-0.9382745520027328,-0.8864155270044011,-0.820001985973903,-0.7401241915785544,-0.6480936519369755,-0.5454214713888396,-0.4337935076260451,-0.31504267969616334,-0.1911188674736163,-0.06405689286260563,0.06405689286260563,0.1911188674736163,0.31504267969616334,0.4337935076260451,0.5454214713888396,0.6480936519369755,0.7401241915785544,0.820001985973903,0.8864155270044011,0.9382745520027328,0.9747285559713095,0.9951872199970214]
    #weights11 = SA[0.01234122979998764,0.028531388628933584,0.0442774388174198,0.05929858491543669,0.07334648141108027,0.08619016153195319,0.0976186521041139,0.10744427011596554,0.11550566805372557,0.12167047292780338,0.12583745634682833,0.1279381953467522,0.1279381953467522,0.12583745634682833,0.12167047292780338,0.11550566805372557,0.10744427011596554,0.0976186521041139,0.08619016153195319,0.07334648141108027,0.05929858491543669,0.0442774388174198,0.028531388628933584,0.01234122979998764]
    
    # adaption to the interval [0,Inf]
    roots = 2*alpha ./(1 .- roots11) .- alpha
    weights = 2*alpha .*weights11./(1 .- roots11).^2
    return roots,weights
end
function element_V_Num2(lmax,nu1,nu2,vint,gamma_dict,buf,weights,roots)
    return weights' * integrand.(roots,lmax,nu1,nu2,vint)*norm(nu1,lmax,gamma_dict)*norm(nu2,lmax,gamma_dict)
end =#

## functions for the full (lower-triangular) matrices. requires already preallocated matrices S,T,V
## S:
function MatrixS2D(S,lmax,nu_arr)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)
            S[nrow,ncol] = element_S2D(lmax,nu_arr[nrow]',nu_arr[ncol])
            #@show([nrow,ncol,S[nrow,ncol]])
        end
    end
end

## T: csm via global factor
function MatrixT2D(T,lmax,nu_arr,hbar,mur,csm_bool,theta_csm)
    MatrixT2D_cr(T,lmax,nu_arr,hbar,mur)
    csm_bool == 1 && (T .*= exp(-2*im*theta_csm*pi/180))
end
function MatrixT2D_cr(T,lmax,nu_arr,hbar,mur)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)
            T[nrow,ncol] = element_T2D(lmax,nu_arr[nrow]',nu_arr[ncol],hbar,mur)
            #@show([nrow,ncol,T[nrow,ncol]])
        end
    end
end

## V:
function MatrixV2D(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt,csm_bool,theta_csm)
    if csm_bool == 0
        MatrixV2D_cr(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt)
    elseif csm_bool == 1
        MatrixV2D_csm(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt,theta_csm)
    end
end

# no csm
function MatrixV2D_cr(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt)
    if gaussopt[1][1] == 0 # all gaussblooleans are the same -> decide by the first one
        #roots,weights = rootsweights(1.0)
        for ncol = 1:lastindex(nu_arr)
            for nrow = ncol:lastindex(nu_arr)  
                V[nrow,ncol] = element_V2D(lmax,nu_arr[nrow]',nu_arr[ncol],vint,gamma_dict,buf)
                #@show([nrow,ncol,V[nrow,ncol]])
                #V[nrow,ncol] = element_V_Num2(lmax,nu_arr[nrow]',nu_arr[ncol],vint,gamma_dict,buf,roots,weights)
            end
        end
    elseif gaussopt[1][1] == 1
        gmax = lastindex(gaussopt) # how many gaussian potentials
        for gi in 1:gmax        
            v0,mu_g = gaussopt[gi][2:3]
            for ncol = 1:lastindex(nu_arr)
                for nrow = ncol:lastindex(nu_arr)  
                    V[nrow,ncol] += element_VGauss2D(lmax,nu_arr[nrow]',nu_arr[ncol],v0,mu_g)
                end
            end            
        end
    end
end

#csm
function MatrixV2D_csm(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt,theta_csm)
    if gaussopt[1][1] == 0 # all gaussblooleans are the same -> decide by the first one
        #roots,weights = rootsweights(1.0)
        for ncol = 1:lastindex(nu_arr)
            for nrow = ncol:lastindex(nu_arr)  
                V[nrow,ncol] = element_V2D(lmax,nu_arr[nrow]',nu_arr[ncol],r->vint_csm(r,vint,theta_csm),gamma_dict,buf)
                #V[nrow,ncol] = element_V_Num2(lmax,nu_arr[nrow]',nu_arr[ncol],r->vint_csm(r,vint,theta_csm),gamma_dict,buf,roots,weights)
            end
        end
    elseif gaussopt[1][1] == 1
        gmax = lastindex(gaussopt) # how many gaussian potentials
        for gi in 1:gmax        
            v0,mu_g = gaussopt[gi][2:3]
            for ncol = 1:lastindex(nu_arr)
                for nrow = ncol:lastindex(nu_arr)  
                    V[nrow,ncol] += element_V2DGauss(lmax,nu_arr[nrow]',nu_arr[ncol],v0,mu_g*exp(2*im*theta_csm*pi/180))
                end
            end            
        end
    end
end


## D: Terms which include derivatives
## THIS IS NOT ADAPTED FOR 2D YET!
function MatrixWD(WD,lmax,nu_arr,WCC,DCC,gamma_dict,buf,csm_bool,theta_csm,diff_bool)
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
    
    if csm_bool == 0
        MatrixWD_cr(WD,lmax,nu_arr,wdfun,gamma_dict,buf)
    elseif csm_bool == 1
        MatrixWD_csm(WD,lmax,nu_arr,wdfun,gamma_dict,buf,theta_csm)
    end
end

# no csm
function MatrixWD_cr(WD,lmax,nu_arr,wdfun,gamma_dict,buf)
    #roots,weights = rootsweights(1.0)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)
            wdfunr(r) = wdfun(r,lmax,nu_arr[ncol])# function definition in each loop iteration... slow?
            WD[nrow,ncol] = element_V2D(lmax,nu_arr[nrow]',nu_arr[ncol],wdfunr,gamma_dict,buf)
            #WD[nrow,ncol] = element_V_Num2(lmax,nu_arr[nrow]',nu_arr[ncol],wdfunr,gamma_dict,buf,roots,weights)
        end
    end
end

#csm
function MatrixWD_csm(WD,lmax,nu_arr,wdfun,gamma_dict,buf,theta_csm)
    #roots,weights = rootsweights(1.0)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)
            wdfunr(r) = wdfun(r,lmax,nu_arr[ncol])
            WD[nrow,ncol] = element_V2D(lmax,nu_arr[nrow]',nu_arr[ncol],r->vint_csm(r,wdfunr,theta_csm),gamma_dict,buf)
            #WD[nrow,ncol] = element_V_Num2(lmax,nu_arr[nrow]',nu_arr[ncol],r->vint_csm(r,wdfunr,theta_csm),gamma_dict,buf,roots,weights)
        end
    end
end