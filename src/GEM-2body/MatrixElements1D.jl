## 1D
# Separate file for the matrix elements. Includes convenience functions for the full matrices T,V,S (only lower-triangular)

## functions to calculate single matrix elements (normalization prefactor is already included):
element_S1D(lmax,nu1,nu2) = (2*(nu1*nu2)^(1/2)/(nu1+nu2))^(lmax+1/2) # norm-overlap
element_T1D(lmax,nu1,nu2,hbar,mur) = -hbar^2/2/mur * element_S1D(lmax,nu1,nu2)*(-2*(1+2*lmax)*nu1+4*nu1^2*(lmax+1/2)/(nu1+nu2) + lmax*(lmax-1)*2*(nu1+nu2)/(2*lmax-1))# kinetic energy
element_VGauss1D(lmax,nu1,nu2,v0,mu_g) = v0*(2*(nu1*nu2)^(1/2)/(nu1+nu2+mu_g))^(lmax+1/2) # gaussian interaction
# single shifted gaussian interaction:
function element_VShiftedGauss1D(lmax,nu1,nu2,v0,mu_g,z0)
    if lmax == 0
        return v0 * sqrt(2)*(nu1*nu2)^(1/4) * exp(-mu_g*(nu1+nu2)*z0^2/(mu_g+nu1+nu2))/(mu_g+nu1+nu2)^(1/2)
    else
        return v0 * (2(nu1*nu2)^(1/2)/(mu_g+nu1+nu2))^(lmax+1/2)*exp(-mu_g*z0^2) * GSL.sf_hyperg_1F1(0.5+lmax,0.5,mu_g^2*z0^2/(mu_g+nu1+nu2))
    end
end
# central interaction via numerical integration:
norm1D(nu,lmax,gamma_dict) = ((2*nu)^(lmax+1/2)/gamma_dict[lmax+1/2])^(1/2);
integrand1D(r,lmax,nu1,nu2,vint) = r^(2*lmax)*exp(-(nu1+nu2)*r^2)*vint(r)
function element_V1D(lmax,nu1,nu2,vint,gamma_dict,buf)
    return quadgk(r -> integrand1D(r,lmax,nu1,nu2,vint),-Inf,0,Inf;segbuf=buf)[1]*norm1D(nu1,lmax,gamma_dict)*norm1D(nu2,lmax,gamma_dict)
end
vint_csm(r,vint,theta_csm) = vint(r*exp(im*theta_csm*pi/180))

## functions for the full (lower-triangular) matrices. requires already preallocated matrices S,T,V
## S: 
function MatrixS1D(S,lmax,nu_arr)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)
            S[nrow,ncol] = element_S1D(lmax,nu_arr[nrow]',nu_arr[ncol])
        end
    end
end

## T: csm via global factor
function MatrixT1D(T,lmax,nu_arr,hbar,mur,csm_bool,theta_csm)
    MatrixT1D_cr(T,lmax,nu_arr,hbar,mur)
    csm_bool == 1 && (T .*= exp(-2*im*theta_csm*pi/180))
end

function MatrixT1D_cr(T,lmax,nu_arr,hbar,mur)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)
            T[nrow,ncol] = element_T1D(lmax,nu_arr[nrow]',nu_arr[ncol],hbar,mur)
        end
    end
end

## V:
function MatrixV1D(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt,shiftgopt,csm_bool,theta_csm)
    if csm_bool == 0
        MatrixV1D_cr(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt,shiftgopt)
    elseif csm_bool == 1
        MatrixV1D_csm(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt,shiftgopt,theta_csm)
    end
end

# no csm
function MatrixV1D_cr(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt,shiftgopt)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)  
            if gaussopt[1] == 1
                v0,mu_g = gaussopt[2:3]
                V[nrow,ncol] = element_VGauss1D(lmax,nu_arr[nrow]',nu_arr[ncol],v0,mu_g)
            elseif shiftgopt[1] == 1
                for ii in keys(shiftgopt[2])
                    v0,mu_g,z0 = shiftgopt[2][ii],shiftgopt[3][ii],shiftgopt[4][ii]
                    V[nrow,ncol] += element_VShiftedGauss1D(lmax,nu_arr[nrow]',nu_arr[ncol],v0,mu_g,z0)
                end
            else gaussopt[1] == 0
                V[nrow,ncol] = element_V1D(lmax,nu_arr[nrow]',nu_arr[ncol],vint,gamma_dict,buf)
            end
        end
    end
end

#csm
function MatrixV1D_csm(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt,shiftgopt,theta_csm)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)  
            if gaussopt[1] == 1
                v0,mu_g = gaussopt[2:3]
                V[nrow,ncol] = element_VGauss1D(lmax,nu_arr[nrow]',nu_arr[ncol],v0,mu_g*exp(2*im*theta_csm*pi/180))
            elseif shiftgopt[1] == 1
                for ii in keys(shiftgopt[2])
                    v0,mu_g,z0 = shiftgopt[2][ii],shiftgopt[3][ii],shiftgopt[4][ii]
                    V[nrow,ncol] += element_VShiftedGauss1D(lmax,nu_arr[nrow]',nu_arr[ncol],v0,mu_g*exp(2*im*theta_csm*pi/180),z0)
                end
            else gaussopt[1] == 0
                V[nrow,ncol] = element_V1D(lmax,nu_arr[nrow]',nu_arr[ncol],r->vint_csm(r,vint,theta_csm),gamma_dict,buf)
            end
        end
    end
end