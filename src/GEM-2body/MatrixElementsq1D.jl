## 1D
# Separate file for the matrix elements. Includes convenience functions for the full matrices T,V,S (only lower-triangular)

## functions to calculate single matrix elements (normalization prefactor is already included):
element_S1D(lmax,nu1,nu2) = (2*(nu1*nu2)^(1/2)/(nu1+nu2))^(lmax+1/2) # norm-overlap
element_T1D(lmax,nu1,nu2,hbar,mur) = -hbar^2/2/mur * element_S1D(lmax,nu1,nu2)*(-2*(1+2*lmax)*nu1+4*nu1^2*(lmax+1/2)/(nu1+nu2) + lmax*(lmax-1)*2*(nu1+nu2)/(2*lmax-1))# kinetic energy
element_VGauss1D(lmax,nu1,nu2,v0,mu_g) = v0*(2*(nu1*nu2)^(1/2)/(nu1+nu2+mu_g))^(lmax+1/2) # gaussian interaction
# central interaction via numerical integration:
norm1D(nu,lmax,gamma_dict) = ((2*nu)^(lmax+1/2)/gamma_dict[lmax+1/2])^(1/2);
integrand1D(r,lmax,nu1,nu2,vint) = r^(2*lmax)*exp(-(nu1+nu2)*r^2)*vint(r)
function element_V1D(lmax,nu1,nu2,vint,gamma_dict,buf)
    return quadgk(r -> integrand1D(r,lmax,nu1,nu2,vint),-Inf,Inf;segbuf=buf)[1]*norm1D(nu1,lmax,gamma_dict)*norm1D(nu2,lmax,gamma_dict)
end


## functions for the full (lower-triangular) matrices. requires already preallocated matrices S,T,V
function MatrixSq1D(S,lmax,nu_arr)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)
            S[nrow,ncol] = element_S1D(lmax,nu_arr[nrow]',nu_arr[ncol])
        end
    end
end

function MatrixTq1D(T,lmax,nu_arr,hbar,mur)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)
            T[nrow,ncol] = element_T1D(lmax,nu_arr[nrow]',nu_arr[ncol],hbar,mur)
        end
    end
end

function MatrixVq1D(V,lmax,nu_arr,vint,gamma_dict,buf,gaussopt)
    for ncol = 1:lastindex(nu_arr)
        for nrow = ncol:lastindex(nu_arr)  
            if gaussopt[1] == 0
                V[nrow,ncol] = element_V1D(lmax,nu_arr[nrow]',nu_arr[ncol],vint,gamma_dict,buf)
            elseif gaussopt[1] == 1
                v0,mu_g = gaussopt[2:3]
                V[nrow,ncol] = element_VGauss1D(lmax,nu_arr[nrow]',nu_arr[ncol],v0,mu_g)
            end
        end
    end
end