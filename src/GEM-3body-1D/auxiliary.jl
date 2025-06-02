# file for auxiliary functions

function make_phys_params3B1D(;hbar = 1.0, mass_arr, svals=["x","y","z"], vint_arr=[[]], parity=+1)
    return (;hbar, mass_arr, svals, vint_arr, parity)
end

function make_num_params3B1D(;lmin=0,Lmin=0, lmax=0, Lmax=0, gem_params=(nmax=5, r1=1.0, rnmax=10.0, Nmax=5, R1=1.0, RNmax=10.0), theta_csm=0.0, omega_cr=0.9, mu0=0.08, c_shoulder=1.6, kmax_interpol=1000, threshold=10^-10)
    return (;lmin, Lmin, lmax, Lmax, gem_params, theta_csm, omega_cr, mu0, c_shoulder, kmax_interpol, threshold)
end

## for gaussopt and csm_bool
function csmgaussopt(gaussopt,csm_bool,num_params)
    sumgb = 0
    for cc in 1:lastindex(gaussopt)
        for vi in 1:lastindex(gaussopt[cc])
            sumgb += abs(gaussopt[cc][vi][1])
        end
    end
    
    if csm_bool == 1 && sumgb != 0
        gaussoptc = gaussopt .*(1.0+0.0*im)
        for cc in 1:lastindex(gaussopt)
            for vi in 1:lastindex(gaussopt[cc])
                gaussoptc[cc][vi][3] *= exp(2*im*num_params.theta_csm*pi/180)
            end
        end
        gaussopt=gaussoptc
    end
    return gaussopt
end