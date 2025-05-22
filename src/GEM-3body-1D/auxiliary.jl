# file for auxiliary functions

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