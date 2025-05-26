# for finding a the value of v0 to achieve the energy target_e2 for the state defined by stateindex (with intermediate optimization of parameters)
function v0GEMOptim(phys_params, num_params, stateindex, target_e2; cr_bool = 0, rtol=10^-4, atol = 10*eps(), g_tol=10^-9,output=false)
    
    # desctructuring:
    (;hbar,mur,vint_arr,lmax,lmin,dim) = phys_params
    vint = vint_arr[1]
    (;gem_params,omega_cr,theta_csm, threshold) = num_params
    (;nmax,r1,rnmax) = gem_params
    
    v0crit = 20.0; # (arbitrary) starting value, should be deep enough
        
    params = zeros(3).+1;
    iter = 0
    while abs(params[3] - target_e2) > max(abs(rtol*target_e2),abs(atol)) #changed to relative tolerance. previously: abs(params[3] - target_e2) > 10^-8
        v0crit = find_v0crit((;hbar,mur,vint_arr,lmax,lmin,dim), num_params, stateindex, target_e2; cr_bool=cr_bool) # wait where is vint_arr updated?! this should be phys_params, no? ah ok no. only for GEMOptim we need to update phys_params!
        function vint_updt(r)
            return v0crit*vint(r)
        end
        phys_params = (;hbar,mur,vint_arr=[vint_updt],lmax,lmin,dim)

        params[1:3] = GEM_Optim_2B(phys_params, num_params, stateindex; cr_bool=cr_bool, g_tol=g_tol)
        num_params = (;gem_params = (nmax, r1 = params[1], rnmax = params[2]), omega_cr, theta_csm, threshold)

        vscale=v0crit
        iter += 1;
        if output == true
            @show([vscale,params])
            println("Number of iterations: ", iter)
        end
    end
    
    return phys_params,num_params,v0crit
end

# for optimizing GEM ranges of specific state defined by stateindex (1 = ground state)
function GEM_Optim_2B(phys_params, num_params, stateindex; cr_bool=0, g_tol=10^-9)
    
    (;gem_params,omega_cr,theta_csm,threshold) = num_params
    (;nmax,r1,rnmax) = gem_params
    
    target(x) = GEM2B.GEM_solve(phys_params, (gem_params = (nmax, r1 = x[1], rnmax = x[2]), omega_cr, theta_csm, threshold); cr_bool)[stateindex]
    
    init = [r1,rnmax] # arbitrary! Change in case the optimization fails
    opt = optimize(target, init, method=NelderMead(), show_trace = false, store_trace=false, extended_trace=false, g_tol=g_tol);
    return [abs.(Optim.minimizer(opt)') Optim.minimum(opt)] # 2 parameters, 1 energy value
end

# for finding the value of v0 to achieve the energy target_e2 for the state defined by stateindex (without intermediate optimization of parameters)
function find_v0crit(phys_params, num_params, stateindex, target_e2; cr_bool=0)
    (;hbar,mur,vint_arr,lmax,lmin,dim) = phys_params
    vint = vint_arr[1]
    
    function fun1(v0)
        function vint_local(r)
            return v0*vint(r)
        end
        
        return real(GEM2B.GEM_solve((hbar,mur,vint_arr=[vint_local],lmax,min,dim), num_params; cr_bool)[stateindex]) - target_e2;
    end
    
    v0crit = find_zeros(fun1,(0, 200))[1];
    return v0crit
end
