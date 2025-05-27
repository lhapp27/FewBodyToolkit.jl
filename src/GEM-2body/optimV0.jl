# for finding a the value of v0 to achieve the energy target_e2 for the state defined by stateindex (with intermediate optimization of parameters)

"""
    v0GEMOptim(phys_params, num_params, stateindex, target_e2; cr_bool=0, rtol=1e-4, atol=10*eps(), g_tol=1e-9, output=false)

Finds a value `v0crit` to globally scale the potential in order to achieve the target energy `target_e2` for the state specified by `stateindex`. Additionally, performs intermediate optimization of Gaussian Expansion Method (GEM) parameters.

# Arguments
- `phys_params::NamedTuple`: Physical parameters including `hbar`, `mur`, `vint_arr`, `lmax`, `lmin`, and `dim`.
- `num_params::NamedTuple`: Numerical parameters including `gem_params`, `omega_cr`, `theta_csm`, and `threshold`.
- `stateindex::Int`: Index of the state to optimize, where `1` indicates the ground state.
- `target_e2::Float64`: Target energy value.

# Keyword Arguments
- `cr_bool::Int`: Use complex rotation (default: 0).
- `rtol::Float64`: Relative tolerance for energy convergence (default: 1e-4).
- `atol::Float64`: Absolute tolerance for energy convergence (default: 10*eps()).
- `g_tol::Float64`: Tolerance for the optimizer (default: 1e-9).
- `output::Bool`: If true, prints intermediate optimization results (default: false).

# Returns
- `phys_params::NamedTuple`: Updated physical parameters with optimized potential.
- `num_params::NamedTuple`: Updated numerical parameters with optimized GEM ranges.
- `v0crit::Float64`: The value to scale the potential with, to achieve the target energy. This is **not** the overall value of `v0` but rather a scaling factor for the potential.

# Examples
```julia
phys_params_scaled,num_params_optimized,scalingfactor = v0GEMOptim(phys_params, num_params, 1, -2.5)
"""


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

"""
    GEM_Optim_2B(phys_params, num_params, stateindex; cr_bool=0, g_tol=1e-9)

Optimize the ranges used in the Gaussian Expansion Method (GEM) for a specific 2-body state.

# Arguments
- `phys_params::NamedTuple`: Physical parameters such as `hbar`, `mur`, `vint_arr`, `lmax`, `lmin`, and `dim`.
- `num_params::NamedTuple`: Numerical parameters such as `gem_params`, `omega_cr`, `theta_csm`, and `threshold`.
- `stateindex::Int`: An integer specifying the index of the state to optimize (e.g., `1` for the ground state).

# Keyword Arguments
- `cr_bool::Int=0`: Indicates whether to use complex rotation.
- `g_tol::Float64=1e-9`: A value specifying the tolerance for the optimizer.

# Returns
- A vector containing:
  - Optimized GEM parameters: `[r1, rnmax]`.
  - The energy value of the optimized state.

# Examples
```julia
r1opt,rnmaxopt,energy = GEM_Optim_2B(phys_params, num_params, 1) # optimize for the ground state
r1opt,rnmaxopt,energy = GEM_Optim_2B(phys_params, num_params, 3; cr_bool = 1) # optimize for the third state (2nd excited) using complex-ranged basis functions
"""

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
