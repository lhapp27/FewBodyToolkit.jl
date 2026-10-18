# Function for solving the generalized eigenvalue problem Hv=ESv

function solveHS(num_params,fill_arrs,result_arrs,return_wavefunctions::Bool,partner=nothing)
    (;threshold) = num_params
    H,S = fill_arrs.T,fill_arrs.S
    # complex ranges without complex scaling: solve in the real basis (see realbasis in eigen2step.jl), unless
    # the conjugate pairing does not hold
    if !isnothing(partner)
        if max(conj_pair_dev(S,partner),conj_pair_dev(H,partner)) < 1e-10
            H,S = realbasis(H,partner),realbasis(S,partner)
        else
            @warn "complex ranges: basis functions not conjugate-paired, the eigenproblem is solved in complex arithmetic"
            partner = nothing
        end
    end
    if return_wavefunctions
        eigen2step_valvec(result_arrs.energies_arr,result_arrs.wavefun_arr,H,S;threshold=threshold)
        isnothing(partner) || realbasis_back!(result_arrs.wavefun_arr,partner)
    elseif !return_wavefunctions
        eigen2step!(result_arrs.energies_arr,H,S;threshold=threshold)
    else
        println("Error in return_wavefunctions: no output produced")
    end
end
