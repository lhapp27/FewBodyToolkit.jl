# Function for solving the generalized eigenvalue problem Hv=ESv

function solveHS(num_params,fill_arrs,result_arrs,return_wavefunctions::Bool,hermitian::Bool=false)
    (;threshold) = num_params
    if return_wavefunctions
        eigen2step_valvec(result_arrs.energies_arr,result_arrs.wavefun_arr,fill_arrs.T,fill_arrs.S;threshold=threshold,hermitian=hermitian)
    elseif !return_wavefunctions
        eigen2step!(result_arrs.energies_arr,fill_arrs.T,fill_arrs.S;threshold=threshold,hermitian=hermitian)
    else
        println("Error in return_wavefunctions: no output produced")
    end
end