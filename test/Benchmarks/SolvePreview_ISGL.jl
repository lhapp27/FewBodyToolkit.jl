using FewBodyToolkit

function run_single_solve_preview(
    pp,
    np;
    obs_params,
    complex_scaling::Bool = false,
    return_wavefunctions::Bool = false,
    debug::Bool = false,
    npreview::Int = 2,
)
    IS = FewBodyToolkit.ISGL
    size_params = IS.size_estimate(pp, np, obs_params, complex_scaling)
    precomp_arrs, temp_arrs, interpol_arrs, fill_arrs, result_arrs = IS.preallocate_data(pp, np, obs_params, size_params, complex_scaling)
    IS.precompute_ISGL(pp, np, size_params, precomp_arrs, temp_arrs)
    IS.interpolNshoulder(pp, np, obs_params, size_params, precomp_arrs, interpol_arrs, return_wavefunctions, complex_scaling)
    IS.fill_TVS(np, size_params, precomp_arrs, interpol_arrs, fill_arrs, complex_scaling, pp.hbar, debug)
    IS.solveHS(np, fill_arrs, result_arrs, return_wavefunctions)

    nener = min(npreview, length(result_arrs.energies_arr))
    return copy(result_arrs.energies_arr[1:nener])
end
