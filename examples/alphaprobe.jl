using Pkg
Pkg.instantiate();
Pkg.activate("D:/Julia/FewBodyToolkit.jl/");

using FewBodyToolkit

pp = make_phys_params3B1D();
np = make_num_params3B1D();
op=(;stateindices=[],centobs_arr=[[],[],[]],R2_arr=[0,0,0]);

complex_scaling = false;

cr_r = false; cr_R = false;
# size params
size_params = FewBodyToolkit.GEM3B1D.size_estimate(pp,np,op,cr_r,cr_R);
#preallocate
precomp_arrs,interpol_arrs,fill_arrs,result_arrs = FewBodyToolkit.GEM3B1D.preallocate_data(pp,np,op,size_params,cr_r,cr_R,complex_scaling)
# precompute
FewBodyToolkit.GEM3B1D.precompute_3B1D(pp,np,size_params,precomp_arrs,cr_r,cr_R)

FewBodyToolkit.GEM3B1D.precompute_alpha_arr(interpol_arrs.alpha_arr,np.gem_params.r1,np.gem_params.rnmax,np.gem_params.R1,np.gem_params.RNmax,precomp_arrs.nu_arr,precomp_arrs.NU_arr,precomp_arrs.jmat)


## jetzt komplex:
cr_r = true; cr_R = true;
# size params
size_params_cr = FewBodyToolkit.GEM3B1D.size_estimate(pp,np,op,cr_r,cr_R);
#preallocate
precomp_arrs_cr,interpol_arrs_cr,fill_arrs_cr,result_arrs_cr = FewBodyToolkit.GEM3B1D.preallocate_data(pp,np,op,size_params_cr,cr_r,cr_R,complex_scaling)
# precompute
FewBodyToolkit.GEM3B1D.precompute_3B1D(pp,np,size_params_cr,precomp_arrs_cr,cr_r,cr_R)

FewBodyToolkit.GEM3B1D.precompute_alpha_arr(interpol_arrs_cr.alpha_arr,np.gem_params.r1,np.gem_params.rnmax,np.gem_params.R1,np.gem_params.RNmax,precomp_arrs_cr.nu_arr,precomp_arrs_cr.NU_arr,precomp_arrs_cr.jmat)



# okay hier wirds wild, weil alpha_arr und alpha_arr_cr gleich sind!!! --> wie kann das trotzdem funktionieren, und warum gibts keinen error wegen extrapolation in fillTVS??