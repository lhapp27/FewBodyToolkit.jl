# functions to preallocate all arrays for the GEM3B1D program

# contains:
# - preallocate_data: for everything

# TR is the element type of the Gaussian ranges, and of everything derived from the ranges alone
# (norms, the overlap matrix S). It is Float64 for the usual real ranges and ComplexF64 as soon as
# complex-ranged basis functions are used in at least one Jacobi coordinate.
struct PrecomputeStruct{TR}
    gamma_dict::Dict{Float64, Float64}
    jmat::Matrix{SMatrix{2, 2, Float64, 4}}
    murR_arr::MMatrix{2, 3, Float64, 6}
    nu_arr::Vector{TR}
    NU_arr::Vector{TR}
    norm_arr::OffsetMatrix{TR, Matrix{TR}}
    NORM_arr::OffsetMatrix{TR, Matrix{TR}}
end

# T  : element type of the tabulated radial integrals
# TA : element type of the interpolation mesh of effective ranges (complex only for complex ranges,
#      where the mesh covers a sector of the complex plane rather than a real interval)
# IT : type of the interpolants. One-dimensional in log(alpha) for real ranges, two-dimensional over
#      (log|alpha|, arg alpha) for complex ones. Inferred rather than spelled out, since the
#      dimensionality is part of the type.
struct InterpolationStruct{T,TA,IT}
    alpha_arr::Vector{Float64}
    theta_arr::Vector{Float64}
    alpha_grid::Matrix{TA}
    v_arr::Array{T, 3}
    w_interpol_arr::OffsetArray{IT, 3, Array{IT, 3}}
    #v_obs_arr, w_obs_arr, w_obs_interpol_arr: observables are currently not supported in 1D
end

# same TR as above; the arguments of a matrix element are typed by the ranges they carry
struct TempArgs1D{TR}
    rowi::Int64
    coli::Int64
    ranges::NamedTuple{(:nua, :nub, :NUa, :NUb), NTuple{4, TR}}
    norm4::TR
    la::Int64
    La::Int64
    lb::Int64
    Lb::Int64
    Lsum::Int64
    avals_new::Vector{Int64}
    bvals_new::Vector{Int64}
    cvals::Vector{Int64}
    factor_ab::Float64
    avals::Vector{Int64}
    bvals::Vector{Int64}
end

# Outer constructor: takes TR from the ranges and lets the inner constructor convert the remaining
# fields (in particular the integer factor_ab). The auto-generated constructor of a *parametric*
# struct matches the declared field types exactly and would not perform that conversion.
function TempArgs1D(rowi,coli,ranges::NamedTuple{(:nua, :nub, :NUa, :NUb),NTuple{4,TR}},norm4,la,La,lb,Lb,Lsum,avals_new,bvals_new,cvals,factor_ab,avals,bvals) where {TR}
    TempArgs1D{TR}(rowi,coli,ranges,norm4,la,La,lb,Lb,Lsum,avals_new,bvals_new,cvals,factor_ab,avals,bvals)
end


# T  : element type of H, V and everything that can pick up a complex-scaling phase
#      (complex as soon as complex_scaling OR complex ranges are active)
# TR : element type of the ranges and of the overlap S (complex only for complex ranges)
struct FillStruct{T,TR}
    wn_interpol_arr::OffsetVector{T, Vector{T}}
    T::Matrix{T}
    V::Matrix{T}
    S::Matrix{TR}
    temp_args_arr::Vector{TempArgs1D{TR}}
    temp_fill_mat::Matrix{T}
    #wn_obs_interpol_arr::OffsetVector{Float64, Vector{Float64}}
end

# TE: element type of the energies. Complex ranges alone leave H and S hermitian, so the
# eigenvalues stay real; only complex scaling makes them genuinely complex.
struct ResultStruct{TE,TW}
    energies_arr::Vector{TE}
    wavefun_arr::Matrix{TW}
    #centobs_output::Array{Float64}
    #R2_output::Matrix{Float64}
end


function preallocate_data(phys_params,num_params,observ_params,size_params,complex_scaling::Bool,complex_ranged_r::Bool=false,complex_ranged_R::Bool=false)
    complex_ranged = complex_ranged_r || complex_ranged_R

    # TT: H, V, and everything that can carry the complex-scaling phase
    # TR: ranges, norms, overlap S (the complex-scaling rotates the potential, not the ranges)
    # TE: energies. With complex ranges alone, H and S stay hermitian -> real eigenvalues.
    TT = (complex_scaling || complex_ranged) ? ComplexF64 : Float64
    TR = complex_ranged ? ComplexF64 : Float64
    TE = complex_scaling ? ComplexF64 : Float64
    
    #Destructing Structs:
    (;gem_params,kmax_interpol,kmax_theta,complex_range_freq) = num_params
    (;nmax,Nmax) = gem_params
    (;nbasis_total,nlL,nl,maxlmax,maxobs,nintmax) = size_params
    (;stateindices,centobs_arr,R2_arr) = observ_params

    # complex-ranged basis functions come in conjugate pairs, so the corresponding coordinate
    # carries twice as many ranges. The two coordinates are independent.
    nmax_eff = complex_ranged_r ? 2*nmax : nmax
    Nmax_eff = complex_ranged_R ? 2*Nmax : Nmax
    
    # Allocate factorials (gamma function)
    gamma_dict = Dict{Float64, Float64}()
    
    # Allocate Jacobi transofmration matrices jmat and reduced masses murR_arr
    jmat = Matrix{SMatrix{2, 2, Float64, 4}}(undef, 3, 3);
    murR_arr = MMatrix{2,3,Float64,6}(zeros(2,3));
    
    # Allocate ranges
    nu_arr = Vector{TR}(undef, nmax_eff)
    NU_arr = Vector{TR}(undef, Nmax_eff)        
    
    # Allocate norms
    norm_arr = OffsetMatrix{TR}(undef, 0:maxlmax, nmax_eff)*zero(TR) # one norm_arr for all l,n combinations. changed to OffsetArray for all values of l=0:maxlmax
    NORM_arr = OffsetMatrix{TR}(undef, 0:maxlmax, Nmax_eff)*zero(TR) # one norm_arr for all L,N combinations
    
    # arrays for range-interpolation method:
    # With complex ranges the interpolation mesh is a sector of the complex alpha-plane: kmax_interpol
    # values of |alpha| times ntheta angular nodes (theta_mesh). Without them ntheta = 1 and the mesh
    # is the plain real one, so everything below reduces to the previous, one-dimensional arrays.
    alpha_arr = zeros(kmax_interpol)
    theta_arr = theta_mesh(complex_ranged,complex_scaling,complex_range_freq,kmax_theta)
    ntheta = lastindex(theta_arr)
    TA = complex_ranged ? ComplexF64 : Float64 # the mesh itself stays real without complex ranges
    alpha_grid = zeros(TA,kmax_interpol,ntheta)
    v_arr = zeros(TT,kmax_interpol,ntheta,2*2*maxlmax+1) # careful, now NOT an OffsetArray! due to problems in A\v # changed to 4*maxlmax+1. 

    x=range(0.1,0.5,3)
    y=x.^2;interpoltype = typeof(cubic_spline_interpolation(x,y))# just for easy inferring the type of interpolation objects # seems quite hacked, but works for now
    yc=x.^2 .+ zero(TT);
    # the interpolants gain a second (angular) dimension for complex ranges; infer the type from a
    # dummy of the right dimensionality rather than spelling it out
    interpoltypeC = complex_ranged ? typeof(cubic_spline_interpolation((x,x),yc.*transpose(y))) : typeof(cubic_spline_interpolation(x,yc))

    w_interpol_arr=OffsetArray{interpoltypeC}(undef,3,nintmax,0:4*maxlmax) # penultimate dimensionality for different Lsum values! removed that dimension for 1D
    wn_interpol_arr = OffsetArray{TT}(undef,0:4*maxlmax)*zero(TT) # for the interpolated wn_values that are actually used

    #for observables (in range-interpolation) # currently not available in 1D
#=     v_obs_arr = zeros(kmax_interpol,2*maxlmax+1)
    w_obs_arr = zeros(3,maxobs,kmax_interpol,2*maxlmax+1,2*maxlmax+1)
    w_obs_interpol_arr = OffsetArray{interpoltype}(undef,3,maxobs,0:2*maxlmax,0:2*maxlmax)
    wn_obs_interpol_arr = OffsetArray{Float64}(undef,0:2*maxlmax)*0.0 =#
    
    # matrices S: for the norm overlap, T: for the kintetic energy (later for the Hamiltonian T+V), V: for the potential energy
    S = zeros(TR,nbasis_total,nbasis_total) # hermitian (not just symmetric) for complex ranges
    T = zeros(TT,nbasis_total,nbasis_total)
    V = zeros(TT,nbasis_total,nbasis_total);
    
    
    # for results:
    energies_arr = zeros(TE,nbasis_total);# Vector{TT}(undef, nbasis_total) # changed to zeros to avoid bad behavior due to undef values not occupied thresholding (eigen2step)
    wavefun_arr = Matrix{TT}(undef, nbasis_total, nbasis_total)
    #centobs_output = Array{Float64}(undef, 3, maxobs, lastindex(stateindices))
    #R2_output = zeros(3, lastindex(stateindices))
    
    # temporary arrays for filling: function arguments and matrix
    # Combined complex ranges and complex scaling leave H neither symmetric nor hermitian,
    # so the full matrix has to be filled instead of only the lower triangle.
    fill_full = complex_ranged && complex_scaling
    ntot = fill_full ? Int64(nbasis_total^2) : Int64((nbasis_total^2+nbasis_total)/2)
    temp_args_arr = Vector{TempArgs1D{TR}}(undef, ntot)
    temp_fill_mat = Matrix{TT}(undef,nbasis_total, nbasis_total)
        
    
    # now constructing Structs for different steps in the program:
    precomp_arrs = PrecomputeStruct(gamma_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr)
    #temp_arrs = TempStruct(temp_clmk,temp_dlmk,temp_S,temp_D1,temp_D2)
    interpol_arrs = InterpolationStruct(alpha_arr,theta_arr,alpha_grid,v_arr,w_interpol_arr) # ,v_obs_arr,w_obs_arr,w_obs_interpol_arr
    fill_arrs = FillStruct(wn_interpol_arr,T,V,S,temp_args_arr,temp_fill_mat) # ,wn_obs_interpol_arr
    result_arrs = ResultStruct(energies_arr,wavefun_arr) #,centobs_output,R2_output)
    
    return precomp_arrs,interpol_arrs,fill_arrs,result_arrs

    ## again, no huge differences between 1D and 3D.
end