# functions to preallocate all arrays for the GEM3B1D program

# contains:
# - preallocate_data: for everything

struct PrecomputeStruct
    gamma_dict::Dict{Float64, Float64}
    jmat::Matrix{SMatrix{2, 2, Float64, 4}}
    murR_arr::MMatrix{2, 3, Float64, 6}
    nu_arr::Vector{Float64}
    NU_arr::Vector{Float64}
    norm_arr::OffsetMatrix{Float64, Matrix{Float64}}
    NORM_arr::OffsetMatrix{Float64, Matrix{Float64}}
end

struct InterpolationStruct{T}
    alpha_arr::Vector{Float64}
    v_arr::Matrix{T}
    A_mat::Matrix{Float64}
    w_arr::Array{T, 4}
    w_interpol_arr::OffsetArray{Interpolations.Extrapolation{T, 1, ScaledInterpolation{T, 1, Interpolations.BSplineInterpolation{T, 1, OffsetVector{T, Vector{T}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}, 3, Array{Interpolations.Extrapolation{T, 1, ScaledInterpolation{T, 1, Interpolations.BSplineInterpolation{T, 1, OffsetVector{T, Vector{T}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}, 3}} # changed 3s to 2s for 1D! back to 3 isntead of 2! #OffsetArray{Interpolations.Extrapolation{Float64, 1, ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, OffsetVector{Float64, Vector{Float64}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}}
    Ainv_arr_kine::OffsetMatrix{Float64, Matrix{Float64}}
    #v_obs_arr::Matrix{Float64}
    #w_obs_arr::Array{Float64, 5}
    #w_obs_interpol_arr::OffsetArray{Interpolations.Extrapolation{Float64, 1, ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, OffsetVector{Float64, Vector{Float64}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}, 4, Array{Interpolations.Extrapolation{Float64, 1, ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, OffsetVector{Float64, Vector{Float64}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}, 4}}
end

struct FillStruct{T}
    fij_arr::Matrix{Float64}
    Kij_arr::Matrix{Float64}
    w_arr_kine::OffsetVector{Float64, Vector{Float64}}
    wn_interpol_arr::OffsetVector{T, Vector{T}}
    kij_arr::OffsetArray{Float64, 3, Array{Float64, 3}}#OffsetArray{Float64}
    gij_arr::OffsetArray{Float64, 3, Array{Float64, 3}}#OffsetArray{Float64}
    T::Matrix{T}
    V::Matrix{T}
    S::Matrix{Float64}
    temp_args_arr::Vector{NamedTuple{(:rowi, :coli, :ranges, :norm4, :la, :La, :lb, :Lb, :Lsum, :avals_new, :bvals_new, :cvals, :factor_ab, :avals, :bvals), Tuple{Int64, Int64, NamedTuple{(:nua, :nub, :NUa, :NUb), NTuple{4, Float64}}, Float64, Int64, Int64, Int64, Int64, Int64, Vector{Int64}, Vector{Int64}, Vector{Int64}, Float64, Vector{Int64}, Vector{Int64}}}} # prob needs fix
    temp_fill_mat::Matrix{T}
    #wn_obs_interpol_arr::OffsetVector{Float64, Vector{Float64}}
end

struct ResultStruct{T}
    energies_arr::Vector{T}
    wavefun_arr::Matrix{T}
    centobs_output::Array{Float64}
    R2_output::Matrix{Float64}
end


function preallocate_data(phys_params,num_params,observ_params,size_params,csm_bool)
    if csm_bool == 0
        TT = Float64 # type variable
    elseif csm_bool == 1
        TT = ComplexF64
    end
    
    #Destructing Structs:
    (;gem_params,kmax_interpol) = num_params
    (;nmax,Nmax) = gem_params
    (;nbasis_total,nlL,nl,maxlmax,maxobs,nintmax) = size_params
    (;stateindices,centobs_arr,R2_arr) = observ_params
    
    # Allocate factorials (gamma function)
    gamma_dict = Dict{Float64, Float64}()
    
    # Allocate jmat and murR_arr
    jmat = Matrix{SMatrix{2, 2, Float64, 4}}(undef, 3, 3);
    murR_arr = MMatrix{2,3,Float64,6}(zeros(2,3)); # reduced masses
    
    # Allocate ranges
    nu_arr = Vector{Float64}(undef, nmax)
    NU_arr = Vector{Float64}(undef, Nmax)        
    
    # Allocate norms
    norm_arr = OffsetMatrix{Float64}(undef, 0:maxlmax, nmax)*0.0 # one norm_arr for all l,n combinations. changed to OffsetArray for all values of l=0:maxlmax
    NORM_arr = OffsetMatrix{Float64}(undef, 0:maxlmax, Nmax)*0.0 # one norm_arr for all L,N combinations
    
    # arrays for range-interpolation method:
    alpha_arr = zeros(kmax_interpol)
    v_arr = zeros(TT,kmax_interpol,2*2*maxlmax+1) # careful, now NOT a OffsetArray! due to problems in A\v # changed to 4*maxlmax+1. 
    A_mat = zeros(2*maxlmax+1,2*maxlmax+1)
    w_arr = zeros(TT,3,kmax_interpol,4*maxlmax+1,4*maxlmax+1) # first dimensionality is 3 for covering interactions for all possible Jacobi-Sets (=3). in principle could be reduced to nboxes? penultimate dimensionality for different Lsum values! #dwnti 1D?

    x=range(0.1,0.5,3)
    y=x.^2;interpoltype = typeof(cubic_spline_interpolation(x,y))# just for easy inferring the type of interpolation objects # seems quite hacked, but works for now
    yc=x.^2 .+ zero(TT);interpoltypeC = typeof(cubic_spline_interpolation(x,yc))# just for easy inferring the type of interpolation objects # for possibly complex arguments!

    w_interpol_arr=OffsetArray{interpoltypeC}(undef,3,nintmax,0:4*maxlmax) # penultimate dimensionality for different Lsum values! removed that dimension for 1D diff13 dwnti 1D?
    wn_interpol_arr = OffsetArray{TT}(undef,0:4*maxlmax)*0.0 # for the interpolated wn_values that are actually used diff13 dwnti 1D?

    Ainv_arr_kine = OffsetArray{Float64}(undef,0:2*maxlmax,0:2*maxlmax)*0.0 # penultimate dimensionality for different Lsum values! # dwnti 1D?
    w_arr_kine = OffsetArray{Float64}(undef,0:2*maxlmax+1)*0.0 # dwnti 1D?
    #for observables (in range-interpolation) # currently not available in 1D
#=     v_obs_arr = zeros(kmax_interpol,2*maxlmax+1)
    w_obs_arr = zeros(3,maxobs,kmax_interpol,2*maxlmax+1,2*maxlmax+1)
    w_obs_interpol_arr = OffsetArray{interpoltype}(undef,3,maxobs,0:2*maxlmax,0:2*maxlmax)
    wn_obs_interpol_arr = OffsetArray{Float64}(undef,0:2*maxlmax)*0.0 =#
    
    # arrays for shoulder method (used in matrix-element calculations) # dwnti 1D?
    fij_arr = zeros(3,4)
    Kij_arr = zeros(3,4)
    kij_arr = OffsetArray{Float64}(undef,3,4,0:2*maxlmax)*0.0
    gij_arr = OffsetArray{Float64}(undef,3,4,0:2*maxlmax)*0.0
    S = zeros(nbasis_total,nbasis_total) #Matrix{Float64}(undef, nbasis_total, nbasis_total)
    T = zeros(TT,nbasis_total,nbasis_total) #Matrix{Float64}(undef, nbasis_total, nbasis_total)
    V = zeros(TT,nbasis_total,nbasis_total);#Matrix{Float64}(undef, nlL * nmax * Nmax, nlL * nmax * Nmax)
    
    
    # for results:
    energies_arr = zeros(TT,nbasis_total);# Vector{TT}(undef, nbasis_total) # changed to zeros to avoid bad behavior due to undef values not occupied thresholding (eigen2step)
    wavefun_arr = Matrix{TT}(undef, nbasis_total, nbasis_total)
    centobs_output = Array{Float64}(undef, 3, maxobs, lastindex(stateindices))
    R2_output = zeros(3, lastindex(stateindices))
    
    # temporary arrays for filling: function arguments and matrix for parallelization
    ntot = Int64((nbasis_total^2+nbasis_total)/2) # only for lower triangular!
    # Define the type of the tuple for the function arguments
    tuple_type = NamedTuple{(:rowi, :coli, :ranges, :norm4, :la, :La, :lb, :Lb, :Lsum, :avals_new, :bvals_new, :cvals, :factor_ab, :avals, :bvals), Tuple{Int64, Int64, NamedTuple{(:nua, :nub, :NUa, :NUb), NTuple{4, Float64}}, Float64, Int64, Int64, Int64, Int64, Int64, Vector{Int64}, Vector{Int64}, Vector{Int64}, Float64, Vector{Int64}, Vector{Int64}}} #diff13
    temp_args_arr = Vector{tuple_type}(undef, ntot)
    temp_fill_mat = Matrix{TT}(undef,nbasis_total, nbasis_total)
        
    
    # now constructing Structs for different steps in the program:
    precomp_arrs = PrecomputeStruct(gamma_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr)
    #temp_arrs = TempStruct(temp_clmk,temp_dlmk,temp_S,temp_D1,temp_D2)
    interpol_arrs = InterpolationStruct(alpha_arr,v_arr,A_mat,w_arr,w_interpol_arr,Ainv_arr_kine) # ,v_obs_arr,w_obs_arr,w_obs_interpol_arr
    fill_arrs = FillStruct(fij_arr,Kij_arr,w_arr_kine,wn_interpol_arr,kij_arr,gij_arr,T,V,S,temp_args_arr,temp_fill_mat) # ,wn_obs_interpol_arr
    result_arrs = ResultStruct(energies_arr,wavefun_arr,centobs_output,R2_output)
    
    return precomp_arrs,interpol_arrs,fill_arrs,result_arrs

    ## again, no huge differences between 1D and 3D.
end