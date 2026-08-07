# functions to preallocate all arrays for the GEM3B1D program

# contains:
# - preallocate_data: for everything

struct PrecomputeStruct{TR, TN}
    gamma_dict::Dict{Float64, Float64}
    jmat::Matrix{SMatrix{2, 2, Float64, 4}}
    murR_arr::MMatrix{2, 3, Float64, 6}
    nu_arr::Vector{TR}
    NU_arr::Vector{TR}
    norm_arr::OffsetMatrix{TN, Matrix{TN}}
    NORM_arr::OffsetMatrix{TN, Matrix{TN}}
end

struct InterpolationStruct{T}
    alpha_arr::Vector{Float64}
    v_arr::Matrix{T}
    w_interpol_arr::OffsetArray{Interpolations.Extrapolation{T, 1, ScaledInterpolation{T, 1, Interpolations.BSplineInterpolation{T, 1, OffsetVector{T, Vector{T}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}, 3, Array{Interpolations.Extrapolation{T, 1, ScaledInterpolation{T, 1, Interpolations.BSplineInterpolation{T, 1, OffsetVector{T, Vector{T}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}, 3}}
    #v_obs_arr::Matrix{Float64}
    #w_obs_arr::Array{Float64, 5}
    #w_obs_interpol_arr::OffsetArray{Interpolations.Extrapolation{Float64, 1, ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, OffsetVector{Float64, Vector{Float64}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}, 4, Array{Interpolations.Extrapolation{Float64, 1, ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, OffsetVector{Float64, Vector{Float64}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}, 4}}
end

struct TempArgs1D{TR, TN}
    rowi::Int64
    coli::Int64
    ranges::NamedTuple{(:nua, :nub, :NUa, :NUb), NTuple{4, TR}}
    norm4::TN
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

struct FillStruct{TTV, TS, TR, TN}
    wn_interpol_arr::OffsetVector{TTV, Vector{TTV}}
    T::Matrix{TTV}
    V::Matrix{TTV}
    S::Matrix{TS}
    temp_args_arr::Vector{TempArgs1D{TR, TN}}
    temp_fill_mat::Matrix{TTV}
    #wn_obs_interpol_arr::OffsetVector{Float64, Vector{Float64}}
end

struct ResultStruct{TE, TW}
    energies_arr::Vector{TE}
    wavefun_arr::Matrix{TW}
    #centobs_output::Array{Float64}
    #R2_output::Matrix{Float64}
end


function preallocate_data(phys_params,num_params,observ_params,size_params,complex_ranged_r::Bool,complex_ranged_R::Bool,complex_scaling::Bool)
    complex_any = complex_ranged_r || complex_ranged_R
    TTV = (complex_any || complex_scaling) ? ComplexF64 : Float64
    TS = complex_any ? ComplexF64 : Float64
    TE = complex_scaling ? ComplexF64 : Float64
    TR = complex_any ? ComplexF64 : Float64
    TN = complex_any ? ComplexF64 : Float64
    
    #Destructing Structs:
    (;gem_params,kmax_interpol) = num_params
    (;nmax,Nmax) = gem_params
    (;nbasis_total,nlL,nl,maxlmax,maxobs,nintmax) = size_params
    (;stateindices,centobs_arr,R2_arr) = observ_params
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
    norm_arr = OffsetMatrix{TN}(undef, 0:maxlmax, nmax_eff)*zero(TN)
    NORM_arr = OffsetMatrix{TN}(undef, 0:maxlmax, Nmax_eff)*zero(TN)
    
    # arrays for range-interpolation method:
    alpha_arr = zeros(kmax_interpol)
    v_arr = zeros(TTV,kmax_interpol,2*2*maxlmax+1) # careful, now NOT an OffsetArray! due to problems in A\v # changed to 4*maxlmax+1. 

    x=range(0.1,0.5,3)
    y=x.^2;interpoltype = typeof(cubic_spline_interpolation(x,y))# just for easy inferring the type of interpolation objects # seems quite hacked, but works for now
    yc=x.^2 .+ zero(TTV);interpoltypeC = typeof(cubic_spline_interpolation(x,yc))# just for easy inferring the type of interpolation objects # for possibly complex arguments!

    w_interpol_arr=OffsetArray{interpoltypeC}(undef,3,nintmax,0:4*maxlmax) # penultimate dimensionality for different Lsum values! removed that dimension for 1D
    wn_interpol_arr = OffsetArray{TTV}(undef,0:4*maxlmax)*zero(TTV) # for the interpolated wn_values that are actually used

    #for observables (in range-interpolation) # currently not available in 1D
#=     v_obs_arr = zeros(kmax_interpol,2*maxlmax+1)
    w_obs_arr = zeros(3,maxobs,kmax_interpol,2*maxlmax+1,2*maxlmax+1)
    w_obs_interpol_arr = OffsetArray{interpoltype}(undef,3,maxobs,0:2*maxlmax,0:2*maxlmax)
    wn_obs_interpol_arr = OffsetArray{Float64}(undef,0:2*maxlmax)*0.0 =#
    
    # matrices S: for the norm overlap, T: for the kintetic energy (later for the Hamiltonian T+V), V: for the potential energy
    S = zeros(TS,nbasis_total,nbasis_total)
    T = zeros(TTV,nbasis_total,nbasis_total)
    V = zeros(TTV,nbasis_total,nbasis_total);
    
    
    # for results:
    energies_arr = zeros(TE,nbasis_total);# Vector{TT}(undef, nbasis_total) # changed to zeros to avoid bad behavior due to undef values not occupied thresholding (eigen2step)
    wavefun_arr = Matrix{TTV}(undef, nbasis_total, nbasis_total)
    #centobs_output = Array{Float64}(undef, 3, maxobs, lastindex(stateindices))
    #R2_output = zeros(3, lastindex(stateindices))
    
    # temporary arrays for filling: function arguments and matrix
    ntot = (complex_any && complex_scaling) ? Int64(nbasis_total^2) : Int64((nbasis_total^2+nbasis_total)/2)
    # Define the type of the tuple for the function arguments
    temp_args_arr = Vector{TempArgs1D{TR, TN}}(undef, ntot)
    temp_fill_mat = Matrix{TTV}(undef,nbasis_total, nbasis_total)
        
    
    # now constructing Structs for different steps in the program:
    precomp_arrs = PrecomputeStruct(gamma_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr)
    #temp_arrs = TempStruct(temp_clmk,temp_dlmk,temp_S,temp_D1,temp_D2)
    interpol_arrs = InterpolationStruct(alpha_arr,v_arr,w_interpol_arr) # ,v_obs_arr,w_obs_arr,w_obs_interpol_arr
    fill_arrs = FillStruct(wn_interpol_arr,T,V,S,temp_args_arr,temp_fill_mat) # ,wn_obs_interpol_arr
    result_arrs = ResultStruct(energies_arr,wavefun_arr) #,centobs_output,R2_output)
    
    return precomp_arrs,interpol_arrs,fill_arrs,result_arrs

    ## again, no huge differences between 1D and 3D.
end