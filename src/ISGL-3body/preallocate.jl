# functions to preallocate all arrays for the ISGL program

# contains:
# - preallocate_data: for everything

struct PrecomputeStruct
    gamma_dict::Dict{Float64, Float64}
    cleb_arr::OffsetArray{Float64, 5, Array{Float64, 5}}
    spintrafo_dict::Dict{Tuple{Int64,Int64,Float64,Float64,Float64},Float64}
    spinoverlap_dict::Dict{Tuple{Int64,Int64,Int64,Float64,Float64,Float64,Float64},Float64}
    global6j_dict::Dict{Tuple{Float64,Float64,Float64,Float64},Float64}
    facsymm_dict::Dict{Tuple{Int64,Int64,Int64,Int64,Float64,Float64},Float64}
    jmat::Matrix{SMatrix{2, 2, Float64, 4}}
    murR_arr::MMatrix{2, 3, Float64, 6}
    nu_arr::Vector{Float64}
    NU_arr::Vector{Float64}
    norm_arr::OffsetMatrix{Float64, Matrix{Float64}}#OffsetMatrix{Float64, Matrix{Float64}} # not sure if the size here is required more accurately?
    NORM_arr::OffsetMatrix{Float64, Matrix{Float64}}#OffsetMatrix{Float64, Matrix{Float64}}
    Clmk_arr::OffsetArray{Float64, 3, Array{Float64, 3}} # 3,4,5 gibt einfach nur die dimensionen an, nicht deren Länge!
    Dlmk_arr::OffsetArray{ComplexF64, 4, Array{ComplexF64, 4}}
    S_arr::OffsetArray{Float64, 6, Array{Float64, 6}}
    SSO_arr::OffsetArray{Float64,8,Array{Float64, 8}}
end

struct TempStruct#{maximax,maxkmax}
    temp_clmk::Vector{Float64}#MVector{maxkmax, Float64}
    temp_dlmk::Matrix{ComplexF64}#MMatrix{maxkmax, 3, ComplexF64}
    temp_S::Vector{Float64}#MVector{maximax, Float64}
    temp_D1::Vector{ComplexF64}#MVector{3, ComplexF64}
    temp_D2::Vector{ComplexF64}#MVector{3, ComplexF64}
end

struct InterpolationStruct{T}
    alpha_arr::Vector{Float64}
    v_arr::Matrix{T}
    A_mat::Matrix{Float64}
    w_arr::Array{T, 4}
    w_interpol_arr::OffsetArray{Interpolations.Extrapolation{T, 1, ScaledInterpolation{T, 1, Interpolations.BSplineInterpolation{T, 1, OffsetVector{T, Vector{T}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}, 4, Array{Interpolations.Extrapolation{T, 1, ScaledInterpolation{T, 1, Interpolations.BSplineInterpolation{T, 1, OffsetVector{T, Vector{T}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}, 4}}#OffsetArray{Interpolations.Extrapolation{Float64, 1, ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, OffsetVector{Float64, Vector{Float64}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}}
    Ainv_arr_kine::OffsetMatrix{Float64, Matrix{Float64}}
    v_obs_arr::Matrix{Float64}
    w_obs_arr::Array{Float64, 5}
    w_obs_interpol_arr::OffsetArray{Interpolations.Extrapolation{Float64, 1, ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, OffsetVector{Float64, Vector{Float64}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}, 4, Array{Interpolations.Extrapolation{Float64, 1, ScaledInterpolation{Float64, 1, Interpolations.BSplineInterpolation{Float64, 1, OffsetVector{Float64, Vector{Float64}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{Base.OneTo{Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Tuple{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}}, BSpline{Cubic{Line{OnGrid}}}, Throw{Nothing}}, 4}}
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
    temp_args_arr::Vector{NamedTuple{(:rowi, :coli, :ranges, :norm4, :mij_arr, :sa, :JsSa, :sb, :JsSb, :JlLa, :JlLb, :la, :La, :lb, :Lb, :Lsum, :avals_new, :bvals_new, :factor_ab, :avals, :bvals),Tuple{Int64,Int64,NamedTuple{(:nua, :nub, :NUa, :NUb),Tuple{Float64,Float64,Float64,Float64}},Float64,Array{SArray{Tuple{6},Int64,1,6},1},Float64,Float64,Float64,Float64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Array{Int64,1},Vector{Int64},Int64,Vector{Int64},Vector{Int64}}}}
    temp_fill_mat::Matrix{T}
    wn_obs_interpol_arr::OffsetVector{Float64, Vector{Float64}}
end

struct ResultStruct{T}
    energies_arr::Vector{T}
    wavefun_arr::Matrix{T}
    centobs_output::Array{Float64}
    R2_output::Matrix{Float64}
end


function preallocate_data(phys_params,num_params,observ_params,size_params,csm_bool)
    if csm_bool == 0
        TT = Float64
    elseif csm_bool == 1
        TT = ComplexF64
    end
    
    #Destructing Structs:
    (;J_tot) = phys_params
    (;gem_params,kmax_interpol) = num_params
    (;nmax,Nmax) = gem_params
    (;nintmax,nbasis_total,nlL,nl,maxlmax,maximax,maxkmax,maxobs,JlL_complete) = size_params
    (;stateindices,centobs_arr,R2_arr) = observ_params
    
    # Allocate facts & cleb
    gamma_dict = Dict{Float64, Float64}() # necessary dimensionality not clear atm
    #cleb_arr = OffsetArray{Float64}(undef, 0:maxlmax, -maxlmax:maxlmax, 0:maxlmax, -maxlmax:maxlmax)*0.0 #
    #arguments: la,ma,La,Ma
    JlL_range = minimum(JlL_complete):maximum(JlL_complete) # OffsetArray needs range as input
    r3 = min(0,minimum(JlL_complete)):max(maxlmax,maximum(JlL_complete)) # for third argument of cleb_arr instead of 0:maxlmax
    m3 = max(maxlmax,maximum(JlL_complete)) # for m-arguments of cleb_arr instead of -maxlmax:maxlmax
    cleb_arr = OffsetArray{Float64}(undef, 0:maxlmax, -m3:m3, r3, -m3:m3, JlL_range)*0.0 #arguments: la,ma,La,Ma,JlL

    # Allocate spintrafo_dict and spinoverlap_dict
    spintrafo_dict = Dict{Tuple{Int64,Int64,Float64,Float64,Float64},Float64}()
    spinoverlap_dict = Dict{Tuple{Int64,Int64,Int64,Float64,Float64,Float64,Float64},Float64}()
    global6j_dict = Dict{Tuple{Float64,Float64,Float64,Float64},Float64}()

    # Dict for facsymm values:
    facsymm_dict = Dict{Tuple{Int64,Int64,Int64,Int64,Float64,Float64},Float64}()
    
    # Allocate jmat and murR_arr
    jmat = Matrix{SMatrix{2, 2, Float64, 4}}(undef, 3, 3);
    murR_arr = MMatrix{2,3,Float64,6}(zeros(2,3)); # reduced masses
    
    # Allocate ranges
    nu_arr = Vector{Float64}(undef, nmax)
    NU_arr = Vector{Float64}(undef, Nmax)        
    #ranges = (nu_arr, NU_arr)
    
    # Allocate norms
    norm_arr = OffsetMatrix{Float64}(undef, 0:maxlmax, nmax)*0.0 # one norm_arr for all l,n combinations. changed to OffsetArray for all values of l=0:maxlmax
    NORM_arr = OffsetMatrix{Float64}(undef, 0:maxlmax, Nmax)*0.0 # one norm_arr for all L,N combinations
    
    # ISGL arrays:
    Clmk_arr = OffsetArray{Float64}(undef, 0:maxlmax, -maxlmax:maxlmax, maxkmax)*0.0
    Dlmk_arr = OffsetArray{ComplexF64}(undef, 0:maxlmax, -maxlmax:maxlmax, maxkmax, 3)*0.0
    S_arr = OffsetArray{Float64}(undef,0:maxlmax,0:maxlmax,0:maxlmax,0:maxlmax,JlL_range,maximax)*0.0 # changed back to array from dict
    SSO_arr = OffsetArray{Float64}(undef,0:maxlmax,0:maxlmax,0:maxlmax,0:maxlmax,JlL_range,JlL_range,1:6,maximax)*0.0 # for spin-orbit interaction
    
    # arrays for range-interpolation method:
    alpha_arr = zeros(kmax_interpol)
    v_arr = zeros(TT,kmax_interpol,2*maxlmax+1) # careful, now NOT a OffsetArray! due to problems in A\v
    A_mat = zeros(2*maxlmax+1,2*maxlmax+1)
    w_arr = zeros(TT,3,kmax_interpol,2*maxlmax+1,2*maxlmax+1) # first dimensionality is 3 for covering interactions for all possible Jacobi-Sets (=3). in principle could be reduced to nboxes? penultimate dimensionality for different Lsum values!
    x=range(0.1,0.5,3)
    y=x.^2;interpoltype = typeof(cubic_spline_interpolation(x,y))# just for easy inferring the type of interpolation objects
    yc=x.^2 .+ zero(TT);interpoltypeC = typeof(cubic_spline_interpolation(x,yc))# just for easy inferring the type of interpolation objects # for possibly complex arguments!
    w_interpol_arr=OffsetArray{interpoltypeC}(undef,3,nintmax,0:2*maxlmax,0:2*maxlmax) # penultimate dimensionality for different Lsum values!; nintmax for (maximum) number of interactions
    #aac=zeros(27);gac=zeros(27);abc=zeros(27);gbc=zeros(27)
    Ainv_arr_kine = OffsetArray{Float64}(undef,0:2*maxlmax,0:2*maxlmax)*0.0 # penultimate dimensionality for different Lsum values!
    w_arr_kine = OffsetArray{Float64}(undef,0:2*maxlmax+1)*0.0
    wn_interpol_arr = OffsetArray{TT}(undef,0:2*maxlmax)*0.0 # for the interpolated wn_values that are actually used
    #for observables (in range-interpolation)
    v_obs_arr = zeros(kmax_interpol,2*maxlmax+1)
    w_obs_arr = zeros(3,maxobs,kmax_interpol,2*maxlmax+1,2*maxlmax+1)
    w_obs_interpol_arr = OffsetArray{interpoltype}(undef,3,maxobs,0:2*maxlmax,0:2*maxlmax)
    wn_obs_interpol_arr = OffsetArray{Float64}(undef,0:2*maxlmax)*0.0
    
    # arrays for shoulder method (used in matrix-element calculations)
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
    
    # temporary arrays for clmk,dlmk,S,D1,D2: temp_clmk,temp_dlmk,temp_S,temp_D1,temp_D2
    temp_clmk=zeros(maxkmax)#MVector{maxkmax, Float64}(zeros(maxkmax))
    temp_dlmk=zeros(maxkmax,3)#MMatrix{maxkmax,3, ComplexF64}(zeros(maxkmax,3))
    temp_S=zeros(maximax)#MVector{maximax, Float64}(zeros(maximax))
    temp_D1=zeros(3)#MVector{3, ComplexF64}(zeros(3))
    temp_D2=zeros(3)#MVector{3, ComplexF64}(zeros(3))
    
    # temporary arrays for filling: function arguments and matrix for parallelization
    ntot = Int64((nbasis_total^2+nbasis_total)/2) # only for lower triangular!
    # Define the type of the tuple for the function arguments
    #tuple_type = NamedTuple{(:rowi, :coli, :ranges, :norm4, :mij_arr, :S_arr, :la, :La, :lb, :Lb, :Lsum, :avals_new, :bvals_new, :cvals, :factor_ab, :avals, :bvals), Tuple{Int64, Int64, NamedTuple{(:nua, :nub, :NUa, :NUb), NTuple{4, Float64}}, Float64, Vector{SVector{6, Int64}}, OffsetArray{Float64, 5, Array{Float64, 5}}, Int64, Int64, Int64, Int64, Int64, Vector{Int64}, Vector{Int64}, Vector{Int64}, Float64, Vector{Int64}, Vector{Int64}}}
    tuple_type = NamedTuple{(:rowi, :coli, :ranges, :norm4, :mij_arr, :sa, :JsSa, :sb, :JsSb, :JlLa, :JlLb, :la, :La, :lb, :Lb, :Lsum, :avals_new, :bvals_new, :factor_ab, :avals, :bvals),Tuple{Int64,Int64,NamedTuple{(:nua, :nub, :NUa, :NUb),Tuple{Float64,Float64,Float64,Float64}},Float64,Array{SArray{Tuple{6},Int64,1,6},1},Float64,Float64,Float64,Float64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Array{Int64,1},Vector{Int64},Int64,Vector{Int64},Vector{Int64}}}
    temp_args_arr = Vector{tuple_type}(undef, ntot)
    temp_fill_mat = zeros(TT,nbasis_total, nbasis_total)
        
    
    # now constructing Structs for different steps in the program:
    precomp_arrs = PrecomputeStruct(gamma_dict,cleb_arr,spintrafo_dict,spinoverlap_dict,global6j_dict,facsymm_dict,jmat,murR_arr,nu_arr,NU_arr,norm_arr,NORM_arr,Clmk_arr,Dlmk_arr,S_arr,SSO_arr)
    temp_arrs = TempStruct(temp_clmk,temp_dlmk,temp_S,temp_D1,temp_D2)
    interpol_arrs = InterpolationStruct(alpha_arr,v_arr,A_mat,w_arr,w_interpol_arr,Ainv_arr_kine,v_obs_arr,w_obs_arr,w_obs_interpol_arr)
    fill_arrs = FillStruct(fij_arr,Kij_arr,w_arr_kine,wn_interpol_arr,kij_arr,gij_arr,T,V,S,temp_args_arr,temp_fill_mat,wn_obs_interpol_arr)
    result_arrs = ResultStruct(energies_arr,wavefun_arr,centobs_output,R2_output)
    
    return precomp_arrs,temp_arrs,interpol_arrs,fill_arrs,result_arrs
end