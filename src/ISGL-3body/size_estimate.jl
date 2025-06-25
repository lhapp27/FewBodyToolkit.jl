# functions to calculate the parameters that determine the necessary sizes of arrays within the ISGL program

# contains:
# - size_estimate that bundles everything:
# - index_interaction_types to determine the indices of the interaction types
# - abc_size for box-size of TVS Matrices
# - sS/lLsS/lL-coupl for coupling of spins and angular momenta
# - boxsize_fun for box numbers and number of basis functions in each box
# - imax_fun for imax in ISGL
# - kmax_fun for kmax in ISGL

struct SizeParams{T<:Number}
    abvals_arr::Vector{Vector{Int64}}
    cvals::Vector{Int64}
    gauss_indices::Vector{Vector{Int}}
    gaussopt_arr::Vector{Vector{Tuple{Float64,T}}}
    central_indices::Vector{Vector{Int}}
    so_indices::Vector{Vector{Int}}
    nint_arr::Vector{Int64}
    nintmax::Int64
    groupindex_arr::Vector{Int64}
    nboxes::Int64
    abI::Int64
    factor_bf::Int64
    box_size_arr::Vector{Int64}
    nbasis_total::Int64
    starts::Vector{Int64}
    ends::Vector{Int64}
    bvalsdiag::Vector{Vector{Int64}} # actually I think it could be simplified to Vector{Int64}
    s_arr::Vector{Vector{Float64}} # added for spin-systems
    JsS_arr::Vector{Vector{Vector{Float64}}} # added for spin-systems
    s_complete::Vector{Float64} # added for spin-systems
    JsS_complete::Vector{Float64} # added for spin-systems
    JlL_arr::Vector{Vector{Vector{Vector{Int64}}}} # added for spin-systems
    JlL_complete::Vector{Int64} # added for spin-systems
    lL_nested::Vector{Vector{Vector{Vector{Vector{Tuple{Int64, Int64}}}}}} # changed for spin-systems
    #lL_arr::Vector{Vector{Tuple{Int64, Int64}}}
    lL_complete::Vector{Tuple{Int64, Int64}}
    l_complete::Vector{Int64}
    nlL::Int64
    nl::Int64
    maxlmax::Int64
    imax_dict::Dict{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}}, Int64}
    imaxSO_dict::Dict{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}, Int64}, Int64}
    maximax::Int64
    mij_arr_dict::Dict{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}},Vector{SVector{6, Int64}}}
    mijSO_arr_dict::Dict{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}, Int64},Vector{SVector{6, Int64}}}
    kmax_dict::Dict{Tuple{Int64, Int64},Int64}
    maxkmax::Int64
    maxobs::Int64
end

function size_estimate(phys_params,num_params,observ_params,csm_bool)
    
    # input interpretation:
    (;mass_arr,svals,vint_arr,J_tot,parity,spin_arr) = phys_params
    (;lmax,Lmax,gem_params,theta_csm,omega_cr,mu0,c_shoulder,kmax_interpol,lmin,Lmin) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;stateindices,centobs_arr,R2_arr) = observ_params
    
    # size of Hamiltonian matrix due to symmetries and interactions
    cvals = findall(isempty.(vint_arr) .==0 ) # consider only values for c where there are interactions (any type)
    
    # number of interactions per Jacobi-set c:
    gauss_indices, gaussopt_arr, central_indices, so_indices, nint_arr, nintmax = index_interaction_types(vint_arr,csm_bool, theta_csm)

    # box sizes, indices and factors for symmetrization
    abvals_arr,groupindex_arr,nboxes,abI,factor_bf = abc_size(cvals,svals)
    
    # coupling of angular momenta
    s_arr,JsS_arr,s_complete,JsS_complete = sScoupl(spin_arr,cvals)
    JlL_arr,JlL_complete = lLsScoupl(J_tot,JsS_arr,s_arr,cvals)
    lL_nested,lL_complete,l_complete = lLcoupl(J_tot,parity,cvals,svals,spin_arr,s_arr,JsS_arr,JlL_arr,lmin,Lmin,lmax,Lmax)

    # box sizes = number of basis functions in each box; starts,ends = indices for boxes within big matrix; bvalsdiag = simplification for identical particles
    box_size_arr,nbasis_total = boxsize_fun(groupindex_arr,abvals_arr,lL_nested,nmax,Nmax)
    starts = [1; cumsum(box_size_arr[1:end-1]) .+ 1]
    ends = cumsum(box_size_arr)
    #bvalsdiag = [[abvals_arr[boxR][1]] for boxR in groupindex_arr]
    bvalsdiag = [Int64[], Int64[], Int64[]]
    for boxC in groupindex_arr
        bvalsdiag[boxC] = [abvals_arr[boxC][1]]
    end

    # imax for ISGL
    imax_dict,mij_arr_dict = imax_fun(lL_complete)
    imaxSO_dict,mijSO_arr_dict = imax_funSO(lL_complete) # for SO
    
    # kmax for ISGL
    kmax_dict = kmax_fun(l_complete)


    ## max-values:
    nlL = lastindex(lL_complete) # former: lLmax'
    nl = lastindex(l_complete) # former lmax'
    maxlmax=max(lmax,Lmax)
    maximax=findmax(imax_dict)[1]
    maxkmax=findmax(kmax_dict)[1] # [1] is the value; [2] would be the key
    maxobs = maximum(lastindex.(centobs_arr)) # max number of observables

    # Constructing Struct (collective data structure size_params with all the size parameters)
    size_params = SizeParams(abvals_arr,cvals,gauss_indices,gaussopt_arr,central_indices,so_indices,nint_arr,nintmax,groupindex_arr,nboxes,abI,factor_bf,box_size_arr,nbasis_total,starts,ends,bvalsdiag,s_arr,JsS_arr,s_complete,JsS_complete,JlL_arr,JlL_complete,lL_nested,lL_complete,l_complete,nlL,nl,maxlmax,imax_dict,imaxSO_dict,maximax,mij_arr_dict,mijSO_arr_dict,kmax_dict,maxkmax,maxobs)
    
    return size_params
end


"""
    csmgaussopt(gaussopt, csm_bool, theta_csm)

If `csm_bool == 1`, returns a new
`Vector{Vector{Tuple{Float64,ComplexF64}}}` where each
`mu` has been scaled by `exp(2im * theta_csm * pi/180)`.
Otherwise returns the original `gaussopt` unchanged.
"""
function csmgaussopt(gaussopt::Vector{Vector{Tuple{Float64,Float64}}},csm_bool::Integer,theta_csm::Real)
    if csm_bool == 1
        csmfac = exp(2im * theta_csm * pi / 180)
        # build a new array with ComplexF64 mu’s
        return [ 
            [ (v0, mu0 * csmfac) 
              for (v0, mu0) in cc ] 
            for cc in gaussopt 
        ]
    else
        return gaussopt
    end
end


function index_interaction_types(vint_arr,csm_bool, theta_csm)
    gauss_indices = [Int[] for _ in 1:3]
    gaussopt_arr = [Tuple{Float64,Float64}[] for _ in 1:3]
    central_indices = [Int[] for _ in 1:3]
    so_indices = [Int[] for _ in 1:3]
    nint_arr = zeros(Int64,3)

    pushindexpotentialtype!(v::Function, central_indices, gauss_indices, so_indices, i) = push!(central_indices, i) # treat function as a central potential
    pushindexpotentialtype!(v::CentralPotential, central_indices, gauss_indices, so_indices, i) = push!(central_indices, i)
    pushindexpotentialtype!(v::SpinOrbitPotential, central_indices, gauss_indices, so_indices, i) = push!(so_indices, i)
    pushindexpotentialtype!(v::GaussianPotential, central_indices, gauss_indices, so_indices, i) = push!(gauss_indices, i)
    

    for c in 1:3
        gauss_indices[c] = Int[]
        central_indices[c] = Int[]
        so_indices[c] = Int[]
        for (i, v) in enumerate(vint_arr[c])
            pushindexpotentialtype!(v, central_indices[c], gauss_indices[c], so_indices[c], i)

            if i in gauss_indices[c] # if this is a Gaussian potential
                push!(gaussopt_arr[c], (v.v0, v.mu_g)) # store the parameters of the Gaussian potential
            else
                push!(gaussopt_arr[c], (NaN, NaN)) # if not a Gaussian potential, store NaN. This is necessary to keep the order of indices.
            end

        end
        nint_arr[c] = lastindex(vint_arr[c])
    end

    nintmax = maximum(nint_arr)

    gaussopt_arrC = csmgaussopt(gaussopt_arr, csm_bool, theta_csm) # adjust to complex values in case of csm_bool == 1

    return gauss_indices, gaussopt_arrC, central_indices, so_indices, nint_arr, nintmax
end


#this function is maybe a bit confusing due to numerous names and definitions. thats why there are many examples and explanations.
function abc_size(cvals,svals)
    
    ## Step 1: Simplify system based on symmetry of particles
    
    uniq = unique(svals) # what types of particles do we have?
    ntypes = lastindex(uniq) # nr. of different particle types = nr of groups = nr of boxes we need (as of now)
    typepos_arr = Vector{Vector{Int}}(undef, ntypes )#Array{Array}(undef, ntypes) # array of (box-)position values for each particle type; determines sum-indices
    #grouplength_arr = Array{Int}(undef, ntypes) # array of number of possible positions within each group; actually irrelevant. typepos_arr contains all information
    
    # example:          3 ident         2+1             3 diff    
    # svals             ["b","b","b"]   ["b","b","x"]   ["x3","x2","x3"]    # input!
    
    # uniq              ["b"]           ["b","x"]       ["x1","x2","x3"]
    # ntypes            1               2               3
    # typepos_arr       [[1,2,3]]       [[1, 2], [3]]   [[1], [2], [3]]
    # grouplength_arr   [3]             [2, 1]          [1, 1, 1]
    
    for i = 1:ntypes
        typepos_arr[i] = findall(svals .== uniq[i])
        #grouplength_arr[i] = lastindex(typepos_arr[i])
    end
        
    
    ## Step 2: now throw away irrelevant Faddeev components, as governed by the input cvals:
    # cvals = [1,2,3] --> all Faddeev components are relevant; output as above
    
    # cvals = [1,2]   --> Faddeev component 3 is NOT relevant --> box-dimensionality and sum-size are reduced
    # example:          3 ident         2+1             3 diff    
    # svals             ["b","b","b"]   ["b","b","x"]   ["x3","x2","x3"]    # input!
    # typepos_arr       [[1,2,3]]       [[1, 2], [3]]   [[1], [2], [3]]
    # abvals_arr        -makesnosense-  [[1, 2], []]    [[1], [2], []]
    # groupindex_arr    -makesnosense-  [1]             [1, 2]
    # nboxes            -makesnosense-  1               2
    
    abvals_arr = similar(typepos_arr) # abvals is the new typepos_arr
    for i = 1:ntypes
        abvals_arr[i] = intersect(typepos_arr[i],cvals)
    end
    # Which group / box index to keep: useful for what again? for referencing index of ab_arr! -> ab_arr[grpind_arr[i]] contains the values of a,b for consideration
    groupindex_arr = findall(.!isempty.(abvals_arr))
    
    # how many boxes -> box-dimensionality
    nboxes = lastindex(groupindex_arr) # nboxes is the updated ntypes.
    
    #return uniq,ntypes,typepos_arr,grouplength_arr," | ",typepos_arr2,groupindex_arr

    # in which box are the identical particles? -> boxI
    # which Faddeev component gets the symmetry factor +-(-1)^l_a/b? -> abI
    if lastindex(uniq) == 2 ## why do we need it only for 2+1 systems?! because for 3 identicals, it is ok to assume l_i = even/odd, for all i=1,2,3. for 2+1 systems, i.e. bbz, l_3 will be even/odd automatically, but l_1,2 not. therefore the additional factor is required
        if in("b",uniq)
            boxI = findfirst(uniq .== "b")
            abI = abvals_arr[boxI][end]
            factor_bf = 1
        elseif in("f",uniq)
            boxI = findfirst(uniq .== "f")
            abI = abvals_arr[boxI][end] # which component gets the symmetry factor +-(-1)^l_a/b ? end: always the second one of the two identicals
            factor_bf = -1
        else
            error("Error in svals = $svals. If two particles are identical choose \"b\" for bosons or \"f\" for fermions.")
        end
    else
        boxI = 0
        abI = 0
        factor_bf = 1
    end

    return abvals_arr,groupindex_arr,nboxes,abI,factor_bf
end

#ok there are 3 aspects:
# 1. J = JlL + JsS coupling
# 2. parity = (-1)^(l+L)
# 3. (anti-)symmetrization (for each identical pair): (-1)^s * (-1)^l must be +1 (-1) for bosons (fermions)
# and also the normal angular momentum addition JjL = l+L, as well as sij=si+sj and JsS = sij+Sk


# function for spin-spin coupling
# this function does not respect whether some values of sS are not allowed due to J=lL+sS coupling!
function sScoupl(spin_arr,cvals)
    z_arr = spin_arr # single-particle spins
    
    s_arr = ([Vector{Float64}() for _ in 1:3]) # for each c (1:3) there is a vector of possible values for s_k = s_ij = z_i + z_j (two-particle spins)
    #S_arr = zeros(Float64,3) # values for S_k = z_k
    JsS_arr = ([Vector{Vector{Float64}}() for _ in 1:3]) # for each c (1:3) there is a vector of possible values for sS_k = s_k + S_k (effective total spin).
    
    for c in cvals
        a = mod(c,3)+1; # c-2 = c+1
        b = mod(c+1,3)+1; # c-1 = c+2
        za = z_arr[a]
        zb = z_arr[b]
        possible_svals = abs(za-zb):1:abs(za+zb) # coupling of two single-particle spins
        s_arr[c] = Vector(possible_svals)
        
        for (is,sc) in enumerate(s_arr[c])
            zc = z_arr[c]
            possible_sSvals = abs(sc-zc):1:abs(sc+zc) # coupling of the single particle spin z_c to the pair-spin s_c
            #@show(c,is,sc,zc,Vector(possible_sSvals))
            push!(JsS_arr[c],Vector(possible_sSvals))
        end
    end
    
    # find complete list of possible values for s_c
    s_complete = unique(reduce(vcat, s_arr))    
    # similar for sS (this should be independent of c!)
    sS_complete = Vector{Float64}();
    for c=1:3
        for vecs in JsS_arr[c]                                 
            for vals in vecs                                    
                push!(sS_complete,vals)
            end
        end
    end
    sS_complete = unique(sS_complete)
    # this complete list is only based on physically available spins, and is not limited by numerical truncation. it is important for the determination of possible lL values via J = lL+sS
    
    return s_arr,JsS_arr,s_complete,sS_complete
end


#function for J = lL + sS coupling
function lLsScoupl(J_tot,JsS_arr,s_arr,cvals)
    # for each c there are several possible values of sc = s_arr[c], and for each of those values there are several possible values of JsS = JsS_arr[c][is], and we want to find the possible values for JlL = J + (-JsS) (coupling with JsS to J).

    # using nested vector of vectors:
    JlL_arr = [Vector{Vector{Vector{Int64}}}() for _ in 1:3] #before: cvals, but makes problems when c is larger then the amount of cvals.
    for c in cvals # for all c-values
        JlL_arr[c] = [Vector{Vector{Int64}}() for _ in s_arr[c]]
        for (is, sc) in enumerate(s_arr[c])
            JlL_arr[c][is] = [Vector{Int64}() for _ in JsS_arr[c][is]]
            for (iss, JsS) in enumerate(JsS_arr[c][is])
                JlL_arr[c][is][iss] = Vector{Int64}(abs(J_tot - JsS):1:abs(J_tot + JsS))
            end
        end
    end

    # Flatten lL_nested to build the complete list of JlL values
    temp_JlL = Vector{Int64}()
    for c in cvals
        for is in eachindex(JlL_arr[c])
            for iss in eachindex(JlL_arr[c][is])
                for ill in eachindex(JlL_arr[c][is][iss])
                    append!(temp_JlL, JlL_arr[c][is][iss][ill])
                end
            end
        end
    end
    # Determine the maximum number of (l,L)-pairs in any innermost vector
    JlL_complete = unique(temp_JlL)

    return JlL_arr,JlL_complete
end


# function to determine list of (l,L)-tuples to consider, based on J_tot,parity,svals. This version considers also the spins!
function lLcoupl(J,parity,cvals,svals,spin_arr,s_arr,JsS_arr,JlL_arr,lmin,Lmin,lmax,Lmax)
    lL_nested = [Vector{Vector{Vector{Vector{Tuple{Int64, Int64}}}}}() for _ in 1:3] # before: cvals, but makes problems when c is larger then the amount of cvals (same as for JlL_arr)

    for c in cvals # for all c-values
        lL_nested[c] = [Vector{Vector{Vector{Tuple{Int64, Int64}}}}() for _ in s_arr[c]]
        for (is, sc) in enumerate(s_arr[c])
            lL_nested[c][is] = [Vector{Vector{Tuple{Int64, Int64}}}() for _ in JsS_arr[c][is]]
            seff = spin_arr[mod(c + 1, 3) + 1] + spin_arr[mod(c, 3) + 1] - sc# effective seff = s[c]+za+zb which indicates the parity of the spin wave function of particle a+b. only necessary for (anti-)symmetrization of identical particles; pi_s = (-1)^seff
            for (iss, JsS) in enumerate(JsS_arr[c][is])
                lL_nested[c][is][iss] = [Vector{Tuple{Int64, Int64}}() for _ in JlL_arr[c][is][iss]]
                # we could calculate JlJ values also within this loop. for now: separate function
                for (ilL, JlL) in enumerate(JlL_arr[c][is][iss])

                    for l in lmin:lmax
                        for L in Lmin:Lmax
                            if abs(l - L) <= JlL <= l + L # allowed J for l, L coupling
                                if (-1)^(l + L) == parity # allowed combination of l, L from global parity
                                    #(anti-)symmetrization: pi_l = (-1)^l; pi_ls = (-1)^(l+seff) and must be +(-)1 for identical bosons (fermions)
                                    if svals[mod(c + 1, 3) + 1] == svals[mod(c, 3) + 1] == "b"  # if Jacobi-Set c contains relative coordinate between two identical bosons -> even parity: 
                                        mod(l+seff, 2) != 0 && continue #skip if l+seff is not even
                                    elseif svals[mod(c + 1, 3) + 1] == svals[mod(c, 3) + 1] == "f" # if Jacobi-Set c contains relative coordinate between two identical fermions -> odd parity:
                                        mod(l+seff, 2) != 1 && continue #skip if l+seff is not odd
                                    end
                                    
                                    push!(lL_nested[c][is][iss][ilL], (l, L))

                                end
                            end
                        end
                    end

                end
            end
        end
    end


    # Flatten lL_nested to build the complete list of (l,L) tuples
    temp_tuples = Tuple{Int64, Int64}[]
    for c in cvals
        for is in eachindex(lL_nested[c])
            for iss in eachindex(lL_nested[c][is])
                for ill in eachindex(lL_nested[c][is][iss])
                    append!(temp_tuples, lL_nested[c][is][iss][ill])
                end
            end
        end
    end

    # Determine the maximum number of (l,L)-pairs in any innermost vector
    lL_complete = unique(temp_tuples)

    # Build the complete list of individual l and L values needed in the program
    l_vals = Int[]
    for (l, L) in temp_tuples
        push!(l_vals, l, L)
    end
    l_complete = unique(l_vals)

    if lastindex(lL_complete) == 0
        error("Error in size_estimate: No valid (l,L)-couplings found.")
        return
    end

    return lL_nested, lL_complete, l_complete
end


## function to determine the size (= number of basis functions) in each box
function boxsize_fun(groupindex_arr,abvals_arr,lL_nested,nmax,Nmax)
    box_size_arr = zeros(Int64,3)
    # alternative via loops:
    for grp in groupindex_arr # a,b,c
        ab = abvals_arr[grp][1] #relevant groupindex. e.g. for 2+1 system is the third jacobi set in group 2!
        nbasis_in_box = 0
        for is in keys(lL_nested[ab])
            for ijss in keys(lL_nested[ab][is])
                for ijll in keys(lL_nested[ab][is][ijss])
                    for ilL in keys(lL_nested[ab][is][ijss][ijll])
                        for n = 1:nmax
                            for N = 1:Nmax
                                nbasis_in_box += 1
                            end
                        end
                    end
                end
            end
        end
        box_size_arr[grp] = nbasis_in_box
    end
    nbasis_total = Int64(sum(box_size_arr))
    
    return box_size_arr,nbasis_total
end

# imax for ISGL:
function imax_fun(lL_complete)
    imax_dict = Dict{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}}, Int64}()
    mij_arr_dict = Dict{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}},Vector{SVector{6, Int64}}}() # this is a dict (for the different (la,La),(lb,Lb) pairs) of vectors (each has length imax) of SVectors (length 6, for m12,m13,...).
    for i in keys(lL_complete)
        for j in keys(lL_complete)
            imax_dict[lL_complete[i],lL_complete[j]], mij_arr_dict[lL_complete[i],lL_complete[j]] = imax(lL_complete[i],lL_complete[j])
        end
    end
    return imax_dict,mij_arr_dict
end

#function to determine imax and mij-tuple used within S-coeff for ISGL
function imax((la,La),(lb,Lb)) # u1=m12,u2=m13,...
    if mod(la+La-lb-Lb,2) != 0
        println("Error in imax2: Parity_a != Parity_b")
        return
    end
    
    Lsum = Int64((la+La+lb+Lb)/2)
    
    mij_arr = Vector{SVector{6, Int64}}();
    i=0
    for m12 = 0:Lsum # u1
        for m13 = 0:(Lsum-m12) # u2
            m14=la-m12-m13 # u3
            m23=Int64((la+La+lb-Lb)/2) - m12-m13 # u4
            m24=Int64((-la+La-lb+Lb)/2) + m13 # u5
            m34=Int64((-la-La+lb+Lb)/2) + m12 # u6
            (m14<0||m23<0||m24<0||m34<0) && continue # additional constraint: all mij >=0 !
            i=i+1
            push!(mij_arr,SA[m12,m13,m14,m23,m24,m34])
        end
    end
    
    #println("\ni=",i)
    return i,mij_arr
end

#same as above but for spin-orbit interactions
function imax_funSO(lL_complete)
    imaxSO_dict = Dict{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}, Int64}, Int64}() # additional Int64 argument for v
    mijSO_arr_dict = Dict{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}, Int64},Vector{SVector{6, Int64}}}() # additional Int64 argument for v
    for i in keys(lL_complete)
        for j in keys(lL_complete)
            for v in 1:6
                imaxSO_dict[lL_complete[i],lL_complete[j],v], mijSO_arr_dict[lL_complete[i],lL_complete[j],v] = imaxSO(lL_complete[i],lL_complete[j],v)
            end
        end
    end
    return imaxSO_dict,mijSO_arr_dict
end

function imaxSO((la,La),(lb,Lb),v) # u1=m12,u2=m13,...
    # no simple parity check

    Lsum = Int64((la+La+lb+Lb)/2)
    jj = Lsum -1 #variable j

    lambdamat = SMatrix{6,4,Int64}(
        [1 1 0 0;
         1 0 1 0;
         1 0 0 1;
         0 1 1 0;
         0 1 0 1;
         0 0 1 1]
    )
    tla = la - lambdamat[v,1] # \tilde{l}_k = l_k - lambda_{v,k}; k \in \{a,A,b,B\}
    tLa = La - lambdamat[v,2]
    tlb = lb - lambdamat[v,3]
    tLb = Lb - lambdamat[v,4]
    
    mijSO_arr = Vector{SVector{6, Int64}}();
    iSO=0
    for m12 = 0:jj # u1
        for m13 = 0:(jj-m12) # u2
            m14=tla-m12-m13 # u3
            m23=Int64((tla+tLa+tlb-tLb)/2) - m12-m13 # u4
            m24=Int64((-tla+tLa-tlb+tLb)/2) + m13 # u5
            m34=Int64((-tla-tLa+tlb+tLb)/2) + m12 # u6
            (m14<0||m23<0||m24<0||m34<0) && continue # additional constraint: all mij >=0 !
            iSO += 1
            push!(mijSO_arr,SA[m12,m13,m14,m23,m24,m34])
            #@show([m12,m13,m14,m23,m24,m34,iSO])
        end
    end
    
    return iSO,mijSO_arr
end


# kmax for ISGL: same for any interaction type.
function kmax_fun(l_complete)
    #kmax_arr = Array{Array}(undef, lastindex(l_complete))
    kmax_dict = Dict{Tuple{Int64, Int64},Int64}()
    
    for i in keys(l_complete)
        l=l_complete[i]
        for m = -l:l
            kmax_dict[l,m] = kmax(l,m)
        end
    end
    
    return kmax_dict
end

function kmax(l,m)
    f = floor(div(l-m,2))
    # somehow this has a problem for m = -l :
    #kmax=Int64(-1/6*(1+f)*(2+f)*(-3*(1+l-m)*(1+m)+f*(3-2*l+6m+3*f)))
    
    # this works:
    kmax = 0
    for j = 0:f
        for s=0:(l-m-2*j)
            for t = 0:j+m
                for u=0:j
                    kmax += 1
                end
            end
        end
    end
    return kmax
end


