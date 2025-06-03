# functions to calculate the parameters that determine the necessary sizes of arrays to be preallocated

struct SizeParams
    abvals_arr::Vector{Vector{Int64}}
    cvals::Vector{Int64}
    central_indices::Vector{Vector{Int}}
    gauss_indices::Vector{Vector{Int}}
    gaussopt_arr::Vector{Vector{Tuple{Float64,Float64}}}
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
    bvalsdiag::Vector{Vector{Int64}}
    lL_arr::Vector{Vector{Tuple{Int64, Int64}}}
    lL_complete::Vector{Tuple{Int64, Int64}}
    l_complete::Vector{Int64}
    nlL::Int64
    nl::Int64
    maxlmax::Int64
    maxobs::Int64
end

function size_estimate(phys_params,num_params,observ_params)
    
    # input interpretation:
    (;svals,vint_arr,parity) = phys_params #13diff
    (;lmax,Lmax,gem_params,theta_csm,omega_cr,kmax_interpol,lmin,Lmin) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;stateindices,centobs_arr,R2_arr) = observ_params # no observables for 1D atm
    
    # size of Hamiltonian matrix due to symmetries and interactions
    cvals = findall(isempty.(vint_arr) .==0) # consider only values for c where there are interactions (any type)
    # number of interactions per Jacobi-set c:
    central_indices, gauss_indices, gaussopt_arr, nint_arr, nintmax = index_interaction_types(vint_arr)
    
    # box sizes, indices and factors for symmetrization
    abvals_arr,groupindex_arr,nboxes,abI,factor_bf = abc_size(cvals,svals)
    
    # coupling of angular momenta
    lL_arr,lL_complete,l_complete = lLcoupl(parity,lmax,Lmax,cvals,svals,lmin,Lmin)

    # box sizes = number of basis functions in each box; starts,ends = indices for boxes within big matrix; bvalsdiag = simplification for identical particles
    box_size_arr,nbasis_total = boxsize_fun(abvals_arr,groupindex_arr,lL_arr,nmax,Nmax)
    starts = [1; cumsum(box_size_arr[1:end-1]) .+ 1]
    ends = cumsum(box_size_arr)
    #bvalsdiag = [[abvals_arr[boxR][1]] for boxR in groupindex_arr]
    bvalsdiag = [Int64[], Int64[], Int64[]]
    for boxC in groupindex_arr
        bvalsdiag[boxC] = [abvals_arr[boxC][1]]
    end

    ## max-values:
    nlL = lastindex(lL_complete) # former: lLmax'
    nl = lastindex(l_complete) # former lmax'
    maxlmax=max(lmax,Lmax)
    maxobs = maximum(lastindex.(centobs_arr)) # max number of observables

    # Constructing Struct (collective data structure size_params with all the size parameters)
    size_params = SizeParams(abvals_arr,cvals,central_indices,gauss_indices,gaussopt_arr,nint_arr,nintmax,groupindex_arr,nboxes,abI,factor_bf,box_size_arr,nbasis_total,starts,ends,bvalsdiag,lL_arr,lL_complete,l_complete,nlL,nl,maxlmax,maxobs) #13diff

    ## overall not so many differences between 1D and 3D. --> possible to reduce code via additional argument dim?
    return size_params
end

# such that in fillTVS we can use later a sum over the indices of the interaction types. Ideally we would be able to dispatch directly in fillTVS, then we would not need this function. Problem here: needs to be updated, whenever a new interaction type is added.
function index_interaction_types(vint_arr)
    central_indices = [Int[] for _ in 1:3]
    gauss_indices = [Int[] for _ in 1:3]

    gaussopt_arr = [Tuple{Float64,Float64}[] for _ in 1:3]
    nint_arr = zeros(Int64,3)

    pushindexpotentialtype!(v::Function, central_indices, gauss_indices, i) = push!(central_indices, i) # treat function as a central potential
    pushindexpotentialtype!(v::CentralPotential, central_indices, gauss_indices, i) = push!(central_indices, i)
    pushindexpotentialtype!(v::GaussianPotential, central_indices, gauss_indices, i) = push!(gauss_indices, i)    

    for c in 1:3
        for (i, v) in enumerate(vint_arr[c])
            pushindexpotentialtype!(v, central_indices[c], gauss_indices[c], i)

            if i in gauss_indices[c] # if this is a Gaussian potential
                push!(gaussopt_arr[c], (v.v0, v.mu_g)) # store the parameters of the Gaussian potential
            else
                push!(gaussopt_arr[c], (NaN, NaN)) # if not a Gaussian potential, store NaN
            end

        end
        nint_arr[c] = lastindex(vint_arr[c])
    end

    nintmax = maximum(nint_arr)

    return central_indices, gauss_indices, gaussopt_arr, nint_arr, nintmax
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
    if lastindex(uniq) == 2 ## what happens for 3 identical fermions?
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


# function to determine list of (l,L)-tuples to consider, based on parity,svals
function lLcoupl(parity,lmax,Lmax,cvals,svals,lmin,Lmin)
    
    lL_arr = ([Vector{Tuple{Int64, Int64}}() for _ in 1:3])#[[],[],[]] # one list of (l,L)-tuples for each Jacobi-set; only the relevant ones will be filled!

    for c in cvals # we only need  (l,L) list for the relevant jacobi sets
        for l = lmin:lmax
            for L = Lmin:Lmax

                if parity == +1 || parity == -1 # select allowed l,L from parity only if parity is +-1
                    (-1)^(l+L) != parity && continue # allowed combination of l,L from global parity
                elseif parity == 0
                    # do nothing, we use 0 as indicator for parity violation
                else
                    error("Error in lLcoupl: parity=$parity must be +1, -1 or 0.")
                end

                if svals[mod(c+1,3)+1] == svals[(mod(c,3)+1)] == "b"  # if Jacobi-Set c contains relative coordinate between two identical bosons -> only even l
                    (mod(l,2) == 0) && push!(lL_arr[c],(l,L))
                elseif svals[(mod(c+1,3)+1)] == svals[(mod(c,3)+1)] == "f"# if Jacobi-Set c contains relative coordinate between two identical bosons -> only odd l
                    (mod(l,2) == 1) && push!(lL_arr[c],(l,L))
                else
                    push!(lL_arr[c],(l,L))
                end

            end
        end
    end
    
    nlL = maximum(lastindex.(lL_arr))
    lL_complete_index = argmax(lastindex.(lL_arr))
    lL_complete = lL_arr[lL_complete_index]
    
    # most complete list of l or L values necessary in the program
    temp_list = Int[]
    for c in cvals
        for i in eachindex(lL_arr[c])
            push!(temp_list,lL_arr[c][i][1],lL_arr[c][i][2])
        end
    end
    
    l_complete = unique(temp_list)
    if lastindex(lL_complete) == 0
        error("Error in size_estimate: No valid (l,L)-couplings found.")
        return
    end
    
    return lL_arr,lL_complete,l_complete
end


## function to determine the size (= number of basis functions) in each box
function boxsize_fun(abvals_arr,groupindex_arr,lL_arr,nmax,Nmax)
    box_size_arr = zeros(Int64,3)
    for grp in groupindex_arr
        box_size_arr[grp] = Int64(lastindex(lL_arr[abvals_arr[grp][1]])*nmax*Nmax) # number of combinations of (l,L) values in each box, times number of basis functions
    end
    nbasis_total = Int64(sum(box_size_arr))
    return box_size_arr,nbasis_total
end



