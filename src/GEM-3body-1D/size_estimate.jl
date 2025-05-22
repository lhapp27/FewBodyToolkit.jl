# functions to calculate the parameters that determine the necessary sizes of arrays within the ISGL program

# contains:
# - size_estimate that bundles everything:
# - abc_size for box-size of TVS Matrices
# - lLcoupl for coupling of angular momenta
# - imax_fun for imax in ISGL
# - kmax_fun for kmax in ISGL
# - alloc_fun for array allocations

struct SizeParams
    abvals_arr::Vector{Vector{Int64}}
    cvals::Vector{Int64}
    groupindex_arr::Vector{Int64}
    nboxes::Int64
    abI::Int64
    factor_bf::Int64
    box_size_arr::Vector{Int64}
    nbasis_total::Int64
    starts::Vector{Int64}
    ends::Vector{Int64}
    bvalsdiag::Vector{Vector{Int64}} # actually I think it could be simplified to Vector{Int64} 
    lL_arr::Vector{Vector{Tuple{Int64, Int64}}}
    lL_complete::Vector{Tuple{Int64, Int64}}
    l_complete::Vector{Int64}
    nlL::Int64
    nl::Int64
    maxlmax::Int64
    maxobs::Int64
end

function size_estimate(phys_params,num_params,observ_params,gaussopt,coulopt)
    
    # input interpretation:
    (;mass_arr,svals,vint_arr,parity) = phys_params
    (;lmax,Lmax,gem_params,theta_csm,omega_cr,mu0,c_shoulder,kmax_interpol,lmin,Lmin) = num_params
    (;nmax,Nmax,r1,rnmax,R1,RNmax) = gem_params
    (;stateindices,centobs_arr,R2_arr) = observ_params # unnecessary for 1D?
    
    # size of Hamiltonian matrix due to symmetries and interactions
    gaussbool_arr = [gaussopt[i][1][1] for i in 1:3]
    coulbool_arr = [coulopt[i][1][1] for i in 1:3]
    cvals = findall(isempty.(vint_arr) .==0 .|| gaussbool_arr .== 1.0 .|| coulbool_arr .== 1.0) # consider only values for c where there are interactions (or Gaussian interaction!, or Coulomb!)
    abvals_arr,groupindex_arr,nboxes,abI,factor_bf = abc_size(cvals,svals)
    
    # coupling of angular momenta
    lL_arr,lL_complete,l_complete = lLcoupl(parity,lmax,Lmax,cvals,svals,lmin,Lmin)
    
    @show(parity,lL_arr,lL_complete,l_complete)

    # box sizes = number of basis functions in each box; starts,ends = indices for boxes within big matrix; bvalsdiag = simplification for identical particles
    box_size_arr,nbasis_total = boxsize_fun(abvals_arr,groupindex_arr,lL_arr,nmax,Nmax)
    starts = [1; cumsum(box_size_arr[1:end-1]) .+ 1]
    ends = cumsum(box_size_arr)
    #bvalsdiag = [[abvals_arr[boxR][1]] for boxR in groupindex_arr]
    bvalsdiag = [Int64[], Int64[], Int64[]]
    for boxC in groupindex_arr
        bvalsdiag[boxC] = [abvals_arr[boxC][1]]
    end

    @show(box_size_arr,nbasis_total,abvals_arr)

    # not needed for 1D:
#=     # imax for ISGL
    imax_dict,mij_arr_dict = imax_fun(lL_complete)
    
    # kmax for ISGL
    kmax_dict = kmax_fun(l_complete) =#


    ## max-values:
    nlL = lastindex(lL_complete) # former: lLmax'
    nl = lastindex(l_complete) # former lmax'
    maxlmax=max(lmax,Lmax)
    #maximax=findmax(imax_dict)[1]
    #maxkmax=findmax(kmax_dict)[1] # [1] is the value; [2] would be the key
    maxobs = maximum(lastindex.(centobs_arr)) # max number of observables

    # Constructing Struct (collective data structure size_params with all the size parameters)
    size_params = SizeParams(abvals_arr,cvals,groupindex_arr,nboxes,abI,factor_bf,box_size_arr,nbasis_total,starts,ends,bvalsdiag,lL_arr,lL_complete,l_complete,nlL,nl,maxlmax,maxobs)
    
    return size_params
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


# function to determine list of (l,L)-tuples to consider, based of parity,svals ## in 1D: parity = früheres Piz?
function lLcoupl(parity,lmax,Lmax,cvals,svals,lmin,Lmin) # 1D: lmax=Lmax=1 immer? nope, see coulomb!
    
    lL_arr = ([Vector{Tuple{Int64, Int64}}() for _ in 1:3])#[[],[],[]] # one list of (l,L)-tuples for each Jacobi-set; only the relevant ones will be filled!
    
    #println("type,size,arr",typeof(lL_arr),", ",size(lL_arr),", ",lL_arr)

    for c in cvals # we only need  (l,L) list for the relevant jacobi sets
        for l = lmin:lmax
            for L = Lmin:Lmax

                if parity == +1 || parity == -1 # select allowed l,L from parity only if parity is +-1
                    (-1)^(l+L) != parity && continue # allowed combination of l,L from global parity
                end

                    #if (-1)^(l+L) == parity # allowed combination of l,L from global parity
                        if svals[mod(c+1,3)+1] == svals[(mod(c,3)+1)] == "b"  # if Jacobi-Set c contains relative coordinate between two identical bosons -> only even l
                            (mod(l,2) == 0) && push!(lL_arr[c],(l,L))
                        elseif svals[(mod(c+1,3)+1)] == svals[(mod(c,3)+1)] == "f"# if Jacobi-Set c contains relative coordinate between two identical bosons -> only odd l
                            (mod(l,2) == 1) && push!(lL_arr[c],(l,L))
                        else
                            push!(lL_arr[c],(l,L))
                        end 
                    #end

            end
        end
    end
    
    nlL = maximum(lastindex.(lL_arr)) # maximum amount of (l,L)-pairs # nope this is always the number of jacobi-sets?
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
    
    return lL_arr,lL_complete,l_complete # maybe provide also maximum(l_complete)? this should be maxlmax?
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


function imax_fun(lL_complete)
    imax_dict = Dict{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}}, Int64}()
    mij_arr_dict = Dict{Tuple{Tuple{Int64, Int64}, Tuple{Int64, Int64}},Vector{SVector{6, Int64}}}()
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


