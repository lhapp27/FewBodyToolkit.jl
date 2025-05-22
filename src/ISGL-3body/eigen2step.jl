# Function for solving the generalized eigenvalue problem H*x = lambda * S*x in two steps
# "Hermitian pencil method"
# Assumptions here: S = real, symmetric matrix -> deleted for CR-basis functions

function eigen2step(e_arr,H, S; threshold::Float64 = 10^-13)

    if typeof(S[1,1]) == Float64 #copy from 2-body version
        dvec,y = eigen(Symmetric(S));
    elseif typeof(S[1,1]) == ComplexF64
        dvec,y = eigen(Hermitian(S));
    end

    dvec_mask,y_mask = cutSmallEV(dvec,y,threshold=threshold)

    #d = diagm(sqrt.(abs.(dvec_mask)));
    #dinv = diagm(1 ./sqrt.(abs.(dvec_mask)));
    #l = y_mask*dinv;
    l = y_mask*diagm(1 ./sqrt.(abs.(dvec_mask)));

    #@show([issymmetric(H),issymmetric(l' * H * l),maximum(l' * H * l .- transpose(l' * H * l))])
    # l'*H*l is almost symmetric... maybe we can use that for speedup?
    
    e3 = eigvals!((l' * H * l)); # Symmetric(l' * H * l) would speed-up quite a bit, but unclear consequences...
    if (typeof(H[1,1]) == Float64 && issymmetric(H)) || (typeof(H[1,1]) == ComplexF64 && ishermitian(H))
        e_arr[1:lastindex(dvec_mask)] .= real.(e3)
    else
        e_arr[1:lastindex(dvec_mask)] .= e3
    end
end



function eigen2step_valvec(e_arr,v_arr,H, S; threshold::Float64 = 10^-13)

    if typeof(S[1,1]) == Float64
        dvec,y = eigen(Symmetric(S));
    elseif typeof(S[1,1]) == ComplexF64
        dvec,y = eigen(Hermitian(S));
    end

    dvec_mask,y_mask = cutSmallEV(dvec,y,threshold=threshold)

    d = diagm(sqrt.(abs.(dvec_mask)));
    #dinv = diagm(1 ./sqrt.(abs.(dvec_mask)));
    #l = y_mask*dinv;
    l = y_mask*diagm(1 ./sqrt.(abs.(dvec_mask)));
    
    e3,v3 = eigen!(l' * H * l); # why was here "Symmetric()" wrapped ?!

    vecs_output = ((y_mask*d)')\v3;
    foreach(normalize!, eachcol(vecs_output));

    if (typeof(H[1,1]) == Float64 && issymmetric(H)) || (typeof(H[1,1]) == ComplexF64 && ishermitian(H))
        e_arr[1:lastindex(dvec_mask)] .= real.(e3)
    else
        e_arr[1:lastindex(dvec_mask)] .= e3
    end
    v_arr[:,1:lastindex(dvec_mask)] .= vecs_output
end

# Cut out too small eigenvalues
function cutSmallEV(dvec,y;threshold = 10^-13)
    mask = dvec/maximum(dvec) .>= threshold
    # apply mask:
    dvec_mask = dvec[mask]
    y_mask = y[:, mask]
    return dvec_mask,y_mask
end



## old definitions!
function eigen2step_old(H::Matrix{Float64}, S::Matrix{Float64}; threshold::Float64 = 10^-13)

    dvec,y = eigen!(Symmetric(S));

    dvec_mask,y_mask = cutSmallEV(dvec,y,threshold=threshold)

    d = diagm(sqrt.(abs.(dvec_mask)));
    dinv = diagm(1 ./sqrt.(abs.(dvec_mask)));
    l = y_mask*dinv;

    #@show([issymmetric(H),issymmetric(l' * H * l),maximum(l' * H * l .- transpose(l' * H * l))])
    # l'*H*l is almost symmetric... maybe we can use that for speedup?
    
    e3 = eigvals!((l' * H * l)::Matrix{Float64});
    #e_arr[1:lastindex(dvec_mask)]
    
    #return e_arr#e3
    return e3
end

function eigen2step_valvec_old(H::Matrix{Float64}, S::Matrix{Float64}; threshold::Float64 = 10^-13)

    dvec,y = eigen!(Symmetric(S));

    dvec_mask,y_mask = cutSmallEV(dvec,y,threshold=threshold)

    d = diagm(sqrt.(abs.(dvec_mask)));
    dinv = diagm(1 ./sqrt.(abs.(dvec_mask)));
    l = y_mask*dinv;
    
    e3,v3 = eigen!((l' * H * l)::Matrix{Float64});
    
    vecs_output = ((y_mask*d)')\v3;

    foreach(normalize!, eachcol(vecs_output));

    #e_arr[1:lastindex(dvec_mask)] .= e3
    #v_arr[:,1:lastindex(dvec_mask)] .= vecs_output
    
    #return e_arr,v_arr#e3, vecs_output
    return e3,vecs_output
end




