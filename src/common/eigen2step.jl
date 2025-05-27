# Function for solving the generalized eigenvalue problem H*x = lambda * S*x in two steps
# "Hermitian pencil method"
# Assumptions here: S = real, symmetric matrix -> deleted for CR-basis functions

function eigen2step(e_arr,H, S; threshold::Float64 = 10^-13)
    
    if typeof(S[1,1]) == Float64
        dvec,y = eigen!(Symmetric(S)); #careful, this overwrites S
    elseif typeof(S[1,1]) == ComplexF64
        dvec,y = eigen!(Hermitian(S));
    end
    
    dvec_mask,y_mask = cutSmallEV(dvec,y,threshold=threshold)
    
    d = diagm(sqrt.(abs.(dvec_mask)));
    dinv = diagm(1 ./sqrt.(abs.(dvec_mask)));
    l = y_mask*dinv;
    
    #@show([issymmetric(H),issymmetric(l' * H * l),maximum(l' * H * l .- transpose(l' * H * l))])
    # l'*H*l is almost symmetric... maybe we can use that for speedup?
    
    e3 = eigvals!(l' * H * l);
    if (typeof(H[1,1]) == Float64 && issymmetric(H)) || (typeof(H[1,1]) == ComplexF64 && ishermitian(H))
        e_arr[1:lastindex(dvec_mask)] .= real.(e3)
    else
        e_arr[1:lastindex(dvec_mask)] .= e3
    end
end


# quite some code-redundancy. possible simplification via multiple dispatch?
function eigen2step_valvec(e_arr,v_arr,H, S; threshold::Float64 = 10^-13)
    
    if typeof(S[1,1]) == Float64
        dvec,y = eigen(Symmetric(S));
    elseif typeof(S[1,1]) == ComplexF64
        dvec,y = eigen(Hermitian(S));
    end
    
    dvec_mask,y_mask = cutSmallEV(dvec,y,threshold=threshold)
    
    d = diagm(sqrt.(abs.(dvec_mask)));
    dinv = diagm(1 ./sqrt.(abs.(dvec_mask)));
    l = y_mask*dinv;
    
    e3,v3 = eigen!(l' * H * l);
    
    C = ((y_mask*d)')\v3;
    
    # normalization of C = vecs_output:
    #foreach(normalize!, eachcol(C)); incorrect normalization for non-orthogonal basis
    M = C' * S * C
    R = cholesky(Hermitian(M)).U
    C .= C * (R \ I)   # since inv(R) = R \ I; now C_norm' * S * C_norm == I (up to numerical precision)
    
    
    if (typeof(H[1,1]) == Float64 && issymmetric(H)) || (typeof(H[1,1]) == ComplexF64 && ishermitian(H))
        e_arr[1:lastindex(dvec_mask)] .= real.(e3)
    else
        e_arr[1:lastindex(dvec_mask)] .= e3
    end
    v_arr[:,1:lastindex(dvec_mask)] .= C
    #@show(lastindex(dvec_mask))
end

# Cut out too small eigenvalues
function cutSmallEV(dvec,y;threshold = 10^-13)
    mask = dvec/maximum(dvec) .>= threshold
    # apply mask:
    dvec_mask = dvec[mask]
    y_mask = y[:, mask]
    return dvec_mask,y_mask
end