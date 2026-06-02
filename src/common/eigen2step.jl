# Function for solving the generalized eigenvalue problem H*x = lambda * S*x in two steps
# "Hermitian pencil method"
# Assumptions here: S = real, symmetric matrix -> deleted for CR-basis functions

function _eigen_S(S)
    if typeof(S[1,1]) == Float64
        return eigen(Symmetric(S))
    elseif typeof(S[1,1]) == ComplexF64
        return eigen(Hermitian(S))
    else
        error("Unsupported matrix element type in S: $(typeof(S[1,1]))")
    end
end

function _assign_evals!(e_arr, evals, H)
    if (typeof(H[1,1]) == Float64 && issymmetric(H)) || (typeof(H[1,1]) == ComplexF64 && ishermitian(H))
        e_arr[1:lastindex(evals)] .= real.(evals)
    else
        e_arr[1:lastindex(evals)] .= evals
    end
end

function eigen2step(e_arr,H, S; threshold::Float64 = 10^-13)
    
    dvec,y = _eigen_S(S)
    
    dvec_mask,y_mask = cutSmallEV(dvec,y,threshold=threshold)
    
    d = diagm(sqrt.(abs.(dvec_mask)));
    dinv = diagm(1 ./sqrt.(abs.(dvec_mask)));
    l = y_mask*dinv;
    
    e3 = eigvals!(l' * H * l);
    _assign_evals!(e_arr, e3, H)
end


function eigen2step_valvec(e_arr,v_arr,H, S; threshold::Float64 = 10^-13)
    
    dvec,y = _eigen_S(S)
    
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
    
    
    _assign_evals!(e_arr, e3, H)
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