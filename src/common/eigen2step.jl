# Function for solving the generalized eigenvalue problem H*x = lambda * S*x in two steps
# "Hermitian pencil method"
# Assumptions here: S = real, symmetric matrix -> deleted for CR-basis functions

# S arrives as a plain Matrix, so it has to be wrapped: for real elements Hermitian is the same as
# Symmetric and takes the identical LAPACK route. Wrapping explicitly (rather than letting eigen
# detect the symmetry itself) pins the symmetric/hermitian path, which guarantees the real, sorted
# eigenvalues and the orthonormal eigenvectors that reduce_basis relies on.
_eigen_S(S) = eigen(Hermitian(S))

function _assign_evals!(e_arr, evals, H)
    if (typeof(H[1,1]) == Float64 && issymmetric(H)) || (typeof(H[1,1]) == ComplexF64 && ishermitian(H))
        e_arr[1:lastindex(evals)] .= real.(evals)
    else
        e_arr[1:lastindex(evals)] .= evals
    end
end

"""
    reduce_basis(S; threshold=1e-13)

Transformation to an orthonormalized, possibly reduced basis: returns the matrix `L` (size `n` x `m`, `m <= n`) with `L' * S * L == I`.

`L` is built from the eigenvectors of the norm-overlap matrix `S`, dropping all eigenvalues smaller than `threshold` times the largest one. This removes the near-linear dependence of the non-orthogonal basis; `m < n` whenever some eigenvalues are cut. Both real-symmetric `S` and hermitian `S` (complex-ranged basis functions) are supported.

Used by [`eigen2step`](@ref) and [`inverse_solve`](@ref); useful on its own to transform a generalized eigenvalue problem `H*x = E*S*x` into the ordinary one `(L'*H*L)*y = E*y`, with `x = L*y`.
"""
function reduce_basis(S; threshold::Float64 = 10^-13)
    dvec,y = _eigen_S(S)
    dvec_mask,y_mask = cutSmallEV(dvec,y,threshold=threshold)
    return y_mask*diagm(1 ./sqrt.(abs.(dvec_mask)))
end

function eigen2step(e_arr,H, S; threshold::Float64 = 10^-13)
    
    l = reduce_basis(S;threshold=threshold)
    
    e3 = eigvals!(l' * H * l);
    _assign_evals!(e_arr, e3, H)
end


function eigen2step_valvec(e_arr,v_arr,H, S; threshold::Float64 = 10^-13)
    
    l = reduce_basis(S;threshold=threshold)
    
    e3,v3 = eigen!(l' * H * l);
    
    C = l*v3; # back-transformation to the original basis
    
    # normalization of C = vecs_output:
    #foreach(normalize!, eachcol(C)); incorrect normalization for non-orthogonal basis
    M = C' * S * C
    R = cholesky(Hermitian(M)).U
    C .= C * (R \ I)   # since inv(R) = R \ I; now C_norm' * S * C_norm == I (up to numerical precision)
    
    
    _assign_evals!(e_arr, e3, H)
    v_arr[:,1:size(C,2)] .= C
end

"""
    inverse_solve(T, V, S; target_energy=0.0, threshold=1e-13, return_vectors=false)

Solves the inverse problem: instead of the energies at a given interaction strength, it returns the interaction strengths `lambda` at which a state of energy `target_energy` exists.

For a Hamiltonian which is linear in a strength parameter, `H = T + lambda*V`, the condition that `target_energy` is an eigenvalue of `H` is the generalized eigenvalue problem

    (T - target_energy*S) x = -lambda * V x

in the non-orthogonal basis with norm-overlap `S`. One matrix build and one eigen solve therefore replace a whole scan over the strength with bracketing of the energy.

# Arguments
- `T`: kinetic energy matrix (*not* the Hamiltonian `T+V`), as returned by `GEM2B_matrices`, `GEM3B1D_matrices`, `ISGL_matrices`.
- `V`: interaction matrix at unit strength. `lambda` is the factor multiplying *all* interactions, so pass the potential at unit strength.
- `S`: norm-overlap matrix.

# Keywords
- `target_energy=0.0`: the energy that should be reached. Must lie below the lowest continuum threshold of the basis, such that `T - target_energy*S` is positive definite; otherwise an error is raised. `target_energy=0.0` gives the critical strengths at which a state becomes bound.
- `threshold=1e-13`: cut-off for the eigenvalues of `S`, see [`reduce_basis`](@ref). Results near `target_energy=0.0` are sensitive to this value: too small a cut-off can produce spurious, nearly linearly dependent solutions, too large a one cuts the diffuse basis functions needed close to threshold.
- `return_vectors=false`: whether to also return the coefficient vectors of the corresponding states.

Real-symmetric and hermitian matrices are supported, i.e. both the usual and the complex-ranged basis functions. The complex scaling method is not: it makes `T` and `V` complex-symmetric rather than hermitian, and is rejected.

# Returns
- `lambdas`: the positive strengths, in ascending order. The first entry is the smallest strength at which the system supports a state of energy `target_energy`. Negative strengths are discarded; if they are of interest, call the function with `-V` instead.
- `vectors`: (optional) the corresponding coefficient vectors, normalized as `x' * S * x == 1`.

# Example
```julia
T,V,S = GEM2B_matrices(phys_params, num_params)
lambdas = inverse_solve(T,V,S; target_energy=0.0, threshold=1e-8) # critical strengths for binding
```
"""
function inverse_solve(T, V, S; target_energy = 0.0, threshold::Float64 = 10^-13, return_vectors::Bool = false)
    
    # complex scaling leaves T and V complex-symmetric instead of hermitian, and the pencil below
    # would then silently use the wrong matrices
    ishermitian(T) && ishermitian(V) || error("inverse_solve: T and V must be real-symmetric or hermitian. The complex scaling method is not supported.")
    
    L = reduce_basis(S;threshold=threshold)
    
    # The pencil must be reduced with S, not with one of its own matrices: the eigenvalues of T span
    # many orders of magnitude (diffuse Gaussians have tiny kinetic energy), so cutting on T instead
    # would discard exactly the diffuse functions that matter near threshold.
    A = Hermitian(L' * (T .- target_energy .* S) * L) # positive definite below the continuum
    B = Hermitian(L' * (-V) * L)                      # may be indefinite; this is why it is not used as the metric
    
    isposdef(A) || error("inverse_solve: T - target_energy*S is not positive definite in the reduced basis. target_energy = $target_energy must lie below the lowest continuum threshold of the basis.")
    
    # solved as B*y = mu*A*y with the positive-definite A as metric (Cholesky); lambda = 1/mu
    if return_vectors
        mu,Y = eigen(B,A)
    else
        mu = eigvals(B,A)
    end
    
    ind = findall(>(0), mu)                 # mu <= 0 corresponds to lambda <= 0 and is discarded
    ind = ind[sortperm(mu[ind], rev=true)]  # decreasing mu = increasing lambda
    lambdas = 1 ./ mu[ind]
    
    if return_vectors
        Y = Y[:,ind]
        foreach(normalize!, eachcol(Y)) # eigen normalizes as Y'*A*Y == I; unit norm gives x'*S*x == 1
        return lambdas, L*Y
    end
    return lambdas
end

# Cut out too small eigenvalues
function cutSmallEV(dvec,y;threshold = 10^-13)
    mask = dvec/maximum(dvec) .>= threshold
    # apply mask:
    dvec_mask = dvec[mask]
    y_mask = y[:, mask]
    return dvec_mask,y_mask
end
