# Function for solving the generalized eigenvalue problem H*x = lambda * S*x in two steps
# "Hermitian pencil method"
# Assumptions here: S = real, symmetric matrix -> deleted for CR-basis functions

# S arrives as a plain Matrix, so it has to be wrapped: for real elements Hermitian is the same as
# Symmetric and takes the identical LAPACK route. Wrapping explicitly (rather than letting eigen
# detect the symmetry itself) pins the symmetric/hermitian path, which guarantees the real, sorted
# eigenvalues and the orthonormal eigenvectors that reduce_basis relies on.
_eigen_S(S) = eigen(Hermitian(S))

# real-symmetric or hermitian H (no complex scaling): the projected matrix is hermitian as well. It is
# wrapped explicitly, since the product l'*H*l is not hermitian bit by bit and would otherwise go to
# the general (non-hermitian) eigensolver.
_ishermH(H) = (typeof(H[1,1]) == Float64 && issymmetric(H)) || (typeof(H[1,1]) == ComplexF64 && ishermitian(H))
_wrapH(Hr, herm::Bool) = herm ? Hermitian(Hr) : Hr

function _assign_evals!(e_arr, evals, herm::Bool)
    if herm
        e_arr[1:lastindex(evals)] .= real.(evals)
    else
        e_arr[1:lastindex(evals)] .= evals
    end
end

# projection l'*H*l onto the reduced basis
_project(H, l) = l' * (H * l)

"""
    reduce_basis(S; threshold=1e-13)

Transformation to an orthonormalized, possibly reduced basis: returns the matrix `L` (size `n` x `m`, `m <= n`) with `L' * S * L == I`.

`L` is built from the eigenvectors of the norm-overlap matrix `S`, dropping all eigenvalues smaller than `threshold` times the largest one. This removes the near-linear dependence of the non-orthogonal basis; `m < n` whenever some eigenvalues are cut. Both real-symmetric `S` and hermitian `S` (complex-ranged basis functions) are supported.

Used by [`eigen2step`](@ref) and [`inverse_solve`](@ref); useful on its own to transform a generalized eigenvalue problem `H*x = E*S*x` into the ordinary one `(L'*H*L)*y = E*y`, with `x = L*y`.
"""
reduce_basis(S; threshold::Float64 = 10^-13) = _reduce_basis(_eigen_S(S)...; threshold=threshold)

# scales the kept eigenvectors of S in place (y_mask is a view into y)
function _reduce_basis(dvec, y; threshold::Float64 = 10^-13)
    dvec_mask,y_mask = cutSmallEV(dvec,y,threshold=threshold)
    return rmul!(y_mask, Diagonal(1 ./sqrt.(abs.(dvec_mask))))
end

"""
    eigen2step(e_arr, H, S; threshold=1e-13)

Solves the generalized eigenvalue problem `H*x = E*S*x` and writes the energies into `e_arr`.

Similar to `eigvals(H,S)`, but the basis is first orthonormalized by [`reduce_basis`](@ref), which cuts the small eigenvalues of `S`. The number of energies can therefore be smaller than the size of `H`.
"""
function eigen2step(e_arr,H, S; threshold::Float64 = 10^-13)
    
    herm = _ishermH(H)
    l = reduce_basis(S;threshold=threshold)
    
    e3 = eigvals!(_wrapH(_project(H,l),herm));
    _assign_evals!(e_arr, e3, herm)
end

# Same as eigen2step, but overwrites S and H (energies-only path of the 3-body solvers, where neither
# is needed afterwards): S holds H*l after its eigendecomposition, H the projected matrix.
function eigen2step!(e_arr,H, S; threshold::Float64 = 10^-13)
    
    herm = _ishermH(H)
    l = _reduce_basis(eigen!(Hermitian(S))...; threshold=threshold)
    m = size(l,2)

    if eltype(S) == eltype(H)
        Hr = view(H,1:m,1:m)
        mul!(Hr, l', mul!(view(S,:,1:m), H, l))
    else # real S with complex scaling: mul! of real l' and complex H*l has no BLAS path
        Hr = _project(H,l)
    end

    e3 = eigvals!(_wrapH(Hr,herm));
    _assign_evals!(e_arr, e3, herm)
end


"""
    eigen2step_valvec(e_arr, v_arr, H, S; threshold=1e-13)

Same as [`eigen2step`](@ref), but also writes the coefficient vectors into the columns of `v_arr`, normalized as `x' * S * x == 1`.
"""
function eigen2step_valvec(e_arr,v_arr,H, S; threshold::Float64 = 10^-13)
    
    herm = _ishermH(H)
    l = reduce_basis(S;threshold=threshold)
    
    e3,v3 = eigen!(_wrapH(_project(H,l),herm));
    
    C = l*v3; # back-transformation to the original basis
    
    # normalization of C = vecs_output:
    #foreach(normalize!, eachcol(C)); incorrect normalization for non-orthogonal basis
    M = C' * S * C
    R = cholesky(Hermitian(M)).U
    C .= C * (R \ I)   # since inv(R) = R \ I; now C_norm' * S * C_norm == I (up to numerical precision)
    
    
    _assign_evals!(e_arr, e3, herm)
    v_arr[:,1:size(C,2)] .= C
end

"""
    inverse_solve(T, V, S; target_energy=0.0, threshold=1e-13, return_vectors=false)

Solves the inverse problem: instead of the energies at a given interaction strength, it returns the interaction strengths `lambda` at which a state of energy `target_energy` exists.

For a Hamiltonian which is linear in a strength parameter, `H = T + lambda*V`, the condition that `target_energy` is an eigenvalue of `H` is the generalized eigenvalue problem

    (T - target_energy*S) x = -lambda * V x

in the non-orthogonal basis with norm-overlap `S`. One matrix build and one eigen solve therefore replace a whole scan over the strength with bracketing of the energy. For scans over `target_energy`, see [`inverse_prepare`](@ref).

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
inverse_solve(T, V, S; target_energy = 0.0, threshold::Float64 = 10^-13, return_vectors::Bool = false) =
    inverse_solve(inverse_prepare(T,V,S;threshold=threshold); target_energy=target_energy, return_vectors=return_vectors)

# the parts of the pencil that do not depend on target_energy
struct InversePrepared{TL,TM}
    L::TL    # reduced basis, L'*S*L == I
    LTL::TM  # L'*T*L
    B::TM    # -(L'*V*L)
end

"""
    inverse_prepare(T, V, S; threshold=1e-13)

First stage of [`inverse_solve`](@ref) for scans over the target energy: reduces the basis and projects `T` and `V` once. Solve with `inverse_solve(P; target_energy=0.0, return_vectors=false)`, which then costs a single small eigen solve per energy.

# Example
```julia
P = inverse_prepare(T,V,S; threshold=1e-8)
lambdas = [inverse_solve(P; target_energy=E)[1] for E in -1.0:0.1:-0.1]
```
"""
function inverse_prepare(T, V, S; threshold::Float64 = 10^-13)
    
    # complex scaling leaves T and V complex-symmetric instead of hermitian, and the pencil below
    # would then silently use the wrong matrices
    ishermitian(T) && ishermitian(V) || error("inverse_solve: T and V must be real-symmetric or hermitian. The complex scaling method is not supported.")
    
    L = reduce_basis(S;threshold=threshold)
    
    # The pencil must be reduced with S, not with one of its own matrices: the eigenvalues of T span
    # many orders of magnitude (diffuse Gaussians have tiny kinetic energy), so cutting on T instead
    # would discard exactly the diffuse functions that matter near threshold.
    LTL = Hermitian(_project(T,L))
    B = Hermitian(-_project(V,L)) # may be indefinite; this is why it is not used as the metric
    return InversePrepared(L,LTL,B)
end

function inverse_solve(P::InversePrepared; target_energy = 0.0, return_vectors::Bool = false)
    
    # A = L'*(T - target_energy*S)*L with L'*S*L == I; positive definite below the continuum
    F = cholesky!(P.LTL - target_energy*I; check=false)
    issuccess(F) || error("inverse_solve: T - target_energy*S is not positive definite in the reduced basis. target_energy = $target_energy must lie below the lowest continuum threshold of the basis.")
    
    # B*y = mu*A*y with the positive-definite A as metric, reduced with its Cholesky factor U:
    # (U'^-1 * B * U^-1) z = mu z, y = U^-1 z; lambda = 1/mu
    C = Hermitian(rdiv!(ldiv!(F.U', Matrix(P.B)), F.U))
    if return_vectors
        mu,Z = eigen!(C)
    else
        mu = eigvals!(C)
    end
    
    ind = findall(>(0), mu)                 # mu <= 0 corresponds to lambda <= 0 and is discarded
    ind = ind[sortperm(mu[ind], rev=true)]  # decreasing mu = increasing lambda
    lambdas = 1 ./ mu[ind]
    
    if return_vectors
        Y = F.U \ Z[:,ind]
        foreach(normalize!, eachcol(Y)) # Y'*A*Y == I; unit norm gives x'*S*x == 1
        return lambdas, P.L*Y
    end
    return lambdas
end

### Real basis for complex ranges (3-body solvers, without complex scaling) ###
# With complex ranges, every basis function phi comes with its partner phi* (conjugate ranges, same
# channel), so P*M*P == conj(M) for H and S, with the permutation P = partner. The combinations
# c = (phi + phi*)/sqrt2 and s = i*(phi - phi*)/sqrt2 are real functions: U'*M*U is real symmetric, and
# the eigensolve runs in real arithmetic. c is put at the position of the first function of a pair
# (lo = min(i, partner[i])), s at the second one (hi).

# partner of each basis function. Order within each (l,L,...) channel: (n,N) -> (n-1)*Nmax_eff + N, where
# the second half of the ranges of a complex-ranged coordinate are the conjugates of the first half.
function cr_partners(nbasis, gem_params, complex_ranged_r::Bool, complex_ranged_R::Bool)
    complex_ranged_r || complex_ranged_R || return nothing
    nm = complex_ranged_r ? 2*gem_params.nmax : gem_params.nmax
    Nm = complex_ranged_R ? 2*gem_params.Nmax : gem_params.Nmax
    flip(n, m, cr) = cr ? mod1(n + div(m,2), m) : n
    return map(0:nbasis-1) do i
        b, k = divrem(i, nm*Nm)
        n, N = divrem(k, Nm)
        b*nm*Nm + (flip(n+1, nm, complex_ranged_r) - 1)*Nm + flip(N+1, Nm, complex_ranged_R)
    end
end

# max |M[p,p] - conj(M)| relative to max |M|
function conj_pair_dev(M, p)
    d = 0.0
    @inbounds for j in axes(M,2), i in axes(M,1)
        d = max(d, abs(M[p[i],p[j]] - conj(M[i,j])))
    end
    return d/maximum(abs, M)
end

# U'*M*U from the rows lo of M only (uses P*M*P == conj(M)); real symmetric by construction
function realbasis(M, p)
    n = size(M,1)
    Mr = Matrix{Float64}(undef, n, n)
    @inbounds for b in 1:n, a in b:n
        A = M[min(a,p[a]), min(b,p[b])]
        B = M[min(a,p[a]), max(b,p[b])]
        Mr[a,b] = a < p[a] ? (b < p[b] ? real(A) + real(B) : imag(B) - imag(A)) :
                             (b < p[b] ? imag(A) + imag(B) : real(A) - real(B))
        Mr[b,a] = Mr[a,b]
    end
    return Mr
end

# coefficients in the original basis, x = U*y (rows lo, hi of each pair), in place
function realbasis_back!(X, p)
    @views for i in eachindex(p)
        i < p[i] || continue
        yc, ys = X[i,:], X[p[i],:]
        tmp = (yc .- im .* ys) ./ sqrt(2)
        yc .= (yc .+ im .* ys) ./ sqrt(2)
        ys .= tmp
    end
    return X
end

# Cut out too small eigenvalues. They come sorted ascending, so the kept ones are a contiguous tail: views instead of copies.
function cutSmallEV(dvec,y;threshold = 10^-13)
    k = searchsortedfirst(dvec/maximum(dvec), threshold)
    return view(dvec,k:lastindex(dvec)), view(y,:,k:lastindex(dvec))
end
