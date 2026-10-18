# function often used in the examples to print a comparison of two arrays with differences and labels s1,s2
function comparison(num_arr,ref_arr,simax;s1="Numerical", s2="Reference", indexlist=1:simax)
    @printf("%-7s %-15s %-15s %-15s\n", "Index",  s1, s2, "Difference")
    for i in indexlist
        @printf("%-7d %-15.6f %-15.6f %-15.6f\n", i, num_arr[i], ref_arr[i], ref_arr[i] - num_arr[i])
    end
end;

"""
    quadgk_scaled(f, domain::Tuple; segbuf=nothing)

`quadgk` over `domain` (e.g. `(-Inf,0,Inf)`), with a fallback for integrals that vanish or nearly cancel
(e.g. odd moments of even potentials, needed for `parity=0` in 1D). With its purely relative default
tolerance, `quadgk` cannot converge on those and runs until `maxevals` (10^7 evaluations); the fallback
measures the error relative to `∫|f|` instead. Integrals that converge normally are computed as before.
"""
function quadgk_scaled(f, domain::Tuple; segbuf=nothing)
    I, E = quadgk(f, domain...; maxevals=2000, segbuf=segbuf)
    E <= sqrt(eps())*abs(I) && return I
    return Base.invokelatest(quadgk_scaled_fallback, f, domain, segbuf) # runtime dispatch: compiled only when needed
end

function quadgk_scaled_fallback(f, domain, segbuf)
    rtol = sqrt(eps())
    scale = quadgk(x -> abs(f(x)), domain...; rtol=1e-2)[1] # coarse estimate of ∫|f|
    return quadgk(f, domain...; rtol=rtol, atol=rtol*scale, segbuf=segbuf)[1]
end


"""
    parse_complex_ranged(complex_ranged) -> (complex_ranged_r, complex_ranged_R)

Translates the user-facing `complex_ranged` option into two booleans, one per Jacobi coordinate.

Accepted values are the symbols `:none`, `:r`, `:R`, `:both`, and a `Bool` (`true` meaning `:both`,
`false` meaning `:none`) for consistency with `GEM2B_solve`, whose basis has only one coordinate.

Used by both three-body modules (`ISGL` and `GEM3B1D`).
"""
function parse_complex_ranged(complex_ranged)
    complex_ranged isa Bool && return (complex_ranged, complex_ranged)
    complex_ranged === :none && return (false, false)
    complex_ranged === :r    && return (true,  false)
    complex_ranged === :R    && return (false, true)
    complex_ranged === :both && return (true,  true)
    error("complex_ranged must be one of :none, :r, :R, :both (or a Bool); got $(repr(complex_ranged))")
end


# Complete a matrix in place from its lower triangle: hermitian (conjugated mirror) or symmetric.
# In the hermitian case the diagonal is made real explicitly: the O(eps) imaginary parts that the
# matrix elements pick up numerically would otherwise make ishermitian(A) false. eigen2step
# keys its "the eigenvalues are real" decision on exactly that predicate. Mathematically the
# diagonal of S, T and V is real anyway.
function hermitian_fill!(A,hermitian::Bool=true)
    @inbounds for j in axes(A,2)
        hermitian && (A[j,j] = real(A[j,j]))
        for i in j+1:lastindex(A,1)
            A[j,i] = hermitian ? conj(A[i,j]) : A[i,j]
        end
    end
    return A
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

# partner for the fill loop: P*M*P == conj(M) means that of the (at most two) related elements of the lower
# triangle only one has to be computed (see pair_position). It holds for real potentials, but not e.g. for
# identical particles with different interactions on their equivalent pairs; so it is checked on the first
# column, whose elements are all computed, against the partner column. Otherwise nothing (full fill).
function fill_partners(element, nbasis, gem_params, complex_ranged_r::Bool, complex_ranged_R::Bool)
    p = cr_partners(nbasis, gem_params, complex_ranged_r, complex_ranged_R)
    isnothing(p) && return nothing
    d, m = zeros(3), zeros(3)
    for i in 1:nbasis
        pr, pc, cj = pair_position(i,1,p)
        a, b = element(i,1), element(pr,pc)
        d .= max.(d, abs.(b .- (cj ? conj.(a) : a))); m .= max.(m, abs.(a))
    end
    return all(d .<= 1e-10 .* m) ? p : nothing
end

# position of the element related to M[rowi,coli] in the lower triangle, and whether it is the conjugate
function pair_position(rowi, coli, p)
    pr, pc = p[rowi], p[coli]
    return pr >= pc ? (pr, pc, true) : (pc, pr, false)
end

# one element of S and T (H = T + V for add_V, else T and V), conjugated for cj
function store_element!(S,T,V,r,c,sel,tel,vel,add_V::Bool,cj::Bool=false)
    f(x) = cj ? conj(x) : x
    S[r,c] = f(sel)
    if add_V
        T[r,c] = f(tel + vel)
    else
        T[r,c] = f(tel)
        V[r,c] = f(vel)
    end
end

# fill_columns(cols) for all columns: serially, or with parallel in nthreads() tasks with interleaved columns
# (balances the triangle). Every element is written once, so the tasks never write to the same position.
function fill_parallel(fill_columns, nbasis, parallel::Bool)
    parallel || return fill_columns(1:nbasis)
    nt = Threads.nthreads()
    @sync for k in 1:nt
        Threads.@spawn fill_columns(k:nt:nbasis)
    end
end


"""
    theta_mesh(complex_ranged, complex_scaling, complex_range_freq, kmax_theta)

Nodes of the angular mesh of the range-interpolation, i.e. the second dimension of the sector over
which the radial integrals are tabulated.

With complex ranges `nu*(1 +/- i*omega)` the effective Gaussian range at which the interpolant is
looked up (the Schur complement `etaprc`/`etaeff` of the range matrix) is complex, but confined to
the sector `Re > 0`, `|arg| <= atan(omega)`. In `(log|alpha|, theta)` that sector is a rectangle,
so the tensor-product cubic spline simply gains a second dimension.

- Real ranges: a single node at `theta = 0`; the interpolant stays one-dimensional in `log(alpha)`.
- Complex ranges: `kmax_theta` nodes over `[0, atan(omega)]`. The radial integral of a real potential
  obeys `v(conj(alpha)) = conj(v(alpha))` (Schwarz reflection), so the negative half is obtained by
  conjugation and need not be tabulated.
- Complex ranges together with complex scaling: the rotated potential is no longer real and the
  reflection is lost, so the mesh covers the full `[-atan(omega), atan(omega)]` with `2*kmax_theta-1`
  nodes, i.e. the same spacing.

`NGUARD_THETA` further nodes are appended at each end, at the same spacing. Both bounds of the
sector are attained exactly by `arg(etaeff)`, and a cubic spline is markedly less accurate at the
edge of its domain than inside it (measured: an order of magnitude, at unchanged node spacing); the
guard nodes turn the whole used interval into interior. They also remove any need to worry about
roundoff pushing `abs(angle(etaeff))` a few eps past `atan(omega)`, where the interpolants, which
throw on extrapolation, would otherwise error.
"""

function theta_mesh(complex_ranged::Bool,complex_scaling::Bool,complex_range_freq,kmax_theta)
    NGUARD_THETA = 2
    !complex_ranged && return [0.0]
    theta_max = atan(complex_range_freq)
    theta_min = complex_scaling ? -theta_max : 0.0
    nin = complex_scaling ? 2*kmax_theta-1 : kmax_theta
    # omega = 0 collapses the sector to the single ray theta = 0. Keep a non-degenerate (if minute)
    # mesh, so that the interpolant stays well-defined; every lookup then lands on the theta = 0 node.
    h = max((theta_max-theta_min)/(nin-1), sqrt(eps()))
    return collect(range(theta_min-NGUARD_THETA*h, theta_max+NGUARD_THETA*h, nin+2*NGUARD_THETA))
end


# Look up an interpolated radial quantity at the effective Gaussian range eta.
# Real ranges: the interpolant is one-dimensional in log(eta).
# Complex ranges: it is two-dimensional over the sector (log|eta|, arg eta). With `reflect` (no
# complex scaling, see theta_mesh) only arg >= 0 is tabulated and the lower half is reached by
# conjugation.
interpol_lookup(w,eta::Real,reflect::Bool) = w(log(eta))
function interpol_lookup(w,eta::Complex,reflect::Bool)
    theta = angle(eta)
    if reflect && theta < 0
        return conj(w(log(abs(eta)),-theta))
    end
    return w(log(abs(eta)),theta)
end


"""
    check_cr_csm_sector(complex_ranged, complex_scaling, complex_range_freq, complex_scaling_angle)

Guards the one combination in which complex ranges and complex scaling do not simply coexist.

Complex scaling is applied here by rotating the integration contour back onto the real axis, so that
a user-supplied potential is only ever called with a real argument. The price is that the rotation
lands on the Gaussian instead: the radial integrals are taken at `alpha*csmfac^2`, with
`csmfac = exp(-i*theta_csm)` and `alpha` anywhere in the sector `|arg(alpha)| <= atan(omega)`. The
factor `exp(-alpha*csmfac^2*r^2)` decays only while `atan(omega) + 2*theta_csm < pi/2`, and the
corner of the sector where this is tightest (`arg(alpha) = -atan(omega)`) is attained exactly.

This is a limit of the representation, not of the matrix element: the same quantity written without
the contour rotation, `int v(r*exp(i*theta_csm)) r^n exp(-alpha*r^2) dr`, converges throughout (both
forms agree to full precision wherever the rotated one is computable). Lifting the limit would mean
evaluating the potential at complex arguments, which not every user-supplied function supports.

The tail of the potential does not rescue the rotated form. For a slowly decaying potential the
integral genuinely diverges; for a fast decaying one it converges mathematically but still fails
numerically, since `v(r)` and `exp(-alpha*csmfac^2*r^2)` are separate factors and underflow times
overflow gives NaN. Either way `quadgk` reports an unrelated-looking `DomainError`, so the
combination is caught here instead.
"""
function check_cr_csm_sector(complex_ranged::Bool,complex_scaling::Bool,complex_range_freq,complex_scaling_angle)
    (complex_ranged && complex_scaling) || return nothing
    total = atan(complex_range_freq) + 2*complex_scaling_angle*pi/180
    if total >= pi/2
        error("Complex ranges and complex scaling combined: the radial integrals of an interpolated "*
              "potential are evaluated on a contour rotated back to the real axis, which requires "*
              "atan(complex_range_freq) + 2*complex_scaling_angle < 90 deg, but "*
              "atan($(complex_range_freq)) + 2*$(complex_scaling_angle) deg = $(round(total*180/pi,digits=1)) deg. "*
              "Reduce complex_range_freq or complex_scaling_angle. (Analytically treated potentials, "*
              "e.g. GaussianPotential and PowerLawPotential, are never integrated and are not affected.)")
    end
    return nothing
end

# Lookup table for keys made of non-negative integers and half-integers (angular momenta, spins, particle
# indices; all modules). Each key component k is stored at index 2k+1 of a dense array, so that t[k...] is
# an array access instead of a Dict lookup (hashing) in the matrix-element loops. The array grows on
# assignment; entries never assigned are NaN.
mutable struct HalfIntTable{N}
    vals::Array{Float64,N}
end
HalfIntTable{N}() where {N} = HalfIntTable{N}(fill(NaN, ntuple(_ -> 0, N)))
_hidx(k) = Int(2*k) + 1
Base.getindex(t::HalfIntTable{N}, k::Vararg{Real,N}) where {N} = t.vals[map(_hidx, k)...]
function Base.setindex!(t::HalfIntTable{N}, v, k::Vararg{Real,N}) where {N}
    i = map(_hidx, k)
    if any(i .> size(t.vals))
        vals = fill(NaN, max.(i, size(t.vals)))
        vals[CartesianIndices(t.vals)] .= t.vals
        t.vals = vals
    end
    t.vals[i...] = v
    return t
end

# gamma function at the half-integers 0.5, 1.0, ..., ende+0.5 (all modules)
function precompute_gamma!(gamma_dict::HalfIntTable{1},ende)
    for x in 0.5:0.5:ende+0.5
        gamma_dict[x] = gamma(x)
    end
    return gamma_dict
end
