# function often used in the examples to print a comparison of two arrays with differences and labels s1,s2
function comparison(num_arr,ref_arr,simax;s1="Numerical", s2="Reference", indexlist=1:simax)
    @printf("%-7s %-15s %-15s %-15s\n", "Index",  s1, s2, "Difference")
    for i in indexlist
        @printf("%-7d %-15.6f %-15.6f %-15.6f\n", i, num_arr[i], ref_arr[i], ref_arr[i] - num_arr[i])
    end
end;

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


# Reconstruct a hermitian matrix from its lower triangle.
# Hermitian(src,:L) leaves the stored diagonal untouched, so the O(eps) imaginary parts that the
# matrix elements pick up numerically would survive and make ishermitian(dest) false. eigen2step
# keys its "the eigenvalues are real" decision on exactly that predicate, so the diagonal is made
# real explicitly here. Mathematically the diagonal of S, T and V is real anyway.
function hermitian_fill!(dest,src)
    dest .= Hermitian(src,:L)
    for i in axes(dest,1)
        dest[i,i] = real(dest[i,i])
    end
    return dest
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
    const NGUARD_THETA = 2
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
