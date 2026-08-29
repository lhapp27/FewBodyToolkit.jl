# types for the different potential functions

"""
    PotentialFunction

Abstract type for potential functions used in few-body calculations.
This type serves as a parent type for more specific potential implementations
"""
abstract type PotentialFunction end

function (v::PotentialFunction)(r)
    return v.f(r)
end


"""
    CentralPotential(f::Function)
A concrete implementation of `PotentialFunction` that represents a central potential.
It takes a function `f` that defines the potential as a function of the radial distance `r`.
"""
struct CentralPotential{F} <: PotentialFunction
    f::F
end


"""
    GaussianPotential(v0::Float64, mu_g::Float64)
A concrete implementation of `PotentialFunction` that represents a Gaussian potential:
```math
V(r) = v_0  e^{-\\mu_g r^2}
```
where `r` is the radial distance.

# Arguments:
- `v0::Float64`: The strength of the potential.
- `mu_g::Float64`: The width parameter of the Gaussian potential.
"""
struct GaussianPotential <: PotentialFunction
    v0::Float64
    mu_g::Float64

    function GaussianPotential(v0::Float64, mu_g::Float64)
        new(v0, mu_g)
    end
end

"""
    function (gp::GaussianPotential)(r::Float64)

Evaluates the Gaussian potential at a given radial distance `r`.
"""
function (gp::GaussianPotential)(r)
    gp.v0 * exp(-gp.mu_g * r^2)
end


"""
    PowerLawPotential(v0::Float64, p::Float64)
A concrete implementation of `PotentialFunction` that represents a power-law potential:
```math
V(r) = v_0  |r|^{p}
```

Note the **absolute value**: the potential is defined via `|r|`, i.e. it is an even
function. This matters only in 1D, where the coordinate can be negative; there,
`PowerLawPotential(v0,p)` describes `v0*|x|^p`. Odd potentials are not supported.
In 2D and 3D, `r >= 0` anyway and `|r| = r`.

Matrix elements are evaluated analytically (closed-form radial integration), so
neither numerical integration nor range-interpolation is needed. Prominent cases
are the Coulomb interaction (`p=-1`) and the harmonic oscillator (`p=2`).

The matrix elements are finite only if `p` is not too negative. The bound depends
on the dimension and on the angular momentum, and is therefore **not** checked
here but in the individual modules:
- `GEM2B`: requires `p > -2*lmax - dim`, e.g. `p > -1` for `lmax=0` in 1D,
  so a pure `1/|x|` is already divergent in 1D.
- `ISGL`: requires `p > -3`.

# Arguments:
- `v0::Float64`: The strength of the potential.
- `p::Float64`: The exponent of the power law.
"""
struct PowerLawPotential <: PotentialFunction
    v0::Float64
    p::Float64
end

# convenience: allow integer exponents, e.g. PowerLawPotential(-1.0,-1)
PowerLawPotential(v0::Real, p::Real) = PowerLawPotential(Float64(v0), Float64(p))

"""
    function (pp::PowerLawPotential)(r)

Evaluates the power-law potential at a given radial distance `r`. Uses `abs(r)`,
see the note on 1D in [`PowerLawPotential`](@ref).
"""
function (pp::PowerLawPotential)(r)
    pp.v0 * abs(r)^pp.p
end



# postponed to future version
#= """
    SpinOrbitPotential(f::Function)
Defines a potential of the type `SpinOrbitPotential`. The function `f(r)` represents the radial part of a spin-orbit interaction
```math
V_{SO}(r) = f(r) \\vec{l} \\cdot \\vec{s}
```
"""
struct SpinOrbitPotential <: PotentialFunction
    f::Function
end =#


"""
    ContactPotential1D(v0::Float64, z0::Float64)
A concrete implementation of `PotentialFunction` that represents a 1D Contact (Dirac) potential:
```math
V(z) = v_0  \\delta(z - z_0)
```
where `z` is the 1D coordinate.

# Arguments:
- `v0::Float64`: The strength of the potential.
- `z0::Float64`: The position of the delta function.
"""
struct ContactPotential1D <: PotentialFunction
    v0::Float64
    z0::Float64

    function ContactPotential1D(v0::Float64, z0::Float64)
        new(v0, z0)
    end
end

"""
    function (gp::ContactPotential1D)(z::Float64)

Evaluates the contact potential at a given 1D coordinate `z`.
"""
function (gp::ContactPotential1D)(z::Float64)
    z == gp.z0 ? Inf : 0.0
end