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
V(r) = v_0  r^{p}
```
where `r` is the radial distance.

Matrix elements are treated analytically (no numerical integration and no
range-interpolation), since the radial integral is closed-form for
any power law. Prominent cases are the Coulomb interaction (`p=-1`) and the
harmonic oscillator (`p=2`).

# Arguments:
- `v0::Float64`: The strength of the potential.
- `p::Float64`: The exponent of the power law. Must satisfy `p > -3`, otherwise the radial integral diverges.
"""
struct PowerLawPotential <: PotentialFunction
    v0::Float64
    p::Float64

    function PowerLawPotential(v0::Float64, p::Float64)
        p <= -3.0 && error("PowerLawPotential: p = $p is too singular; the radial integral diverges (requires p > -3).")
        new(v0, p)
    end
end

# convenience: allow integer exponents, e.g. PowerLawPotential(-1.0,-1)
PowerLawPotential(v0::Real, p::Real) = PowerLawPotential(Float64(v0), Float64(p))

"""
    function (pp::PowerLawPotential)(r)

Evaluates the power-law potential at a given radial distance `r`.
"""
function (pp::PowerLawPotential)(r)
    pp.v0 * r^pp.p
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