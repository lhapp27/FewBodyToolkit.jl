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
struct CentralPotential <: PotentialFunction
    f::Function
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
    SpinOrbitPotential(f::Function)
Defines a potential of the type `SpinOrbitPotential`. The function `f(r)` represents the radial part of a spin-orbit interaction
```math
V_{SO}(r) = f(r) \\vec{l} \\cdot \\vec{s}
```
"""
struct SpinOrbitPotential <: PotentialFunction
    f::Function
end


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