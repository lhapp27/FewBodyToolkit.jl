# types for the different potential functions
abstract type PotentialFunction end
struct CentralPotential <: PotentialFunction
    f::Function
end

struct GaussianPotential <: PotentialFunction
    v0::Float64
    mu_g::Float64

    function GaussianPotential(v0::Float64, mu_g::Float64)
        new(v0, mu_g)
    end
end

# callable function for GaussianPotential
function (gp::GaussianPotential)(r)
    gp.v0 * exp(-gp.mu_g * r^2)
end


struct SpinOrbitPotential <: PotentialFunction
    f::Function
end

# for evaluating these functions, e.g. in quadgk; used also for observables
function (v::PotentialFunction)(r)
    return v.f(r)
end