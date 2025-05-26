# types for the different potential functions
abstract type PotentialFunction end
struct CentralPotential <: PotentialFunction
    f::Function
end
struct GaussianPotential <: PotentialFunction
    f::Function
end
struct SpinOrbitPotential <: PotentialFunction
    f::Function
end

# for evaluating these functions, e.g. in quadgk; used also for observables
function (v::PotentialFunction)(r)
    return v.f(r)
end