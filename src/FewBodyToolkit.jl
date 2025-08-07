module FewBodyToolkit

using LinearAlgebra, SpecialFunctions, QuadGK, Optim, Roots, StaticArrays, Interpolations, PartialWaveFunctions, WignerSymbols, OffsetArrays, Printf, Reexport

### common types and functions
include("common/potentialtypes.jl")
include("common/eigen2step.jl")
include("common/auxiliary.jl")
export PotentialFunction, CentralPotential, GaussianPotential, ContactPotential1D, SpinOrbitPotential, comparison

### GEM-2body
include("GEM-2body/GEM2B.jl")
@reexport using .GEM2B

### GEM-3body-1D
include("GEM-3body-1D/GEM3B1D.jl")
@reexport using .GEM3B1D

### ISGL-3body
include("ISGL-3body/ISGL.jl")
@reexport using .ISGL

end # module FewBodyToolkit
