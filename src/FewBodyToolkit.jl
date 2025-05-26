module FewBodyToolkit

using LinearAlgebra, SpecialFunctions, QuadGK, Optim, Roots, StaticArrays, Interpolations, HypergeometricFunctions, PartialWaveFunctions, WignerSymbols, OffsetArrays, GSL, Printf

### GEM-2body
include("GEM-2body/GEM2B.jl")
using .GEM2B

# q1D currently unsupported
#include("GEM-2body/GEMq1D.jl")
#using .GEMq1D

### GEM-3body-1D
include("GEM-3body-1D/GEM3B1D.jl")
using .GEM3B1D

### ISGL-3body
include("ISGL-3body/ISGL.jl")
using .ISGL

end # module FewBodyToolkit
