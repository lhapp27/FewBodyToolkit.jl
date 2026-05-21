using FewBodyToolkit
using Test

# GEM2B:
@time include("test2B1D.jl")
@time include("test2B2D.jl")
@time include("test2B3D.jl")
@time include("testMatrixElements.jl")
@time include("testWavefunctions.jl")

# GEM3B1D:
@time include("test3B1D.jl")

# ISGL:
@time include("testISGL.jl")