using Base.Test 
using Base.Threads


set_zero_subnormals(true)

include("parallel.jl")
include("particles.jl")
include("tree.jl")
include("gravity.jl")
include("testing.jl")