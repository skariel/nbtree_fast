using Base.Threads
using Base.Test

set_zero_subnormals(true)

include("constants.jl")
include("particles.jl")
include("tree.jl")
include("dtree.jl")
include("gravity.jl")
include("testing.jl")