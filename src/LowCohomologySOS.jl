module LowCohomologySOS

using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = 4
LinearAlgebra.BLAS.set_num_threads(2)
using Kronecker
using IntervalArithmetic

using StarAlgebras
using Groups

import JuMP

include("group_rings.jl")
include("fox_derivatives.jl")
include("positive_approx.jl")
include("certification.jl")

end
