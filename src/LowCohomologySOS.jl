module LowCohomologySOS

using LinearAlgebra
using IntervalArithmetic

using StarAlgebras
using GroupsCore
using Groups

import JuMP
import JuMP.MOI

include("tensors.jl")
include("group_rings.jl")
include("fox_derivatives.jl")
include("positive_approx.jl")
include("certification.jl")
include("star_conjugation.jl")

end
