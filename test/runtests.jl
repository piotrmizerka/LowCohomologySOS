using StarAlgebras
using Groups
using Test
using LowCohomologySOS
using AbstractAlgebra
using Suppressor

@testset "LowCohomologySOS" begin
   include("group_rings_tests.jl")
   include("fox_derivatives_tests.jl")
   include("positive_approx_tests.jl")
end