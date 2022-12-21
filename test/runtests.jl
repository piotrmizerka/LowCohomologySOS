using StarAlgebras
using Groups
using Test
using LowCohomologySOS
using IntervalArithmetic
using SymbolicWedderburn
using PermutationGroups

import Logging
import JuMP
import JuMP.MOI

include(joinpath(@__DIR__, "..", "scripts", "optimizers.jl"))

function cyclic_group(n::Integer)
    A = Alphabet([:a, :A], [2, 1])
    F = FreeGroup(A)
    a, = Groups.gens(F)
    e = one(F)
    Cₙ = FPGroup(F, [a^n => e])

    return Cₙ
end

function test_homomorphism(hom)
    F = hom.source
    @test isone(hom(one(F)))
    @test all(inv(hom(g)) == hom(inv(g)) for g in gens(F))
    @test all(isone(hom(g)*hom(inv(g))) for g in gens(F))
    @test all(hom(g*h) == hom(g)*hom(h) for g in gens(F) for h in gens(F))
end

@testset "LowCohomologySOS" begin
    include("group_rings_tests.jl")
    include("fox_derivatives_tests.jl")
    include("wedderburn_tests.jl")
    include("positive_approx_tests.jl")
    include("positive_approx_symmetrized_tests.jl")
    include("certification_tests.jl")
    include("integration_tests.jl")

    include("Δ1_SL3Z.jl")
    include("Klein_group_script.jl")
end
