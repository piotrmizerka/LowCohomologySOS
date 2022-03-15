using StarAlgebras
using Groups
using Test
using LowCohomologySOS
using IntervalArithmetic
using PropertyT_new

import AbstractAlgebra
import Logging
import JuMP
import SCS
import PropertyT_new.MOI

function cyclic_group(n::Integer)
    A = Alphabet([:a, :A], [2, 1])
    F = FreeGroup(A)
    a, = Groups.gens(F)
    e = one(F)
    Cₙ = FPGroup(F, [a^n => e])

    return Cₙ
end

function scs_opt(;
    eps = 1e-5,
    acceleration_lookback = 0,
    max_iters = 20_000,
    verbose = true,
)
    return JuMP.optimizer_with_attributes(
        SCS.Optimizer,
        "eps" => eps,
        "acceleration_lookback" => acceleration_lookback,
        "max_iters" => max_iters,
        "verbose" => verbose,
    )
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
    include("positive_approx_tests.jl")
    include("certification_tests.jl")
    include("integration_tests.jl")
end
