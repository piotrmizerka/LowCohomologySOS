using StarAlgebras
using Groups
using Test
using LowCohomologySOS
using IntervalArithmetic
using PropertyT_new

import Logging
import JuMP
import SCS
import JuMP.MOI

function cyclic_group(n::Integer)
    A = Alphabet([:a, :A], [2, 1])
    F = FreeGroup(A)
    a, = Groups.gens(F)
    e = one(F)
    Cₙ = FPGroup(F, [a^n => e])

    return Cₙ
end

function scs_opt(;
    accel = 10,
    alpha = 1.5,
    eps = 1e-9,
    max_iters = 10_000,
    verbose = true,
)
    return JuMP.optimizer_with_attributes(
        SCS.Optimizer,
        "acceleration_lookback" => accel,
        "acceleration_interval" => max(abs(accel), 1),
        "alpha" => alpha,
        "eps_abs" => eps,
        "eps_rel" => eps,
        "linear_solver" => SCS.DirectSolver,
        "max_iters" => max_iters,
        "warm_start" => true,
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

    include("Klein_group_script.jl" )
end
