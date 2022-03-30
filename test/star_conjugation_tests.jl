@testset "associated_elements" begin
    F₂ = FreeGroup(2)
    a, b = Groups.gens(F₂)
    B₁ = [one(F₂), a, b, a^(-1), b^(-1)]

    gsⱼ⁻¹_F₂, sᵢg_F₂, sᵢgsⱼ⁻¹_F₂ = LowCohomologySOS.associated_elements(F₂, B₁, B₁)

    @test gsⱼ⁻¹_F₂ == [[(1,4),(2,5)], [(1,1)], [(2,1)], [], []]
    @test sᵢg_F₂ == [[(1,2),(2,3)], [], [], [(1,1)], [(2,1)]]
    @test sᵢgsⱼ⁻¹_F₂ == [[(1,1,1),(2,2,1)], [(1,1,2),(2,1,3)], [(1,2,2),(2,2,3)],
                        [(1,1,4),(1,2,5)], [(2,1,4),(2,2,5)]]

    ℤ = FreeGroup(1)
    x, = Groups.gens(ℤ)

    B₃ = [one(ℤ), x, x^(-1), x^2, x^(-2), x^3, x^(-3)]
    gsⱼ⁻¹_ℤ, sᵢg_ℤ, sᵢgsⱼ⁻¹_ℤ = LowCohomologySOS.associated_elements(ℤ, B₃, B₃)

    @test gsⱼ⁻¹_ℤ == [[(1,3)], [(1,1)], [(1,5)], [(1,2)], [(1,7)], [(1,4)], []]
    @test sᵢg_ℤ == [[(1,2)], [(1,4)], [(1,1)], [(1,6)], [(1,3)], [], [(1,5)]]
    @test sᵢgsⱼ⁻¹_ℤ == [[(1,1,1)], [(1,1,2)], [(1,1,3)], [(1,1,4)],
                        [(1,1,5)], [(1,1,6)], [(1,1,7)]]

    B₂ = [one(ℤ), x, x^(-1), x^2, x^(-2)]
    gsⱼ⁻¹_ℤ_smaller_ball, sᵢg_ℤ_smaller_ball, sᵢgsⱼ⁻¹_ℤ_smaller_ball = LowCohomologySOS.associated_elements(ℤ, B₃, B₂)
    @test gsⱼ⁻¹_ℤ_smaller_ball == [[(1,3)], [(1,1)], [(1,5)], [(1,2)], [], [(1,4)], []]
    @test sᵢg_ℤ_smaller_ball == [[(1,2)], [(1,4)], [(1,1)], [], [(1,3)], [], [(1,5)]]
    @test sᵢgsⱼ⁻¹_ℤ_smaller_ball == [[(1,1,1)], [(1,1,2)], [(1,1,3)], [(1,1,4)],
                                     [(1,1,5)], [], []]
end

@testset "sos_problem_conjugated" begin
    ℤ = FreeGroup(1)
    x, = Groups.gens(ℤ)

    B₁ = [one(ℤ), x, x^(-1)]
    RℤB₂ = LowCohomologySOS.group_ring(ℤ, B₁, star_multiplication = true)
    B₄ = [one(ℤ), x, x^(-1), x^(-2), x^2, x^3, x^(-3), x^4, x^(-4)] # the order of the basis of RℤB₂ is: one(ℤ), x, x^(-1), x^(-2), x^2
    RℤB₄ = StarAlgebra(ℤ, StarAlgebras.Basis{Int}(B₄))

    Δ = 2*one(RℤB₄)-RℤB₄(x^(-1))-RℤB₄(x)

    sos_problem_computed = LowCohomologySOS.sos_problem_conjugated(Δ^2, Δ, G = ℤ, smaller_group_ring = RℤB₂)
    output = sprint(print, sos_problem_computed)

    @test occursin("Max λ", output)
    @test occursin("Subject to", output)
    @test occursin("-2 P[1,1] + 2 P[1,2] - 2 P[2,2] + 2 P[1,3] - 2 P[3,3] - 2 λ = -6.0", output)
    @test occursin("P[1,1] - 2 P[1,2] + P[2,2] - 2 P[1,3] + P[2,3] + P[3,3] + λ = 4.0", output)
    @test occursin("P[1,2] + P[1,3] - 2 P[2,3] = -1.0", output)
    @test occursin("P[2,3] = 0.0", output)
    @test !occursin("P[4,4]", output)

    inclusion_sign = Sys.iswindows() ? "in" : "∈"
    @test occursin(inclusion_sign*" JuMP.PSDCone()", sprint(print, sos_problem_computed[:sdp]))
end

# TODO: tests for property_t_conjugated_approx