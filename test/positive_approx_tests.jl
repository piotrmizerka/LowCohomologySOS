@testset "constraints" begin
    A = reshape([i for i in 1:9], 3, 3)
    cnstrs_proper = [Vector{Int}() for _ in 1:9]
    cnstrs_proper[1] = [1]
    cnstrs_proper[2] = [2]
    cnstrs_proper[3] = [3]
    cnstrs_proper[4] = [4]
    cnstrs_proper[5] = [5]
    cnstrs_proper[6] = [6]
    cnstrs_proper[7] = [7]
    cnstrs_proper[8] = [8]
    cnstrs_proper[9] = [9]
    @test LowCohomologySOS.constraints(A) == cnstrs_proper

    B = reshape([1 for i in 1:9], 3, 3)
    cnstrs_proper_B = [Vector{Int}()]
    cnstrs_proper_B[1] = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    @test LowCohomologySOS.constraints(B) == cnstrs_proper_B

    pm = reshape([1 for i in 1:9], 3, 3)
    pm[1, :] = [1, 2, 3]
    pm[2, :] = [2, 3, 1]
    pm[3, :] = [3, 1, 2]
    cnstrs_proper_pm = [Vector{Int}() for _ in 1:3]
    cnstrs_proper_pm[1] = [1, 6, 8]
    cnstrs_proper_pm[2] = [2, 4, 9]
    cnstrs_proper_pm[3] = [3, 5, 7]
    @test LowCohomologySOS.constraints(pm) == cnstrs_proper_pm
end

@testset "entry_constraint" begin
    pm = reshape([1 for i in 1:9], 3, 3)
    pm[1, :] = [1, 2, 3]
    pm[2, :] = [2, 3, 1]
    pm[3, :] = [3, 1, 2]
    cnstrs_proper = LowCohomologySOS.constraints(pm)

    @test all(
        LowCohomologySOS.entry_constraint(cnstrs_proper, 1, 1, k, 3, 1) ==
        cnstrs_proper[k] for k in 1:3
    )

    proper_1 = [
        100 * 3^2 * 234 + 18 * 3 + 1,
        100 * 3^2 * 234 + 3 * 234 + 18 * 3 + 3,
        100 * 3^2 * 234 + 2 * 3 * 234 + 18 * 3 + 2,
    ]
    @test LowCohomologySOS.entry_constraint(
        cnstrs_proper,
        19,
        101,
        1,
        3,
        234,
    ) == proper_1
    proper_2 = [
        100 * 3^2 * 234 + 18 * 3 + 2,
        100 * 3^2 * 234 + 3 * 234 + 18 * 3 + 1,
        100 * 3^2 * 234 + 2 * 3 * 234 + 18 * 3 + 3,
    ]
    @test LowCohomologySOS.entry_constraint(
        cnstrs_proper,
        19,
        101,
        2,
        3,
        234,
    ) == proper_2
    proper_3 = [
        100 * 3^2 * 234 + 18 * 3 + 3,
        100 * 3^2 * 234 + 3 * 234 + 18 * 3 + 2,
        100 * 3^2 * 234 + 2 * 3 * 234 + 18 * 3 + 1,
    ]
    @test LowCohomologySOS.entry_constraint(
        cnstrs_proper,
        19,
        101,
        3,
        3,
        234,
    ) == proper_3
end

@testset "sos_problem_matrix" begin
    C₃ = cyclic_group(3)
    a, = Groups.gens(C₃)
    RC₃ = LowCohomologySOS.group_ring(C₃, [one(C₃), a, a^2])
    M = reshape([RC₃(a)], 1, 1)
    order_unit = reshape([one(RC₃)], 1, 1)
    m = LowCohomologySOS.sos_problem_matrix(M, order_unit)
    output = sprint(print, m)

    @test occursin("Max λ", output)
    @test occursin("P[1,1]", sprint(print, m[:P]))
    @test occursin("P[3,3]", output) && !occursin("P[4,4]", output)

    inclusion_sign = Sys.iswindows() ? "in" : "∈"
    @test occursin(inclusion_sign*" JuMP.PSDCone()", sprint(print, m[:sdp]))
end

@testset "sos_problem_solution" begin
    C₃ = cyclic_group(3)
    a, = Groups.gens(C₃)
    RC₃_star = LowCohomologySOS.group_ring(C₃, [one(C₃), a, a^2], star_multiplication = true)

    M_1 = reshape([RC₃_star(a)], 1, 1)
    order_unit_1 = reshape([one(RC₃_star)], 1, 1)
    sos_problem_infeasible =
        LowCohomologySOS.sos_problem_matrix(M_1, order_unit_1)

    λ_1, P_1, termination_status_1 = LowCohomologySOS.sos_problem_solution(
        sos_problem_infeasible,
        optimizer = scs_opt(verbose = false),
    )
    @test termination_status_1 == MOI.INFEASIBLE

    A = Alphabet([:x, :X], [2, 1])
    Z = FreeGroup(A)
    x, = Groups.gens(Z)
    RZ = LowCohomologySOS.group_ring(Z, 2)
    RZ_star = LowCohomologySOS.group_ring(Z, 2, true)

    Δ = reshape([2 * one(RZ) - RZ(x^(-1)) - RZ(x)], 1, 1)
    M_2 = LowCohomologySOS.embed.(identity, Δ^2, Ref(RZ_star))
    order_unit_2 = LowCohomologySOS.embed.(identity, Δ, Ref(RZ_star))
    sos_problem_infeasible_2 = let elt = M_2, unit = order_unit_2
        m = LowCohomologySOS.sos_problem_matrix(M_2, order_unit_2)
        JuMP.@constraint(m, m[:λ] >= 0.2)
        m
    end

    λ_2, P_2, termination_status_2 = LowCohomologySOS.sos_problem_solution(
        sos_problem_infeasible_2,
        optimizer = scs_opt(verbose = false),
    )

    @test termination_status_2 == MOI.INFEASIBLE

    M_3 = [
        4*one(RZ_star) zero(RZ_star) zero(RZ_star)
        zero(RZ_star) 3*one(RZ_star) zero(RZ_star)
        zero(RZ_star) zero(RZ_star) 5*one(RZ_star)
    ]
    order_unit_3 = [
        one(RZ_star) zero(RZ_star) zero(RZ_star)
        zero(RZ_star) one(RZ_star) zero(RZ_star)
        zero(RZ_star) zero(RZ_star) one(RZ_star)
    ]
    sos_problem_3 = LowCohomologySOS.sos_problem_matrix(M_3, order_unit_3)

    λ_3, P_3, termination_status_3 = LowCohomologySOS.sos_problem_solution(
        sos_problem_3,
        optimizer = scs_opt(verbose = false),
    )

    @test termination_status_3 == MOI.OPTIMAL
    @test λ_3 ≈ 3 rtol = 4.0e-3
end

@testset "spectral_gaps_approximated" begin
    A = Alphabet([:x, :X], [2, 1])
    F = FreeGroup(A)
    x, = Groups.gens(F)
    C₃ = FPGroup(F, [x^3 => one(F)])
    xx, = Groups.gens(C₃)
    function quotient_hom(u::FPGroupElement)
        result = one(C₃)
        for i in 1:length(word(u))
            if F(word(u)[i:i]) == x
                result *= xx
            elseif F(word(u)[i:i]) == inv(x)
                result *= inv(xx)
            end
        end
        return result
    end
    relations = [x^3]
    half_basis, sizes = Groups.wlmetric_ball([xx, xx^(-1)], radius = 1)

    λ_1, P_1, termination_status_1, RG_ball_star, Δ₁_1, I_1 =
        LowCohomologySOS.spectral_gaps_approximated(
            quotient_hom,
            relations,
            half_basis,
            optimizer = scs_opt(verbose = false),
        )

    Δ₁_1_proper = reshape(
        [5 * one(RG_ball_star) + 2 * RG_ball_star(xx) + 2 * RG_ball_star(xx^2)],
        1,
        1,
    )
    I_1_proper = reshape([one(RG_ball_star)], 1, 1)
    @test termination_status_1 == MOI.OPTIMAL
    @test λ_1 ≈ 3 rtol = 1e-3
    @test Δ₁_1 == Δ₁_1_proper
    @test I_1 == I_1_proper
end
