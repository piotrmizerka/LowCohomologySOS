@testset "constraints" begin
    A = reshape([i for i in 1:9], 3, 3)
    cnstrs = [
        (x = zeros(Int, 3, 3); x[i] = 1; x)
        for i in 1:9]

    @test all(LowCohomologySOS.constraints(A) .== cnstrs)

    B = reshape([2 for i in 1:9], 3, 3)
    cnstrs_proper_B = [zeros(Int, 3,3), ones(Int, 3,3)]
    @test LowCohomologySOS.constraints(B) == cnstrs_proper_B

    pm = [1 2 3;
          2 3 1;
          3 1 5]

    cnstrs = [zeros(Int, 3, 3) for _ in 1:5]
    cnstrs[1][[1, 6, 8]] .= 1
    cnstrs[2][[2, 4]] .= 1
    cnstrs[3][[3, 5, 7]] .= 1
    cnstrs[5][[9]] .= 1

    @test LowCohomologySOS.constraints(pm) == cnstrs
end

@testset "Tensors" begin
    @test LowCohomologySOS.KroneckerDelta{3}(3,2) isa AbstractMatrix{Int}
    @test LowCohomologySOS.KroneckerDelta{3}(3,2) == [0 0 0; 0 0 0; 0 1 0]

    a = LowCohomologySOS.KroneckerDelta{3}(3,2)
    b = LowCohomologySOS.KroneckerDelta{2}(1,2)
    @test b == [0 1; 0 0]

    @test LowCohomologySOS.Tensor(a, b) == kron(a, b)

    A = rand(3,2)
    @test LowCohomologySOS.Tensor(a, A) == kron(a, A)

    a = LowCohomologySOS.KroneckerDelta{3}(3,2)
    x = rand(3,3)
    @test LowCohomologySOS.Tensor(a, x)[a] == x
    @test kron(a, x)[a] == x

    @test @view(LowCohomologySOS.Tensor(a, x)[a]) == x
end

@testset "indexing by KroneckerDelta" begin
    pm = reshape([1 for i in 1:9], 3, 3)
    pm[1, :] = [1, 2, 3]
    pm[2, :] = [2, 3, 1]
    pm[3, :] = [3, 1, 2]
    cnstrs = LowCohomologySOS.constraints(pm)

    # TODO: finish this


end

@testset "sos_problem (non symmetrized version)" begin
    C₃ = cyclic_group(3)
    a, = Groups.gens(C₃)
    RC₃ = LowCohomologySOS.group_ring(C₃, [one(C₃), a, a^2])
    M = reshape([RC₃(a)], 1, 1)
    order_unit = reshape([one(RC₃)], 1, 1)
    m = LowCohomologySOS.sos_problem(M, order_unit)
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
        LowCohomologySOS.sos_problem(M_1, order_unit_1)

    let model = sos_problem_infeasible
        JuMP.set_optimizer(model, scs_opt(verbose = false))
        JuMP.optimize!(model)
        # λ, Q = LowCohomologySOS.get_solution(model)
        @test JuMP.termination_status(model) == MOI.INFEASIBLE
    end

    A = Alphabet([:x, :X], [2, 1])
    Z = FreeGroup(A)
    x, = Groups.gens(Z)
    RZ = LowCohomologySOS.group_ring(Z, 2)
    RZ_star = LowCohomologySOS.group_ring(Z, 2, true)

    Δ = reshape([2 * one(RZ) - RZ(x^(-1)) - RZ(x)], 1, 1)
    M_2 = LowCohomologySOS.embed.(identity, Δ^2, Ref(RZ_star))
    order_unit_2 = LowCohomologySOS.embed.(identity, Δ, Ref(RZ_star))
    sos_problem_infeasible_2 = let elt = M_2, unit = order_unit_2
        m = LowCohomologySOS.sos_problem(M_2, order_unit_2)
        JuMP.@constraint(m, m[:λ] >= 0.2)
        m
    end

    let model = sos_problem_infeasible_2
        JuMP.set_optimizer(model, scs_opt(verbose = false))
        JuMP.optimize!(model)
        # λ, Q = LowCohomologySOS.get_solution(model)
        @test JuMP.termination_status(model) == MOI.INFEASIBLE
    end

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
    sos_problem_3 = LowCohomologySOS.sos_problem(M_3, order_unit_3)

    λ_3, Q_3, ts_3 = let model = sos_problem_3
        JuMP.set_optimizer(model, scs_opt(verbose = false))
        JuMP.optimize!(model)
        λ, Q = LowCohomologySOS.get_solution(model)
        λ, Q, JuMP.termination_status(model)
    end

    @test ts_3 == MOI.OPTIMAL
    @test λ_3 ≈ 3 rtol = 4.0e-10
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

    solution = LowCohomologySOS.spectral_gaps_approximated(
        quotient_hom,
        relations,
        half_basis,
        optimizer = scs_opt(verbose = false)
    )

    RG = parent(first(solution.laplacian))

    Δ₁ = reshape([5 * one(RG) + 2 * RG(xx) + 2 * RG(xx^2)], 1, 1)
    unit = reshape([one(RG)], 1, 1)
    @test solution.termination_status == MOI.OPTIMAL
    @test solution.λ ≈ 3 rtol = 1e-3
    @test solution.laplacian == Δ₁
    @test solution.unit == unit

    # TODO: test this with symmetric S
end
