@testset "constraints" begin
    A = reshape([i for i in 1:9], 3, 3)
    cnstrs_proper = [Vector{Int}() for _ in 1:9]
    cnstrs_proper[1] = [1]; cnstrs_proper[2] = [2]; cnstrs_proper[3] = [3]
    cnstrs_proper[4] = [4]; cnstrs_proper[5] = [5]; cnstrs_proper[6] = [6]
    cnstrs_proper[7] = [7]; cnstrs_proper[8] = [8]; cnstrs_proper[9] = [9]
    @test LowCohomologySOS.constraints(A) == cnstrs_proper

    B = reshape([1 for i in 1:9], 3, 3)
    cnstrs_proper_B = [Vector{Int}()]
    cnstrs_proper_B[1] = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    @test LowCohomologySOS.constraints(B) == cnstrs_proper_B

    pm = reshape([1 for i in 1:9], 3, 3)
    pm[1,:] = [1, 2, 3]
    pm[2,:] = [2, 3, 1] 
    pm[3,:] = [3, 1, 2]
    cnstrs_proper_pm = [Vector{Int}() for _ in 1:3]
    cnstrs_proper_pm[1] = [1, 6, 8] 
    cnstrs_proper_pm[2] = [2, 4, 9] 
    cnstrs_proper_pm[3] = [3, 5, 7]
    @test LowCohomologySOS.constraints(pm) == cnstrs_proper_pm
end

@testset "entry_constraint" begin
    pm = reshape([1 for i in 1:9], 3, 3)
    pm[1,:] = [1, 2, 3]
    pm[2,:] = [2, 3, 1] 
    pm[3,:] = [3, 1, 2]
    cnstrs_proper = LowCohomologySOS.constraints(pm)

    @test all(LowCohomologySOS.entry_constraint(cnstrs_proper, 1, 1, k, 3, 1) == cnstrs_proper[k] for k in 1:3)

    proper_1 = [100*3^2*234+18*3+1, 100*3^2*234+3*234+18*3+3, 100*3^2*234+2*3*234+18*3+2] 
    @test LowCohomologySOS.entry_constraint(cnstrs_proper, 19, 101, 1, 3, 234) == proper_1
    proper_2 = [100*3^2*234+18*3+2, 100*3^2*234+3*234+18*3+1, 100*3^2*234+2*3*234+18*3+3] 
    @test LowCohomologySOS.entry_constraint(cnstrs_proper, 19, 101, 2, 3, 234) == proper_2
    proper_3 = [100*3^2*234+18*3+3, 100*3^2*234+3*234+18*3+2, 100*3^2*234+2*3*234+18*3+1] 
    @test LowCohomologySOS.entry_constraint(cnstrs_proper, 19, 101, 3, 3, 234) == proper_3
end

@testset "sos_problem_matrix" begin
    C₃ = cyclic_group(3)
    a, = Groups.gens(C₃)

    RC₃ = LowCohomologySOS.group_ring(C₃, [one(C₃), a, a^2], false)
    M_1 = [RC₃(a)]
    order_unit_1 = [one(RC₃)]
    @test (@capture_out print(LowCohomologySOS.sos_problem_matrix(M_1, order_unit_1))) == 
    "Max λ\n"*
    "Subject to\n"*
    " -P[1,1] - P[3,2] - P[2,3] - λ = 0.0\n"*
    " -P[2,1] - P[1,2] - P[3,3] = -1.0\n"*
    " -P[3,1] - P[2,2] - P[1,3] = 0.0\n"*
    " sdp : [P[1,1]  P[1,2]  P[1,3];\n"*
    "  P[2,1]  P[2,2]  P[2,3]; ∈ JuMP.PSDCone()\n"*
    "  P[3,1]  P[3,2]  P[3,3]]\n"

    RC₃_star = LowCohomologySOS.group_ring(C₃, [one(C₃), a, a^2], true)
    M_2 = [RC₃_star(a)]
    order_unit_2 = [one(RC₃_star)]
    @test (@capture_out print(LowCohomologySOS.sos_problem_matrix(M_2, order_unit_2))) == 
    "Max λ\n"*
    "Subject to\n"*
    " -P[1,1] - P[2,2] - P[3,3] - λ = 0.0\n"*
    " -P[3,1] - P[1,2] - P[2,3] = -1.0\n"*
    " -P[2,1] - P[3,2] - P[1,3] = 0.0\n"*
    " sdp : [P[1,1]  P[1,2]  P[1,3];\n"*
    "  P[2,1]  P[2,2]  P[2,3]; ∈ JuMP.PSDCone()\n"*
    "  P[3,1]  P[3,2]  P[3,3]]\n"

    C₂ = LowCohomologySOS.cyclic_group(2)
    x, = Groups.gens(C₂)
    RC₂ = LowCohomologySOS.group_ring(C₂, [one(C₂), x])
    M_3 = [RC₂(x) zero(RC₂); zero(RC₂) RC₂(x)]
    order_unit_3 = [one(RC₂) zero(RC₂); zero(RC₂) one(RC₂)]
    
    @test (@capture_out print(LowCohomologySOS.sos_problem_matrix(M_3, order_unit_3))) == 
    "Max λ\n"*
    "Subject to\n"*
    " -P[1,1] - P[2,2] - λ = 0.0\n"*
    " -P[2,1] - P[1,2] = -1.0\n"*
    " -P[1,3] - P[2,4] = 0.0\n"*
    " -P[2,3] - P[1,4] = 0.0\n"*
    " -P[3,1] - P[4,2] = 0.0\n"*
    " -P[4,1] - P[3,2] = 0.0\n"*
    " -P[3,3] - P[4,4] - λ = 0.0\n"*
    " -P[4,3] - P[3,4] = -1.0\n"*
    " sdp : [P[1,1]  P[1,2]  P[1,3]  P[1,4];\n"*
    "  P[2,1]  P[2,2]  P[2,3]  P[2,4];\n"*
    "  P[3,1]  P[3,2]  P[3,3]  P[3,4]; ∈ JuMP.PSDCone()\n"*
    "  P[4,1]  P[4,2]  P[4,3]  P[4,4]]\n"

    @test (@capture_out print(LowCohomologySOS.sos_problem_matrix(M_3, order_unit_3, 0.5))) == 
    "Max λ\n"*
    "Subject to\n"*
    " -P[1,1] - P[2,2] - λ = 0.0\n"*
    " -P[2,1] - P[1,2] = -1.0\n"*
    " -P[1,3] - P[2,4] = 0.0\n"*
    " -P[2,3] - P[1,4] = 0.0\n"*
    " -P[3,1] - P[4,2] = 0.0\n"*
    " -P[4,1] - P[3,2] = 0.0\n"*
    " -P[3,3] - P[4,4] - λ = 0.0\n"*
    " -P[4,3] - P[3,4] = -1.0\n"*
    " sdp : [P[1,1]  P[1,2]  P[1,3]  P[1,4];\n"*
    "  P[2,1]  P[2,2]  P[2,3]  P[2,4];\n"*
    "  P[3,1]  P[3,2]  P[3,3]  P[3,4]; ∈ JuMP.PSDCone()\n"*
    "  P[4,1]  P[4,2]  P[4,3]  P[4,4]]\n"*
    " λ ≤ 0.5\n"
end

@testset "sos_problem_solution_scs" begin
    C₃ = cyclic_group(3)
    a, = Groups.gens(C₃)
    RC₃_star = LowCohomologySOS.group_ring(C₃, [one(C₃), a, a^2], true)
    M_1 = reshape([RC₃_star(a)], 1, 1)
    order_unit_1 = reshape([one(RC₃_star)], 1, 1)
    sos_problem_infeasible = LowCohomologySOS.sos_problem_matrix(M_1, order_unit_1)
    λ_1, P_1 = 1, 1
    @suppress begin
        λ_1, P_1 = LowCohomologySOS.sos_problem_solution_scs(sos_problem_infeasible, true)
    end
    @test isnan(λ_1) && isnan.(P_1) == [1 1 1; 1 1 1; 1 1 1]

    A = Alphabet([:x, :X], [2, 1])
    Z = FreeGroup(A)
    x, = Groups.gens(Z)
    RZ = LowCohomologySOS.group_ring(Z, 2)
    RZ_star = LowCohomologySOS.group_ring(Z, 2, true)
    Δ = reshape([2*one(RZ)-RZ(x^(-1))-RZ(x)], 1, 1)
    M_2 = LowCohomologySOS.embed_to_group_ring.(Δ^2, Ref(RZ_star), identity)
    order_unit_2 = LowCohomologySOS.embed_to_group_ring.(Δ, Ref(RZ_star), identity)
    sos_problem_infeasible_2 = LowCohomologySOS.sos_problem_matrix(M_2, order_unit_2)
    λ_2, P_2 = 1, 1
    @suppress begin
        λ_2, P_2 = LowCohomologySOS.sos_problem_solution_scs(sos_problem_infeasible_2, true)
    end
    @test λ_2+1 ≈ 1 rtol=4.0e-3 # something strange happens when comparing values close to 0

    M_3 = [4*one(RZ_star) zero(RZ_star) zero(RZ_star);
           zero(RZ_star) 3*one(RZ_star) zero(RZ_star);
           zero(RZ_star) zero(RZ_star) 5*one(RZ_star)]
    order_unit_3 = [one(RZ_star) zero(RZ_star) zero(RZ_star);
                    zero(RZ_star) one(RZ_star) zero(RZ_star);
                    zero(RZ_star) zero(RZ_star) one(RZ_star)]
    sos_problem_3 = LowCohomologySOS.sos_problem_matrix(M_3, order_unit_3)
    λ_3 = 0
    P_3 = reshape([0 for i in 1:15^2], 15, 15)
    @suppress begin
        λ_3, P_3 = LowCohomologySOS.sos_problem_solution_scs(sos_problem_3, true)
    end
    @test λ_3 ≈ 3 rtol=4.0e-3
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
    λ_1, P_1, RG_ball_star, Δ₁_1, I_1 = 1, 1, 1, 1, 1
    @suppress begin
        λ_1, P_1, RG_ball_star, Δ₁_1, I_1 = LowCohomologySOS.spectral_gaps_approximated(quotient_hom, relations, half_basis, true)
    end
    Δ₁_1_proper = reshape([5*one(RG_ball_star)+2*RG_ball_star(xx)+2*RG_ball_star(xx^2)], 1, 1)
    I_1_proper = reshape([one(RG_ball_star)], 1, 1)
    @test λ_1 ≈ 3 rtol=1e-3
    @test Δ₁_1 == Δ₁_1_proper
    @test I_1 == I_1_proper
end