using Pkg

Pkg.activate(@__DIR__)
using StarAlgebras
using Groups
using AbstractAlgebra
using LowCohomologySOS
using JuMP
using SCS
using MathOptInterface

A = Alphabet([:x, :X, :y, :Y], [2, 1, 4, 3])
    F = FreeGroup(A)
    RF = LowCohomologySOS.group_ring(F, 1)
    x, y = Groups.gens(F)
    G = FPGroup(F, [x * y => y * x])
    RG = LowCohomologySOS.group_ring(G, 1)

testt = let 
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
    A = Alphabet([:x, :X], [2, 1])
    Z = FreeGroup(A)
    x, = Groups.gens(Z)
    RZ = LowCohomologySOS.group_ring(Z, 2)
    RZ_star = LowCohomologySOS.group_ring(Z, 2, true)

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

    @info termination_status_3
    termination_status_3
end

SL₃ƵSpectralGaps = let halfRadius = 2
    A = Alphabet([:e12, :E12, :e21, :E21, :e13, :E13, :e31, :E31, :e23, :E23, :e32, :E32], [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11])
    F = FreeGroup(A)
    e12, e21, e13, e31, e23, e32 = Groups.gens(F)
 
    N = 3
    SL₃Ƶ = MatrixAlgebra(zz, N)
    E(SL₃Ƶ, i,j) = (e_ij = one(SL₃Ƶ); e_ij[i,j] = 1; e_ij)
    S = [E(SL₃Ƶ, i,j) for i in 1:N for j in 1:N if i≠j]
    S = unique([S; inv.(S)])

    function h(u::FPGroupElement)
        result = one(SL₃Ƶ)

        for i in 1:length(word(u))
            if F(word(u)[i:i]) == e12
                result *= S[1]
            elseif F(word(u)[i:i]) == inv(e12)
                result *= inv(S[1])
            elseif F(word(u)[i:i]) == e21
                result *= S[3]
            elseif F(word(u)[i:i]) == inv(e21)
                result *= inv(S[3])
            elseif F(word(u)[i:i]) == e13
                result *= S[2]
            elseif F(word(u)[i:i]) == inv(e13)
                result *= inv(S[2])
            elseif F(word(u)[i:i]) == e31
                result *= S[5]
            elseif F(word(u)[i:i]) == inv(e31)
                result *= inv(S[5])
            elseif F(word(u)[i:i]) == e23
                result *= S[4]
            elseif F(word(u)[i:i]) == inv(e23)
                result *= inv(S[4])
            elseif F(word(u)[i:i]) == e32
                result *= S[6]
            elseif F(word(u)[i:i]) == inv(e32)
                result *= inv(S[6])
            end
        end

        return result
    end

    relations = [e12*e13*e12^(-1)*e13^(-1), e12*e32*e12^(-1)*e32^(-1), e13*e23*e13^(-1)*e23^(-1), 
                 e23*e21*e23^(-1)*e21^(-1), e21*e31*e21^(-1)*e31^(-1), e31*e32*e31^(-1)*e32^(-1),
                 e12*e23*e12^(-1)*e23^(-1)*e13^(-1), e13*e32*e13^(-1)*e32^(-1)*e12^(-1), 
                 e21*e13*e21^(-1)*e13^(-1)*e23^(-1), e23*e31*e23^(-1)*e31^(-1)*e21^(-1), 
                 e31*e12*e31^(-1)*e12^(-1)*e32^(-1), e32*e21*e32^(-1)*e21^(-1)*e31^(-1)]
    
    halfBasis, sizes = Groups.wlmetric_ball(S, radius = halfRadius)
    LowCohomologySOS.spectral_gaps_certification(h, relations, halfBasis, optimizer=scs_opt(eps = 1e-8, verbose = false))
end
