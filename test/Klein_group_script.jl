include(joinpath(@__DIR__, "..", "scripts", "utils.jl"))

C₂²_spectral_gaps2 = let half_radius = 1
    A = Alphabet(
            [
                :a,
                :A,
                :b,
                :B,
            ],
            [2, 1, 4, 3],
        )
    F₂ = FreeGroup(A)
    a, b = Groups.gens(F₂)
    C₂² = FPGroup(F₂, [a^2 => one(F₂), b^2 => one(F₂), a * b => b * a])

    quotient_hom = let source = F₂, target = C₂²
        function f(i, source, target)
            if alphabet(source) == alphabet(target)
                Groups.word_type(target)([i])
            else
                throw("Unsupported")
            end
        end
        Groups.Homomorphism(f, source, target)
    end

    relations = [a^2, b^2, a*b*a^(-1)*b^(-1)]

    gensx = gens(C₂²)
    S = let s = gensx
        [s; inv.(s)]
    end
    half_basis, sizes = Groups.wlmetric_ball(S, radius = half_radius)

    Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.spectral_gap_elements(
        quotient_hom,
        relations,
        half_basis
    )

    Δ₁_sos_problem = LowCohomologySOS.sos_problem(
        Δ₁,
        Iₙ
    )

    C₂²_data = let
        (
            M = Δ₁,
            order_unit = Iₙ,
            half_basis = half_basis,
            RG = parent(first(Δ₁)),
        )
    end

    solve_in_loop(
        Δ₁_sos_problem,
        logdir = "./test_logs",
        optimizer = scs_opt(eps = 1e-5, max_iters = 100_000),
        data = C₂²_data
    )
end

rm("./test_logs", recursive = true)
