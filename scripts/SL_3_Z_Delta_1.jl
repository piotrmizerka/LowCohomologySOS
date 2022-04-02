using Pkg

Pkg.activate(joinpath(@__DIR__, ".."))

using Revise

include(joinpath(@__DIR__, "..", "scripts", "utils.jl"))

SL₃ℤ_spectral_gaps = let half_radius = 2
    SL(n, R) = PropertyT_new.SpecialLinearGroup{n}(R)
    SL₃ℤ = SL(3, Int8)

    gensx = gens(SL₃ℤ)
    S = let s = gensx
        [s; inv.(s)]
    end

    half_basis, sizes = Groups.wlmetric_ball(S, radius = half_radius)

    A = Alphabet(
            [
                :e12,
                :E12,
                :e21,
                :E21,
                :e13,
                :E13,
                :e31,
                :E31,
                :e23,
                :E23,
                :e32,
                :E32,
            ],
            [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11],
        )
    F_sl_3_z = FreeGroup(A)
    e12, e21, e13, e31, e23, e32 = Groups.gens(F_sl_3_z)

    quotient_hom = let source = F_sl_3_z, target = SL₃ℤ
        function hom(i, source, target::PropertyT_new.SpecialLinearGroup{N, T}) where {N, T}
            As = alphabet(source)
            s = String(As[i])
            (i,j,isinv) = let m = match(r"^(E|e)(\d)(\d)$", s)
                @assert !isnothing(m)
                isinv = isuppercase(Char(first(m[1])))
                i = parse(Int, m[2])
                j = parse(Int, m[3])
                i,j, isinv
            end

            At = alphabet(target)
            eij = PropertyT_new.MatrixGroups.ElementaryMatrix{N}(i, j, T(isinv ? -1 : 1))
            return Groups.word_type(target)([At[eij]])
        end
        PropertyT_new.Homomorphism(hom, source, target)
    end

    relations = [e12*e13*e12^(-1)*e13^(-1), e12*e32*e12^(-1)*e32^(-1), e13*e23*e13^(-1)*e23^(-1),
                 e23*e21*e23^(-1)*e21^(-1), e21*e31*e21^(-1)*e31^(-1), e31*e32*e31^(-1)*e32^(-1),
                 e12*e23*e12^(-1)*e23^(-1)*e13^(-1), e13*e32*e13^(-1)*e32^(-1)*e12^(-1),
                 e21*e13*e21^(-1)*e13^(-1)*e23^(-1), e23*e31*e23^(-1)*e31^(-1)*e21^(-1),
                 e31*e12*e31^(-1)*e12^(-1)*e32^(-1), e32*e21*e32^(-1)*e21^(-1)*e31^(-1)]

    Δ₁_sos_problem, Δ₁, Iₙ, half_basis, RG_ball_star = LowCohomologySOS.sos_problem_delta_1(
        quotient_hom, 
        relations, 
        half_basis
    )

    SL₃ℤ_data = let 
        (
            M = Δ₁,
            order_unit = Iₙ,
            half_basis = half_basis,
            RG = RG_ball_star,
        )
    end

    solve_in_loop(
        Δ₁_sos_problem, 
        logdir = "./logs", 
        optimizer = scs_opt(eps = 1e-5, max_iters = 100_000),
        data = SL₃ℤ_data
    )
end
