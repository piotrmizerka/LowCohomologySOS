using LowCohomologySOS
using Groups
using PropertyT_new
using JuMP
using SCS

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

SAut_F₄SpectralGaps = let halfRadius = 2
    SAut_F₄ = SpecialAutomorphismGroup(FreeGroup(4))
    gensx = gens(SAut_F₄)
    S = let s = gensx
        [s; inv.(s)]
    end

    halfBasis, sizes = Groups.wlmetric_ball(S, radius = halfRadius)

    F24 = FreeGroup(alphabet(SAut_F₄))
    quotient_hom = let source = F24, target = SAut_F₄
        function f(i, source, target)
            if alphabet(source) == alphabet(target)
                Groups.word_type(target)([i])
            else
                throw("Unsupported")
            end
        end
        PropertyT_new.Homomorphism(f, source, target)
    end

    relations_word_notation = filter((r->(length(r.first)+length(r.second)<12)), Groups.relations(SAut_F₄)) # removing the longest relations
    relations = [F24(r.first)*inv(F24(r.second)) for r in relations_word_notation]

    LowCohomologySOS.spectral_gaps_certification(quotient_hom, relations, halfBasis, optimizer = scs_opt(eps = 1e-5, max_iters = 100_000))
end
