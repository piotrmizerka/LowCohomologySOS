using LowCohomologySOS
using Groups
using PropertyT_new
using JuMP
using SCS
using Revise

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
