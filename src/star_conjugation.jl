"""
Different strategy of proving property (T). 

We are motivated by the following equality:

    Δ²-λΔ = d₀*(Δ₁-λIₙ)d₀

for an fpp group G=⟨s₁,…,sₙ|…⟩ and d₀ = [1-sᵢ]ₙₓ₁ ∈ Mₙ,₁(ℝG). 

If we manage to show that there exists an SOS matrix M ∈ Mₙ(ℝG) such that

    Δ²-λΔ = d₀*Md₀,

this would prove property (T) for G as *-conjugating preserves SOS (that's why d₀*Md₀ would be an SOS).
This approach has an advantage that we can look for smaller balls: in order to access the soultion
supported on a ball of radius r, we only have to consider ball of radius r-1. Obviously, this
approach won't cover the whole ball of radius r but the sole fact of the access to some of its elements
may lead to new results.

Technicalities:

Let r be the radius of the ball on which we would like to find an SOS (that is we look for
the support of the SOS to be contained in the ball of radius r).
Denote this ball by Bᵣ and let Bᵣ = {g₁,…,gₘ}. For proving propery (T) suffices then that
M is an SOS - that is, there exists a semi-definite positive matrix P ∈ Mₘₙ(ℝ) such that, 
for any g ∈ B₂ᵣ₊₂, one has 

    (Δ²-λΔ)(g) = (d₀*Md₀)(g) ⇔ (Δ²-λΔ)(g) = (d₀*𝕩*P𝕩d₀)(g)

where 𝕩 = Iₙ⊗[gᵢ]ₘₓ₁ ∈ Mₘₙ,ₙ(ℝG). This boils down to the following equations
to hold for any g ∈ B₂ᵣ₊₂:

    (Δ²-λΔ)(g) = ∑ᵢ,ⱼ(Mᵢ,ⱼ(g)-Mᵢ,ⱼ(gsⱼ⁻¹)-Mᵢ,ⱼ(sᵢg)+Mᵢ,ⱼ(sᵢgsⱼ⁻¹)) (*)

For any h ∈ B₂ᵣ we know how to compute the constraints on P corresponding to Mᵢ,ⱼ(h) (this is 
implemented in "positive_approx.jl" file) and on B₂ᵣ₊₂∖B₂ᵣ, the equations (*) tautologically hold
as both sides are not supported on B₂ᵣ₊₂∖B₂ᵣ. Therefore, we have to know which of the elements
of the form gsⱼ⁻¹, sᵢg, and sᵢgsⱼ⁻¹ are contained in B₂ᵣ (we handle this with "associated_elements"
function). 
"""

function associated_elements(
    G,
    bigger_ball,
    smaller_ball
)
    ran = 1:length(bigger_ball)
    gsⱼ⁻¹ = [NTuple{2, Int}[] for _ in ran]
    sᵢg = [NTuple{2, Int}[] for _ in ran]
    sᵢgsⱼ⁻¹ = [NTuple{3, Int}[] for _ in ran]

    for (k, g) in pairs(bigger_ball)
        for (i, s) in pairs(gens(G))
            idx = findfirst(x -> x==g*inv(s), smaller_ball)
            if idx ≠ nothing
                push!(gsⱼ⁻¹[k], (i, idx))
            end

            idx = findfirst(x -> x==s*g, smaller_ball)
            if idx ≠ nothing
                push!(sᵢg[k], (i, idx))
            end

            for (j, t) in pairs(gens(G))
                idx = findfirst(x -> x==s*g*inv(t), smaller_ball)
                if idx ≠ nothing
                    push!(sᵢgsⱼ⁻¹[k], (i, j, idx))
                end
            end
        end
    end

    return gsⱼ⁻¹, sᵢg, sᵢgsⱼ⁻¹
end

function sos_problem_conjugated(
    ξ::AlgebraElement,
    order_unit::AlgebraElement,
    upper_bound::Float64 = Inf;
    G::Group,
    smaller_group_ring::StarAlgebra,
)
    @assert parent(ξ) === parent(order_unit)

    m = size(smaller_group_ring.mstructure, 1)
    n = ngens(G)
    mn = m * n

    result = JuMP.Model()

    JuMP.@variable(result, P[1:mn, 1:mn], Symmetric)
    JuMP.@constraint(result, sdp, P in JuMP.PSDCone())

    λ = JuMP.@variable(result, λ)
    JuMP.@objective(result, Max, λ)

    if isfinite(upper_bound)
        JuMP.@constraint(result, λ <= upper_bound)
    end

    cnstrs = constraints(smaller_group_ring.mstructure)

    gsⱼ⁻¹, sᵢg, sᵢgsⱼ⁻¹ =
        associated_elements(G, basis(parent(ξ)), basis(smaller_group_ring))

    for (k, g) in pairs(basis(parent(ξ)))
        rhs = if k <= length(cnstrs)
            sum(
                sum(P[entry_constraint(cnstrs, i, j, k, m, n)]) for i = 1:n
                for j = 1:n
            )
        else
            zero(JuMP.AffExpr)
        end

        if !isempty(gsⱼ⁻¹[k])
            rhs -= sum(
                P[p] for i = 1:n for (j, k) in gsⱼ⁻¹[k] for
                p in entry_constraint(cnstrs, i, j, k, m, n)
            )
        end

        if !isempty(sᵢg[k])
            rhs -= sum(
                P[p] for (i, k) in sᵢg[k] for j = 1:n for
                p in entry_constraint(cnstrs, i, j, k, m, n)
            )
        end

        if !isempty(sᵢgsⱼ⁻¹[k])
            rhs += sum(
                P[p] for (i, j, k) in sᵢgsⱼ⁻¹[k] for
                p in entry_constraint(cnstrs, i, j, k, m, n)
            )
        end

        JuMP.@constraint(result, ξ(g) - λ * order_unit(g) == rhs)
    end

    return result
end

function property_t_conjugated_approx(
    G,
    half_radius;
    optimizer
)
    S = gens(G)
    S = unique([S; inv.(S)])
    bigger_basis, sizes = Groups.wlmetric_ball(S, one(G), radius = 2*half_radius+2)
    ℝG_bigger = group_ring(G, bigger_basis, star_multiplication = true, additive_only = true)
    ℝG_smaller = group_ring(G, bigger_basis[1:sizes[half_radius]], star_multiplication = true)

    Δ = length(S)*one(ℝG_bigger)-sum(ℝG_bigger(s) for s in S)

    sos_problem = sos_problem_conjugated(Δ^2, Δ, G = G, smaller_group_ring = ℝG_smaller)

    λ, P, termination_status = sos_problem_solution(sos_problem; optimizer)

    return λ, P, termination_status
end

# TODO: certification in this setting
