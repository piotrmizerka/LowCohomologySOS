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
    gsⱼ⁻¹ = [[] for k in 1:length(bigger_ball)]
    sᵢg = [[] for k in 1:length(bigger_ball)]
    sᵢgsⱼ⁻¹ = [[] for k in 1:length(bigger_ball)]

    for k in 1:length(bigger_ball)
        for i in 1:length(gens(G))
            gsⱼ⁻¹_elt = bigger_ball[k]*gens(G,i)^(-1)
            gsⱼ⁻¹_id = findall(x -> x==gsⱼ⁻¹_elt, smaller_ball)
            if length(gsⱼ⁻¹_id) == 1 
                append!(gsⱼ⁻¹[k], [[i, gsⱼ⁻¹_id[1]]])
            end

            sᵢg_elt = gens(G,i)*bigger_ball[k]
            sᵢg_id = findall(x -> x==sᵢg_elt, smaller_ball)
            if length(sᵢg_id) == 1 
                append!(sᵢg[k], [[i, sᵢg_id[1]]])
            end

            for j in 1:length(gens(G))
                sᵢgsⱼ⁻¹_elt = gens(G,i)*bigger_ball[k]*gens(G,j)^(-1)
                sᵢgsⱼ⁻¹_id = findall(x -> x==sᵢgsⱼ⁻¹_elt, smaller_ball)
                if length(sᵢgsⱼ⁻¹_id) == 1 
                    append!(sᵢgsⱼ⁻¹[k], [[i, j, sᵢgsⱼ⁻¹_id[1]]])
                end
            end
        end
    end
    
    return gsⱼ⁻¹, sᵢg, sᵢgsⱼ⁻¹
end

function sos_problem_conjugated(
    ξ,
    order_unit,
    upper_bound::Float64 = Inf;
    G,
    smaller_group_ring
)

    m = size(smaller_group_ring.mstructure, 1)
    n = length(gens(G))
    mn = m * n

    result = JuMP.Model()

    JuMP.@variable(result, P[1:mn, 1:mn], Symmetric)
    JuMP.@constraint(result, sdp, P in JuMP.PSDCone())

    if upper_bound < Inf
        λ = JuMP.@variable(result, λ <= upper_bound)
    else
        λ = JuMP.@variable(result, λ)
    end

    cnstrs = constraints(smaller_group_ring.mstructure)

    bigger_group_ring = parent(ξ)
    gsⱼ⁻¹, sᵢg, sᵢgsⱼ⁻¹ = associated_elements(G, bigger_group_ring.basis, smaller_group_ring.basis)

    u = StarAlgebras.coeffs(order_unit)
    @assert length(StarAlgebras.coeffs(ξ)) == length(u) == length(bigger_group_ring.basis)
    @assert length(cnstrs) == length(smaller_group_ring.basis)
    JuMP.@constraint(
        result,
        [k = 1:length(u)],
        ξ[k] - λ * u[k] ==
        ((k <= length(cnstrs)) ? sum(sum(P[entry_constraint(cnstrs, i, j, k, m, n)]) for i in 1:n for j in 1:n) : 0) -
        sum(sum(P[entry_constraint(cnstrs, i, jk[1], jk[2], m, n)]) for jk in gsⱼ⁻¹[k] for i in 1:n)-
        sum(sum(P[entry_constraint(cnstrs, ik[1], j, ik[2], m, n)]) for ik in sᵢg[k] for j in 1:n)+
        sum(sum(P[entry_constraint(cnstrs, ijk[1], ijk[2], ijk[3], m, n)]) for ijk in sᵢgsⱼ⁻¹[k])
    )

    JuMP.@objective(result, Max, λ)

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
