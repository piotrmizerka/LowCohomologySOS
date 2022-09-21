using SymbolicWedderburn
using PermutationGroups

struct MatrixAction <: SymbolicWedderburn.ByPermutations end
struct TensorSupportElement
    row_generator
    column_generator
    entry
end
Base.:(==)(s::TensorSupportElement, t::TensorSupportElement) =
    s.generator == t.generator && s.entry == t.entry
Base.hash(se::TensorSupportElement, h::UInt = UInt(0)) = hash(se.generator_id, hash(se.entry, h))

# function permute_generator(
#     s, 
#     p::PermutationGroups.AbstractPerm
# )
#     return s^p
# end

# Our action of Σ on the basis {sᵢ⊗eᵢ|sᵢ∈S, eᵢ∈E} is given by the tensor representation:
# Σ→GLₙₘ(ℝ), g↦ϕ(g)⊗σ(g), where σ:Σ→GLₘ(ℝ) is the permutation representation defining the action
# of Σ on the set E.
function SymbolicWedderburn.action(
    ::MatrixAction,
    p::PermutationGroups.AbstractPerm,
    tensor_support_element::TensorSupportElement,
)
    # permuted_row_generator = permute_generator(tensor_support_element.row_generator, p) # TODO
    # permuted_column_generator = permute_generator(tensor_support_element.column_generator, p)
    # entry_after_action = act_on_entry(tensor_support_element.generator, p) # TODO
    return TensorSupportElement(
        permuted_row_generator^p,
        permuted_column_generator^p,
        entry_after_action^p
    )
end

function matrix_bases(basis, half_basis)
    constraints_basis = TensorSupportElement[]
    psd_basis = TensorSupportElement[]
    G = parent(first(support))
    for i in length(gens(G))
        for j in length(gens(G))
            for e in basis
                push!(constraints_basis, TensorSupportElement(gens(G,i), gens(G,j), e))
            end
            if i == j
                for e in half_basis
                    push!(psd_basis, TensorSupportElement(gens(G,i), gens(G,j), e)) # half_basis on diagonal
                end
            end
        end
    end
    return constraints_basis, psd_basis
end

function wedderburn_decomposition_matrix(
    Σ,
    G,
    basis,
    half_basis
)
    action = MatrixAction()
    constraints_basis, psd_basis = matrix_bases(basis, half_basis)
    return SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end


################# Example ###############################
using Groups

A = Alphabet([:a, :A, :b, :B], [2, 1, 4, 3])
F₂ = FreeGroup(A)
a, b = Groups.gens(F₂)

const half_radius = 1

# Π₀₁₋₁₋₁ ##############################################################################################
G = D₆ = S₃ = FPGroup(F₂, [a^3 => one(F₂), b^2 => one(F₂), b*a*b => b^(-1)])

S = let s = gens(G)
    [s; inv.(s)]
end
half_basis, sizes = Groups.wlmetric_ball(S, radius = half_radius)


