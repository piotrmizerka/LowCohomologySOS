function _conj(
    t::Groups.Transvection,
    σ::PermutationGroups.AbstractPerm,
)
    return Groups.Transvection(t.id, t.i^σ, t.j^σ, t.inv)
end

function _conj(
    t::Groups.Transvection,
    x::Groups.Constructions.WreathProductElement,
)
    tσ = _conj(t, inv(x.p))
    dual_id = ifelse(t.id == :ϱ, :λ, :ϱ)
    dual_inv = ifelse(t.inv, false, true)
    new_id = isone(x.n.elts[t.i]) ? t.id : dual_id
    new_inv = isone(x.n.elts[t.i]*x.n.elts[t.j]) ? t.inv : dual_inv

    return Groups.Transvection(new_id, tσ.i, tσ.j, new_inv)
end

struct AlphabetPermutation{GEl,I} <: SymbolicWedderburn.ByPermutations
    perms::Dict{GEl,PermutationGroups.Perm{I}}
end

function AlphabetPermutation(A::Alphabet, G, op)
    return AlphabetPermutation(
        Dict(
            g => PermutationGroups.Perm([A[op(l, g)] for l in A.letters]) for
            g in G
        ),
    )
end

function Base.:^(
    w::Groups.AbstractWord,
    p::PermutationGroups.AbstractPerm,
)
    return typeof(w)([l^p for l in w])
end

function SymbolicWedderburn.action(
    act::AlphabetPermutation,
    g::Groups.GroupElement,
    gel,
)
    return parent(gel)(word(gel)^(act.perms[g]))
end

function subset_permutation(
    subset,
    g::Groups.GroupElement,
    act::AlphabetPermutation
)
    subset_idies = Dict(subset[i] => i for i in UInt32.(eachindex(subset)))

    return PermutationGroups.Perm([subset_idies[SymbolicWedderburn.action(act, g, subset[i])] for i in UInt32.(eachindex(subset))])
end

struct WedderburnActions{AP,CEH1,CEH2} <: SymbolicWedderburn.ByPermutations
    alphabet_perm::AP
    S_action::CEH1
    basis_action::CEH2
end

SymbolicWedderburn._int_type(::Type{<:InducedActionHomomorphism{<:WedderburnActions}}) = UInt32

function WedderburnActions(A::Alphabet, G, op, S, basis)
    alphabet_perm = AlphabetPermutation(
        Dict(
            g => PermutationGroups.Perm([A[op(l, g)] for l in A.letters]) for
            g in G
        ),
    )

    S_action = SymbolicWedderburn.CachedExtensionHomomorphism(
        G,
        alphabet_perm,
        S
    )

    basis_action = SymbolicWedderburn.CachedExtensionHomomorphism(
        G,
        alphabet_perm,
        basis
    )

    return WedderburnActions(alphabet_perm, S_action, basis_action)
end

# action on psd_basis ###########################
struct PSDBasisElement{GEl}
    generator::GEl # generator
    basis_elt::GEl # elt of the half-basis
end

Base.:(==)(a::PSDBasisElement, b::PSDBasisElement) =
    a.generator == b.generator && a.basis_elt == b.basis_elt
Base.hash(pbe::PSDBasisElement, h::UInt) =
    hash(pbe.generator, hash(pbe.basis_elt, h))

function SymbolicWedderburn.action(
    act::WedderburnActions,
    g::Groups.GroupElement,
    pbe::PSDBasisElement,
)
    return SymbolicWedderburn.action(act.alphabet_perm, g, pbe)
end

function SymbolicWedderburn.action(
    act::AlphabetPermutation,
    g::Groups.GroupElement,
    pbe::PSDBasisElement,
)
    s = SymbolicWedderburn.action(act, g, pbe.generator)
    g = SymbolicWedderburn.action(act, g, pbe.basis_elt)

    return PSDBasisElement(s, g)
end
#################################################

function SymbolicWedderburn.action(
    act::WedderburnActions,
    g::Groups.GroupElement,
    idx::LinIdx{T}
) where {T}
    N = length(basis(act.basis_action))
    n = length(basis(act.S_action))

    C_indicies = CartesianIndices((N, n, n))
    L_indices = LinearIndices(C_indicies)
    k, j, i = Tuple(C_indicies[idx.id])

    ig = i^SymbolicWedderburn.induce(act.S_action, g)
    jg = j^SymbolicWedderburn.induce(act.S_action, g)
    kg = k^SymbolicWedderburn.induce(act.basis_action, g)

    return LinIdx{T}(L_indices[kg, jg, ig])
end
#########################################################

struct LinIdx{T}
    id::T
end

function matrix_bases(basis, half_basis, S)
    n = length(S)^2 * length(basis)
    constr_basis = [LinIdx{UInt32}(i) for i in 1:n]
    psd_basis = [PSDBasisElement(s, g) for s in S for g in half_basis]

    return constr_basis, psd_basis
end

# struct TensorSupportElement{GEl}
#     i::GEl
#     j::GEl
#     k::GEl
# end

# Base.:(==)(s::TensorSupportElement, t::TensorSupportElement) =
#     s.i == t.i && s.j == t.j && s.k == t.k
# Base.hash(se::TensorSupportElement, h::UInt = UInt(0)) = hash(se.i, hash(se.j, hash(se.k, h)))

# function SymbolicWedderburn.action(
#     act::AlphabetPermutation,
#     g::Groups.GroupElement,
#     tse::TensorSupportElement,
# )
#     s = SymbolicWedderburn.action(act, g, tse.i)
#     t = SymbolicWedderburn.action(act, g, tse.j)
#     g = SymbolicWedderburn.action(act, g, tse.k)

#     return TensorSupportElement(s, t, g)
# end
