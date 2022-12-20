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

function preprocess_actions(
    S, 
    basis, 
    G,
    act::AlphabetPermutation
)
    S_action = Dict(g => subset_permutation(S, g, act) for g in G)
    basis_action = Dict(g => subset_permutation(basis, g, act) for g in G)

    return S_action, basis_action
end

struct WedderburnActions <: SymbolicWedderburn.ByPermutations
    alphabet_perm::AlphabetPermutation
    S_action::Dict
    basis_action::Dict
    S
    basis
    S_size::Integer
    basis_size::Integer
end

function WedderburnActions(A::Alphabet, G, op, S, basis)
    alphabet_perm = AlphabetPermutation(
        Dict(
            g => PermutationGroups.Perm([A[op(l, g)] for l in A.letters]) for
            g in G
        ),
    )
    S_action, basis_action = preprocess_actions(S, basis, G, alphabet_perm)

    return WedderburnActions(alphabet_perm, S_action, basis_action, S, basis, UInt32(length(S)), UInt32(length(basis)))
end

# action on psd_basis ###########################
struct PSDBasisElement{GEl}
    s::GEl
    g::GEl
end

Base.:(==)(a::PSDBasisElement, b::PSDBasisElement) =
    a.s == b.s && a.g == b.g
Base.hash(pbe::PSDBasisElement, h::UInt = UInt(0)) = hash(pbe.s, hash(pbe.g, h))

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

# action on constraints_basis ###################
function id_2_triple(
    id::Integer,
    bs::Integer, # stands for basis_size
    n::Integer # stands for generators' number
)
    _1 = UInt32(1)
    return div(id-_1,bs*n)+_1, div((id-_1)%(bs*n),bs)+_1, (id-_1)%bs+_1
end

function triple_2_id(
    i::Integer,
    j::Integer,
    k::Integer,
    bs::Integer,
    n::Integer
)
    _1 = UInt32(1)
    return (i-_1)*n*bs+(j-_1)*bs+k
end

function SymbolicWedderburn.action(
    act::WedderburnActions,
    g::Groups.GroupElement,
    id::Integer
)
    i, j, k = id_2_triple(id, act.basis_size, act.S_size)
    iσ, jσ, kσ = act.S_action[g][i], act.S_action[g][j], act.basis_action[g][k]

    return triple_2_id(iσ, jσ, kσ, act.basis_size, act.S_size)
end
#########################################################

function matrix_bases(basis, half_basis, S)
    constr_basis = [i for i in UInt32(1):UInt32(length(S))^2*UInt32(length(basis))]
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
