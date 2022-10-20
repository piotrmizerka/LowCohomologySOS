struct AlphabetPermutation{GEl,I} <: SymbolicWedderburn.ByPermutations
    perms::Dict{GEl,PermutationGroups.Perm{I}}
end

function AlphabetPermutation(A::Alphabet, G::Groups.Constructions.WreathProduct, op)
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

struct TensorSupportElement{GEl}
    i::GEl
    j::GEl
    k::GEl
end

Base.:(==)(s::TensorSupportElement, t::TensorSupportElement) =
    s.i == t.i && s.j == t.j && s.k == t.k
Base.hash(se::TensorSupportElement, h::UInt = UInt(0)) = hash(se.i, hash(se.j, hash(se.k, h)))

function SymbolicWedderburn.action(
    act::AlphabetPermutation,
    g::Groups.GroupElement,
    gel,
)
    return parent(gel)(word(gel)^(act.perms[g]))
end

function SymbolicWedderburn.action(
    act::AlphabetPermutation,
    g::Groups.GroupElement,
    tse::TensorSupportElement,
)
    s = SymbolicWedderburn.action(act, g, tse.i)
    t = SymbolicWedderburn.action(act, g, tse.j)
    g = SymbolicWedderburn.action(act, g, tse.k)

    return TensorSupportElement(s, t, g)
end

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
    tσ = _conj(t, x.p)
    dual_id = ifelse(t.id == :ϱ, :λ, :ϱ)
    dual_inv = ifelse(t.inv, false, true)
    new_id = isone(x.n.elts[tσ.i]) ? t.id : dual_id
    new_inv = isone(x.n.elts[tσ.i]*x.n.elts[tσ.j]) ? t.inv : dual_inv

    return Groups.Transvection(new_id, tσ.i, tσ.j, new_inv)
end

function matrix_bases(basis, half_basis, S)
    constr_basis = [TensorSupportElement(s, t, g) for s in S for t in S for g in basis]
    psd_basis = [TensorSupportElement(s, s, g) for s in S for g in half_basis]
    
   return constr_basis, psd_basis
end
