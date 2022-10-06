struct AlphabetPermutation{GEl,I} <: SymbolicWedderburn.ByPermutations
    perms::Dict{GEl,PermutationGroups.Perm{I}}
end

function AlphabetPermutation(
    A::Alphabet,
    Γ::PermutationGroups.AbstractPermutationGroup,
    op,
)
    return AlphabetPermutation(
        Dict(γ => inv(PermutationGroups.Perm([A[op(l, γ)] for l in A])) for γ in Γ),
    )
end

function AlphabetPermutation(A::Alphabet, W::Groups.Constructions.WreathProduct, op)
    return AlphabetPermutation(
        Dict(
            # w => inv(PermutationGroups.Perm([A[op(op(l, w.n), w.p)] for l in A.letters])) for
            w => inv(PermutationGroups.Perm([A[op(l, w)] for l in A.letters])) for # I think this is the right way to do - no need to compose op twice as in the line above
            w in W
        ),
    )
end

function Base.:^(
    w::Groups.AbstractWord,
    p::PermutationGroups.AbstractPerm,
)
    return typeof(w)([l^p for l in w])
end

struct TensorSupportElement
    row_generator::Groups.GroupElement
    column_generator::Groups.GroupElement
    entry::Groups.GroupElement
end
Base.:(==)(s::TensorSupportElement, t::TensorSupportElement) =
    s.row_generator == t.row_generator && s.column_generator == t.column_generator && s.entry == t.entry
Base.hash(se::TensorSupportElement, h::UInt = UInt(0)) = hash(se.row_generator, hash(se.column_generator, hash(se.entry, h)))

function SymbolicWedderburn.action(
    act::AlphabetPermutation,
    γ::Groups.Constructions.WreathProductElement,
    tensor_support_element::TensorSupportElement,
)
    G = parent(tensor_support_element.column_generator)
    w_row = word(tensor_support_element.row_generator)^(act.perms[γ])
    w_col = word(tensor_support_element.column_generator)^(act.perms[γ])
    w_entry = word(tensor_support_element.entry)^(act.perms[γ])
    
    return TensorSupportElement(
        G(w_row),
        G(w_col),
        G(w_entry)
    )
end
