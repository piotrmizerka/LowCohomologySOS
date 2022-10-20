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

function matrix_bases(
    basis, 
    half_basis,
    S # stands for the generating set we choose
)
    constraints_basis = LowCohomologySOS.TensorSupportElement[]
    psd_basis = LowCohomologySOS.TensorSupportElement[]

    for i in eachindex(S)
        for j in eachindex(S)
            for e in basis
                push!(constraints_basis, LowCohomologySOS.TensorSupportElement(S[i], S[j], e))
            end
            if i == j
                for e in half_basis
                    push!(psd_basis, LowCohomologySOS.TensorSupportElement(S[i], S[j], e)) # half_basis on diagonal
                end
            end
        end
    end

    return constraints_basis, psd_basis
end

function wedderburn_decomposition_matrix(
    Σ,
    constraints_basis,
    psd_basis,
    S
)
    action = LowCohomologySOS.AlphabetPermutation(alphabet(parent(first(S))), Σ, _conj)
    # constraints_basis, psd_basis = matrix_bases(basis, half_basis, S)

    return SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end
