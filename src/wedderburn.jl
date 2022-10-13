using Revise
using SymbolicWedderburn
using PermutationGroups
using Groups
using LowCohomologySOS

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
    basis,
    half_basis,
    S # stands for the generating set we choose
)
    action = LowCohomologySOS.AlphabetPermutation(alphabet(parent(first(S))), Σ, _conj)
    constraints_basis, psd_basis = matrix_bases(basis, half_basis, S)

    return SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end


const n = 2

SAutFn = Groups.SpecialAutomorphismGroup(FreeGroup(n))
Σ = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(n))
S = let s = Groups.gens(SAutFn)
    [s; inv.(s)]
end
S = unique!(S)

const half_radius = 1

basis, sizes = Groups.wlmetric_ball(S, radius = 2*half_radius)
half_basis = basis[1:sizes[half_radius]]

w_dec_matrix = wedderburn_decomposition_matrix(Σ, basis, half_basis, S)


# check if the convention for the group operation in SAutFns is to compose the automorphisms 
# first right then left (the classical way) or the other way round.
S[1]
S[2]
elt = S[1]*S[2]
s1, s2 = Groups.domain(S[1])
evaluate(elt)
evaluate(elt) == (s1*s2, s2*s1*s2) # first right then left
evaluate(elt) == (s1*s2*s1, s2*s1) # first left then right
# It turns out that the convention is the classical one: compose automorphisms first
# from the right then from the left.