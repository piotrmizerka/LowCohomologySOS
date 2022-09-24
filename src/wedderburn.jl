using Revise
using SymbolicWedderburn
using PermutationGroups
using Groups

struct MatrixAction <: SymbolicWedderburn.ByPermutations end
struct TensorSupportElement
    row_generator::Groups.GroupElement
    column_generator::Groups.GroupElement
    entry::Groups.GroupElement
end
Base.:(==)(s::TensorSupportElement, t::TensorSupportElement) =
    s.row_generator == t.row_generator && s.column_generator == t.column_generator && s.entry == t.entry
Base.hash(se::TensorSupportElement, h::UInt = UInt(0)) = hash(se.row_generator, hash(se.column_generator, hash(se.entry, h)))

function determine_transvection(
    g::Groups.GroupElement # has to be a genrating transvection or its inverse
)
    SAutFn = parent(g)
    Fn = SAutFn.group
    Fn_gens = Groups.gens(Fn)
    g_id, g_i, g_j, g_inv = 0, 0, 0, 0
    for k in eachindex(Fn_gens)
        if g(Fn_gens[k]) != Fn_gens[k]
            g_i = k
            for l in eachindex(Fn_gens)
                if g(Fn_gens[k]) == Fn_gens[k]*Fn_gens[l]
                    g_id = :ϱ
                    g_j = l
                    g_inv = false
                    break
                elseif g(Fn_gens[k]) == Fn_gens[k]*Fn_gens[l]^(-1)
                    g_id = :ϱ
                    g_j = l
                    g_inv = true
                    break
                elseif g(Fn_gens[k]) == Fn_gens[l]*Fn_gens[k]
                    g_id = :λ
                    g_j = l
                    g_inv = false
                    break
                elseif g(Fn_gens[k]) == Fn_gens[l]^(-1)*Fn_gens[k]
                    g_id = :λ
                    g_j = l
                    g_inv = true
                    break
                end
            end
            break    
        end
    end

    return Groups.Transvection(g_id, g_i, g_j, g_inv)
end

id_transvection_dict = Dict{Integer,Groups.Transvection}()
tranvection_id_dict = Dict{Groups.Transvection,Integer}()

################# SAut(F₂) example ###############################
SAutF₂ = Groups.SpecialAutomorphismGroup(FreeGroup(2))
S = let s = Groups.gens(SAutF₂)
    [s; inv.(s)]
end
S = unique!(S)

tranvection_id_dict = Dict([(determine_transvection(s), s.word[1]) for s in S])

# Caution: S and tranvection_id_dict are global variables - TO FIX!!

function act_on_transvection(
    σ::Groups.Constructions.WreathProductElement,
    g::Groups.GroupElement # has to be a genrating transvection or its inverse
)
    g_transvection = determine_transvection(g)
    g_id, g_i, g_j, g_inv = g_transvection.id, g_transvection.i, g_transvection.j, g_transvection.inv

    a_sigma_j = (σ.n.elts[σ.p[g_j]] == parent(first(σ.n.elts))() ? 0 : 1)

    if g_id == :ϱ
        if σ.n.elts[σ.p[g_i]] == parent(first(σ.n.elts))()
            result_transvection = Groups.Transvection(:ϱ, σ.p[g_i], σ.p[g_j], a_sigma_j)
        else
            result_transvection = Groups.Transvection(:λ, σ.p[g_i], σ.p[g_j], (a_sigma_j+1)%2)
        end
    else
        if σ.n.elts[σ.p[g_i]] == parent(first(σ.n.elts))()
            result_transvection = Groups.Transvection(:λ, σ.p[g_i], σ.p[g_j], a_sigma_j)
        else
            result_transvection = Groups.Transvection(:ϱ, σ.p[g_i], σ.p[g_j], (a_sigma_j+1)%2)
        end
    end
    
    t_inv = (g_inv == true ? -1 : 1)
    result_transvection = (t_inv == -1 ? result_transvection^(-1) : result_transvection)

    return S[tranvection_id_dict[result_transvection]] # check if indicies are consistent with S!!
end

function wreath_conjugation(
    σ::Groups.Constructions.WreathProductElement,
    g::Groups.GroupElement
)
    SAutFn = parent(g)
    result = one(SAutFn)
    word = g.word
    for i in 1:length(word)
        result = result*act_on_transvection(σ, SAutFn(word[i:i]))
    end
    return result
end

# function AutFG_emb(A::AutomorphismGroup, g::Groups.Constructions.WreathProductElement)
#     isa(A.group, FreeGroup) || throw("Not an Aut(Fₙ)")
#     parent(g).P.n == length(A.group.gens) || throw("No natural embedding of $(parent(g)) into $A")
#     powers = [(elt == parent(elt)() ? 0 : 1) for elt in g.n.elts]
#     elt = reduce(*, [A(Groups.flip_AutSymbol(i))^pow for (i,pow) in enumerate(powers)])
#     Groups.r_multiply!(elt, [Groups.perm_autsymbol(g.p)])
#     return elt
# end

# function wreath_conjugation(
#     g::Groups.Constructions.WreathProductElement,
#     a::Groups.GroupElement
# )
#     g = AutFG_emb(parent(a),g)
#     return g*a*g^-1
# end

# Our action of Σ on the basis {sᵢ⊗eᵢ|sᵢ∈S, eᵢ∈E} is given by the tensor representation:
# Σ→GLₙₘ(ℝ), g↦ϕ(g)⊗σ(g), where σ:Σ→GLₘ(ℝ) is the permutation representation defining the action
# of Σ on the set E.
function SymbolicWedderburn.action(
    ::MatrixAction,
    g::Groups.Constructions.WreathProductElement,
    tensor_support_element::TensorSupportElement,
)
    return TensorSupportElement(
        wreath_conjugation(g, tensor_support_element.row_generator),
        wreath_conjugation(g, tensor_support_element.column_generator),
        wreath_conjugation(g, tensor_support_element.entry)
    )
end

function matrix_bases(
    basis, 
    half_basis,
    S # stands for the generating set we choose
)
    constraints_basis = TensorSupportElement[]
    psd_basis = TensorSupportElement[]

    for i in eachindex(S)
        for j in eachindex(S)
            for e in basis
                push!(constraints_basis, TensorSupportElement(S[i], S[j], e))
            end
            if i == j
                for e in half_basis
                    push!(psd_basis, TensorSupportElement(S[i], S[j], e)) # half_basis on diagonal
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
    action = MatrixAction()
    constraints_basis, psd_basis = matrix_bases(basis, half_basis, S)

    # @info length(constraints_basis)
    # @info length(psd_basis)

    return SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end


Σ = let G = PermutationGroups.SymmetricGroup(2),
    P = PermutationGroups.SymmetricGroup(2)
    Groups.Constructions.WreathProduct(G, P)
end

const half_radius = 1

basis, sizes = Groups.wlmetric_ball(S, radius = 2*half_radius)
half_basis = basis[1:sizes[half_radius]]

w_dec_matrix = wedderburn_decomposition_matrix(Σ, basis, half_basis, S)
