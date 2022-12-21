using Revise
using SymbolicWedderburn
using PermutationGroups
using PropertyT

struct Word{T}
    alphabet::Vector{T}
    letters::Vector{Int}

    function Word(a::AbstractVector{T}, l::AbstractVector{<:Integer}) where {T}
        all(i -> 1 <= i <= length(a), l) ||
            throw(ArgumentError("Invalid word over alphabet $a: $w"))
        return new{T}(a, l)
    end
end
Base.show(io::IO, w::Word) = join(io, w.alphabet[w.letters], "·")
Base.:(==)(w::Word, v::Word) =
    w.alphabet == v.alphabet && w.letters == v.letters
Base.hash(w::Word, h::UInt = UInt(0)) =
    hash(w.alphabet, hash(w.letters, hash(Word, h)))

# defining the symmetry group Σ
Σ = S₃ = PermGroup(perm"(1,2,3)", perm"(1,2)")

# defining the set on which Σ acts
E = words = let A = [:x, :y, :z], radius = 1
    words = [Word(A, [1]), Word(A, [2]), Word(A, [3])]
    for r in 2:radius
        append!(
            words,
            [
                Word(A, collect(w)) for
                w in Iterators.product(fill(1:3, r)...)
            ],
        )
    end
    words
end

#####################################################################################################

# classical, non-matrix example - compare with the notes from 2022_07_16

# defining the action (by letter permutations in this case)
struct OnLetters <: SymbolicWedderburn.ByPermutations end
function SymbolicWedderburn.action(
    ::OnLetters,
    p::PermutationGroups.AbstractPerm,
    w::Word,
)
    return Word(w.alphabet, [l^p for l in w.letters])
end
action = OnLetters()

sa = SymbolicWedderburn.symmetry_adapted_basis(Rational{Int}, Σ, action, E)

mπs = SymbolicWedderburn.multiplicity.(sa)

wdec = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, E, E)

Uπs = direct_summands(wdec)

# let's compute the invariant matrix from sample block ones (see the notes pp. 9 and 10)

# Uπs - just to compare with the notes (the may be rescaled), we don't need them further, however, as the whole computation
# is enclosed in the "reconstruct" function from the PropertyT-new package
U_trivial = Uπs[2]    
U₂ = Uπs[1]

# first reconstruction example
P_trivial = [1]
P₂ = [1]
P_blocks = [P₂, P_trivial] # I noticed the order of direct summands is non-standard here

P_invariant = round.(PropertyT.reconstruct(P_blocks, wdec), digits=3)

# second reconstruction example
P_trivial = [1]
P₂ = [2]
P_blocks = [P₂, P_trivial] # I noticed the order of direct summands is non-standard here

P_invariant = round.(PropertyT.reconstruct(P_blocks, wdec), digits=3)

# third reconstruction example
P_trivial = [10]
P₂ = [1]
P_blocks = [P₂, P_trivial] # I noticed the order of direct summands is non-standard here

P_invariant = round.(PropertyT.reconstruct(P_blocks, wdec), digits=3)


# Matrix example - in the similar manner, compare with the notes "2022_07_28_matrix_Wedderburn"

# defining the action (by matrix entries permutations tensored with the action induced by letter permutations)
struct MatrixAction <: SymbolicWedderburn.ByPermutations end
struct TensorSupportElement
    generator_id::Integer
    entry::Word
end
Base.:(==)(s::TensorSupportElement, t::TensorSupportElement) =
    s.generator_id == t.generator_id && s.entry == t.entry
Base.hash(se::TensorSupportElement, h::UInt = UInt(0)) = hash(se.generator_id, hash(se.entry, h))

# denoting Σ={1,a,a²,b,ba,ba²}, we define the representation ϕ:Σ→GL₂(ℝ) by a↦[1 0; 0 1], b↦[0 1; 1 0]
function permute_C₂(i::Integer, p::PermutationGroups.AbstractPerm)::Integer 
    id_perm = p^6
    if i∈[1,2] && p!=id_perm && p^2==id_perm
        return i%2+1
    end
    return i
end

# Our action of Σ on the basis {1⊗x,1⊗y,1⊗z,2⊗x,2⊗y,2⊗z} is given by the tensor representation:
# Σ→GL₆(ℝ), g↦ϕ(g)⊗σ(g), where σ:Σ→GL₃(ℝ) is the permutation representation defining the action
# of Σ on the set E={x,y,z}.
function SymbolicWedderburn.action(
    ::MatrixAction,
    p::PermutationGroups.AbstractPerm,
    tensor_support_element::TensorSupportElement,
)
    permuted_generator_id = permute_C₂(tensor_support_element.generator_id, p)
    return TensorSupportElement(
        permuted_generator_id,
        Word(tensor_support_element.entry.alphabet, [l^p for l in tensor_support_element.entry.letters])
    )
end
matrix_action = MatrixAction()

matrix_basis = let 
    matrix_basis = TensorSupportElement[]
    for i in [1,2]
        for e in E
            push!(matrix_basis, TensorSupportElement(i, e))
        end
    end
    matrix_basis
end

sa_matrix = SymbolicWedderburn.symmetry_adapted_basis(Rational{Int}, Σ, matrix_action, matrix_basis)

mπs_matrix = SymbolicWedderburn.multiplicity.(sa_matrix)

wdec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, matrix_action, matrix_basis, matrix_basis)

Uπs_matrix = direct_summands(wdec_matrix)
U_trivial_matrix = Uπs_matrix[2]
U₁_matrix = Uπs_matrix[3]
U₂_matrix = Uπs_matrix[1]

P_trivial_matrix = [1]
P₁_matrix = [10]
P₂_matrix = [1 0;-1 2]
P_blocks_matrix = [P₂_matrix, P_trivial_matrix, P₁_matrix]

P_invariant_matrix = round.(PropertyT.reconstruct(P_blocks_matrix, wdec_matrix), digits=2)