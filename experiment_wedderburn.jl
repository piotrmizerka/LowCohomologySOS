# This is the code replicating the example from my sheet notes from 2022_07_16

using Revise
using SymbolicWedderburn
using PermutationGroups

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
using PropertyT_new

# Uπs - just to compare with the notes (the may be rescaled), we don't need them further, however, as the whole computation
# is enclosed in the "reconstruct" function from the PropertyT-new package
U_trivial = Uπs[2]    
U₂ = Uπs[1]

# first reconstruction example
P_trivial = [1]
P₂ = [1]
P_blocks = [P₂, P_trivial] # I noticed the order of direct summands is non-standard here

P_invariant = round.(PropertyT_new.reconstruct(P_blocks, wdec), digits=3)

# second reconstruction example
P_trivial = [1]
P₂ = [2]
P_blocks = [P₂, P_trivial] # I noticed the order of direct summands is non-standard here

P_invariant = round.(PropertyT_new.reconstruct(P_blocks, wdec), digits=3)

# third reconstruction example
P_trivial = [10]
P₂ = [1]
P_blocks = [P₂, P_trivial] # I noticed the order of direct summands is non-standard here

P_invariant = round.(PropertyT_new.reconstruct(P_blocks, wdec), digits=3)


# Matrix example

# defining the action (by matrix entries permutations tensored with the action induced by letter permutations)
struct MatrixAction <: SymbolicWedderburn.ByPermutations end
struct SingleEntry
    row_id::Integer
    col_id::Integer
    entry::Word
end
Base.:(==)(s::SingleEntry, t::SingleEntry) =
    s.row_id == t.row_id && s.col_id == t.col_id && s.entry == t.entry
Base.hash(se::SingleEntry, h::UInt = UInt(0)) =
    hash(se.row_id, hash(se.col_id, hash(se.entry, h)))


function permute_C₂(i::Integer, p::PermutationGroups.AbstractPerm)::Integer
    if (1^p == 1 && 2^p == 2 && 3^p == 3) || (1^p == 2 && 2^p == 1 && 3^p == 3) # check if inside C₂ subgroup
        return i^p
    end
    return i
end
function SymbolicWedderburn.action(
    ::MatrixAction,
    p::PermutationGroups.AbstractPerm,
    single_entry::SingleEntry,
)
    permuted_row_id = permute_C₂(single_entry.row_id, p)
    perumuted_col_id = permute_C₂(single_entry.col_id, p)
    return SingleEntry(
        # single_entry.row_id^p, # this shall be good for the desired action, not for this example
        # single_entry.col_id^p,
        permuted_row_id,
        perumuted_col_id,
        Word(single_entry.entry.alphabet, [l^p for l in single_entry.entry.letters])
    )
end
matrix_action = MatrixAction()

matrix_basis = let 
    matrix_basis = SingleEntry[]
    for i in [1,2]
        for j in [1,2]
            for e in E
                push!(matrix_basis, SingleEntry(i, j, e))
            end
        end
    end
    matrix_basis
end

sa_matrix = SymbolicWedderburn.symmetry_adapted_basis(Rational{Int}, Σ, matrix_action, matrix_basis)
