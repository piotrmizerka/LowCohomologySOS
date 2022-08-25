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

# defining the action (by permutations in this case)
struct OnLetters <: SymbolicWedderburn.ByPermutations end
function SymbolicWedderburn.action(
    ::OnLetters,
    p::PermutationGroups.AbstractPerm,
    w::Word,
)
    return Word(w.alphabet, [l^p for l in w.letters])
end
action = OnLetters()

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

sa = SymbolicWedderburn.symmetry_adapted_basis(Rational{Int}, Σ, action, E)

mπs = SymbolicWedderburn.multiplicity.(sa)

wdec = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, E, E)

Uπs = direct_summands(wdec)

pairs(Uπs)

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

P_invariant = PropertyT_new.reconstruct(P_blocks, wdec)

# second reconstruction example
P_trivial = [1]
P₂ = [2]
P_blocks = [P₂, P_trivial] # I noticed the order of direct summands is non-standard here

P_invariant = PropertyT_new.reconstruct(P_blocks, wdec)

# third reconstruction example
P_trivial = [10]
P₂ = [1]
P_blocks = [P₂, P_trivial] # I noticed the order of direct summands is non-standard here

P_invariant = PropertyT_new.reconstruct(P_blocks, wdec)