struct KroneckerDelta{N,T} <: AbstractMatrix{T}
    i::Int
    j::Int
    val::T

    function KroneckerDelta{N,T}(i, j, val) where {N,T}
        @assert 1 ≤ i ≤ N
        @assert 1 ≤ j ≤ N
        return new{N,T}(i, j, val)
    end
end

KroneckerDelta{N}(i::Integer, j::Integer) where {N} =
    KroneckerDelta{N,Int}(i, j, 1)

Base.size(δ::KroneckerDelta{N}) where {N} = (N, N)
Base.@propagate_inbounds function Base.getindex(
    δ::KroneckerDelta{N},
    i::Integer,
    j::Integer,
) where {N}
    @boundscheck checkbounds(δ, i, j)
    return ifelse(i == δ.i && j == δ.j, δ.val, zero(δ.val))
end

function _getindex_indices(a::AbstractMatrix, δ::KroneckerDelta{N}) where N
    m1,r1 = divrem(size(a, 1), N)
    m2,r2 = divrem(size(a, 2), N)

    @assert r1 == 0
    @assert r2 == 0

    rb = (δ.i-1)*m1 + 1
    re = rb + m1 - 1
    cb = (δ.j-1)*m2 + 1
    ce = cb + m2 - 1

    @boundscheck begin
        checkbounds(a, rb, cb)
        checkbounds(a, re, ce)
    end

    return rb:re, cb:ce
end

Base.@propagate_inbounds function Base.getindex(a::AbstractMatrix, δ::KroneckerDelta{N}) where N
    rows, cols = _getindex_indices(a, δ)
    return δ.val * a[rows, cols]
end

function Base.view(a::AbstractMatrix, δ::KroneckerDelta{N}) where N
    rows, cols = _getindex_indices(a, δ)

    isone(δ.val) || error("View with KroneckerDelta is defined only for value `1`, got $(δ.val)")

    return @view a[rows, cols]
end

struct Tensor{T,A,B} <: AbstractMatrix{T}
    a::A
    b::B
end

Tensor(a::AbstractMatrix{T}, b::AbstractMatrix{S}) where {T,S} =
    Tensor{promote_type(T, S),typeof(a),typeof(b)}(a, b)
Base.size(t::Tensor) = size(t.a) .* size(t.b)

function _tensor_index(i, dim)
    d, r = divrem(i, dim)
    if 1 ≤ r < dim
        d = d+1
    else
        r = dim
    end
    return d, r
end

Base.@propagate_inbounds function Base.getindex(
    t::Tensor,
    i::Integer,
    j::Integer,
)
    @boundscheck checkbounds(t, i, j)

    d1, r1 = _tensor_index(i, size(t.b, 1))
    d2, r2 = _tensor_index(j, size(t.b, 2))

    return t.a[d1, d2]*t.b[r1, r2]
end

struct BinaryMatrix{T} <: AbstractMatrix{T}
    nzeros::Vector{Int} # list of indices
    n::Int # nrows
    m::Int # ncols
    val::T

    function BinaryMatrix(nzeros, n, m, val::T; sorted::Bool = false) where T
        @assert n ≥ 1
        @assert m ≥ 1

        if !sorted && !issorted(nzeros)
            sort!(nzeros)
        end
        !isempty(nzeros) && @assert first(nzeros) ≥ 1 && last(nzeros) ≤ n*m

        return new{T}(nzeros, n, m, val)
    end
end

Base.size(bm::BinaryMatrix) = (bm.n, bm.m)
Base.@propagate_inbounds function Base.getindex(
    bm::BinaryMatrix,
    i::Integer,
    j::Integer,
)
    li = LinearIndices(bm)
    idx = li[i,j]

    return idx ∈ bm.nzeros ? bm.val : zero(bm.val)
end

LinearAlgebra.dot(bm::BinaryMatrix, m::AbstractMatrix) =
    sum(bm.val*m[idx] for idx in bm.nzeros)
