basis(A::StarAlgebras.StarAlgebra) = A.basis
basis(w::SymbolicWedderburn.WedderburnDecomposition) = w.basis

function coeffs(
    M::AbstractMatrix{<:AlgebraElement},
    n::Integer # generators number
)
    # basis = basis(parent(first(M))) # this ain't working somehow....
    basis = parent(first(M)).basis

    return [M[i,j](g) for i in 1:n for j in 1:n for g in basis]
end

function dot_fast(A::SparseMatrixCSC, B::SparseMatrixCSC)
    result = zero(eltype(B))

    for idx in findall(!iszero,A)
        JuMP.add_to_expression!(result, A[idx], B[idx])
    end

    return result
end

function invariant_constraint_matrix(
    v_inv::AbstractVector,
    A_gs_cart,
    m::Integer, # stands for half_basis size
    Σ_order::Integer
)
    bs = length(A_gs_cart)
    n = isqrt(div(length(v_inv), bs))
    result = spzeros(Int, m*n, m*n)
    dropzeros!(result)
    for it in SparseArrays.nonzeroinds(v_inv)
        i, j, k = id_2_triple(it, bs, n)
        out_left = (j-1)*m^2*n
        out_up = (i-1)*m
        for ci in A_gs_cart[k]
            e, f = ci[1], ci[2]
            in_left_f = (f-1)*m*n
            result[out_left+in_left_f+out_up+e] += 1
        end
    end

    return result/Σ_order
end

function sos_problem(
    M::AbstractMatrix{<:AlgebraElement},
    order_unit::AbstractMatrix{<:AlgebraElement},
    w_dec_matrix::SymbolicWedderburn.WedderburnDecomposition,
    Σ_order::Integer,
    upper_bound::Float64 = Inf,
)
    @assert size(M) == size(order_unit)
    @assert !isempty(M)

    A = parent(first(M))
    @assert all(x -> parent(x) === A, M)
    @assert all(x -> parent(x) === A, order_unit)

    result = JuMP.Model()

    P = map(SymbolicWedderburn.direct_summands(w_dec_matrix)) do ds
        dim = size(ds, 1)
        P = JuMP.@variable(result, [1:dim, 1:dim], Symmetric)
        JuMP.@constraint(result, P in JuMP.PSDCone())
        P
    end
    
    JuMP.@variable(result, λ)
    JuMP.@objective(result, Max, λ)

    if upper_bound < Inf
        JuMP.@constraint(result, λ <= upper_bound)
    end

    cnstrs = dropzeros!.(sparse.(constraints(A.mstructure)))
    A_gs_cart = [findall(!iszero,c) for c in cnstrs]

    @assert length(cnstrs) == length(basis(A))

    m = size(first(cnstrs))[1] # this is equal to the half_basis' size
    n = isqrt(div(length(basis(w_dec_matrix)), length(basis(A))))
    M_c = convert(Vector{eltype(w_dec_matrix)}, coeffs(M,n))
    U_c = convert(Vector{eltype(w_dec_matrix)}, coeffs(order_unit,n))
    
    it = 1
    Uπs_Pπs_degrees = [(sparse(SymbolicWedderburn.image_basis(ds)), sparse(P[π]), degree(ds)) 
                            for (π, ds) in pairs(SymbolicWedderburn.direct_summands(w_dec_matrix))]
    for v_inv in w_dec_matrix.invariants
        m_c = dot(M_c, v_inv)
        u_c = dot(U_c, v_inv)
        lhs = m_c-λ*u_c
        rhs = zero(typeof(lhs))
        matrix_inv = invariant_constraint_matrix(v_inv, A_gs_cart, m, Σ_order)
        for (Uπ, Pπ, deg) in Uπs_Pπs_degrees
            JuMP.add_to_expression!(rhs, deg, dot_fast(Uπ*matrix_inv*Uπ', Pπ))
        end
        JuMP.@constraint(result, lhs == rhs)

        if it%100_000 == 0
            @info it
        end
        it += 1
    end
    
    return result, P
end

function get_solution(
    m::JuMP.Model,
    P, # positive definite blocks
    w_dec_matrix::SymbolicWedderburn.WedderburnDecomposition
)
    λ = JuMP.value(m[:λ])
    Q = let 
        P_blocks = [JuMP.value.(P[π]) for (π, ds) in pairs(SymbolicWedderburn.direct_summands(w_dec_matrix))]
        P = PropertyT.reconstruct(P_blocks, w_dec_matrix)
        if any(isnan, P) || any(isinf, P)
            @error "obtained solution contains NaNs or ±Inf"
            P
        else
            real(sqrt(Symmetric((P + P')./2)))
        end
    end
    return λ, Q
end