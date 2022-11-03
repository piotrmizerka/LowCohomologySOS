using JuMP # without this won't precompile (why???? - it is already present in the definition of LowCohomologySOS module, though through "import")

basis(A::StarAlgebras.StarAlgebra) = A.basis
basis(w::WedderburnDecomposition) = w.basis

function coeffs(
    M::AbstractMatrix{<:AlgebraElement},
    n::Integer # generators number
)
    # basis = basis(parent(first(M))) # this ain't working somehow....
    basis = parent(first(M)).basis

    return [M[i,j](g) for i in 1:n for j in 1:n for g in basis]
end

function invariant_constraint_matrix(
    v_inv::AbstractVector,
    A_gs,
    m::Integer # stands for half_basis size
)
    bs = length(A_gs)
    n = floor(Int, sqrt(div(length(v_inv), bs)))
    result = zeros(eltype(first(A_gs)[2]), m*n, m*n)
    for it in SparseArrays.nonzeroinds(v_inv)
        i, j, k = div(it-1,bs*n)+1, div((it-1)%(bs*n),bs)+1, (it-1)%bs+1
        result += v_inv[it]*KroneckerDelta{n}(i, j)⊗A_gs[k]
    end

    return result
end

function sos_problem(
    M::AbstractMatrix{<:AlgebraElement},
    order_unit::AbstractMatrix{<:AlgebraElement},
    w_dec_matrix::SymbolicWedderburn.WedderburnDecomposition,
    upper_bound::Float64 = Inf
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
        @constraint(result, P in JuMP.PSDCone())
        P
    end
    
    JuMP.@variable(result, λ)
    JuMP.@objective(result, Max, λ)

    if upper_bound < Inf
        JuMP.@constraint(result, λ <= upper_bound)
    end

    cnstrs = constraints(A.mstructure)
    @assert length(cnstrs) == length(basis(A))

    m = size(first(cnstrs))[1] # this is equal to the half_basis' size
    n = floor(Int, sqrt(div(length(basis(w_dec_matrix)), length(basis(A)))))
    M_c = convert(Vector{eltype(w_dec_matrix)}, coeffs(M,n))
    U_c = convert(Vector{eltype(w_dec_matrix)}, coeffs(order_unit,n))
    
    it = 1
    for v_inv in w_dec_matrix.invariants
        m_c = dot(M_c, v_inv)
        u_c = dot(U_c, v_inv)
        lhs = m_c-λ*u_c
        rhs = zero(typeof(lhs))
        matrix_inv = invariant_constraint_matrix(v_inv, cnstrs, m)
        # for (π, ds) in pairs(SymbolicWedderburn.direct_summands(w_dec_matrix))
        #     Uπ = SymbolicWedderburn.image_basis(ds)
        #     rhs += degree(ds)*dot(Uπ*matrix_inv*Uπ', P[π])
        # end
        JuMP.@constraint(result, lhs == rhs)

        @info it
        it += 1
    end
    
    return result
end

function get_solution_symmetrized(
    m::JuMP.Model,
    w_dec_matrix::SymbolicWedderburn.WedderburnDecomposition
)
    λ = JuMP.value(m[:λ])
    Q = let 
        P_diag = JuMP.value.(m[:P])
        P_blocks = [P_diag[π] for (π, ds) in pairs(SymbolicWedderburn.direct_summands(w_dec_matrix))]
        P = PropertyT_new.reconstruct(P_blocks, w_dec_matrix)
        if any(isnan, P) || any(isinf, P)
            @error "obtained solution contains NaNs or ±Inf"
            P
        else
            real(sqrt(Symmetric((P + P')./2)))
        end
    end
    return λ, Q
end