using JuMP # without this won't precompile (why???? - it is already present in the definition of LowCohomologySOS module, though through "import")

function invariant_constraint_matrix(
    orbit_representative::TensorSupportElement,
    Σ::Groups.Constructions.WreathProduct,
    action,
    psd_basis,
    A_gs::Dict,
    generators_number::Integer
)
    dim = length(psd_basis)
    result = zeros(eltype(first(A_gs)[2]), dim, dim)
    for σ ∈ Σ
        σ_orbit_representative = SymbolicWedderburn.action(action, σ, orbit_representative)
        σi, σj = word(σ_orbit_representative.i)[1], word(σ_orbit_representative.j)[1]
        δ_σi_σj = KroneckerDelta{generators_number}(σi, σj)
        A_σg = A_gs[σ_orbit_representative.k]
        result += δ_σi_σj⊗A_σg
    end
    result /= length(collect(Σ))

    return result
end

function sos_problem_symmetrized(
    M::AbstractMatrix{<:AlgebraElement},
    order_unit::AbstractMatrix{<:AlgebraElement},
    w_dec_matrix::SymbolicWedderburn.WedderburnDecomposition,
    Σ::Groups.Constructions.WreathProduct,
    psd_basis, # try to deal without this parameter!
    S, # try to deal without this parameter!
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
        λ = JuMP.@constraint(result, λ <= upper_bound)
    end

    cnstrs = constraints(A.mstructure)
    @assert length(cnstrs) == length(basis(A))

    constraints_basis = w_dec_matrix.basis
    considered_constraints = Dict{eltype(constraints_basis), Bool}()
    for constraint in constraints_basis
        considered_constraints[constraint] = false
    end 

    action = AlphabetPermutation(alphabet(parent(first(S))), Σ, _conj)
    A_gs = Dict(g => A_g for (A_g, g) in zip(cnstrs, basis(A)))

        i, j = word(constraint.i)[1], word(constraint.j)[1]
        k = basis_idies[constraint.k]
        if considered_constraints[i,j,k] == false
            i, j = word(constraint.i)[1], word(constraint.j)[1]
            lhs = M[i,j](constraint.k)
                Uπ = SymbolicWedderburn.image_basis(ds)
                rhs += dot(degree(ds)*Uπ*δ_i_j_g_invariant*Uπ', P[π])
            end
            JuMP.@constraint(result, lhs == rhs)
            for σ ∈ Σ
                σ_i, σ_j = word(σ_constraint.i)[1], word(σ_constraint.j)[1]
                σ_k = basis_idies[σ_constraint.k]
                considered_constraints[σ_i,σ_j,σ_k] = true
            end
        end
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