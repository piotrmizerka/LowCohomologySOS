##### Certification - see certify_SOS_decomposition from Marek's code (Property T) - in 1712.07167.jl file

function sos_from_matrix(Q, support, RG::StarAlgebra)
    mn = size(Q,1)

    # Changing Q to the corresponding interval-entry matrix
    Q_interval = reshape([@interval(0) for i in 1:(mn*mn)], mn, mn)
    for i in 1:mn
        for j in i:mn
            Q_interval[i,j] = @interval(Q[i,j])
            Q_interval[j,i] = Q_interval[i,j]
        end
    end

    P_interval = Q_interval^2
    P_interval_RG = reshape([@interval(0)*zero(RG) for i in 1:(mn*mn)], mn, mn)
    for i in 1:mn
        for j in 1:mn
            P_interval_RG[i,j] = P_interval[i,j]*one(RG)
        end
    end

    m = length(support)
    n = floor(Int, mn/m)
    Iₙ = reshape([zero(RG) for i in 1:(n*n)], n, n)
    for i in 1:n
        Iₙ[i,i] = one(RG)
    end

    x = reshape([RG(support[i]) for i in 1:length(support)], m, 1)
    xx = collect(Iₙ⊗x)
    result = xx'*P_interval_RG*xx

    return result
end

function certify_sos_decomposition(X, order_unit, λ::Number, Q::AbstractMatrix, support, RG::StarAlgebra)
    λ_interval = @interval(λ)
    eoi = X - λ_interval*order_unit

    residual = eoi - sos_from_matrix(Q, support, RG)
    l1_norm = 0
    mn = size(X,1)
    for i in 1:mn
        for j in 1:mn
            l1_norm += norm(residual[i,j],1)
        end
    end

    @info "l₁ norm of the error in interval arithmetic:"
    @info l1_norm

    result = λ_interval-l1_norm

    return result
end

function spectral_gaps_certification(h::Function, relations, half_basis)
    λₐₚ, Pₐₚ, RG, Δ₁, Iₙ = spectral_gaps_approximated(h, relations, half_basis)
    Qₐₚ = real(sqrt(Symmetric( (Pₐₚ.+ Pₐₚ')./2 )))

    @info "Approximated λ:"
    @info λₐₚ
    
    result = certify_sos_decomposition(Δ₁, Iₙ, λₐₚ, Qₐₚ, half_basis, RG)

    @info "Certified λ (interval atithmetic):"
    @info result

    return result
end
