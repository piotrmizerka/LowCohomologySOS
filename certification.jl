using Pkg
Pkg.activate(@__DIR__)

using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = 4
LinearAlgebra.BLAS.set_num_threads(2)
using Kronecker
using IntervalArithmetic

using JuMP
using SCS
using ProxSDP

include("SOSDecompositions.jl")

##### CERTIFICATION - see certify_SOS_decomposition from Marek's code (Property T) - in 1712.07167.jl file
# Base.adjoint(X::AlgebraElement) = StarAlgebra.star(X)
function SOSFromMatrix(Q, support, RG, G)
    mn = size(Q)[1]

    # Changing Q to the corresponding interval-entry matrix
    QInterval = reshape([@interval(0) for i in 1:(mn*mn)], mn, mn)
    for i in 1:mn
        for j in 1:mn
            QInterval[i,j] = @interval(Q[i,j])
        end
    end

    QIntervalᵀ = QInterval' 
    PInterval = QIntervalᵀ*QInterval
    # PInterval = QInterval^2 # we don't need Cholesky decomposition - we can take normal squares insted
    PIntervalRG = reshape([@interval(0)*RG(0) for i in 1:(mn*mn)], mn, mn)
    for i in 1:mn
        for j in 1:mn
            PIntervalRG[i,j] = PInterval[i,j]*RG(one(G))
        end
    end

    m = length(support)
    n = floor(Int, mn/m)

    Iₙ = reshape([RG(0) for i in 1:(n*n)], n, n)
    for i in 1:n
        Iₙ[i,i] = RG(one(G))
    end

    x = reshape([RG(support[i]) for i in 1:length(support)], m, 1)
    xx = collect(Iₙ⊗x)
    xxᵀ = starOfMatrixOverGroupRing(xx)

    result = xxᵀ*PIntervalRG*xx

    return result
end

function certifySOSDecomposition(X, orderunit, λ::Number, Q::AbstractMatrix, support, RG, G)
    λInterval = @interval(λ)
    eoi = X - λInterval*orderunit

    @info "Element to be certified in interval arithmetic:"
    @info eoi

    residual = eoi - SOSFromMatrix(Q, support, RG, G)
    l1Norm = 0
    mn = size(X)[1]
    for i in 1:mn
        for j in 1:mn
            l1Norm += norm(residual[i,j],1)
        end
    end

    @info "l₁ norm in interval arithmetic:"
    @info l1Norm

    result = λInterval-l1Norm

    return result
end

function spectralGapsCertification(G, supportSize) # change support to halfbasis
    λₐₚ, Pₐₚ, RG, Δ₁, Iₙ = spectralGapsApproximated(G, supportSize)
    # Qₐₚ = real(sqrt(Symmetric( (P.+ P')./2 ))) # we do not need Cholesky decomposition - we can take the normal square to ensure that Q^2 will be >= using interval arithmetic
    PSymm = Symmetric((Pₐₚ.+ Pₐₚ')./2)
    minEigenvalue = minimum(eigvals(PSymm))
    if minEigenvalue < 0
        subtrahend = prevfloat(prevfloat(minEigenvalue))
        for i in 1:size(PSymm)[1]
            PSymm[i,i] -= subtrahend
        end
    end
    @info eigvals(PSymm)
    Ch = cholesky(PSymm)
    Qₐₚ = Ch.U

    @info "Approximated λ:"
    @info λₐₚ
    @info "Approximated P:"
    @info Pₐₚ
    @info "Approximated Q such that P = QᵀQ:"
    @info Qₐₚ
    
    result = certifySOSDecomposition(Δ₁, Iₙ, λₐₚ, Qₐₚ, RG.basis, RG, G)

    @info "Certified λ (interval atithmetic):"
    @info result

    return result
end
