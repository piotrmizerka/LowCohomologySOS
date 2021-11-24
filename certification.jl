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
function SOSFromMatrix(Q, support, RG)
    mn = size(Q)[1]

    # Changing Q to the corresponding interval-entry matrix
    QInterval = reshape([@interval(0) for i in 1:(mn*mn)], mn, mn)
    for i in 1:mn
        for j in i:mn
            QInterval[i,j] = @interval(Q[i,j])
            QInterval[j,i] = QInterval[i,j]
        end
    end

    PInterval = QInterval^2
    PIntervalRG = reshape([@interval(0)*RG(0) for i in 1:(mn*mn)], mn, mn)
    for i in 1:mn
        for j in 1:mn
            PIntervalRG[i,j] = PInterval[i,j]*one(RG)
        end
    end

    m = length(support)
    n = floor(Int, mn/m)

    Iₙ = reshape([RG(0) for i in 1:(n*n)], n, n)
    for i in 1:n
        Iₙ[i,i] = one(RG)
    end

    x = reshape([RG(support[i]) for i in 1:length(support)], m, 1)
    xx = collect(Iₙ⊗x)
    xxᵀ = starOfMatrixOverGroupRing(xx)

    result = xxᵀ*PIntervalRG*xx

    return result
end

function certifySOSDecomposition(X, orderunit, λ::Number, Q::AbstractMatrix, support, RG)
    λInterval = @interval(λ)
    eoi = X - λInterval*orderunit

    residual = eoi - SOSFromMatrix(Q, support, RG)
    l1Norm = 0
    mn = size(X)[1]
    for i in 1:mn
        for j in 1:mn
            l1Norm += norm(residual[i,j],1)
        end
    end

    @info "l₁ norm of the error in interval arithmetic:"
    @info l1Norm

    result = λInterval-l1Norm

    return result
end

function spectralGapsCertification(h::Function, relations, halfBasis) # change support to halfbasis
    λₐₚ, Pₐₚ, RG, Δ₁, Iₙ = spectralGapsApproximated(h, relations, halfBasis)
    Qₐₚ = real(sqrt(Symmetric( (Pₐₚ.+ Pₐₚ')./2 )))

    @info "Approximated λ:"
    @info λₐₚ
    
    result = certifySOSDecomposition(Δ₁, Iₙ, λₐₚ, Qₐₚ, halfBasis, RG)

    @info "Certified λ (interval atithmetic):"
    @info result

    return result
end
