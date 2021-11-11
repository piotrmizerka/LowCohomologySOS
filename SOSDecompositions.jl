using Pkg
Pkg.activate(@__DIR__)

using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = 4
LinearAlgebra.BLAS.set_num_threads(2)
using Kronecker

using JuMP
using SCS
using ProxSDP

# include("starAlgebras.jl")
include("FoxDerivatives.jl")
 
function constraints(pm::AbstractMatrix{<:Integer}, total_length=maximum(pm)) # this function has to be customized for matrices
   cnstrs = [Vector{Int}() for _ in 1:total_length]
   li = LinearIndices(CartesianIndices(size(pm)))
   for i in eachindex(pm)
      push!(cnstrs[pm[i]], li[i])
   end
   return cnstrs
end

# order: column snake (from Szczecin to Rzeszow)
function entryConstraint(cnstrs, i, j, k, m, n)
   B = (j-1)*m^2*n+(i-1)*m
   result = copy(cnstrs[k])
   for l in 1:length(cnstrs[k])
      summand = cnstrs[k][l]%m
      if summand == 0
         summand = m
      end
      factor = cnstrs[k][l]-summand
      result[l] = B+factor*n+summand
   end

   return result
end

function SOSProblemMatrix(M, orderUnit, upper_bound::Float64=Inf)
   underlyingGroupRing = parent(M[1,1])
   m = size(underlyingGroupRing.mstructure, 1)
   n = size(M)[1]
   mn = m*n
   result = JuMP.Model();

   JuMP.@variable(result, P[1:mn, 1:mn])
   JuMP.@SDconstraint(result, sdp, P >= 0)

   if upper_bound < Inf
      λ = JuMP.@variable(result, λ <= upper_bound)
   else
      λ = JuMP.@variable(result, λ)
   end

   cnstrs = constraints(underlyingGroupRing.mstructure)
   
   for i in 1:n
      for j in 1:n
         mij, u = M[i,j].coeffs, orderUnit[i,j].coeffs
         # @assert length(cnstrs) == length(mij) == length(u)
         JuMP.@constraint(result, [k=1:length(cnstrs)], mij[k] - λ*u[k] == sum(P[entryConstraint(cnstrs, i, j, k, m, n)]))
      end
   end

   JuMP.@objective(result, Max, λ)

   return result
end

# Show SOS summands - has to be customized for the matrix case
function SOSSummands(SOSProblem, groupRing)
   Q = real.(sqrt(value.(SOSProblem[:P])))
   Q = [round.(r;digits=1) for r in eachrow(Q)]
   @info Q
   @info groupRing.basis
end

Base.adjoint(X::AlgebraElement) = StarAlgebra.star(X)
function SOSFromMatrix(P, support, RG, G)
   # @info "Support:"
   # @info support

   mn = size(P)[1]
   m = length(support)
   n = floor(Int, mn/m)

   Iₙ = reshape([RG(0) for i in 1:(n*n)], n, n)
   for i in 1:n
      Iₙ[i,i] = RG(one(G))
   end

   # @info "Iₙ:"
   # @info Iₙ

   PRG = reshape([RG(0) for i in 1:(mn*mn)], mn, mn)
   for i in 1:mn
      for j in 1:mn
         PRG[i,j] = P[i,j]*RG(one(G))
      end
   end

   # @info "PRG:"
   # @info PRG
   # @info PRG*PRG

   x = reshape([RG(support[i]) for i in 1:length(support)], m, 1)
   xx = collect(Iₙ⊗x)
   
   # @info "xx:"
   # @info xx

   # xxᵀ = xx' # not working!
   xxᵀ = starOfMatrixOverGroupRing(xx)

   @info "xxᵀ:"
   @info xxᵀ

   result = xxᵀ*PRG*xx
   # result = xxᵀ*xx

   @info "result:"
   @info result

   return result
end

let n = 3
   Cₙ = cyclicGroup(n)
   RCₙ = groupRing(Cₙ, n)
   P = [1 1 1 1 1 1;1 1 1 1 1 1;1 1 1 1 1 1;1 1 1 1 1 1;1 1 1 1 1 1;1 1 1 1 1 1]
   M = SOSFromMatrix(P, RCₙ.basis, RCₙ, Cₙ) # M shall be equal to []
end

function SOSProblemSolutionSCS(SOSProblem)
   with_scs = with_optimizer(SCS.Optimizer, eps=1e-2, acceleration_lookback=0)
   set_optimizer(SOSProblem, with_scs)
   optimize!(SOSProblem)
   λ, P = value(SOSProblem[:λ]), value.(SOSProblem[:P])
   # Q = real.(sqrt(P))

   @info "λ:"
   @info λ
   @info "Size of P:"
   @info size(P)
   # @info "P:"
   # @info P
   # @info "Q s.t. Q^TQ = P:"
   # @info Q
end

function spectralGapsApproximated(G, supportSize)
   generators = gens(G)
   jacobianMatrixEncodedx = jacobianMatrixEncoded(G.relations, generators)
   RGDifferentials = suitableGroupRing(G, generators, jacobianMatrixEncodedx)

   Δ₁⁺x = Δ₁⁺(G, jacobianMatrixEncodedx, generators, RGDifferentials)
   Δ₁⁻x = Δ₁⁻(G, generators, RGDifferentials)
   Δ₁x = Δ₁(G, jacobianMatrixEncodedx, generators, RGDifferentials)

   RGBallStar = groupRing(G, supportSize, true)
   # RGBallStar = RGDifferentials # we may try to find a solution with a prescibed support , e.g. the same as for computing the differentials - advantage: goes fast, con: may not find a solution

   # @info "Basis of RGDifferentials:"
   # @info RGDifferentials.basis
   # @info "Basis of RGBallStar:"
   # @info RGBallStar.basis
   @info "Size of basis of RGDifferentials and dimension of its multiplication table:"
   @info [length(RGDifferentials.basis) size(RGDifferentials.mstructure)[1]]
   @info "Size of basis of RGBallStar and dimension of its multiplication table:"
   @info [length(RGBallStar.basis) size(RGBallStar.mstructure)[1]]

   Δ₁⁺xx = changeUnderlyingGroupRing(Δ₁⁺x, RGDifferentials, RGBallStar, G)
   Δ₁⁻xx = changeUnderlyingGroupRing(Δ₁⁻x, RGDifferentials, RGBallStar, G)
   Δ₁xx = changeUnderlyingGroupRing(Δ₁x, RGDifferentials, RGBallStar, G)

   Iₙ = [RGBallStar(0) for i in 1:length(Δ₁⁻xx)]
   Iₙ = reshape(Iₙ, size(Δ₁⁻xx)[1], size(Δ₁⁻xx)[2])
   for i in 1:size(Δ₁⁻xx)[1]
      Iₙ[i,i] = RGBallStar(one(G))
   end

   # Δ₁⁺SOSProblem = SOSProblemMatrix(Δ₁⁺xx^2, Δ₁⁺xx) # CAUTION: may require potentially twice the basis as Δ₁
   # Δ₁⁻SOSProblem = SOSProblemMatrix(Δ₁⁻xx^2, Δ₁⁻xx) # as above
   Δ₁SOSProblem = SOSProblemMatrix(Δ₁xx, Iₙ)

   # @info "Solution for (Δ₁⁺)²-λΔ₁⁺ = SOS:"
   # SOSProblemSolutionSCS(Δ₁⁺SOSProblem) # CAUTION: may require potentially twice the basis as Δ₁
   # @info "Solution for (Δ₁⁻)²-λΔ₁⁻ = SOS:"
   # SOSProblemSolutionSCS(Δ₁⁻SOSProblem) # CAUTION: may require potentially twice the basis as Δ₁
   @info "Solution for Δ₁-λIₙ = SOS:"
   SOSProblemSolutionSCS(Δ₁SOSProblem)
end
 
SL₃ƵShorterPresentationSpectralGaps = let maxRules = 1000, supportSize = 3
   A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
   F = FreeGroup(A)
   x,y,z = Groups.gens(F)
   ε = one(F);
   SL₃ƵShorterPresentation = FPGroup(F, [x^3 => ε, y^3 => ε, z^2 => ε, (x*z)^3 => ε, (y*z)^3 => ε, 
                   (x^(-1)*z*x*y)^2 => ε, (y^(-1)*z*y*x)^2 => ε, (x*y)^6 => ε], maxrules = maxRules )
   
   # spectralGapsApproximated(SL₃ƵShorterPresentation, supportSize)
   differentials(SL₃ƵShorterPresentation)
end

SL₃ƵElementaryMatrixPresentationSpectralGaps = let maxRules = 5000, supportSize = 2
   A = Alphabet([:e12, :E12, :e21, :E21, :e13, :E13, :e31, :E31, :e23, :E23, :e32, :E32], [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11])
   F = FreeGroup(A)
   e12, e21, e13, e31, e23, e32 = Groups.gens(F)
   ε = one(F);
   SL₃ƵElementaryMatrixPresentation = FPGroup(F, [e12*e13 => e13*e12, e12*e32 => e32*e12, e13*e23 => e23*e13, e23*e21 => e21*e23, e21*e31 => e31*e21, e31*e32 => e32*e31,
                   e12*e23 => e13*e23*e12, e13*e32 => e12*e32*e13, e21*e13 => e23*e13*e21, e23*e31 => e21*e31*e23, 
                   e31*e12 => e32*e12*e31, e32*e21 => e31*e21*e32,
                   (e12*e21^(-1)*e12)^4 => ε], maxrules = maxRules )

   spectralGapsApproximated(SL₃ƵElementaryMatrixPresentation, supportSize)
end

SL₃ƵSteinbergGroupSpectralGaps = let maxRules = 5000, supportSize = 2
   A = Alphabet([:e12, :E12, :e21, :E21, :e13, :E13, :e31, :E31, :e23, :E23, :e32, :E32], [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11])
   F = FreeGroup(A)
   e12, e21, e13, e31, e23, e32 = Groups.gens(F)
   ε = one(F);
   SteinbergGroupPresentation = FPGroup(F, [e12*e13 => e13*e12, e12*e32 => e32*e12, e13*e23 => e23*e13, e23*e21 => e21*e23, e21*e31 => e31*e21, e31*e32 => e32*e31,
                   e12*e23 => e13*e23*e12, e13*e32 => e12*e32*e13, e21*e13 => e23*e13*e21, e23*e31 => e21*e31*e23, 
                   e31*e12 => e32*e12*e31, e32*e21 => e31*e21*e32], maxrules = maxRules )

   # differentials(SteinbergGroupPresentation)
   spectralGapsApproximated(SteinbergGroupPresentation, supportSize)
end