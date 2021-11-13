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

function SOSProblemSolutionSCS(SOSProblem)
   with_scs = with_optimizer(SCS.Optimizer, eps=1e-8, acceleration_lookback=0)
   set_optimizer(SOSProblem, with_scs)
   optimize!(SOSProblem)
   λ, P = value(SOSProblem[:λ]), value.(SOSProblem[:P])
   # Q = real.(sqrt(P))

   # @info "λ:"
   # @info λ
   # @info "Size of P:"
   # @info size(P)
   # @info "P:"
   # @info P
   # @info "Q s.t. Q^TQ = P:"
   # @info Q

   return λ, P
end

function spectralGapsApproximated(G, supportSize)
   generators = gens(G)
   jacobianMatrixEncodedx = jacobianMatrixEncoded(G.relations, generators)
   RGDifferentials = suitableGroupRing(G, generators, jacobianMatrixEncodedx)

   Δ₁⁺x = Δ₁⁺(G, jacobianMatrixEncodedx, generators, RGDifferentials)
   Δ₁⁻x = Δ₁⁻(G, generators, RGDifferentials)
   Δ₁x = Δ₁(G, jacobianMatrixEncodedx, generators, RGDifferentials)

   @info "Δ₁:"
   printMatrix(Δ₁x)

   RGBallStar = groupRing(G, supportSize, true)
   # RGBallStar = RGDifferentials # we may try to find a solution with a prescibed support , e.g. the same as for computing the differentials - advantage: goes fast, con: may not find a solution

   # @info "Basis of RGDifferentials:"
   # @info RGDifferentials.basis
   # @info "Basis of RGBallStar:"
   # @info RGBallStar.basis
   # @info "Size of basis of RGDifferentials and dimension of its multiplication table:"
   # @info [length(RGDifferentials.basis) size(RGDifferentials.mstructure)[1]]
   # @info "Size of basis of RGBallStar and dimension of its multiplication table:"
   # @info [length(RGBallStar.basis) size(RGBallStar.mstructure)[1]]

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
   # @info "Solution for Δ₁-λIₙ = SOS:"
   # SOSProblemSolutionSCS(Δ₁SOSProblem)

   λ, P = SOSProblemSolutionSCS(Δ₁SOSProblem)

   result = λ, P, RGBallStar, Δ₁xx, Iₙ

   return result
end
 

