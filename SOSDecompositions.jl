# include("starAlgebras.jl")
include("FoxDerivatives.jl")
 
function constraints(pm::AbstractMatrix{<:Integer}, total_length=maximum(pm))
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
   with_scs = with_optimizer(SCS.Optimizer, eps=1e-5, acceleration_lookback=0, max_iters = 5000000)
   set_optimizer(SOSProblem, with_scs)
   optimize!(SOSProblem)
   λ, P = value(SOSProblem[:λ]), value.(SOSProblem[:P])

   return λ, P
end

# h:Free group --> our group G
function spectralGapsApproximated(h::Function, relations, halfBasis)
   F = parent(relations[1])
   G = parent(h(relations[1]))

   jacobianFreeGroup = jacobianMatrix(relations)
   RGDifferentials = suitableGroupRing(jacobianFreeGroup, h)

   D₀x = D₀(G, RGDifferentials, [h(x) for x in Groups.gens(F)])
   D₁ = jacobianMatrix(jacobianFreeGroup, h, RGDifferentials)

   Δ₁⁺ = starOfMatrixOverGroupRing(D₁)*D₁
   Δ₁⁻ = D₀x*starOfMatrixOverGroupRing(D₀x)
   Δ₁ = Δ₁⁺+Δ₁⁻

   RGBallStar = groupRing(G, halfBasis, true)

   Δ₁⁺x = changeUnderlyingGroupRing(Δ₁⁺, RGDifferentials, RGBallStar, G)
   Δ₁⁻x = changeUnderlyingGroupRing(Δ₁⁻, RGDifferentials, RGBallStar, G)
   Δ₁x = changeUnderlyingGroupRing(Δ₁, RGDifferentials, RGBallStar, G)

   Iₙ = [RGBallStar(0) for i in 1:length(Δ₁⁻x)]
   Iₙ = reshape(Iₙ, size(Δ₁⁻x)[1], size(Δ₁⁻x)[2])
   for i in 1:size(Δ₁⁻x)[1]
      Iₙ[i,i] = RGBallStar(one(G))
   end

   Δ₁SOSProblem = SOSProblemMatrix(Δ₁x, Iₙ)
   λ, P = SOSProblemSolutionSCS(Δ₁SOSProblem)

   result = λ, P, RGBallStar, Δ₁x, Iₙ

   return result
end
 

