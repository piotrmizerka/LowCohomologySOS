using Pkg
Pkg.activate(@__DIR__)

using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = 4
LinearAlgebra.BLAS.set_num_threads(2)

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

   # @info "Optimization problem:"
   # @info result

   return result
end

# Show SOS summands - has to be customized for the matrix case
function SOS(SOSProblem, groupRing)
   Q = real.(sqrt(value.(SOSProblem[:P])))
   Q = [round.(r;digits=1) for r in eachrow(Q)]
   @info Q
   @info groupRing.basis
end


## cyclic group example ###################################################3
cyclicGroupOptimizationProblem = let n = 3
   Cₙ = cyclicGroup(n)
   ID = one(Cₙ)
   RCₙ = groupRing(Cₙ, n, true)
   S = collect(RCₙ.basis)

   a = S[2]
   X = 2*RCₙ(ID)-RCₙ(a)-RCₙ(inv(a))+n*sum(RCₙ(s) for s in S)
   
   M = [X RCₙ(0);RCₙ(0) X]

   # @info "Matrix to be certified:"
   # @info M
 
   orderUnit = [RCₙ(ID) RCₙ(0); RCₙ(0) RCₙ(ID)]

   # @info "Order unit matrix:"
   # @info orderUnit

   SOSProblemMatrix(M, orderUnit)
end

function SOSProblemSolutionSCS(SOSProblem)
   with_scs = with_optimizer(SCS.Optimizer, eps=1e-8)
   set_optimizer(SOSProblem, with_scs)
   optimize!(SOSProblem)
   λ, P = value(SOSProblem[:λ]), value.(SOSProblem[:P])
   # Q = real.(sqrt(P))

   @info "λ:"
   @info λ
   @info "P:"
   @info P
   # @info "Q s.t. Q^TQ = P:"
   # @info Q
end

# When running further commands after completing the one below, unexpected problems occur in VSCode
λ, cyclicGroupSolutionProxSDP = let SOS_problem = cyclicGroupOptimizationProblem
   with_ProxSDP = with_optimizer(ProxSDP.Optimizer, log_verbose=true, tol_gap=1e-4, tol_feasibility=1e-4)
   set_optimizer(SOS_problem, with_ProxSDP)
   optimize!(SOS_problem)
   λ, P_Cₙ = value(SOS_problem[:λ]), value.(SOS_problem[:P])
   Q = real.(sqrt(P_Cₙ))
   λ, Q, svdvals(Q), svdvals(P_Cₙ)
end

λ, cyclicGroupSolutionSCS = let SOS_problem = cyclicGroupOptimizationProblem
   with_scs = with_optimizer(SCS.Optimizer, eps=1e-8)
   set_optimizer(SOS_problem, with_scs)
   optimize!(SOS_problem)
   λ, P_Cₙ = value(SOS_problem[:λ]), value.(SOS_problem[:P])
   Q = real.(sqrt(P_Cₙ))
   λ, Q, svdvals(Q), svdvals(P_Cₙ)
end
##############################################

## another example of the definition of an SDP matrix problem #################3
let
   A = Alphabet([:x, :X, :y, :Y], [2, 1, 4, 3])
   F = FreeGroup(A)
   x,y = Groups.gens(F)
   ε = one(F)
   # G = FPGroup(F, [x^2 => ε, y^2 => ε, x*y => y*x] )
   G = FPGroup(F, [x^2 => ε, y^2 => ε] )

   RG = groupRing(G, 1, true)
   S = collect(RG.basis)

   @info "Basis:"
   @info S

   M = [RG(one(G)) RG(one(G)); RG(one(G)) RG(one(G))]

   @info "Matrix to be certified:"
   @info M

   SOSProblemMatrix(M, M)
end
##############################################################################

function spectralGapsApproximated(G, supportSize)
   generators = gens(G)
   jacobianMatrixEncodedx = jacobianMatrixEncoded(G.relations, generators)
   RGDifferentials = suitableGroupRing(G, generators, jacobianMatrixEncodedx)

   Δ₁⁺x = Δ₁⁺(G, jacobianMatrixEncodedx, generators, RGDifferentials)
   Δ₁⁻x = Δ₁⁻(G, generators, RGDifferentials)
   Δ₁x = Δ₁(G, jacobianMatrixEncodedx, generators, RGDifferentials)

   RGBallStar = groupRing(G, supportSize, true)

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

   # @info "Δ₁⁺ in both rings - visibly should be the same:"
   # printMatrix(Δ₁⁺x)
   # printMatrix(Δ₁⁺xx)
   # @info "Δ₁⁻ in both rings - visibly should be the same:"
   # printMatrix(Δ₁⁻x)
   # printMatrix(Δ₁⁻xx)
   # @info "Δ₁ in both rings - visibly should be the same:"
   # printMatrix(Δ₁x)
   # printMatrix(Δ₁xx)

   Iₙ = [RGBallStar(0) for i in 1:length(Δ₁xx)]
   Iₙ = reshape(Iₙ, size(Δ₁xx)[1], size(Δ₁xx)[2])
   for i in 1:size(Δ₁xx)[1]
      Iₙ[i,i] = RGBallStar(one(G))
   end

   Δ₁⁺SOSProblem = SOSProblemMatrix(Δ₁⁺xx^2, Δ₁⁺xx) # CAUTION: may require potentially twice the basis as Δ₁
   Δ₁⁻SOSProblem = SOSProblemMatrix(Δ₁⁻xx^2, Δ₁⁻xx) # as above
   Δ₁SOSProblem = SOSProblemMatrix(Δ₁xx, Iₙ)

   @info "Solution for (Δ₁⁺)²-λΔ₁⁺ = SOS:"
   SOSProblemSolutionSCS(Δ₁⁺SOSProblem) # CAUTION: may require potentially twice the basis as Δ₁
   @info "Solution for (Δ₁⁻)²-λΔ₁⁻ = SOS:"
   SOSProblemSolutionSCS(Δ₁⁻SOSProblem) # CAUTION: may require potentially twice the basis as Δ₁
   @info "Solution for Δ₁-λIₙ = SOS:"
   SOSProblemSolutionSCS(Δ₁SOSProblem)
end

let n = 3
   Cₙ = cyclicGroup(n)
   spectralGapsApproximated(Cₙ, n)
end
 
SL₃ƵShorterPresentationSpectralGaps = let maxRules = 1000, supportSize = 5
   A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
   F = FreeGroup(A)
   x,y,z = Groups.gens(F)
   ε = one(F);
   SL₃ƵShorterPresentation = FPGroup(F, [x^3 => ε, y^3 => ε, z^2 => ε, (x*z)^3 => ε, (y*z)^3 => ε, 
                   (x^(-1)*z*x*y)^2 => ε, (y^(-1)*z*y*x)^2 => ε, (x*y)^6 => ε], maxrules = maxRules )
   
   spectralGapsApproximated(SL₃ƵShorterPresentation, supportSize)
end

SL₃ƵElementaryMatrixPresentationSpectralGaps = let maxRules = 1000, supportSize = 2
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