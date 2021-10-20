using Pkg
Pkg.activate(@__DIR__)

using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = 4
LinearAlgebra.BLAS.set_num_threads(2)

using JuMP
using SCS
using ProxSDP

include("starAlgebras.jl")

 
function constraints(pm::AbstractMatrix{<:Integer}, total_length=maximum(pm)) # this function has to be customized for matrices
   cnstrs = [Vector{Int}() for _ in 1:total_length]
   li = LinearIndices(CartesianIndices(size(pm)))
   for i in eachindex(pm)
      push!(cnstrs[pm[i]], li[i])
   end
   return cnstrs
end

function SOS_problem_primal(X::AlgebraElement, orderunit::AlgebraElement;
   upper_bound::Float64=Inf) # can stay

   @assert parent(X) === parent(orderunit) # can stay
   Al = parent(X) # can stay
   
   @info "Multiplicative structure", Al.mstructure

   N = size(Al.mstructure, 1) # ??
   m = JuMP.Model(); # can stay

   JuMP.@variable(m, P[1:N, 1:N]) # N has to be replaced with N*M; N = matrix size, M = supp size
   JuMP.@SDconstraint(m, sdp, P >= 0) # can stay

   if iszero(StarAlgebras.aug(X)) && iszero(StarAlgebras.aug(orderunit)) # to replace
      JuMP.@constraint(m, augmentation, sum(P) == 0) #to replace
   end
   if upper_bound < Inf # can stay
      λ = JuMP.@variable(m, λ <= upper_bound) # can stay
   else
      λ = JuMP.@variable(m, λ) # can stay
   end

   cnstrs = constraints(Al.mstructure) # can stay provided "contraints" function will be customized appropriately
   
   @info "Constraints", cnstrs

   @assert length(cnstrs) == length(X.coeffs) == length(orderunit.coeffs) # as above
   x, u = X.coeffs, orderunit.coeffs # as above

   @info "Coefficients of the group ring element to be certified", x
   @info "Coefficients of the order unit", u
   @info "Length of contraints", length(cnstrs)

   JuMP.@constraint(m, lincnstr[i=1:length(cnstrs)], x[i] - λ*u[i] == sum(P[cnstrs[i]])) # as above
   JuMP.@objective(m, Max, λ) # can stay

   @info m

   return m
end

function SOS(SOSProblem, groupRing) # has to be customized for the matrix case
   Q = real.(sqrt(value.(SOSProblem[:P])))
   Q = [round.(r;digits=1) for r in eachrow(Q)]
   @info Q
   @info groupRing.basis
end

cyclicGroupOptimizationProblem, RCₙ = let n = 3
   RCₙ, ID = cyclicGroupRing(n)
   S = collect(RCₙ.basis)
   a = S[2]
   X = 2*RCₙ(ID)-RCₙ(a)-RCₙ(inv(a))+n*sum(RCₙ(s) for s in S)
   
   @info "Group ring element to be certified:"
   @info X

   SOS_problem_primal(X, 1*RCₙ(ID)), RCₙ
end;

λ, cyclicGroupSolutionProxSDP = let SOS_problem = cyclicGroupOptimizationProblem
   with_ProxSDP = with_optimizer(ProxSDP.Optimizer, log_verbose=true, tol_gap=1e-4, tol_feasibility=1e-4)
   set_optimizer(SOS_problem, with_ProxSDP)
   optimize!(SOS_problem)
   λ, P_Cₙ = value(SOS_problem[:λ]), value.(SOS_problem[:P])
   Q = real.(sqrt(P_Cₙ))
   λ, Q, svdvals(Q), svdvals(P_Cₙ)
end;

λ, cyclicGroupSolutionSCS = let SOS_problem = cyclicGroupOptimizationProblem
   with_scs = with_optimizer(SCS.Optimizer, eps=1e-8)
   set_optimizer(SOS_problem, with_scs)
   optimize!(SOS_problem)
   λ, P_Cₙ = value(SOS_problem[:λ]), value.(SOS_problem[:P])
   Q = real.(sqrt(P_Cₙ))
   λ, Q, svdvals(Q), svdvals(P_Cₙ)
end

G1SOSOptimizationProblem, RG₁ = let # change X for matrix and define appropriately
   RG1, ID = G1GroupRing(2)
   S = let s = gens(RG1.object)
      unique!([s; inv.(s)])
   end
   a = S[1]
   b = S[2]
   c = S[3]
   A = S[4]
   B = S[5]
   C = S[6]
   X = 3*RG1(ID)+RG1(a)+RG1(A)+RG1(b)+RG1(B)+RG1(c)+RG1(C)
   SOS_problem_primal(X^2, X), RG1
end; # semicolon - the let block does not return anything - just computes the let block and there is no possibility to see the result unless the @info command is present somewhere

λ, G1SOSSolution = let SOS_problem = G1SOSOptimizationProblem # can stay
   with_scs = with_optimizer(SCS.Optimizer, eps=1e-8)
   set_optimizer(SOS_problem, with_scs)
   optimize!(SOS_problem)
   status = termination_status(SOS_problem)
   λ, P_G1 = value(SOS_problem[:λ]), value.(G1SOSOptimizationProblem[:P])
   @info status λ
   @info SOS(SOS_problem,RG₁)
   λ, P_G1
end;

SnSOSOptimizationProblem = let n = 4
   RSn = symmetricGroupRing(n)
   S = collect(RSn.object)
   Δ = RSn(length(S)) - sum(RSn(s) for s in S)
   SOS_problem_primal(Δ^2, Δ)
end;

λ, SnSOSSolution = let SOS_problem = SnSOSOptimizationProblem
   with_scs = with_optimizer(SCS.Optimizer, eps=1e-8)
   set_optimizer(SOS_problem, with_scs)
   optimize!(SOS_problem)
   status = termination_status(SOS_problem)
   λ, P_G1 = value(SOS_problem[:λ]), value.(SnSOSOptimizationProblem[:P])
   @info status λ
   λ, P_G1
end;

function foxDerivative(relatorWord, generator, underlyingGroupRing, underlyingGroup)
   relatorWordLength = length(relatorWord)
   result = underlyingGroupRing(0)
   multiplier = underlyingGroupRing(one(underlyingGroup))
   generators = gens(underlyingGroup)

   for i in 1:relatorWordLength
      if relatorWord[i] == generator.word[1] # generator
         result += multiplier
      elseif relatorWord[i] == (generator.word[1]+1) # inverse of a generator
         result -= (multiplier*underlyingGroupRing(inv(generator)))
      end

      if relatorWord[i]%2 == 1
         multiplier *= underlyingGroupRing(generators[floor(Int, (relatorWord[i]+1)/2)])
      else
         multiplier *= underlyingGroupRing(inv(generators[floor(Int, relatorWord[i]/2)]))
      end
   end

   return result
end

function jacobianMatrix(G, sufficientBallRadius)
   RG = groupRing(G, sufficientBallRadius)

   relations = G.relations
   generators = filter(x -> x != one(G), gens(G))
   relationsNumber = length(relations)
   generatorsNumber = length(generators)


   result = [RG(0) for i in collect(1:relationsNumber*generatorsNumber)]
   result = reshape(result, relationsNumber, generatorsNumber)

   for i in 1:relationsNumber
      for j in 1:generatorsNumber
         relator = relations[i].first*inv(relations[i].second)
         result[i,j] = foxDerivative(relator.word, generators[j], RG, G)
      end
   end

   return result
end

CₙJacobianMatrix = let n = 5
   Cₙ = cyclicGroup(n)
   jacobianMatrixx = jacobianMatrix(Cₙ, n)

   @info jacobianMatrixx
end

SL₃ƵShorterPresentation = let maxRules = 50000
   A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
   F = FreeGroup(A)
   x,y,z = Groups.gens(F)
   ε = one(F);
   G = FPGroup(F, [x^3 => ε, y^3 => ε, z^2 => ε, (x*z)^3 => ε, (y*z)^3 => ε, 
                   (x^(-1)*z*x*y)^2 => ε, (y^(-1)*z*y*x)^2 => ε, (x*y)^6 => ε], maxrules = maxRules )
   
   G
end

SL₃ƵJacobianShorterPresentation = let supportSize = 4
   jacobianMatrixx = jacobianMatrix(SL₃ƵShorterPresentation, supportSize)

   @info jacobianMatrixx
end

SL₃ƵElementaryMatrixPresentation = let maxRules = 50000
   A = Alphabet([:e12, :E12, :e21, :E21, :e13, :E13, :e31, :E31, :e23, :E23, :e32, :E32], [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11])
   F = FreeGroup(A)
   e12, e21, e13, e31, e23, e32 = Groups.gens(F)
   ε = one(F);
   G = FPGroup(F, [e12*e13 => e13*e12, e12*e32 => e32*e12, e13*e23 => e23*e13, e23*e21 => e21*e23, e21*e31 => e31*e21, e31*e32 => e32*e31,
                   e12*e23 => e13*e23*e12, e13*e32 => e12*e32*e13, e21*e13 => e23*e13*e21, e23*e31 => e21*e31*e23, 
                   e31*e12 => e32*e12*e31, e32*e21 => e31*e21*e32,
                   (e12*e21^(-1)*e12)^4 => ε], maxrules = maxRules )
   
   G
end

SL₃ƵJacobianElementaryMatrixPresentation = let supportSize = 4
   jacobianMatrixx = jacobianMatrix(SL₃ƵElementaryMatrixPresentation, supportSize)

   @info jacobianMatrixx
end