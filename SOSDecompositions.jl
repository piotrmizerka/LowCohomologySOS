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

   N = size(Al.mstructure, 1) # mstructure has to be replaced with something else for matrices
   m = JuMP.Model(); # can stay

   JuMP.@variable(m, P[1:N, 1:N]) # can stay
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
   @assert length(cnstrs) == length(X.coeffs) == length(orderunit.coeffs) # as above
   x, u = X.coeffs, orderunit.coeffs # as above
   JuMP.@constraint(m, lincnstr[i=1:length(cnstrs)], x[i] - λ*u[i] == sum(P[cnstrs[i]])) # as above
   JuMP.@objective(m, Max, λ) # can stay

   return m
end

function SOS(SOSProblem, groupRing) # has to be customized for the matrixc case
   Q = real.(sqrt(value.(SOSProblem[:P])))
   @info size(Q)
   Q = [round.(sort(r);digits=2) for r in eachrow(Q)]
   # Q = [round.(r;digits=1) for r in eachrow(Q)]
   @info Q
   @info groupRing.basis
   # return sum([AlgebraElement(collect(c),groupRing)^2 for c in eachrow(Q)])
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

cyclicGroupOptimizationProblem, RCₙ = let n = 3
   RCₙ, ID = cyclicGroupRing(n)
   S = collect(RCₙ.basis)
   a = S[2]
   X = 2*RCₙ(ID)-RCₙ(a)-RCₙ(inv(a))+n*sum(RCₙ(s) for s in S)
   # X = 2*RCₙ(ID)-RCₙ(a)-RCₙ(inv(a))+0.01*sum(RCₙ(s) for s in S)
   
   @info X

   SOS_problem_primal(X, 1*RCₙ(ID)), RCₙ
end;

λ, cyclicGroupSolution = let SOS_problem = cyclicGroupOptimizationProblem
   # with_scs = with_optimizer(SCS.Optimizer, eps=1e-8)
   with_ProxSDP = with_optimizer(ProxSDP.Optimizer, log_verbose=true, tol_gap=1e-4, tol_feasibility=1e-4)
   # set_optimizer(SOS_problem, with_scs)
   set_optimizer(SOS_problem, with_ProxSDP)
   optimize!(SOS_problem)
   # status = termination_status(SOS_problem)
   λ, P_Cₙ = value(SOS_problem[:λ]), value.(SOS_problem[:P])
   # @info status λ
   # @info SOS(SOS_problem,RCₙ)
   # Q = real.(sqrt(P_Cₙ))
   # λ, Q, svdvals(Q)
end
