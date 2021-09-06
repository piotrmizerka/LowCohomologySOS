using Pkg
Pkg.activate(@__DIR__)

using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = 4
LinearAlgebra.BLAS.set_num_threads(2)

using JuMP
using SCS

include("starAlgebras.jl")

 
function constraints(pm::AbstractMatrix{<:Integer}, total_length=maximum(pm))
   cnstrs = [Vector{Int}() for _ in 1:total_length]
   li = LinearIndices(CartesianIndices(size(pm)))
   for i in eachindex(pm)
      push!(cnstrs[pm[i]], li[i])
   end
   return cnstrs
end;

function SOS_problem_primal(X::AlgebraElement, orderunit::AlgebraElement;
   upper_bound::Float64=Inf)

   @assert parent(X) === parent(orderunit)
   Al = parent(X)

   N = size(Al.mstructure, 1)
   m = JuMP.Model();

   JuMP.@variable(m, P[1:N, 1:N])
   JuMP.@SDconstraint(m, sdp, P >= 0)

   if iszero(StarAlgebras.aug(X)) && iszero(StarAlgebras.aug(orderunit))
      JuMP.@constraint(m, augmentation, sum(P) == 0)
   end
   if upper_bound < Inf
      λ = JuMP.@variable(m, λ <= upper_bound)
   else
      λ = JuMP.@variable(m, λ)
   end

   cnstrs = constraints(Al.mstructure)
   @assert length(cnstrs) == length(X.coeffs) == length(orderunit.coeffs)
   x, u = X.coeffs, orderunit.coeffs
   JuMP.@constraint(m, lincnstr[i=1:length(cnstrs)], x[i] - λ*u[i] == sum(P[cnstrs[i]]))
   JuMP.@objective(m, Max, λ)

   return m
end;

G1SOSOptimizationProblem = let 
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
   SOS_problem_primal(X^2, X)
end;

λ, G1SOSSolution = let SOS_problem = G1SOSOptimizationProblem
   with_scs = with_optimizer(SCS.Optimizer, eps=1e-8)
   set_optimizer(SOS_problem, with_scs)
   optimize!(SOS_problem)
   status = termination_status(SOS_problem)
   λ, P_G1 = value(SOS_problem[:λ]), value.(G1SOSOptimizationProblem[:P])
   @info status λ
   λ, P_G1
end;

S3SOSOptimizationProblem = let n = 4
   RSn = symmetricGroupRing(n)
   S = collect(RSn.object)
   Δ = RSn(length(S)) - sum(RSn(s) for s in S)
   SOS_problem_primal(Δ^2, Δ)
end;

λ, S3SOSSolution = let SOS_problem = S3SOSOptimizationProblem
   with_scs = with_optimizer(SCS.Optimizer, eps=1e-8)
   set_optimizer(SOS_problem, with_scs)
   optimize!(SOS_problem)
   status = termination_status(SOS_problem)
   λ, P_G1 = value(SOS_problem[:λ]), value.(G1SOSOptimizationProblem[:P])
   @info status λ
   λ, P_G1
end;
