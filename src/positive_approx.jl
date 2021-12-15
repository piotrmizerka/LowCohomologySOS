# The function below concerns definition of constraints defined for the semi-positive definite matrix P to be computed.
function constraints(pm::AbstractMatrix{<:Integer}, total_length=maximum(pm))
   cnstrs = [Vector{Int}() for _ in 1:total_length]
   li = LinearIndices(CartesianIndices(size(pm)))

   for i in eachindex(pm)
      push!(cnstrs[pm[i]], li[i])
   end

   return cnstrs
end

# As the constraints function, this function also concerns defining the constraints arising from matrix P.
# More precisely, it defines constraints arising from the (row_id,column_id)-entry of the matrix equation for our problem.
# Me can apply the constraints function written by M. Kaluba to define the constraints arising from each entry.
# Order of linear indices for matrices which has to be applied: column snake (from Seattle to Miami).
function entry_constraint(cnstrs, row_id, column_id, constraind_id, relations_number, generators_number)
   B = (column_id-1)*relations_number^2*generators_number+(row_id-1)*relations_number
   result = copy(cnstrs[constraind_id])
   for l in 1:length(cnstrs[constraind_id])
      summand = (cnstrs[constraind_id][l]%relations_number != 0) ? cnstrs[constraind_id][l]%relations_number : relations_number
      factor = cnstrs[constraind_id][l]-summand
      result[l] = B+factor*generators_number+summand
   end

   return result
end

function sos_problem_matrix(M, order_unit, upper_bound::Float64=Inf)
   underlying_group_ring = parent(rand(M))
   m = size(underlying_group_ring.mstructure, 1)
   n = size(M,1)
   mn = m*n
   result = JuMP.Model();

   JuMP.@variable(result, P[1:mn, 1:mn])
   JuMP.@SDconstraint(result, sdp, P >= 0)

   if upper_bound < Inf
      λ = JuMP.@variable(result, λ <= upper_bound)
   else
      λ = JuMP.@variable(result, λ)
   end

   cnstrs = constraints(underlying_group_ring.mstructure)
   
   for i in 1:n
      for j in 1:n
         mij, u = M[i,j].coeffs, order_unit[i,j].coeffs
         @assert length(cnstrs) == length(mij) == length(u)
         JuMP.@constraint(result, [k=1:length(cnstrs)], mij[k] - λ*u[k] == sum(P[entry_constraint(cnstrs, i, j, k, m, n)]))
      end
   end

   JuMP.@objective(result, Max, λ)

   return result
end

function sos_problem_solution_scs(sos_problem)
   with_scs = with_optimizer(SCS.Optimizer, eps=1e-5, acceleration_lookback=0, max_iters = 5000000)
   set_optimizer(sos_problem, with_scs)
   optimize!(sos_problem)
   λ, P = value(sos_problem[:λ]), value.(sos_problem[:P])

   return λ, P
end

# h:Free group --> our group G
function spectral_gaps_approximated(h::Function, relations, half_basis)
   F = parent(rand(relations))
   G = parent(h(rand(relations)))

   d₁ = jacobian_matrix(relations)
   d₀x = d₀(parent(rand(d₁)), Groups.gens(F))
   Δ₁⁺ = d₁'*d₁
   Δ₁⁻ = d₀x*d₀x'
   Δ₁ = Δ₁⁺+Δ₁⁻

   RG_ball_star = group_ring(G, half_basis, true)

   Δ₁x = embed_to_group_ring.(Δ₁, Ref(RG_ball_star), h)

   n = length(Groups.gens(F))
   @assert size(Δ₁x,1) === size(Δ₁x,2) === n
   Iₙ = reshape([zero(RG_ball_star) for i in 1:(n*n)], n, n)
   for i in 1:n
      Iₙ[i,i] = one(RG_ball_star)
   end

   Δ₁_sos_problem = sos_problem_matrix(Δ₁x, Iₙ)
   λ, P = sos_problem_solution_scs(Δ₁_sos_problem)

   return λ, P, RG_ball_star, Δ₁x, Iₙ
end
 