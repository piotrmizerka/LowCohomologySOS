using JuMP

include("starAlgebras.jl")

mutable struct GroupRingElem{T, A<:AbstractVector, GR<:StarAlgebra} <: NCRingElem
    coeffs::A
    parent::GR
 
    function GroupRingElem{T, A, GR}(c::AbstractVector{T}, RG::GR, check=true) where {T, A, GR}
       if check
          if hasbasis(RG)
             length(c) == length(RG.basis) || throw(
             "Can't create GroupRingElem -- lengths differ: length(c) =
             $(length(c)) != $(length(RG.basis)) = length(RG.basis)")
          else
             # @warn("Basis of the GroupRing is not defined.")
          end
       end
       return new{T, A, GR}(c, RG)
    end
 end

 mutable struct GroupRingElem2
    coeffs::AbstractVector
    parent::StarAlgebra
 
    function GroupRingElem2(c::AbstractVector, RG::StarAlgebra, check=true)
       if check
          if isdefined(RG,:basis)
             length(c) == length(RG.basis) || throw(
             "Can't create GroupRingElem -- lengths differ: length(c) =
             $(length(c)) != $(length(RG.basis)) = length(RG.basis)")
          else
             # @warn("Basis of the GroupRing is not defined.")
          end
       end
       return new(c, RG)
    end
 end

 aug(X::GroupRingElem2) = sum(X.coeffs)

 function constraints(pm::StarAlgebras.MTable{UInt32, true, Matrix{UInt32}}, total_length=maximum(pm)) where {I<:Integer}
   cnstrs = [Vector{CartesianIndex{2}}() for _ in 1:total_length]
   for i in eachindex(pm)
       push!(cnstrs[pm[i]], i)
   end
   return cnstrs
end

function SOS_problem_primal(X::GroupRingElem2, orderunit::GroupRingElem2;
    upper_bound::Float64=Inf)

   #  N = size(X.parent.pm, 1)
    N = size(X.parent.mstructure, 1)
    m = JuMP.Model();

    JuMP.@variable(m, P[1:N, 1:N])
    # SP = Symmetric(P)
    JuMP.@SDconstraint(m, sdp, P >= 0)

    if iszero(aug(X)) && iszero(aug(orderunit))
        JuMP.@constraint(m, augmentation, sum(P) == 0)
    end

    if upper_bound < Inf
        λ = JuMP.@variable(m, λ <= upper_bound)
    else
        λ = JuMP.@variable(m, λ)
    end

   #  cnstrs = constraints(X.parent.pm)
    cnstrs = constraints(X.parent.mstructure)
    @assert length(cnstrs) == length(X.coeffs) == length(orderunit.coeffs)

    x, u = X.coeffs, orderunit.coeffs
    JuMP.@constraint(m, lincnstr[i=1:length(cnstrs)],
        x[i] - λ*u[i] == sum(P[cnstrs[i]]))

    JuMP.@objective(m, Max, λ)

    return m
end

function SOS_G1()
   RG1, RADIUS, G1 = G1GroupRing()
   x = GroupRingElem2([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],RG1)
   orderunit = GroupRingElem2([0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],RG1)
   m = SOS_problem_primal(x,orderunit)
   println(m)
end

function SOS_S3()
   RS3 = S3GroupRing()
   x = GroupRingElem2([1,1,1,1,1,1],RS3)
   orderunit = GroupRingElem2([0,1,1,1,1,1],RS3)
   m = SOS_problem_primal(x,orderunit)
   println(m)
end

SOS_G1()
# SOS_S3() problems with types - need to be resolved yet