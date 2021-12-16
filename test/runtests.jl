using StarAlgebras
using Groups
using Test
using LowCohomologySOS


@testset "cyclic_group" begin n = 3
   Cₙ = LowCohomologySOS.cyclic_group(n)

   a, = Groups.gens(Cₙ)

   @test a^n == one(Cₙ)
   @test all(a^i*a^j == a^(i+j) for i in 2:n for j in 2:n)
end

@testset "group_ring" begin n = 15
   Cₙ = LowCohomologySOS.cyclic_group(n)
   a, = Groups.gens(Cₙ)

   RCₙ = LowCohomologySOS.group_ring(Cₙ, n)
   @test RCₙ.object == Cₙ
   @test Set(collect(RCₙ.basis)) == Set([[a];[a^i for i in 2:n]])
   mstr = RCₙ.mstructure
   b = RCₙ.basis
   @test mstr[getindex(b, a),getindex(b, a)] == getindex(b, a^2)
   @test all(mstr[getindex(b, a),getindex(b, a^j)] == getindex(b, a^(j+1)) for j in 2:n)
   @test all(mstr[getindex(b, a^i),getindex(b, a^j)] == getindex(b, a^(i+j)) for i in 2:n for j in 2:n)

   RCₙ_star = LowCohomologySOS.group_ring(Cₙ, n, true)
   @test RCₙ_star.object == Cₙ
   @test Set(collect(RCₙ_star.basis)) == Set([[a];[a^i for i in 2:n]])
   mstr_star = RCₙ_star.mstructure
   b_star = RCₙ_star.basis
   @test mstr_star[getindex(b_star, a),getindex(b_star, a)] == getindex(b_star, one(Cₙ))
   @test all(mstr_star[getindex(b_star, a),getindex(b_star, a^j)] == getindex(b_star, a^(-1+j)) for j in 2:n if -1+j != 1)
   @test all(mstr_star[getindex(b_star, a^i),getindex(b_star, a^j)] == getindex(b_star, a^(-i+j)) for i in 2:n for j in 2:n if -i+j != 1 && i-j != 1)
end

@testset "fox_derivative" begin
   A = Alphabet([:x, :X, :y, :Y], [2, 1, 4, 3])
   F = FreeGroup(A)
   RF = LowCohomologySOS.group_ring(F, 2)

   @test LowCohomologySOS.fox_derivative(RF, one(F), 1) == zero(RF)
   @test LowCohomologySOS.fox_derivative(RF, one(F), 2) == zero(RF)

   x, y = Groups.gens(F)

   @test LowCohomologySOS.fox_derivative(RF, x, 1) == one(RF)
   @test LowCohomologySOS.fox_derivative(RF, x, 2) == zero(RF)
   @test LowCohomologySOS.fox_derivative(RF, y, 1) == zero(RF)
   @test LowCohomologySOS.fox_derivative(RF, y, 2) == one(RF)

   @test LowCohomologySOS.fox_derivative(RF, inv(x), 1) == -RF(inv(x))
   @test LowCohomologySOS.fox_derivative(RF, inv(y), 2) == -RF(inv(y))

   @test LowCohomologySOS.fox_derivative(RF, x*y, 1) == one(RF)
   @test LowCohomologySOS.fox_derivative(RF, x*y, 2) == RF(x)

   @test LowCohomologySOS.fox_derivative(RF, x*y*x*y, 1) == one(RF)+RF(x*y)
   @test LowCohomologySOS.fox_derivative(RF, x*y*x*y, 2) == RF(x)+RF(x*y*x)
end
