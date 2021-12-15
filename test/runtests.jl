using StarAlgebras
using Groups
using Test
using LowCohomologySOS

@testset "Fox derivatives" begin
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
