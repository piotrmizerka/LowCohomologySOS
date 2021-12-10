using StarAlgebras
using AbstractAlgebra
using Groups
using Test

include("../src/SOSStarAlgebras.jl")

function differentialsExamples()
    A = Alphabet([:x, :X, :y, :Y], [2, 1, 4, 3])
    F = FreeGroup(A)
    x,y = Groups.gens(F)
    ε = one(F)
    # G = FPGroup(F, [x*y => ε] ) # proper matrix: [1 x]
    # G = FPGroup(F, [x*y => ε, y*x => ε] ) # proper matrix: [1 x; X 1]
    # G = FPGroup(F, [x*y*x^(-1)*y^(-1) => ε] ) # proper matrix: [1-y x-1]
    G = FPGroup(F, [x*y => y*x] ) # proper matrix: as above
 
    # A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
    # F = FreeGroup(A)
    # x,y,z = Groups.gens(F)
    # ε = one(F)
    # G = FPGroup(F, [x*y*z => y*x] ) # proper matrix: [1-y x-1 xy]
    # G = FPGroup(F, [x*y*z => y*x, z => ε] ) # proper matrix: [1-y x-1 xy; 0 0 1]
 
    # TODO: do the examples above and add cyclic group
end

 
## Defining optimization problems ###################################################3
cyclicGroupOptimizationProblem = let n = 3
    Cₙ = cyclic_group(n)
    ID = one(Cₙ)
    RCₙ, halfBasis = group_ring(Cₙ, n, true)
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
 
 ## Solving optimization problems (in particular establishing approx. spectral gaps) #################3
let
   A = Alphabet([:x, :X, :y, :Y], [2, 1, 4, 3])
   F = FreeGroup(A)
   x,y = Groups.gens(F)
   ε = one(F)
   # G = FPGroup(F, [x^2 => ε, y^2 => ε, x*y => y*x] )
   # G = FPGroup(F, [x^2 => ε, y^2 => ε] )
   # G = FPGroup(F, [x*y => y*x] ) # G = Z²
   G = FPGroup(F, [y => ε] ) # G = Z

   @info "Generators of G:"
   @info gens(G)

   RG = group_ring(G, 2, true)
   S = collect(StarAlgebras.basis(RG))

   @info "Basis:"
   @info S

   # M = [RG(one(G)) RG(one(G)); RG(one(G)) RG(one(G))]
   # x, y = gens(G)
   # Δ = 4*RG(one(G))-RG(x)-RG(inv(x))-RG(y)-RG(inv(y))
   x = gens(G)[1]
   Δ = 1*RG(one(G))-RG(x)-RG(inv(x))
   M = Δ^2

   @info "Matrix to be certified:"
   @info M

   # SOSProblemMatrixx = SOSProblemMatrix([M], [Δ])
   SOSProblemMatrixx = SOSProblemMatrix([Δ], [RG(one(G))])
   SOSProblemSolutionSCS(SOSProblemMatrixx)
end
##############################################################################

###### test examples - certification ###################################
let n = 3
   Cₙ = cyclic_group(n)
   # TODO spectralGapsCertification(Cₙ, n)
end

SL₃ƵSpectralGaps = let halfRadius = 2
   A = Alphabet([:e12, :E12, :e21, :E21, :e13, :E13, :e31, :E31, :e23, :E23, :e32, :E32], [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11])
   F = FreeGroup(A)
   e12, e21, e13, e31, e23, e32 = Groups.gens(F)

   N = 3
   SL₃Ƶ = MatrixAlgebra(zz, N)
   E(SL₃Ƶ, i,j) = (e_ij = one(SL₃Ƶ); e_ij[i,j] = 1; e_ij)
   S = [E(SL₃Ƶ, i,j) for i in 1:N for j in 1:N if i≠j]
   S = unique([S; inv.(S)])

   function h(w::FPGroupElement)
      result = one(SL₃Ƶ)
      for i in 1:length(w.word)
         if w.word[i]%2 ==1
            if Groups.gens(F,floor(Int,w.word[i]/2)+1) == e12
               result *= S[1]
            elseif Groups.gens(F,floor(Int,w.word[i]/2)+1) == e21
               result *= S[3]
            elseif Groups.gens(F,floor(Int,w.word[i]/2)+1) == e13
               result *= S[2]
            elseif Groups.gens(F,floor(Int,w.word[i]/2)+1) == e31
               result *= S[5]
            elseif Groups.gens(F,floor(Int,w.word[i]/2)+1) == e23
               result *= S[4]
            else
               result *= S[6]
            end
         else
            if Groups.gens(F, Int(w.word[i]/2)) == e12
               result *= inv(S[1])
            elseif Groups.gens(F, Int(w.word[i]/2)) == e21
               result *= inv(S[3])
            elseif Groups.gens(F, Int(w.word[i]/2)) == e13
               result *= inv(S[2])
            elseif Groups.gens(F, Int(w.word[i]/2)) == e31
               result *= inv(S[5])
            elseif Groups.gens(F, Int(w.word[i]/2)) == e23
               result *= inv(S[4])
            else
               result *= inv(S[6])
            end
         end
      end

      return result
   end

   relations = [e12*e13*e12^(-1)*e13^(-1), e12*e32*e12^(-1)*e32^(-1), e13*e23*e13^(-1)*e23^(-1), 
                e23*e21*e23^(-1)*e21^(-1), e21*e31*e21^(-1)*e31^(-1), e31*e32*e31^(-1)*e32^(-1),
                e12*e23*e12^(-1)*e23^(-1)*e13^(-1), e13*e32*e13^(-1)*e32^(-1)*e12^(-1), 
                e21*e13*e21^(-1)*e13^(-1)*e23^(-1), e23*e31*e23^(-1)*e31^(-1)*e21^(-1), 
                e31*e12*e31^(-1)*e12^(-1)*e32^(-1), e32*e21*e32^(-1)*e21^(-1)*e31^(-1)]
   
   halfBasis, sizes = Groups.wlmetric_ball(S, radius = halfRadius)
   spectralGapsCertification(h, relations, halfBasis)
end
#########################################################################
