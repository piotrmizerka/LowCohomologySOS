include("certification.jl")

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
 
    differentials(G)
 end

 function CₙDifferentials(n)
    differentials(cyclicGroup(n))
 end
 
## Defining optimization problems ###################################################3
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
 
    RG = groupRing(G, 2, true)
    S = collect(RG.basis)
 
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

 ## stuff for tests ############################
let n = 3
    Cₙ = cyclicGroup(n)
    spectralGapsApproximated(Cₙ, n)
 end
 ###############################################

 SL₃ƵShorterPresentationSpectralGaps = let maxRules = 1000, supportSize = 3
   A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
   F = FreeGroup(A)
   x,y,z = Groups.gens(F)
   ε = one(F);
   SL₃ƵShorterPresentation = FPGroup(F, [x^3 => ε, y^3 => ε, z^2 => ε, (x*z)^3 => ε, (y*z)^3 => ε, 
                   (x^(-1)*z*x*y)^2 => ε, (y^(-1)*z*y*x)^2 => ε, (x*y)^6 => ε], maxrules = maxRules )
   
   # spectralGapsApproximated(SL₃ƵShorterPresentation, supportSize)
   differentials(SL₃ƵShorterPresentation)
end

SL₃ƵElementaryMatrixPresentationSpectralGaps = let maxRules = 5000, supportSize = 2
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

SL₃ƵSteinbergGroupSpectralGaps = let maxRules = 5000, supportSize = 2
   A = Alphabet([:e12, :E12, :e21, :E21, :e13, :E13, :e31, :E31, :e23, :E23, :e32, :E32], [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11])
   F = FreeGroup(A)
   e12, e21, e13, e31, e23, e32 = Groups.gens(F)
   ε = one(F);
   SteinbergGroupPresentation = FPGroup(F, [e12*e13 => e13*e12, e12*e32 => e32*e12, e13*e23 => e23*e13, e23*e21 => e21*e23, e21*e31 => e31*e21, e31*e32 => e32*e31,
                   e12*e23 => e13*e23*e12, e13*e32 => e12*e32*e13, e21*e13 => e23*e13*e21, e23*e31 => e21*e31*e23, 
                   e31*e12 => e32*e12*e31, e32*e21 => e31*e21*e32], maxrules = maxRules )

   # differentials(SteinbergGroupPresentation)
   spectralGapsApproximated(SteinbergGroupPresentation, supportSize)
end

###### test examples - certification ###################################
let n = 3
   Cₙ = cyclicGroup(n)
   spectralGapsCertification(Cₙ, n)
end


#########################################################################