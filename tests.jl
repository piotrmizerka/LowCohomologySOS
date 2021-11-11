include("SOSDecompositions.jl")

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
 
## stuff for tests - cyclic group example ###################################################3
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
 
 # stuff for tests - When running further commands after completing the one below, unexpected problems occur in VSCode
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
 
 ## stuff for tests - another example of an SDP matrix problem #################3
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