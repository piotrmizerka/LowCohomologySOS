using Pkg

Pkg.activate(@__DIR__)
# using Revise # CHECK THIS
# includet("src/LowCohomologySOS.jl") # CHECK THIS
include("src/LowCohomologySOS.jl")


SL₃ƵSpectralGaps = let halfRadius = 2
    A = Alphabet([:e12, :E12, :e21, :E21, :e13, :E13, :e31, :E31, :e23, :E23, :e32, :E32], [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11])
    F = FreeGroup(A)
    e12, e21, e13, e31, e23, e32 = Groups.gens(F)
 
    N = 3
    SL₃Ƶ = MatrixAlgebra(zz, N)
    E(SL₃Ƶ, i,j) = (e_ij = one(SL₃Ƶ); e_ij[i,j] = 1; e_ij)
    S = [E(SL₃Ƶ, i,j) for i in 1:N for j in 1:N if i≠j]
    S = unique([S; inv.(S)])

    function h(u::FPGroupElement)
        result = one(SL₃Ƶ)

        for i in 1:length(word(u))
            if F(word(u)[i:i]) == e12
                result *= S[1]
            elseif F(word(u)[i:i]) == inv(e12)
                result *= inv(S[1])
            elseif F(word(u)[i:i]) == e21
                result *= S[3]
            elseif F(word(u)[i:i]) == inv(e21)
                result *= inv(S[3])
            elseif F(word(u)[i:i]) == e13
                result *= S[2]
            elseif F(word(u)[i:i]) == inv(e13)
                result *= inv(S[2])
            elseif F(word(u)[i:i]) == e31
                result *= S[5]
            elseif F(word(u)[i:i]) == inv(e31)
                result *= inv(S[5])
            elseif F(word(u)[i:i]) == e23
                result *= S[4]
            elseif F(word(u)[i:i]) == inv(e23)
                result *= inv(S[4])
            elseif F(word(u)[i:i]) == e32
                result *= S[6]
            elseif F(word(u)[i:i]) == inv(e32)
                result *= inv(S[6])
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
    spectral_gaps_certification(h, relations, halfBasis)
end
