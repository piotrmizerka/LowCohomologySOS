using Pkg

Pkg.activate(@__DIR__)
using Revise # I DON'T UNDERSTAND HOW IT LOADS THE FILES - GOT ISSUES WITH THIS
# includet("src/SOSStarAlgebras.jl") # I DON'T UNDERSTAND HOW IT LOADS THE FILES - GOT ISSUES WITH THIS
include("src/SOSStarAlgebras.jl")

RG, = let n = 3
    G = cyclicGroup(n)
    RG, halfBasis = groupRing(G,n)
end

function cyclicGroup2(n)
    A = Alphabet([:a, :A], [2,1])
    F = FreeGroup(A)
    a, = Groups.gens(F)
    e = one(F)
    Cₙ = FPGroup(F, [a^n => e])

    return Cₙ
end

let
    

    x = [0 1 0;0 0 1;1 0 0]
    z = [0 1 0;1 0 0;-1 -1 -1]
    id = [1 0 0;0 1 0;0 0 1]

    function check_y(y)
        # println(Int(det(y)))

        if Int(det(y)) != 1
            return false
        end

        if y^3 != id 
            return false
        end
        if (y*z)^3 != id
            return false
        end
        if (x^(-1)*z*x*y)^2 != id
            return false
        end
        if (y^(-1)*z*y*x)^2 != id
            return false
        end
        if (x*y)^6 != id
            return false
        end

        return true
    end

    # println("dsfsff")

    S = [-1,0,1]
    for a11 in S
        for a12 in S
            for a13 in S
                for a21 in S
                    for a22 in S
                        for a23 in S
                            for a31 in S
                                for a32 in S
                                    for a33 in S
                                        y = [a11 a12 a13;a21 a22 a23;a31 a32 a33]

                                        if check_y(y)
                                            println(y)
                                            return
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
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

SL₃ƵShorterPresentationSpectralGaps = let halfRadius = 5
    A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
    F = FreeGroup(A)
    x,y,z = Groups.gens(F)

    relations = [x^3, y^3, z^2]
    # relations = [x^3, y^3, z^2, (x*z)^3, (y*z)^3, (x^(-1)*z*x*y)^2, (y^(-1)*z*y*x)^2, (x*y)^6]

    N = 3
    SL₃Ƶ = MatrixAlgebra(zz, N)
    E(SL₃Ƶ, i,j) = (e_ij = one(SL₃Ƶ); e_ij[i,j] = 1; e_ij)
    S = [E(SL₃Ƶ, i,j) for i in 1:N for j in 1:N if i≠j]
    S = unique([S; inv.(S)])
    halfBasis, sizes = Groups.wlmetric_ball(S, radius = halfRadius)
    hx = one(SL₃Ƶ)
    hx[1,1] = 0
    hx[1,2] = 1
    hx[1,3] = 0
    hx[2,1] = 0
    hx[2,2] = 0
    hx[2,3] = 1
    hx[3,1] = 1
    hx[3,2] = 0
    hx[3,3] = 0
    hy = one(SL₃Ƶ)
    hy[1,1] = 1
    hy[1,2] = 0
    hy[1,3] = 1
    hy[2,1] = 0
    hy[2,2] = -1
    hy[2,3] = -1
    hy[3,1] = 0
    hy[3,2] = 1
    hy[3,3] = 0
    hz = one(SL₃Ƶ)
    hz[1,1] = 0
    hz[1,2] = 1
    hz[1,3] = 0
    hz[2,1] = 1
    hz[2,2] = 0
    hz[2,3] = 0
    hz[3,1] = -1
    hz[3,2] = -1
    hz[3,3] = -1
    # append!(S,[hx,inv(hx),hy,inv(hy),hz,inv(hz)])
    # S = unique(S)
    S = unique([hx,inv(hx),hy,inv(hy),hz,inv(hz)])

    function h(w::FPGroupElement)
        result = one(SL₃Ƶ)
        for i in 1:length(w.word)
            if parent(w)(word(w)[i:i]) == x
                result *= hx
            elseif parent(w)(word(w)[i:i]) == inv(x)
                result *= inv(hx)
            elseif parent(w)(word(w)[i:i]) == y
                result *= hy
            elseif parent(w)(word(w)[i:i]) == inv(y)
                result *= inv(hy)
            elseif parent(w)(word(w)[i:i]) == z
                result *= hz
            elseif parent(w)(word(w)[i:i]) == inv(z)
                result *= inv(hz)
            end
        end
  
        return result
     end

     halfBasis, sizes = Groups.wlmetric_ball(S, radius = halfRadius)

    #  @info length(halfBasis)

    spectralGapsCertification(h, relations, halfBasis)
end