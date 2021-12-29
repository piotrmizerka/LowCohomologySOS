function matrix_in_interval(M, M_interval)
    for k in 1:length(M)
        for g in supp(M[k]) 
            if M[k](g) ∉ M_interval[k](g)
                return false
            end
        end
    end

    return true
end

C₃ = LowCohomologySOS.cyclic_group(3)
a, = Groups.gens(C₃)
half_basis_C_3 = [one(C₃), a, a^2]
RC₃ = LowCohomologySOS.group_ring(C₃, 3, true)
Q_1 = [2 0 0;
       0 3 0;
       0 0 4]
Q_2 = [1 0 1;
       0 1 0;
       1 0 1]

C₂ = LowCohomologySOS.cyclic_group(2)
x, = Groups.gens(C₂)
half_basis_C_2 = [one(C₂), x]
RC₂ = LowCohomologySOS.group_ring(C₂, 2, true)
Q_3 = [1 0 0 1;
       0 1 0 0;
       0 0 1 0;
       1 0 0 1]

@testset "sos_from_matrix" begin
    proper_sos_1 = reshape([29*one(RC₃)], 1, 1)
    @test matrix_in_interval(proper_sos_1, LowCohomologySOS.sos_from_matrix(Q_1, half_basis_C_3, RC₃))

    proper_sos_2 = reshape([5*one(RC₃)+2*RC₃(a)+2*RC₃(a^2)], 1, 1)
    @test matrix_in_interval(proper_sos_2, LowCohomologySOS.sos_from_matrix(Q_2, half_basis_C_3, RC₃))

    proper_sos_3 = [3*one(RC₂) 2*RC₂(x);
                    2*RC₂(x) 3*one(RC₂)]
    @test matrix_in_interval(proper_sos_3, LowCohomologySOS.sos_from_matrix(Q_3, half_basis_C_2, RC₂))
end

@testset "certify_sos_decomposition" begin
    order_unit_1 = reshape([one(RC₃)], 1, 1)

    X_1a = reshape([33*one(RC₃)], 1, 1)
    certified_interval_1a = @interval(0)
    @suppress begin
        certified_interval_1a = LowCohomologySOS.certify_sos_decomposition(X_1a, order_unit_1, 4, Q_1, half_basis_C_3, RC₃)
    end
    @test 4 ∈ certified_interval_1a

    X_1b = reshape([28*one(RC₃)], 1, 1)
    certified_interval_1b = @interval(0)
    @suppress begin
        certified_interval_1b = LowCohomologySOS.certify_sos_decomposition(X_1b, order_unit_1, 1, Q_1, half_basis_C_3, RC₃)
    end
    @test certified_interval_1b.hi < 0

    order_unit_2 = [zero(RC₂) 2*RC₂(x);
                    2*RC₂(x) zero(RC₂)] # actually I did not verify if this is an order unit but this is irrelevant for the sake of the test
    X_2 = [3*one(RC₂) 4*RC₂(x);
           4*RC₂(x) 3*one(RC₂)]
    certified_interval_2 = @interval(0)
    @suppress begin
        certified_interval_2 = LowCohomologySOS.certify_sos_decomposition(X_2, order_unit_2, 1, Q_3, half_basis_C_2, RC₂)
    end
    @test 1 ∈ certified_interval_2
end

@testset "spectral_gaps_certification" begin
    A = Alphabet([:x, :X], [2, 1])
    F = FreeGroup(A)
    ax, = Groups.gens(F)
    function quotient_hom_1(u::FPGroupElement)
        result = one(C₃)
        for i in 1:length(word(u))
            if F(word(u)[i:i]) == ax
                result *= a
            elseif F(word(u)[i:i]) == inv(ax)
                result *= inv(a)
            end
        end
        return result
    end
    relations = [ax^3]
    certified_interval = @interval(0)
    @suppress begin
        certified_interval = LowCohomologySOS.spectral_gaps_certification(quotient_hom_1, relations, half_basis_C_3, true)
    end
    @test 2.999998988664523 ∈ certified_interval
end