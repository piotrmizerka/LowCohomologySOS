function matrix_in_interval(M, M_interval)
    @assert size(M) == size(M_interval)
    return all(zip(M, M_interval)) do (x, Ix)
        S = supp(x)
        S == supp(Ix) || false
        return all(x(s) == Ix(s) for s in S)
    end
end

@testset "sos_from_matrix" begin

    @testset "Cyclic(3)" begin
        C₃ = cyclic_group(3)
        a, = Groups.gens(C₃)
        half_basis_C_3 = [one(C₃), a, a^2]
        RC₃ = LowCohomologySOS.group_ring(C₃, 3, true)
        Q_1 = [
            2 0 0
            0 3 0
            0 0 4
        ]
        Q_2 = [
            1 0 1
            0 1 0
            1 0 1
        ]

        proper_sos_1 = reshape([29 * one(RC₃)], 1, 1) 
        @test matrix_in_interval(
            proper_sos_1,
            LowCohomologySOS.sos_from_matrix(RC₃, Q_1, half_basis_C_3),
        )

        proper_sos_2 = reshape([5 * one(RC₃) + 2 * RC₃(a) + 2 * RC₃(a^2)], 1, 1)
        @test matrix_in_interval(
            proper_sos_2,
            LowCohomologySOS.sos_from_matrix(RC₃, Q_2, half_basis_C_3),
        )

        @testset "certify" begin
            order_unit_1 = reshape([one(RC₃)], 1, 1)

            X_1a = reshape([33 * one(RC₃)], 1, 1)

            certified_interval_1a = Logging.with_logger(Logging.NullLogger()) do
                return LowCohomologySOS.certify_sos_decomposition(
                    X_1a,
                    order_unit_1,
                    4,
                    Q_1,
                    half_basis_C_3,
                )
            end

            @test 4 ∈ certified_interval_1a

            X_1b = reshape([28 * one(RC₃)], 1, 1)

            certified, certified_interval_1b = Logging.with_logger(Logging.NullLogger()) do
                return LowCohomologySOS.certify_sos_decomposition(
                    X_1b,
                    order_unit_1,
                    1,
                    Q_1,
                    half_basis_C_3,
                )
            end
            @test certified == false
            @test sup(certified_interval_1b) < 0
        end
    end

    @testset "Cyclic(2)" begin

        C₂ = cyclic_group(2) # this causes an error (at least in Windows 10 and online, on Ubuntu 20.4 works :))
        x, = Groups.gens(C₂)
        half_basis_C_2 = [one(C₂), x]
        RC₂ = LowCohomologySOS.group_ring(C₂, 2, true)
        Q_3 = [
            1 0 0 1
            0 1 0 0
            0 0 1 0
            1 0 0 1
        ]

        proper_sos_3 = [
            3*one(RC₂) 2*RC₂(x)
            2*RC₂(x) 3*one(RC₂)
        ]
        @test matrix_in_interval(
            proper_sos_3,
            LowCohomologySOS.sos_from_matrix(RC₂, Q_3, half_basis_C_2),
        )

        @testset "certify" begin
            # actually I did not verify if this is an order unit but this is irrelevant for the sake of the test
            order_unit_2 = [
                zero(RC₂) 2*RC₂(x)
                2*RC₂(x) zero(RC₂)
            ]
            X_2 = [
                3*one(RC₂) 4*RC₂(x)
                4*RC₂(x) 3*one(RC₂)
            ]

            certified, certified_interval_2 = Logging.with_logger(Logging.NullLogger()) do
                return LowCohomologySOS.certify_sos_decomposition(
                    X_2,
                    order_unit_2,
                    1,
                    Q_3,
                    half_basis_C_2,
                )
            end
            @test certified
            @test 1 ∈ certified_interval_2
        end
    end
end
