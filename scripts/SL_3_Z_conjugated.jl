using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using LowCohomologySOS
using Groups
import Groups.MatrixGroups
using JuMP

include(joinpath(@__DIR__, "optimizers.jl"))

SL₃ℤ_spectral_gaps = let half_radius = 1
    SL(n, R) = MatrixGroups.SpecialLinearGroup{n}(R)
    SL₃ℤ = SL(3, Int8)

    λ, P, termination_status = LowCohomologySOS.property_t_conjugated_approx(SL₃ℤ, half_radius, optimizer = scs_opt(eps = 1e-5, max_iters = 100_000))
end
