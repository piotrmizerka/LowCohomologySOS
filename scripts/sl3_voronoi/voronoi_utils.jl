using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)
using Groups
using LowCohomologySOS

include(joinpath(@__DIR__, "../optimizers.jl"))
include(joinpath(@__DIR__, "../utils.jl"))

sl3 = MatrixGroups.SpecialLinearGroup{3}(Int8)
S = gens(sl3)
S_inv =[S; inv.(S)]
ball8, sizes = Groups.wlmetric_ball(S_inv, radius = 8)

function gelt_from_matrix(M::AbstractMatrix, group_subset)
    for g in group_subset
        if MatrixGroups.matrix_repr(g) == M
            return g
        end
    end
end

function saturate(S, stabiliser)
    @assert typeof(first(S)) == eltype(stabiliser)
    saturated_list = [g^(-1)*s*g for s in S for g in stabiliser]
    return collect(Set(saturated_list))
end

function averaged_rep(matrix_list, half_basis, RG)
    sum_ = sum(RG(gelt_from_matrix(M, half_basis)) for M in matrix_list)
    stab_order = length(matrix_list)
    return sum_//stab_order
end

function rg_elt(differential_data)
    return sum(dd[1]*dd[2] for dd in differential_data)
end

m2_arrays = [
    [1 0 0; 0 1 0; 0 0 1], [-1 0 0; 0 -1 0; 0 0 1], [-1 0 0; 0 1 0; 0 0 -1], 
    [1 0 0; 0 -1 0; 0 0 -1], [-1 0 0; 0 0 1; 0 1 0], [1 0 0; 0 0 -1; 0 1 0], 
    [1 0 0; 0 0 1; 0 -1 0], [-1 0 0; 0 0 -1; 0 -1 0], [0 -1 0; 1 0 0; 0 0 1], 
    [0 1 0; -1 0 0; 0 0 1], [0 1 0; 1 0 0; 0 0 -1], [0 -1 0; -1 0 0; 0 0 -1], 
    [0 1 0; 0 0 1; 1 0 0], [0 -1 0; 0 0 -1; 1 0 0], [0 -1 0; 0 0 1; -1 0 0], 
    [0 1 0; 0 0 -1; -1 0 0], [0 0 1; 1 0 0; 0 1 0], [0 0 -1; -1 0 0; 0 1 0], 
    [0 0 -1; 1 0 0; 0 -1 0], [0 0 1; -1 0 0; 0 -1 0], [0 0 -1; 0 1 0; 1 0 0], 
    [0 0 1; 0 -1 0; 1 0 0], [0 0 1; 0 1 0; -1 0 0], [0 0 -1; 0 -1 0; -1 0 0]
]
m31_arrays = [
    [1 0 0; 0 1 0; 0 0 1], [-1 0 0; 0 -1 0; 0 0 1], [-1 0 0; 1 1 0; 0 0 -1], 
    [1 0 0; -1 -1 0; 0 0 -1], [0 1 0; 1 0 0; 0 0 -1], [0 -1 0; -1 0 0; 0 0 -1], 
    [0 -1 0; 1 1 0; 0 0 1], [0 1 0; -1 -1 0; 0 0 1], [-1 -1 0; 1 0 0; 0 0 1], 
    [1 1 0; -1 0 0; 0 0 1], [-1 -1 0; 0 1 0; 0 0 -1], [1 1 0; 0 -1 0; 0 0 -1]
]
m32_arrays = [
    [1 0 0; 0 1 0; 0 0 1], [-1 0 0; 0 -1 0; 1 1 1], [-1 0 0; 0 0 -1; 0 -1 0], 
    [1 0 0; 0 0 1; -1 -1 -1], [1 0 0; -1 -1 -1; 0 1 0], [-1 0 0; 1 1 1; 0 0 -1], 
    [0 -1 0; -1 0 0; 0 0 -1], [0 1 0; 1 0 0; -1 -1 -1], [0 1 0; 0 0 1; 1 0 0], 
    [0 -1 0; 0 0 -1; 1 1 1], [0 -1 0; 1 1 1; -1 0 0], [0 1 0; -1 -1 -1; 0 0 1], 
    [0 0 1; 1 0 0; 0 1 0], [0 0 -1; -1 0 0; 1 1 1], [0 0 -1; 0 -1 0; -1 0 0], 
    [0 0 1; 0 1 0; -1 -1 -1], [0 0 1; -1 -1 -1; 1 0 0], [0 0 -1; 1 1 1; 0 -1 0], 
    [1 1 1; -1 0 0; 0 -1 0], [-1 -1 -1; 1 0 0; 0 0 1], [-1 -1 -1; 0 1 0; 1 0 0], 
    [1 1 1; 0 -1 0; 0 0 -1], [1 1 1; 0 0 -1; -1 0 0], [-1 -1 -1; 0 0 1; 0 1 0]
]
m4_arrays = [
    [1 0 0; 0 1 0; 0 0 1], [-1 0 0; 0 -1 0; 1 0 1], [-1 0 0; 0 0 -1; 0 -1 0], [1 0 0; 0 0 1; -1 -1 0], 
    [-1 0 0; 1 1 0; 0 0 -1], [1 0 0; -1 -1 0; -1 0 -1], [1 0 0; -1 0 -1; 0 1 0], [-1 0 0; 1 0 1; 1 1 0]
]
m2_stab = [gelt_from_matrix(M,ball8) for M in m2_arrays]
m31_stab = [gelt_from_matrix(M,ball8) for M in m31_arrays]
m32_stab = [gelt_from_matrix(M,ball8) for M in m32_arrays]
m4_stab = [gelt_from_matrix(M,ball8) for M in m4_arrays]