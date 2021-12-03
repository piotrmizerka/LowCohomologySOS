using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = 4
LinearAlgebra.BLAS.set_num_threads(2)
using Kronecker
using IntervalArithmetic

using StarAlgebras
using AbstractAlgebra
using Groups

using JuMP
using SCS
using ProxSDP

include("starAlgebras.jl")
include("FoxDerivatives.jl")
include("SOSDecompositions.jl")
include("certification.jl")