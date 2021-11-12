module EscapeSimulator
##
using Random
using Distributions
using SpecialFunctions
using DataFrames
using StatsBase
using BitBasis
using HypergeometricFunctions






include("count_likelihood.jl") # analyzes count data using mle and bayesian methods
include("sensitivity_analysis.jl") # Minimum disparity functions
include("WrightSampler.jl")
include("virus_gillespie.jl")
include("virus_trace.jl")
include("equilibrium.jl")


end # module
