using Revise
using EscapeSimulator
using HDF5
using Statistics
using BenchmarkTools
using BitBasis
using StatsBase
using Random

include("file_access.jl")

# trialsimulations_statistics.h5 is

h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_statistics.h5", "r") do fid
    global reservoir_factor = read(fid["best_fit_multplier"])
end
##


## Combinations

using Combinatorics

function rebound_dict(ab_list, n_antibodies; quartile_estimator_samples = 40) 
    out = Dict{Array{String,1},Array{Float64,1}}() 
    for combo in combinations( ab_list, n_antibodies)
        list = Float64[]
        for _ in 1:quartile_estimator_samples
			p = trial_rebound_prob(get_start_theta("all") .* reservoir_factor, 
                ab_profile_list_bayes(combo); n_samples = 50)
			push!(list, p)
        end
        out[combo] = list
    end
    return out
end

function write_rebound_h5(outdict)
    ii = length(first(keys(outdict)))
    h5open(mainfolder*"therapy_ranking.h5","r+") do fid
        loc = create_group(fid, "$(ii)_antibodies")
        for k in keys(outdict)
            kk = create_group(loc, (*)(k...))
            kk["therapy"] = k
            kk["rebounds"] = outdict[k]
        end
    end
end

h5open(mainfolder*"therapy_ranking.h5","w")

##
out1 = rebound_dict(ablist, 1 ) 
write_rebound_h5(out1)
##
out2 = rebound_dict(ablist, 2 ) 
write_rebound_h5(out2)
##
out3 = rebound_dict(ablist, 3 ) 
write_rebound_h5(out3)
##
