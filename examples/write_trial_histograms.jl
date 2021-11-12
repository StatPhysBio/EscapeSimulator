include("file_access.jl")
include("minimum_disparity.jl")
simulate_trial(trial; kwds...) = simulate_trial_bayes(trial,trial_antibodies(trial);kwds...)
outpath = mainfolder*"histogramdata.h5"

## Write the transformed histograms out
h5open(outpath,"w") do fid
    for trial in trial_list
        fid[trial] = simulate_trial(trial;diversity_multiplier = 2.07)
    end
end
# write the untransformed hisotgrams out
h5open(outpath,"r+") do fid
    for trial in trial_list
        fid[trial*"ref"] = simulate_trial(trial;diversity_multiplier = 1.0)
    end
end

##

# h5open(outpath,"w") do fid
#     for trial in trial_list
#         fid[trial] = simulate_trial(trial;diversity_multiplier = 2.45)
#     end
# end