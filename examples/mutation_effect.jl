include("file_access.jl")
include("minimum_disparity.jl")

## mutation effect on rebound probability for the three trials
rebound_prob_1074 = [
    simulate_trial_mle_prob("10-1074",["10-1074"]; n_samples= 2000, mutations = true),
    simulate_trial_mle_prob("10-1074",["10-1074"]; n_samples= 2000, mutations = false)
];

rebound_prob_3BNC = [
    simulate_trial_mle_prob("3BNC",["3BNC117"]; n_samples= 2000, mutations = true),
    simulate_trial_mle_prob("3BNC",["3BNC117"]; n_samples= 2000, mutations = false)
];

rebound_prob_combo = [
    simulate_trial_mle_prob("combo",["10-1074","3BNC117"]; n_samples= 2000, mutations = true),
    simulate_trial_mle_prob("combo",["10-1074","3BNC117"]; n_samples= 2000, mutations = false)
];

fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_mutvsnomut.h5", "w")
fid["10-1074"] =  rebound_prob_1074
fid["3BNC117"] = rebound_prob_3BNC
fid["combo"] =  rebound_prob_combo
close(fid)