
include("file_access.jl")


prof = ab_profile_list_mle(["10-1074","3BNC117"])
vp = initialize_viral_population(.4, prof, antibody = false)
times  = 0:10:500;

out = ensemble_evolution(.1, prof, times; mut_per_gen = 10^(-5.), n_ensemble = 1000, λ = 1.0)
##


using Plots

plot(transpose(dropdims(mean(out,dims = 3),dims=3)))


##
h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/equilibrium", "w") do fid
    fid["sims"] = out
end