include("file_access.jl")
include("minimum_disparity.jl")

div_mult = sort(vcat(collect(2 .^(.8:.1:1.8)), collect(10 .^(-1:.25:1))))

##
# generate a div_mult x trial matrix of mean divergences
mult_out = 
	reduce(hcat,
	reduce(vcat, 
		concordance_bayes(trial; trefs = [get_observed_rebound_times(trial)], selection_multiplier = 1.0, diversity_multiplier = d, n_pars = 20)
			for d in div_mult) 
			for trial in ["10-1074","3BNC","combo"])
##
bestind = argmin(ind -> sum(mult_out, dims = 2)[ind], eachindex(div_mult))
minind = argmin(ii -> log(div_mult[ii])^2, eachindex(div_mult))
stat_diversity = sum(mult_out, dims = 2)[minind] - sum(mult_out, dims = 2)[bestind] # difference between min and optimal
best_diversity = div_mult[bestind] # actually optimal diversity
##
get_Δr(row) = (row[end] - minimum(row))

stat_dist = reduce((x,y) -> cat(x,y,dims=3), simulate_concordance_bayes(trial, 
    reference_diversity_multiplier = 1.0, # define the null hypothesis
    diversity_vec = collect(2 .^(-1.2:.2:1.2))
    )
    for trial in ["10-1074","3BNC","combo"])

simulated_Δr = get_Δr.(eachrow(sum(stat_dist, dims=3)[:,:,1]))
p_val = mean(stat_diversity .< simulated_Δr)

##
fid = h5open(mainfolder*"disparity_data.h5", "cw")
create_group(fid,"diversity")
divh5 = fid["diversity"]

divh5["theta_mult"] = div_mult
divh5["theta_mult_disparity"] = mult_out
divh5["observed_DeltaR"] = stat_diversity
divh5["simulated_DeltaR"] = simulated_Δr
divh5["p_val"] = p_val

close(fid)

##

sel_mult = sort(vcat(collect(2 .^(-1.8:.1:-.2)), collect(10 .^(-1:.25:1))))

sel_mult_out = 
	reduce(hcat,
	reduce(vcat, 
		concordance_bayes(trial; trefs = [get_observed_rebound_times(trial)], selection_multiplier = s, diversity_multiplier = 2.07, n_pars = 20)
			for s in sel_mult) 
			for trial in ["10-1074","3BNC","combo"])

sel_stat_dist = reduce((x,y) -> cat(x,y,dims=3), simulate_concordance_bayes(trial, 
    reference_diversity_multiplier = 2.07, # define the null hypothesis
    selection_vec = collect(2 .^(-2:.25:2))
    ) for trial in ["10-1074","3BNC","combo"])


simulated_Δrs = get_Δr.(eachrow(sum(sel_stat_dist, dims=3)[:,:,1]))
# bestind = argmin(ind -> sum(mle_mult_out, dims = 2)[ind], eachindex(div_mult))
# minind = argmin(ii -> log(div_mult[ii])^2, eachindex(div_mult))
# stat_diversity = sum(mle_mult_out, dims = 2)[minind] - sum(mle_mult_out, dims = 2)[bestind] # difference between min and optimal

fid = h5open(mainfolder*"disparity_data.h5", "cw")
create_group(fid,"selection")
selh5 = fid["selection"]
selh5["simulated_DeltaR"] = simulated_Δrs
close(fid)