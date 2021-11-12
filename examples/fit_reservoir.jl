using HDF5: eachindex




include("file_access.jl")
include("minimum_disparity.jl")
##

## Collect data at base parameter values 
#----------------------------------------------------------------#
# Statistical analysis required for reservoir stuff
#----------------------------------------------------------------#
mult = 2 .^(-1:.1:2.5)
sel_mult = 2 .^(-2.5:.2:1)

## Test diversity
@time rebs1074 = mapreduce(
	m->vec(simulate_trial_mle("10-1074", ["10-1074"]; n_samples = 240, diversity_multiplier=m)),
	hcat,  
	mult
	)

@time rebs3BNC =  mapreduce(
	m->vec(simulate_trial_mle("3BNC", ["3BNC117"]; n_samples = 240, diversity_multiplier=m)),
	hcat,  
	mult
	)
@time rebscombo =  mapreduce(
	m->vec(simulate_trial_mle("combo", ["10-1074","3BNC117"]; n_samples = 240, diversity_multiplier=m)),
	hcat,  
	mult
	)
fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_diversity.h5", "w")
fid["10-1074"] = rebs1074
fid["3BNC117"] = rebs3BNC
fid["combo"] = rebscombo
close(fid)
##


tempfun = generate_statistic_closure(
	mult, # vector of multipliers
	rebs1074, rebs3BNC, rebscombo
	)

tempfun2 = generate_concord_vector_closure(
		mult, # vector of multipliers
		rebs1074, rebs3BNC, rebscombo
		)
## for generating the null distribution
mat1 = simulate_trial_mle("10-1074",["10-1074"], n_samples = 200)
mat21 = simulate_trial_mle("3BNC",["3BNC117"], n_samples = 200)
mat22 = simulate_trial_mle("3BNC",["3BNC117"], n_samples = 200)
mat2 = mapreduce(
	(x,y)->vcat(x,y)[sample(1:18,14,replace = false)], hcat,
	eachcol(mat21), eachcol(mat22)
	)
mat3 = simulate_trial_mle("combo",["10-1074","3BNC117"], n_samples = 200)
##
@time dist = map(x->tempfun(x...),
	zip(eachcol(mat1), eachcol(mat2), eachcol(mat3)))
dist2 = map(dist) do x 
    x[1] 
    end
(observed, best_fit_multplier) = tempfun(get_observed_rebound_times("10-1074"),get_observed_rebound_times("3BNC"),get_observed_rebound_times("combo"))
concord = tempfun2(get_observed_rebound_times("10-1074"),get_observed_rebound_times("3BNC"),get_observed_rebound_times("combo"))



##
mtest = best_fit_multplier
newλ = 1.0
tmat1 = simulate_trial_mle("10-1074",["10-1074"];  diversity_multiplier= mtest, n_samples = 200, λ =newλ)
tmat21 = simulate_trial_mle("3BNC",["3BNC117"]; diversity_multiplier= mtest,  n_samples = 200, λ = newλ)
tmat22 = simulate_trial_mle("3BNC",["3BNC117"]; diversity_multiplier= mtest, n_samples = 200, λ = newλ)
tmat2 = mapreduce(
	(x,y)->vcat(x,y)[sample(1:18,14,replace = false)], hcat,
	eachcol(tmat21), eachcol(tmat22)
	)
tmat3 = simulate_trial_mle("combo",["10-1074","3BNC117"]; diversity_multiplier= mtest, n_samples = 200,  λ = newλ, mut_per_gen = μtest)
##

h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_statistics.h5", "r+") do fid
    # fid["multiplier"] = mult
    # fid["statistic"] = dist2
    # fid["observed_stat"] = observed
    # fid["best_fit_multplier"] = best_fit_multplier
    # fid["untransformed/10-1074"] = mat1
    # fid["untransformed/3BNC117"] = mat2
    # fid["untransformed/combo"] = mat3
    # fid["lowmut/10-1074"] = tmat1
    # fid["lowmut/3BNC117"] = tmat2
    # fid["lowmut/combo"] = tmat3
    fid["himut/10-1074"] = tmat1
    fid["himut/3BNC117"] = tmat2
    fid["himut/combo"] = tmat3
    # fid["transformedl/10-1074"] = tmat1
    # fid["transformedl/3BNC117"] = tmat2
    # fid["transformedl/combo"] = tmat3
    # fid["concordvec"] = concord
end


## Bayesian analysis of diveristy

tmat1 = simulate_trial_bayes("10-1074",["10-1074"]; diversity_multiplier= best_fit_multplier)
tmat3 =simulate_trial_bayes("combo",["10-1074","3BNC117"], diversity_multiplier= best_fit_multplier)
tmat21 = simulate_trial_bayes("3BNC", ["3BNC117"]; diversity_multiplier= best_fit_multplier)
tmat22 = simulate_trial_bayes("3BNC", ["3BNC117"]; diversity_multiplier= best_fit_multplier)

##
smooshy(x) = reshape(x, size(x)[1], :) # get everything into a 
tmat1 = smooshy(tmat1)
tmat21 = smooshy(tmat21)
tmat22 = smooshy(tmat22)
tmat3 = smooshy(tmat3)
tmat2 = mapreduce(
	(x,y)->vcat(x,y)[sample(1:18,14,replace = false)], hcat,
	eachcol(tmat21), eachcol(tmat22)
	)
##
fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_bayesian.h5", "w")
fid["transformed/10-1074"] = tmat1
fid["transformed/3BNC117"] = tmat2
fid["transformed/combo"] = tmat3
close(fid)

##

pvalue = 1-ecdf(dist2)(observed)

## Collect data at base parameter values 
@time rebs1074 = simulate_trial_rebounds("10-1074", ["10-1074"]; samples_per_patient = 300)
@time rebs3BNC = simulate_trial_rebounds("3BNC", ["3BNC117"]; samples_per_patient = 300)
@time rebscombo = simulate_trial_rebounds("combo", ["3BNC117","10-1074"]; samples_per_patient = 300)
fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_279.h5", "w")
fid["10-1074"] = rebs1074
fid["3BNC117"] = rebs3BNC
fid["combo"] = rebscombo
close(fid)



## Generate the distribution of rebound probabilities

out = rebound_dict(ablist, 2)

## Test selection
@time rebs1074 = mapreduce(
	m->vec(simulate_trial_mle("10-1074", ["10-1074"]; n_samples = 240, diversity_multiplier=m)),
	hcat,  
	mult
	)

@time rebs3BNC =  mapreduce(
	m->vec(simulate_trial_mle("3BNC", ["3BNC117"]; n_samples = 240, diversity_multiplier=m)),
	hcat,  
	mult
	)
@time rebscombo =  mapreduce(
	m->vec(simulate_trial_mle("combo", ["10-1074","3BNC117"]; n_samples = 240, diversity_multiplier=m)),
	hcat,  
	mult
	)
fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_diversity.h5", "w")
fid["10-1074"] = rebs1074
fid["3BNC117"] = rebs3BNC
fid["combo"] = rebscombo
close(fid)

# generate rebound times for different trials, averaging over parameters and runs, comparing to t the true values...
##
test1 = simulate_trial_bayes("3BNC",["3BNC117"], diversity_multiplier = 2.07)

##
fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_diversity.h5", "w")
fid["10-1074"] = rebs1074
fid["3BNC117"] = rebs3BNC
fid["combo"] = rebscombo
close(fid)

## 
div_mult = 2 .^(-1:.1:2.5)
sel_mult = 10 .^(-1:.1:1)

## Current best data for reservoir
mult_out = 
	reduce(hcat,
	reduce(vcat, 
		concordance_bayes(trial; trefs = [get_observed_rebound_times(trial)], selection_multiplier = 1.0, diversity_multiplier = d, n_pars = 100)
			for d in div_mult) 
			for trial in ["10-1074","3BNC","combo"])
##
bestind = argmin(ind -> sum(mult_out, dims = 2)[ind], eachindex(div_mult))
minind = argmin(ii -> log(div_mult[ii])^2, eachindex(div_mult))
stat_diversity = sum(mult_out, dims = 2)[minind] - sum(mult_out, dims = 2)[bestind] # difference between min and optimal
best_diversity = div_mult[bestind] # actually optimal diversity
##
sel_out = 
	reduce(hcat,
	reduce(vcat, 
		concordance_bayes(trial; trefs = [get_observed_rebound_times(trial)], selection_multiplier = s, diversity_multiplier = best_diversity, n_pars = 100)
			for s in sel_mult) 
			for trial in ["10-1074","3BNC","combo"])
##


##

h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_selection.h5", "w") do fid
	fid["10-1074"] = out[:,1]
	fid["3BNC117"] = out[:,2]
	fid["combo"] = out[:,3]
	fid["sel_mult"] = sel_mult
end
