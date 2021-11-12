

function concordance_vector(
	t_rebounds, # unadjusted simulations n_patients list of rebound times
	simulation_matrix # (trebounds x sims)
	)
	(_, sims) = size(simulation_matrix)
	out = zeros(sims)
	Threads.@threads for ii in 1:sims
		out[ii] = length(t_rebounds) * observed_discrepancy(t_rebounds,simulation_matrix[:,ii])
	end
	#concordance_vector = map(x->length(t_rebounds) * Tomoko.observed_discrepancy(t_rebounds,x),
	#	 eachcol(simulation_matrix))
	return out # returns a vector as long as the number of columns of simulation matrix
end

function generate_statistic_closure(
	mult, # vector of multipliers
	rebsmats... # a list as long as the number of trials
	)
	let null = argmin(ii -> log(mult[ii])^2,  eachindex(mult)); # find the multplier closest to 1
	function stat_fun(t_rebs...) # returns a function of the true rebounds
		total_concord = sum(
			concordance_vector(t_reb,reb_mat) for (t_reb, reb_mat) in Iterators.zip(t_rebs,rebsmats)
			)
		min = argmin(ii -> total_concord[ii], eachindex(total_concord)) 
		(total_concord[null] - total_concord[min],mult[min]) # that function returns the minimial disparity and the observed mulitplier at the minimum
	end
	return stat_fun
	end
end

function generate_concord_vector_closure(
	mult, # vector of multipliers
	rebsmats... # a tuple as long as the number of trials being analyzed together
	# contains distributional information about the 
	)
	let null = argmin(ii -> log(mult[ii])^2,  eachindex(mult)); # determine null hypothesis
	function stat_fun(t_rebs...)
		total_concord = sum(
			concordance_vector(t_reb,reb_mat) for (t_reb, reb_mat) in Iterators.zip(t_rebs,rebsmats))
		return total_concord
	end
	return stat_fun
	end
end


function concordance_mle(trial; trefs = [get_observed_rebound_times(trial)], kwds...)
	out = simulate_trial_mle(trial, trial_antibodies(trial); kwds...)
	[length(get_observed_rebound_times(trial))*observed_discrepancy(tref, vec(out)) for tref in trefs]
end



function concordance_bayes(trial; trefs = [get_observed_rebound_times(trial)], kwds...)
	out = simulate_trial_bayes(trial, trial_antibodies(trial); kwds...) # patients x replicates x parameters
	[length(get_observed_rebound_times(trial))*mean(observed_discrepancy(tref, vec(parsim)) # take the mean over all parameter values
		for parsim in eachslice(out, dims = 3)) for tref in trefs]
end

function smooshy(array)
	[ii for ii in eachcol(reshape(array, first(size(array)), :))]
end

function simulate_concordance_bayes(trial; reference_diversity_multiplier = 1.0, 
	reference_selection_multiplier = 1.0, diversity_vec = [], selection_vec = [], kwds...)
	refs = smooshy(simulate_trial_bayes(trial, trial_antibodies(trial); 
		selection_multiplier = reference_selection_multiplier, 
		diversity_multiplier = reference_diversity_multiplier,
		n_pars = 20,
		n_runs = 10
		))
	if !isempty(diversity_vec)
		return reduce(hcat, concordance_bayes(trial; trefs = refs, 
			diversity_multiplier = d, 
			selection_multiplier = reference_selection_multiplier) 
				for d in vcat(diversity_vec, reference_diversity_multiplier))
	elseif !isempty(selection_vec)
		return reduce(hcat, concordance_bayes(trial; trefs = refs, n_pars = 20, 
		selection_multiplier = d, 
		diversity_multiplier = reference_diversity_multiplier) 
			for d in vcat(selection_vec, reference_selection_multiplier))
	end
end
# returns a vector of the null hypothesis disparity and then the disparity under 
# all the alternate hypotheses. Number of simulations is size of refs.

#= 
Example usage
optimal_diversity = 
	reduce(hcat,
	reduce(vcat, 
		concordance_bayes(trial; trefs = [get_observed_rebound_times(trial)], selection_multiplier = s, diversity_multiplier = 2.07, n_pars = 100)
			for s in sel_mult) 
			for trial in ["10-1074","3BNC","combo"])



out = 
	reduce(hcat,
	reduce(vcat, 
		concordance_bayes(trial; trefs = [get_observed_rebound_times(trial)], selection_multiplier = s, diversity_multiplier = 2.07, n_pars = 100)
			for s in sel_mult) 
			for trial in ["10-1074","3BNC","combo"])
generates a list of concordance matrices

stat_dist = 
	simulate_concordance_bayes("10-1074", reference_diversity_multiplier = 2.07, selection_vec = [.8,1.1])
=#


# function discrepancy_stochastic_rebound(trial, ab_list ; 
#     data_size=16, min= 1, max=56, # parameters of the data
#     samples_per_patient = 10^1,
# 	selection_multiplier = 1.0,
# 	diversity_multiplier = 1.0) # * length of θ_vec = generated samples
# 	t_rebounds = get_observed_rebound_times(trial)
# 	t_simulated = simulate_trial_rebounds(trial, ab_list; 
# 		samples_per_patient, selection_multiplier, diversity_multiplier)
# 	Tomoko.observed_discrepancy(t_rebounds, t_simulated)
# end