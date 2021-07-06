function concordance_vector(
	t_rebounds, # unadjusted simulations n_patients list of rebound times
	simulation_matrix # (sims vs. trebounds)
	)
	(_, sims) = size(simulation_matrix)
	out = zeros(sims)
	Threads.@threads for ii in 1:sims
		out[ii] = length(t_rebounds) * Tomoko.observed_discrepancy(t_rebounds,simulation_matrix[:,ii])
	end
	#concordance_vector = map(x->length(t_rebounds) * Tomoko.observed_discrepancy(t_rebounds,x),
	#	 eachcol(simulation_matrix))
	return out
end

function generate_statistic_closure(
	mult, # vector of multipliers
	rebsmats... # a list as long as the number of trials
	)
	let null = argmin(ii -> log(mult[ii])^2,  eachindex(mult));
	function stat_fun(t_rebs...)
		total_concord = sum(
			concordance_vector(t_reb,reb_mat) for (t_reb, reb_mat) in Iterators.zip(t_rebs,rebsmats)
			)
		min = argmin(ii -> total_concord[ii], eachindex(total_concord)) 
		(total_concord[null] - total_concord[min],mult[min])
	end
	return stat_fun
	end
end

function generate_concord_vector_closure(
	mult, # vector of multipliers
	rebsmats... # a tuple as long as the number of trials being analyzed together
	# contains distributional information about the 
	)
	let null = argmin(ii -> log(mult[ii])^2,  eachindex(mult));
	function stat_fun(t_rebs...)
		total_concord = sum(
			concordance_vector(t_reb,reb_mat) for (t_reb, reb_mat) in Iterators.zip(t_rebs,rebsmats)
			)
	end
	return stat_fun
	end
end



function discrepancy_stochastic_rebound(trial, ab_list ; 
    data_size=16, min= 1, max=56, # parameters of the data
    samples_per_patient = 10^1,
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0) # * length of θ_vec = generated samples
	t_rebounds = get_observed_rebound_times(trial)
	t_simulated = simulate_trial_rebounds(trial, ab_list; 
		samples_per_patient, selection_multiplier, diversity_multiplier)
	Tomoko.observed_discrepancy(t_rebounds, t_simulated)
end