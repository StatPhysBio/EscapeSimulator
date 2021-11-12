

using Combinatorics

# Functions for running simulations of trials

include("file_access.jl")

#= main function returns a matrix of samples organized by patient
	This allows us to simulate trial outcomes when choosing without replacement in fit analysis.
	In the fit analysis, we assume no multiplicity in the θ resampling
	But also gets flattened for raw histogram construction.
	no_mutations keyword turns off mutations after treatment.
	multiplier keywords are used in the sensitivity analysis. =#


##
#-------------------------------
# Actual gillespie code 
#-------------------------------


##
# --------------------------------------------------------------
# Funcitons for generating traces and measuring quantities
# --------------------------------------------------------------



# function simulate_trial_rebounds(trial, ab_list; samples_per_patient = 10)
# 	theta_vec = get_start_theta(trial)
# 	rebound_times = Float64[]
# 	for θ in theta_vec
# 		for _ in  1:samples_per_patient
# 			vp = initialize_viral_population(θ,ab_list)
# 			(t_rebound,_) =  simulate_viral_rebound_time(vp)
# 			push!(rebound_times,t_rebound)
# 		end
# 	end
# 	return rebound_times
# end




function simulate_trial_rebounds(trial, ab_list; samples_per_patient = 10, 
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0
	)
	rebound_times = vec(simulate_trial_rebounds_matrix(trial, ab_list;
		samples_per_patient, selection_multiplier, diversity_multiplier))
	return rebound_times
end

function simulate_trial_rebound_prob(trial, ab_list; samples_per_patient = 10, 
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, 
	no_mutations = false, breakpoint = .05
	)
	rebound_times = vec(simulate_trial_rebounds_matrix(trial, ab_list; diversity_multiplier, selection_multiplier,
		samples_per_patient, no_mutations, breakpoint))
	return mean(x < 56.0 for x in rebound_times)
end

function simulate_trial_rebound_prob_bayes(ab_list; samples_per_patient = 10, 
	trial = "all",
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, 
	no_mutations = false, breakpoint = 0.01
	)
	rebound_times = vec(simulate_trial_rebounds_bayes(trial, ab_list; diversity_multiplier, selection_multiplier,
		samples_per_patient, no_mutations, breakpoint))
	return mean(x < 56.0 for x in rebound_times)
end



function simulate_trial_rebounds_bayes(trial, ab_list; samples_per_patient = 10, 
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, 
	no_mutations = false, 
	breakpoint = .8
	)
#= main function returns a matrix of samples organized by patient
	This allows us to simulate trial outcomes when choosing without replacement in fit analysis.
	But also gets flattened for raw histogram construction.
	no_mutations keyword turns off mutations after treatment.
	multiplier keywords are used in the sensitivity analysis. =#
	theta_vec = get_start_theta(trial)
	rebound_times = ones(length(theta_vec), samples_per_patient )
	ws_fun = wright_sampler_array_bayes_closure(ab_list; selection_multiplier)
	# creates a random function of θ to generate correlated frequencies
	for (ii,θ_0) in enumerate(theta_vec)
		θ = θ_0 * diversity_multiplier
		ws = ws_fun(θ)
		xxvec = [Tomoko.sample!.(ws) for jj in 1:samples_per_patient ]
		vp = initialize_viral_population(θ,ab_list; selection_multiplier)
		if no_mutations
			vp.mutation_operator = x -> nothing # if mutations are turned off set mutations to do nothing
		end
		Threads.@threads for jj in  1:samples_per_patient
			xx = xxvec[jj]
			new = restart_viral_population(vp, xx)
			(t_rebound,_) =  simulate_viral_rebound_time(new; breakpoint)
			rebound_times[ii,jj] = t_rebound
		end
	end
	return rebound_times
end







## 
# --------------------------------
# Statistical analysis
# --------------------------------



###

##
function concordance_vector(
	t_rebounds, # unadjusted simulations n_patients list of rebound times
	simulation_matrix # (sims vs. trebounds)
	)
	(_, sims) = size(simulation_matrix)
	out = zeros(sims)
	Threads.@threads for ii in 1:sims
		out[ii] = length(t_rebounds) * observed_discrepancy(t_rebounds,simulation_matrix[:,ii])
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

#----------------------------------------------------------------#
# Rebound probability 
#----------------------------------------------------------------#
function rebound_dict(ab_list, n_antibodies; 
		diversity_multiplier = 2.07, # best fit multiplier
		samples_per_patient = 8, 
		quartile_estimator_samples = 20) 
    out = Dict{Array{String,1},Array{Float64,1}}()
    for combo in combinations( ab_list, n_antibodies)
        list = Float64[]
        for ii in 1:quartile_estimator_samples
			p = simulate_trial_rebound_prob_bayes(combo; samples_per_patient, diversity_multiplier)
			push!(list, p)
        end
        out[combo] = list
    end
    return out
end




# future function for making bayesian inference
# function run_patient_bayes(trial, ab_list;
# 	trial_runs = 40,
# 	samples_per_patient = 100, 
# 	selection_multiplier = 1.0,
# 	diversity_multiplier = 1.0, 
# 	no_mutations = false, 
# 	breakpoint = .8
# 	)
# #=
# 	Returns an array with the first dimension runs, the second dimension thetas 
# 	and the third dimension bayesian realizations.
# =#
# 	theta_vec = get_start_theta(trial)
# 	ab_profile_iterator = Iterators.flatten(repeat(
# 		ab_profile_list_bayes(ablist; selection_multiplier),
# 		length(theta_vec)) for ii in 1:trial_runs)
# 	theta_iterator = Iterators.flatten(theta_vec for ii in 1:trial_runs)
# 	rebound_times = ones(length(theta_vec), trial_runs, samples_per_patient)
# 	for (ii,θ_0) in enumerate(zip(theta_iterator, ab_profile_iterator))
# 		θ = θ_0 * diversity_multiplier
# 		ws = wright_sampler_array(θ, ab_list; selection_multiplier)
# 		xxvec = [Tomoko.sample!.(ws) for jj in 1:samples_per_patient ]
# 		vp = initialize_viral_population(θ,ab_list; selection_multiplier)
# 		if no_mutations
# 			vp.mutation_operator = x -> nothing # if mutations are turned off set mutations to do nothing
# 		end
# 		Threads.@threads for jj in  1:samples_per_patient
# 			xx = xxvec[jj]
# 			new = restart_viral_population(vp, xx)
# 			(t_rebound,_) =  simulate_viral_rebound_time(new; breakpoint)
# 			rebound_times[ii,jj] = t_rebound
# 		end
# 	end
# 	return rebound_times
# end