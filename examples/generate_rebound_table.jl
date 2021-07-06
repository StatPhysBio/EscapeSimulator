

## set of helper functions and path definitions for reading in data






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

function trial_traces(trial, ab_list; samples_per_patient = 10, 
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, 
	no_mutations = false, 
	breakpoint = 2,
	timepoints = 0:.5:56
	)
#= returns a 3 d array of viral rebound traces [:,1,:] being wt and [:,2,:] being mutant  =#
	theta_vec = get_start_theta(trial)
	rebound_times = ones(length(theta_vec), samples_per_patient )
	trace = []
	out_mat = permutedims(hcat(collect(timepoints),collect(timepoints)))
	for (ii,θ_0) in enumerate(theta_vec)
		θ = θ_0 * diversity_multiplier
		ws = wright_sampler_array(θ, ab_list; selection_multiplier)
		xxvec = [Tomoko.sample!.(ws) for jj in 1:samples_per_patient ]
		vp = initialize_viral_population(θ,ab_list; selection_multiplier)
		if no_mutations
			vp.mutation_operator = x -> nothing # if mutations are turned off set mutations to do nothing
		end

		for jj in  1:samples_per_patient
			xx = xxvec[jj]
			new = restart_viral_population(vp, xx)
			trace =  virus_time_trace(new, timepoints; breakpoint = 5.0) ./ new.capacity
			out_mat = cat(out_mat, trace, dims = 3)
		end
	end
	return out_mat
end



function simulate_trial_rebounds_matrix(trial, ab_list; samples_per_patient = 10, 
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
	for (ii,θ_0) in enumerate(theta_vec)
		θ = θ_0 * diversity_multiplier
		ws = wright_sampler_array(θ, ab_list; selection_multiplier)
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
	rebound_times = vec(simulate_trial_rebounds_matrix(trial, ab_list;
		samples_per_patient, no_mutations, breakpoint))
	return mean(x < 56.0 for x in rebound_times)
end

function simulate_trial_rebound_prob_bayes(ab_list; samples_per_patient = 10, 
	trial = "all",
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, 
	no_mutations = false, breakpoint = 0.01
	)
	rebound_times = vec(simulate_trial_rebounds_bayes(trial, ab_list; diversity_multiplier,
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

#----------------------------------------------------------------#
# Rebound probability and
#----------------------------------------------------------------#
using Combinatorics

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





#----------------------------------------------------------------#
# Make some plots about mutations
#----------------------------------------------------------------#
##
fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_traces.h5", "w")
fid["10-1074"] = trial_traces("10-1074",["10-1074"], samples_per_patient = 100)
fid["3BNC117"] = trial_traces("3BNC",["3BNC117"], samples_per_patient = 100)
fid["combo"] = trial_traces("combo",["10-1074","3BNC117"], samples_per_patient = 100)
close(fid)
##





using CairoMakie

lines(mult,val1 .- mean(val1))
CairoMakie.lines!(mult,val2 .- mean(val2), color = :red)
CairoMakie.lines!(mult,val3 .- mean(val3), color = :blue)
current_figure()

lines(mult,16 .* val1 .+ 14 .* val2 .+ 7 .* val3)

## Collect data at base parameter values 
#----------------------------------------------------------------#
# Statistical analysis required for reservoir stuff
#----------------------------------------------------------------#
mult = 2 .^(-2.5:.05:2.5)


@time rebs1074 = mapreduce(
	m->simulate_trial_rebounds("10-1074", ["10-1074"]; samples_per_patient = 2400, diversity_multiplier=m),
	hcat,  
	mult
	)
@time rebs3BNC = mapreduce(
		m->simulate_trial_rebounds("3BNC", ["3BNC117"]; samples_per_patient = 2400, diversity_multiplier=m),
		hcat,  
		mult
	)
@time rebscombo = mapreduce(
		m->simulate_trial_rebounds("combo", ["3BNC117","10-1074"]; samples_per_patient = 2400, diversity_multiplier=m),
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
##
mat1 = simulate_trial_rebounds_matrix("10-1074",["10-1074"],samples_per_patient = 2000)
mat21 = simulate_trial_rebounds_matrix("3BNC",["3BNC117"],samples_per_patient = 2000)
mat22 = simulate_trial_rebounds_matrix("3BNC",["3BNC117"],samples_per_patient = 2000)
mat2 = mapreduce(
	(x,y)->vcat(x,y)[sample(1:18,14,replace = false)], hcat,
	eachcol(mat21), eachcol(mat22)
	)
mat3 = simulate_trial_rebounds_matrix("combo",["10-1074","3BNC117"],samples_per_patient = 2000)

@time dist = map(tempfun,
	eachcol(mat1),eachcol(mat2),eachcol(mat3))

(observed, best_fit_multplier) = tempfun(get_observed_rebound_times("10-1074"),get_observed_rebound_times("3BNC"),get_observed_rebound_times("combo"))
concord = tempfun2(get_observed_rebound_times("10-1074"),get_observed_rebound_times("3BNC"),get_observed_rebound_times("combo"))



##
tmat1 = simulate_trial_rebounds_matrix("10-1074",["10-1074"],samples_per_patient = 2000, diversity_multiplier = best_fit_multplier)
tmat21 = simulate_trial_rebounds_matrix("3BNC",["3BNC117"],samples_per_patient = 2000, diversity_multiplier = best_fit_multplier)
tmat22 = simulate_trial_rebounds_matrix("3BNC",["3BNC117"],samples_per_patient = 2000, diversity_multiplier = best_fit_multplier)
tmat2 = mapreduce(
	(x,y)->vcat(x,y)[sample(1:18,14,replace = false)], hcat,
	eachcol(tmat21), eachcol(tmat22)
	)
tmat3 = simulate_trial_rebounds_matrix("combo",["10-1074","3BNC117"],samples_per_patient = 2000, diversity_multiplier = best_fit_multplier)
##

fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_statistics.h5", "w")
fid["multiplier"] = mult
fid["statistic"] = dist
fid["observed_stat"] = observed
fid["best_fit_multplier"] = best_fit_multplier
fid["untransformed/10-1074"] = mat1
fid["untransformed/3BNC117"] = mat2
fid["untransformed/combo"] = mat3
fid["transformed/10-1074"] = tmat1
fid["transformed/3BNC117"] = tmat2
fid["transformed/combo"] = tmat3
fid["concordvec"] = concord
close(fid)


pvalue = 1-ecdf(dist)(observed)

## Collect data at base parameter values 
@time rebs1074 = simulate_trial_rebounds("10-1074", ["10-1074"]; samples_per_patient = 300)
@time rebs3BNC = simulate_trial_rebounds("3BNC", ["3BNC117"]; samples_per_patient = 300)
@time rebscombo = simulate_trial_rebounds("combo", ["3BNC117","10-1074"]; samples_per_patient = 300)
fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_279.h5", "w")
fid["10-1074"] = rebs1074
fid["3BNC117"] = rebs3BNC
fid["combo"] = rebscombo
close(fid)

## mutation effect on rebound probability for the three trials
 rebound_prob_1074 = [simulate_trial_rebound_prob("10-1074",["10-1074"]; samples_per_patient = 10000, no_mutations = false),
 simulate_trial_rebound_prob("10-1074",["10-1074"]; samples_per_patient = 10000, no_mutations = true)];

 rebound_prob_3BNC = [simulate_trial_rebound_prob("3BNC",["3BNC117"]; samples_per_patient = 10000, no_mutations = false),
 simulate_trial_rebound_prob("3BNC",["3BNC117"]; samples_per_patient = 10000, no_mutations = true)];

 rebound_prob_combo = [simulate_trial_rebound_prob("combo",["10-1074","3BNC117"]; samples_per_patient = 20000, no_mutations = false),
 simulate_trial_rebound_prob("combo",["10-1074","3BNC117"]; samples_per_patient = 20000, no_mutations = true)];

 fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_mutvsnomut.h5", "w")
 fid["10-1074"] =  rebound_prob_1074
 fid["3BNC117"] = rebound_prob_3BNC
 fid["combo"] =  rebound_prob_combo
 close(fid)

## Generate the distribution of rebound probabilities

out = rebound_dict(ablist, 2)



 ## Tests and scratchwork








ab_list = ["10-1074"]
##
vp = initialize_viral_population(.02,ab_list)
@time (t_rebound,trace) =  simulate_viral_rebound_time(vp::ViralPopulation)
##
@btime vp = initialize_viral_population(.1,ab_list)

@btime vp2 = restart_viral_population(vp, [.5,.7,.8])

trace = virus_time_trace(vp, 0:56)

function packbits2(vec::BitVector)
	x::Int64 = 0
	for ii in Iterators.reverse(vec)
		x += ii
		x <<= 1
	end
	x >>=1
	return x
end

function random_int2(frequency_vector::Vector{Float64})
	packbits(rand(length(frequency_vector)) .> frequency_vector)
end

##
lines(0:56, trace[1,:] .+ tracey[2,:], color = :blue)
lines!(0:56, trace[1,:], color = :red)
current_figure()

##