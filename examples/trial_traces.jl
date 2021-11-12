include("file_access.jl")


function trial_traces(trial, ab_list; samples_per_patient = 10, 
	# save full traces for plotting trajectory swarm
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, 
	no_mutations = false, 
	breakpoint = 2,
	timepoints = 0:.2:56
	)
	theta_vec = get_start_theta(trial)
	trace = []
	out_mat = permutedims(hcat(collect(timepoints),collect(timepoints)))
	for (ii,θ_0) in enumerate(theta_vec .* diversity_multiplier)
		ab_profile_list = ab_profile_list_bayes(ab_list; selection_multiplier)
		vp = initialize_viral_population(θ_0, ab_profile_list; λ = 20.0)
		if no_mutations
			vp.mutation_operator = x -> nothing # if mutations are turned off set mutations to do nothing
		end
		for jj in  1:samples_per_patient
			restart_viral_population!(vp)
			tr = virus_time_trace(vp, timepoints; breakpoint)
			trace =  hcat(map(tr) do x 
				[x...] ./ vp.capacity
			end...)
			out_mat = cat(out_mat, trace, dims = 3)
		end
	end
	return out_mat
end

##
tt = trial_traces("10-1074",["10-1074"], samples_per_patient = 100, diversity_multiplier=2.07)

##
#----------------------------------------------------------------#
# Make some plots about mutations
#----------------------------------------------------------------#
##
fid = h5open("/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialsimulations_traces.h5", "w")
fid["10-1074"] = tt
close(fid)
##

