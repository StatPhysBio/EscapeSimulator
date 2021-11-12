# defines the functions which actually run simulations

export rebound_time, # outputs a single rebound time
    viral_rebound_times, # rebound vector as a funciton of a single viral population
    trial_rebound_times, # matrix of rebound times
    trial_rebound_prob, # function of a theta vector 
	virus_time_trace # full trace for drawing rebound trajectories


function wt_esc_counts(vp)
	# note that this name is exactly backward
	# returns a tuple 
	esc::Int = 0
	for virus in vp.pop
		if virus.escape
			esc += 1
		end
	end
	(esc, length(vp.pop) - esc)
end

function virus_time_trace(vp::ViralPopulation, time_points; breakpoint = 5.0)
	# records the viral state (wt,esc) counts as a tuple
	# breakpoint sets a threshold at which to end the simulation and save on time
	# because escape has already been acieved
	# Default traces are set way above capacity, assuring that the trace has
	# the same length as time_points.
	out = Tuple{Int64, Int64}[]
	for time in time_points
		run_viruses_stop(vp, time)
		(esc, wt) = wt_esc_counts(vp)
		push!(out, (esc,wt))
		if esc/vp.capacity > breakpoint
			break
		end
	end
	return out
end

function rebound_time(vp::ViralPopulation; breakpoint = 0.8)
	# function to get a single rebound time and trace
	# breakpoint is a keyword that determines when a simulation should be treated deterministically
	# it is set very conservatively here. The main purpose is to prevent running simulations 
	# at full capacity all the way to the endpoint time (56 days)
	# returns a tuple of rebound time and trace if required.
	(esc, _) = wt_esc_counts(vp)
	trace = virus_time_trace(vp, 0:56; breakpoint)
	if esc/vp.capacity > breakpoint
		t_rebound = 0.0
	elseif trace[end][1]/vp.capacity < breakpoint
		t_rebound = 57
	else
		ind = findfirst(col->col[1]/vp.capacity > breakpoint, trace)
		t_cross = (0:56)[ind]
		q_cross = min(trace[ind][1]/vp.capacity,1)
		t_rebound = 3*log(1 + exp(t_cross / 3) * (1-q_cross)/q_cross)
	end 
	return (t_rebound,trace)
end

function viral_rebound_times(vp::ViralPopulation;
    breakpoint = 0.8,
    n_samples = 10)
	# Simulation with reused viral population structs
	# Single threaded loop used to rerun a population multiple times
    rebounds = Float64[]
    for _ in 1:n_samples
        (t_rebound, _) = rebound_time(vp::ViralPopulation; breakpoint)
        push!(rebounds, t_rebound)
        restart_viral_population!(vp)
    end
    return rebounds
end


function trial_rebound_times(theta_vec, ab_profile_list; 
	# ab_profile_list is a list of lists of antibody mutation and selection properties
	# the lengths of the sublists determines the escape function in the dynamics.
    n_samples = 100, # number of samples per bayesian realization * n_runs
	n_runs = 10, # parameter for better thread efficiency
    breakpoint = 0.8, pars...)
    times = zeros(length(theta_vec)*n_runs, n_samples) # float matrix of rebound times 
	Threads.@threads for (ii, θ) in collect(enumerate(repeat(theta_vec, n_runs)))
        vp = initialize_viral_population(θ, ab_profile_list; pars...) 
        times[ii,:] .= viral_rebound_times(vp; n_samples,  breakpoint)
    end
    return reshape(times,(length(theta_vec), n_runs*n_samples))
end

function trial_rebound_prob(theta_vec, ab_profile_list; 
    breakpoint = 0.1, kwds...)
    mean(vec(trial_rebound_times(theta_vec, ab_profile_list; breakpoint, kwds...)) .< 56.0)
end