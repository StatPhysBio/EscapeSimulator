# defines the functions which actually run simulations

export rebound_time, # outputs a single rebound time
    viral_rebound_times, # rebound vector as a funciton of a single viral population
    trial_rebound_times, # matrix of rebound times
    trial_rebound_prob # function of a theta vector 


function wt_esc_counts(vp)
	esc::Int = 0
	for virus in vp.pop
		if virus.escape
			esc += 1
		end
	end
	(esc, length(vp.pop) - esc)
end

function virus_time_trace(vp::ViralPopulation, time_points; breakpoint = 5.0)
	# breakpoint sets a threshold at which to end the simulation and save on time
	# because escape has already been acieved
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
	# breakpoint is a keyword that determines when a simulation should be treated deterministically
	# it is set very conservatively here. The main purpose is to prevent running simulations 
	# at full capacity all the way to the endpoint time (56 days)
	(esc, _) = wt_esc_counts(vp)
	trace = virus_time_trace(vp, 0:56; breakpoint)
	if esc/vp.capacity > breakpoint
		t_rebound = 0.0
	elseif trace[end][1]/vp.capacity < breakpoint
		t_rebound = 57
	else
		ind = findfirst(col->col[1]/vp.capacity > breakpoint, trace)
		t_cross = (0:56)[ind]
		q_cross = trace[ind][1]/vp.capacity
		t_rebound = 3*log(1 + exp(t_cross / 3) * (1-q_cross)/q_cross)
	end 
	return (t_rebound,trace)
end

function viral_rebound_times(vp::ViralPopulation;
    breakpoint = 0.8,
    n_samples = 10)
	# ensemble simulation with reused viral population structs
    rebounds = Float64[]
    for _ in 1:n_samples
        (t_rebound, _) = rebound_time(vp::ViralPopulation; breakpoint)
        push!(rebounds, t_rebound)
        restart_viral_population!(vp)
    end
    return rebounds
end


function trial_rebound_times(theta_vec, ab_profile_list; 
    n_samples = 1000,
    breakpoint = 0.8,
    mutations = true,
    mut_per_gen = 3*10^(-6), # the mutation rate measured in growth rate 10^-6
    decayrate = 0.4,
    f = 1.0/3)
    times = zeros(length(theta_vec),n_samples)
	Threads.@threads for (ii,θ) in collect(enumerate(theta_vec))
        vp = initialize_viral_population(θ, ab_profile_list;
		    mut_per_gen, # the mutation rate measured in growth rate 10^-6
		    decayrate, f) 
        if !mutations
            vp.mutation_operator = x -> nothing
        end
        times[ii,:] .= viral_rebound_times(vp; n_samples,  breakpoint)
    end
    return times
end

function trial_rebound_prob(theta_vec, ab_profile_list; 
    n_samples = 1000,
    breakpoint = 0.8,
    mutations = true,
    mut_per_gen = 3*10^(-6), # the mutation rate measured in growth rate 10^-6
    decayrate = 0.4,
    f = 1.0/3)
    mean(vec(trial_rebound_times(theta_vec, ab_profile_list; 
        n_samples, breakpoint, mutations, mut_per_gen, decayrate, f)) .< 56.0)
end