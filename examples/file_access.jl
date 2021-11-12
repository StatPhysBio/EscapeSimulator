using EscapeSimulator
using HDF5
using Statistics
using BenchmarkTools
using StatsBase
using Random

# Basic parameters of the simulation and their default values.
# mut_per_gen = 3*1.1*10^(-5), # the base (avg'd) transition rate measured in growth rate
# decayrate = 0.31,
# λ = 2.0,
# f = 1.0/3,
# mutations = true

## This folder sets up functions that read in data from .h5 files.
# it allows us to read in the neccesary data for generating a viral population from the bayesian distribtuions

# These is the main function for generating populations
# It returns a list of lists of site profiles.
# This is the main parameter for the initialization of a viral population
# vp = initialize_viral_population(θ::Float64, ab_profile_list; kwds)
ab_profile_list_bayes(ablist; selection_multiplier = 1.0) = map(ab->get_antibody_profile_bayes(ab; selection_multiplier), ablist) 
ab_profile_list_mle(ablist; selection_multiplier = 1.0) = map(ab->get_antibody_profile_mle(ab; selection_multiplier), ablist) 


function simulate_trial_mle(trial, ab_list; 
	n_samples = 10, # *10 from n_runs
	selection_multiplier = 1.0, diversity_multiplier = 1.0, 
    kwds...)
	theta_vec = get_start_theta(trial) .* diversity_multiplier
	ab_profile_list = ab_profile_list_mle(ab_list; selection_multiplier)
	trial_rebound_times(theta_vec, ab_profile_list; n_samples, kwds...)
end

function repeat_minimal_noise(vec, new_length)
    # a function to make the length of a vector longer without doing multinomial draws
    reps = fld(new_length,length(vec));
    vcat(repeat(vec,reps), sample(vec, mod(new_length,length(vec)), replace=false))
end

function simulate_trial_bayes(trial, ab_list; # trial determines θ distribution
    # trial one of ["10-1074","3BNC","combo","all"]
	n_samples = 5, # *10 from n_runs = number of simulations for a given parameter choice
    n_patients = length(get_observed_rebound_times(trial)),
    n_pars = 100, # number of realizations of parameter values
	selection_multiplier = 1.0, diversity_multiplier = 1.0, 
	kwds...)
	theta_vec = 
        repeat_minimal_noise(get_start_theta(trial) .* diversity_multiplier, n_patients) 
        # The purpose of this auxiliary function is to make sure that the theta distribuiton matches the number of patients
    reduce((x,y) -> cat(x,y; dims = 3), 
    trial_rebound_times(theta_vec, ab_profile_list_bayes(ab_list; selection_multiplier); n_samples, kwds...) 
        for ii in 1:n_pars) # generates a n_patients x 10 n_samples x n_pars array
end

function simulate_trial_mle_prob(trial, ab_list; 
	n_samples = 100, # length(theta_vec)*10 from n_runs multithreading
	selection_multiplier = 1.0,
	diversity_multiplier = 1.0, kwds...
	)
	theta_vec = get_start_theta(trial) .* diversity_multiplier
    ab_profile_list = ab_profile_list_mle(ab_list; selection_multiplier)
    trial_rebound_prob(theta_vec, ab_profile_list; n_samples, kwds...)
end



# get_start_theta can be used to access the theta information from each of the three trials
# in trial list ["10-1074","3BNC","combo"]
# "all" uses the pooled diversity estimates.

mainfolder = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/"
snpanalysis = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/snpanalysis.h5"
trialpath(x) = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialanalysis$x.h5"
trial_list = ["10-1074","3BNC","combo"]

function trial_antibodies(trial)
    Dict(
        "10-1074" => ["10-1074"], 
        "3BNC" => ["3BNC117"],
        "combo" => ["10-1074","3BNC117"]
    )[trial]
end

function get_start_theta(trial)
	if trial != "all"
    h5trial = h5open(trialpath(trial), "r")
    thetas = Float64[]
    for patient in h5trial["genetic"]
        for day in patient
            if (read(day["day"]) < 2) & haskey(day, "theta")
                push!(thetas,read(day["theta"]))
            end
        end
    end
    close(h5trial)
    return thetas
	else 
	return mapreduce(get_start_theta, vcat, trial_list)
	end
end


h5open(snpanalysis, "r") do fid
	global ablist = keys(fid)
	ablist = [ab for ab in ablist if (&)(ab != "101074", ab != "VRC01_dms")]
end 


function get_start_theta(trial)
	if trial != "all"
    h5trial = h5open(trialpath(trial), "r")
    thetas = Float64[]
    for patient in h5trial["genetic"]
        for day in patient
            if (read(day["day"]) < 2) & haskey(day, "theta")
                push!(thetas,read(day["theta"]))
            end
        end
    end
    close(h5trial)
    return thetas
	else 
	return mapreduce(get_start_theta, vcat, trial_list)
	end
end


function get_antibody_profile_mle(ab; selection_multiplier = 1.0)
    # Takes the antibody name as a string
    # return the mle parameters for the antibody
    sites = Array{Float64,1}[]
    h5open(snpanalysis, "r") do fid
            for site in fid[ab]
                rs = read(site["rsel"])[1] * selection_multiplier
                mut = reverse(read(site["mut"])) # [backmut, forwardmut]
                if rs > 0
                    push!(sites,vcat(mut...,rs))
                end
            end
    end
    return sites
end

function get_antibody_profile_bayes(ab; selection_multiplier = 1.0)
    # return bayes sampled parameters for the antibody
    sites = Array{Float64,1}[]
    h5open(snpanalysis, "r") do fid
            for site in fid[ab]
                rs = rand(read(site["bayes_avg1.0"])) * selection_multiplier
                mut = reverse(read(site["mut"]))  # [backmut, forwardmut]
                if rs > 0
                    push!(sites,vcat(mut...,rs))
                end
            end
    end
    return sites
end

function get_observed_rebound_times(trial)
    # return the mle parameters for the antibody
	t_rebounds = Float64[]
    h5trial = h5open(trialpath(trial), "r")
	growth_rate = read(h5trial["viraemic/growth"])
    for patient_key in keys(h5trial["viraemic/patients"])
		patient = h5trial["viraemic/patients/$patient_key"]
		if (|)(trial != "3BNC", !(in(["2A1","2A3","2A4"])(patient_key))) 
			push!(t_rebounds,-log(read(patient["x"])) / growth_rate)
		end
    end
    close(h5trial)
    return t_rebounds
end
