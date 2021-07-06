using EscapeSimulator
using HDF5
using Statistics
using BenchmarkTools
using StatsBase
using Random

mainfolder = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/"
snpanalysis = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/snpanalysis.h5"
trialpath(x) = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/trialanalysis$x.h5"
trial_list = ["10-1074","3BNC","combo"]

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
    # return the mle parameters for the antibody
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
