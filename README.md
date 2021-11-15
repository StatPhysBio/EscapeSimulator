# EscapeSimulator

EscapeSimulator is a pop-gen simulator for Julia. This is a derivative of the Tomoko.jl package but with differences:
	1. Genotypes are represented as unsigned integers (`UInt64`) instead of `BitVectors` which have unlimited length. This means that genotypes have a maximum length of 64, thus the maximal number of escape sites in the simulation is 64. The benefit is that UInt64 is a static memory object which leads to significant speed-ups.
	2. The mutation model is much more flexible in that mutations can occur at different forward and backward rates µ ≠ µ_dagger, and they can differ based on the site of interest.
	3. Other functions related to the data analysis and discrepancy are also included

# Basic Operations

## Fitting fitness ratios
See the file: bayes_posterior.jl in HIVTreatmentOptimization/JuliaScripts for example implementations. At each time-point for each patient in the longitudinal-deep-sequencing patients we construct an object, a `LikelihoodSample`, which contains the importance-sampling estimator for the likelihood at different values for the observable fitness difference ratio `rs = σ/θ = (f_wt - f_mut)/μ`.

Afterward we can combine likelihood estimates using different averaging procedures 

```julia
posterior_sample = baysian_pop_rs(listof_listof_LikelihoodSamples; 
	samples = 2*10^3, burn_in = 10^2, avg = avg)
```

`avg` is a keyword between 0 and 1 that toggles between independent `avg = 0` and patient averaged `avg = 1.0` samples.

## Running Simulations
To run simulations and calculate and estimate a rebound time, you must first create a viral popualtion with `initialize_viral_population` which has the following signiture
```julia
vp = initialize_viral_population(θ, ab_profile_list;
		# θ is the diversity of the viral population at treatment initiation, and sets the effective carrying capacity via θ = 2 N_e μ
		# ab_profile_list consists of a list of ab_profiles. Each profile is a list of site profiles, which itself is a list cosisting of 
		# [forward mutation rate, backward mutation rate, fitness cost]
		mut_per_gen = 3*1.1*10^(-5), # the base (avg'd) transition rate measured in growth rate
		decayrate = 0.31, # the rate of removal 
        λ = 2.0, # total noise rate of the populaiton
		f = 1.0/3, # absolute growth rate
        mutations = true, # whether or not to include post-treatment mutations
        antibody = true)
```
We then can evolve a population forward in time or get a rebound time using several functions with increasing sophistication

	- `virus_time_trace` generates a time trace from 0-56 days recording mutant and wildtype fractions
	- `rebound_time` generates a truncated trace up until `breakpoint = 0.8`, a threshold for the mutant fraction to acieve afterwhich the dynamics are deterministic. It returns `(reboundtime,trace)`
	- `viral_rebound_times` efficiently samples multiple rebound times with the same parameters by restarting the viral popuplation after every breakpoint
	- `trial_rebound_times` multi-threaded version which simulates the rebound times of a trial with a constant antibody profile, but which includes diversity variation.


## Calculating Discrepancies

This includes one all purpose function which calculates the discrepancy between the reference data and the simulated data.
```julia
	observed_discrepancy(reference_data, simulated_data; min = 0.001, max = 56)`
```
The keys `min` and `max` indicate truncation points for the data categories (in this case the minimum and maximum rebound times)