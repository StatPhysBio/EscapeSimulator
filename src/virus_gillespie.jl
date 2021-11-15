# Defines the basic dynamics of the virus
# Depends on WrightSampler.jl
# Gets used by virus_trace.jl

export initialize_viral_population, restart_viral_population!

struct Virus
	genotype::Int64 # this is the genotype of bits.
	fitness::Float64 # the additive fitness of the virus
	escape::Bool # whether or not the virus is escaped
end

mutable struct ViralPopulation{F1,F2,F3} #function types must be explicitly seen by compiler
	time::Float64 # current time
    fitness_total::Float64 # current fitness. Divided by cap gives φ
    pop::Vector{Virus} # State of population length is the total population size
    capacity::Int 	# carrying capacity
    r::Float64 	# rate of decay of neutralized viruses
    λ::Float64 	# noise of system
	fitness_function::F1
	mutation_operator::F2
	escape_function::F3
    equilibrium_sampler::Vector{WrightSampler} # to generate new starting conditions for doing mulitple runs
end


function bitmask_toggle(b::Int64, vec0::Vector{Float64}, vec1::Vector{Float64})
    out = Vector{Float64}(undef,length(vec0))
    for i in 1:bit_length(b)
        if readbit(b, i) == 0
            out[i] = vec0[i]
        else 
            out[i] = vec1[i]
        end
    end
    return out
end


function bitmask_toggle_choose(b::Int64, vec0::Vector{Float64}, vec1::Vector{Float64}, tot::Float64)
    limit = rand() * tot
    @inbounds for i in 1:length(vec0)
        if readbit(b, i) == 0
            limit -= vec1[i] # first (0,1), testing (1,0)
        else
            limit -= vec0[i]
        end
        if limit < 0 
            return i # terminate early
        end
    end
    return 0
end


function mutation_operator_closure(ab_profile_list, ν)
    let mut0 = mapreduce(y -> map(x->x[1],y), vcat, ab_profile_list),
        mut1 = mapreduce(y -> map(x->x[2],y), vcat, ab_profile_list), 
		tot = sum(mut0) + sum(mut1)
		# this is a let-block to capture a fixed
		function mutate!(vp::ViralPopulation) # attempt to mutate a virus at index ind
            if rand() < tot*ν # if there is a shot at mutation
                ind = rand(1:length(vp.pop)) # choose a random virus
				@inbounds v = vp.pop[ind]
				#mutvec = mut[ (1:2:2*l) .+ bitarray(v.genotype,l)] # get the mutational parameters
                ii = bitmask_toggle_choose(v.genotype, mut0, mut1, tot)
				if ii > 0
					# flippy the dippy and make a new virus
					new_genotype = flip(v.genotype, 2^(ii-1))
					new_fitness = vp.fitness_function(new_genotype)
                    new_escape = vp.escape_function(new_genotype)
                    vp.fitness_total += new_fitness - v.fitness
                    vp.pop[ind] = Virus(new_genotype,new_fitness,new_escape)
				end
			end
        end
        return mutate!
	end
end

function escape_function_closure(ab_profile_list)
	# true or false on being escaped or not
	# determines binding in simulations
	ab_site_pattern = length.(ab_profile_list)
	let bitpatterns = Int64[packbits(vcat([
						    if (n == m) 
								ones(ab_site_pattern[m]) else
								zeros(ab_site_pattern[m]) 
							end
		for m in 1:length(ab_site_pattern)]...))
				for n in 1:length(ab_site_pattern)]
		# construct the bit masks for each antibody
		function is_escaped(genotype::Int64)
            mapreduce( x->anyone(genotype,x), (&), bitpatterns)
        end
        return is_escaped
	end
end

function bitmask_sum(b::Int64, vec::Vector{Float64})
    s = 0.0
    for i in 1:bit_length(b)
        if readbit(b, i) == 1
            s += vec[i]
        end
    end
    return s
end


function fitness_function_closure(ab_profile_list, f, μ)
	# absolute fitness and mutation rate
    let deltavec::Vector{Float64} =  mapreduce(y -> map(x->x[3],y), vcat, ab_profile_list) .* μ
		# this is a let-block to capture a fixed Δ-vector
		function fitness(b::Int64)
			f - bitmask_sum(b::Int64, deltavec)
        end
        return fitness
	end
end

function random_genotype(frequency_vector::Vector{Float64})
	x::Int64 = 0
	for ii in Iterators.reverse(frequency_vector) #start from the most significant bit
		x += rand() > ii
		x <<= 1 # bitshift to a larger power of 2
	end
	x >>=1 #undo the last bitshift
	return x
end

function equilibrium_sampler(θ, ab_profile_list)
	# absolute fitness and mutation rate
    mapreduce( y -> map(x-> WrightSampler((x .* θ)...), y), vcat, ab_profile_list)
end

function initialize_viral_population(θ, ab_profile_list;
    # θ is the diversity of the viral population at treatment initiation, and sets the effective carrying capacity via θ = 2 N_e μ
    # ab_profile_list consists of a list of ab_profiles. Each profile is a list of site profiles, which itself is a list cosisting of 
    # [forward mutation rate, backward mutation rate, fitness cost]
        mut_per_gen = 3*1.1*10^(-5), # the base (avg'd) transition rate measured in growth rate
        decayrate = 0.31, # the rate of removal 
        λ = 2.0, # total noise rate of the populaiton
        f = 1.0/3, # absolute growth rate
        mutations = true, # whether or not to include post-treatment mutations
        antibody = true)
	μ = mut_per_gen * f # absolute transition mutation rate in equilirbrium
	capacity = ceil(Int64, λ * θ / (2 * μ))
    ν = θ / capacity / 2 # per-birth-event mutation rate; ν = μ/λ = θ / 2 capacity
	if λ<f 
		error("pop_size too small to support growth rate") 
	end
    sampler_vec = equilibrium_sampler(θ, ab_profile_list)
	fitness_function = fitness_function_closure(ab_profile_list, f, μ)
	mutation_operator = mutations ? mutation_operator_closure(ab_profile_list, ν) : x -> nothing
	escape_function = antibody ? escape_function_closure(ab_profile_list) : x -> true
	pop = Virus[]
	fitness_total = 0.0
	vp = ViralPopulation(0.0, 
        fitness_total, 
        pop, capacity, decayrate, λ, fitness_function, mutation_operator, escape_function, sampler_vec)
    restart_viral_population!(vp)
    return vp
    return equilibrium_sampler
end

function restart_viral_population!(vp::ViralPopulation)
	# this looks like it is still based on the static struct form
	# really only need to update time and pop state
	# have to do some profiling so see if this is really a problem
	# fair amount of memory allocation,
    vp.pop = Vector{Virus}(undef,vp.capacity)
    start_freq = sample!.(vp.equilibrium_sampler)
    for ii in 1:vp.capacity
        genotype = random_genotype(start_freq)
        fitness = vp.fitness_function(genotype)
        escape = vp.escape_function(genotype)
        @inbounds vp.pop[ii] = Virus(genotype, fitness, escape)
    end
    vp.fitness_total = sum(virus.fitness for virus in vp.pop)
    vp.time = 0.0
end

function random_swap!(array::Vector) # bring an element of a vector to the end
    # My thinking is that this will be good for memory management
    ii = rand(1:length(array))
    @inbounds (array[end],array[ii]) = (array[ii],array[end]);
end


function gillespie_step!(vp::ViralPopulation)
    N = length(vp.pop)
    vp.time += randexp() / (N*vp.λ) # step the population forward a little bit
    if N > 0
        vp.mutation_operator(vp) # see if you can mutate something
        random_swap!(vp.pop)
        virus = pop!(vp.pop) # get the virus stack-style
        if virus.escape
            if rand() < ((vp.λ - vp.fitness_total/vp.capacity + virus.fitness) / vp.λ / 2 )
                push!(vp.pop,virus)
                push!(vp.pop,virus)
                vp.fitness_total += virus.fitness
            else
                vp.fitness_total -= virus.fitness
            end # else get rid of that virus, it is dead
        else # if the virus is neutralized
            if rand() >  vp.r /vp.λ # check and see if you are actually getting rid of it at the proper rate
                push!(vp.pop,virus) # put that virus back if you tried to kill it too quickly
            else
                vp.fitness_total -= virus.fitness
            end # else get rid of that shit, it is dead
        end
    end
end

function run_viruses_stop(vp::ViralPopulation, end_time)
	while vp.time < end_time
		gillespie_step!(vp)
	end
end
