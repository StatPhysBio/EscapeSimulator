export mut_frequency_vec, time_sample, ensemble_evolution

function bitaccumulate!(vec::Vector,b::Int64)
    for i in 1:bit_length(b)
        if readbit(b, i) == 1
            vec[i] += 1 
        end
    end
end

function mut_frequency_vec(vp)
    n_gene = length(vp.equilibrium_sampler)
    vec = zeros(Int64,n_gene);
    for virus in vp.pop
        bitaccumulate!(vec,virus.genotype)
    end
    return vec ./ length(vp.pop)
end

function time_sample(vp, time_points)
	reduce(hcat,
        begin
        run_viruses_stop(vp, time)
        mut_frequency_vec(vp) 
        end
        for time in time_points)
end


function ensemble_evolution(θ, ab_profile, time_points; n_ensemble = 100, antibody = false, kwds...)
    reduce((x,y)->cat(x,y;dims = 3),
        time_sample(initialize_viral_population(θ, ab_profile; antibody, kwds...), time_points) 
        for _ in 1:n_ensemble)
end
