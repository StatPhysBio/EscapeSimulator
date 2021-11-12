export WrightSampler

#=
This document implements a WrightSampler mutable struct type.
This object is used to generate Gibbs-sampler output on Wright equilibria
Also included are methods of the current state to generate stochastic gradients
And fisher information estimates
 =#


# This lists the neccesary parameters 
# for sampling from the Wright-equilibrium distribution. The partition function is the most costly to compute and so it cahched
# θ_mut, θ_wt, σ are forward and backward diversity parameters respectively. At the moment, there is only code for  σ>0.
mutable struct WrightSampler
    θ_wt::Float64
    θ_mut::Float64
    σ::Float64
    k::Int64   # auxiliary variable for gibbs sampling
    x::Float64 # The freqeuncy of wt
end


function WrightSampler(θ_wt, θ_mut, σ)
#    z = gamma(θ_wt) * HypergeometricFunctions.pFq([θ_wt],[θ_mut + θ_wt],σ) /gamma(θ_wt +θ_mut)
    ws = WrightSampler(θ_wt, θ_mut, σ, 0, .5)
    for ii in 1:10
    incrementk!(ws)
    incrementx!(ws)
    end
    return ws
end

# This sampler is robuset for small-exponent gamma distributions
# Code was copied from somewhere unknown (citation needed)
function randlogGamma(a;scale = 1)
    if a > .2
        return log.(rand(Gamma(a)))
    else
        L = 1 / a - 1
        w = a / MathConstants.e / (1 - a)
        ww = 1 / (1 + w)
        η = z -> ( z >= 0 ? exp(-z) : w * L * exp(L * z))
        h = z -> exp(-z - exp(-z / a))
        function rh(a)
            while true
                U = rand()
                z = (U <= ww) ? -log(U / ww) : log(rand()) / L
                h(z) / η(z) > rand() && return(z)
            end
        end
        return log(scale) - rh(a) / a
    end
end

function sigmoid(x) # numerically stable sigmoid
    if x >= 0
        1 / ( 1 + exp(-x) )
    else
        exp(x) / ( 1 + exp(x) )
    end
end

function incrementk!(ws::WrightSampler)
    if ws.σ >= 0 
        ws.k = rand(Poisson(abs(ws.σ * ws.x)))
    else
        ws.k = rand(Poisson(abs(ws.σ * (1-ws.x))))
    end
end

function incrementx!(ws::WrightSampler)
    if ws.σ >= 0 
        a = randlogGamma(ws.θ_wt + ws.k;scale = 1)
        b = randlogGamma(ws.θ_mut;scale = 1)
    else
        a = randlogGamma(ws.θ_wt;scale = 1)
        b = randlogGamma(ws.θ_mut + ws.k;scale = 1)
    end
    ws.x = sigmoid(a-b)
end

function sample!(ws::WrightSampler)
    for ii in 1:5
    incrementk!(ws)
    incrementx!(ws)
    end
    return ws.x
end


function rejection_sample_k(θ_wt, θ_mut,σ)
    k = 0 
    while true
        k = rand(Poisson(σ))
        k == 0 ? break : nothing
        p = exp( loggamma(θ_wt + θ_mut) + loggamma(k + θ_wt) - loggamma(k + θ_wt + θ_mut) -loggamma(θ_wt))
        rand(Bernoulli(p)) ? break : nothing
    end
    return k
end