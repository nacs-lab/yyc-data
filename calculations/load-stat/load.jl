#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import Distributions: Binomial
import NaCsCalc.Utils: binomial_estimate
import NaCsPlot
using PyPlot

single_sample(n, p) = binomial_estimate(rand(Binomial(n, p)), n)[1]

function collect_samples(n, p, N)
    r = 0.0
    r2 = 0.0
    for i in 1:N
        s = single_sample(n, p)
        r += s
        r2 += s^2
    end
    avg = r / N
    avg2 = r2 / N
    return avg, sqrt(avg2 - avg^2)
end

function plot_samples(p, ns, N=100000; kw...)
    avgs = Vector{Float64}(length(ns))
    stds = Vector{Float64}(length(ns))
    for i in 1:length(ns)
        avgs[i], stds[i] = collect_samples(ns[i], p, N)
    end
    errorbar(ns, avgs, stds; kw...)
end

plot_samples(0.01, 10:10:200, label="1.0%")
plot_samples(0.015, (10:10:200) .+ 2, label="1.5%")
plot_samples(0.02, (10:10:200) .+ 4, label="2.0%")
grid()
ylim([0, 0.1])
legend()
show()
