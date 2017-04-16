#!/usr/bin/julia

include("script.jl")

using PyPlot

function plot_f12(rates, tscale)
    init = Float64[1, 0, 0, 0, 0, 0, 0, 0]
    init_f1 = Float64[0, 0, 0, 0, 0, 1, 0, 0]
    r = sum(rates * init)
    r_f1 = sum(rates * init_f1)
    τ = 1 / r
    τ_f1 = 1 / r_f1
    τ_max = max(τ, τ_f1)
    ts = linspace(0, τ_max * tscale, 100)
    f1 = gen_model(rates, init).(ts)
    A = rates_to_A(rates)
    full_f = t->expm(A * t) * init
    full_phi = full_f.(ts)
    plot(ts, f1)
    plot(ts, full_phi)
    ylim(0, 1)
    maxt = maximum(ts)
    xlim(0, maxt)
    axvline(τ_f1, color="blue")
    axvline(τ, color="red")
end

figure()
title("F2 Counter OP")
plot_f12(rates_f2_counterop, 100)
figure()
title("F2 Co-prop")
plot_f12(rates_f2_coprop, 100)
figure()
title("F1 Up/Down")
plot_f12(rates_f1_up, 100)
figure()
title("F1 Co-prop")
plot_f12(rates_f1_coprop, 100)
show()
