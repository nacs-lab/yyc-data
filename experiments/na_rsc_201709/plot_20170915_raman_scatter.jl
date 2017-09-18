#!/usr/bin/julia

include("load_20170915_raman_scatter.jl")

const prefix = joinpath(@__DIR__, "imgs", "data_20170915_rabi_scatter")

using NaCsPlot
using PyPlot

figure()
NaCsPlot.plot_survival_data(f1_down, fmt="o-")
grid()
ylim([0, ylim()[2]])
title("F1 down")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f1_down")

figure()
NaCsPlot.plot_survival_data(f1_up, fmt="o-")
grid()
ylim([0, ylim()[2]])
title("F1 up")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f1_up")

figure()
NaCsPlot.plot_survival_data(f1_diagonal, fmt="o-")
grid()
ylim([0, ylim()[2]])
title("F1 diagonal")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f1_diagonal")

figure()
NaCsPlot.plot_survival_data(f2_counterop, fmt="o-")
grid()
ylim([0, ylim()[2]])
title("F2 counterop")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f2_counterop")

figure()
NaCsPlot.plot_survival_data(f2_diagonal, fmt="o-")
grid()
ylim([0, ylim()[2]])
title("F2 diagonal")
xlabel("Time (ms)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_f2_diagonal")

NaCsPlot.maybe_show()
