#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
# using DataStructures
# using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180108.csv")
const data = readdlm(iname_a, ',', Float64, skipstart=1)

const prefix = joinpath(@__DIR__, "imgs", "data_20180108")

figure()
errorbar(data[:, 1], data[:, 2], data[:, 3], label="First")
errorbar(data[:, 1], data[:, 6], data[:, 7], label="Second")
title("Signal")
grid()
xlabel("Phase (ns)")
ylabel("Count")
legend()
NaCsPlot.maybe_save("$(prefix)_signal")

figure()
errorbar(data[:, 1], data[:, 4], data[:, 5], label="First")
errorbar(data[:, 1], data[:, 8], data[:, 9], label="Second")
title("Load")
grid()
xlabel("Phase (ns)")
ylabel("Survival")
legend()
NaCsPlot.maybe_save("$(prefix)_load")

NaCsPlot.maybe_show()
