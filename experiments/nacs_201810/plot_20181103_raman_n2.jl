#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181103_165031.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20181104_020816.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
# const iname_c = joinpath(@__DIR__, "data", "data_20181104_134413.mat")
# const params_c, logicals_c = NaCsData.load_striped_mat(iname_c)

data_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2), (3, 4)))
data_b = NaCsData.select_count(params_b, logicals_b, NaCsData.select_single((1, 2), (3, 4)))
# data_c = NaCsData.select_count(params_c, logicals_c, NaCsData.select_single((1, 2), (3, 4)))

const spec_a = (81:0.06:90) .* 2 .- 60
const spec_b = (79:0.06:88) .* 2 .- 60
const spec_c = ((81.9:0.06:85.5) .* 2 .- 60,
                (82.1:0.06:86) .* 2 .- 60,
                (83:0.06:89) .* 2 .- 60)

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)
# const split_c = NaCsData.split_data(data_c, spec_c)

data_100 = split_a
data_60 = split_b
# data_73 = split_c[1]
# data_86 = split_c[2]
# data_120 = split_c[3]

const prefix = joinpath(@__DIR__, "imgs", "data_20181103_raman_n2")

figure()
NaCsPlot.plot_survival_data(data_60, fmt="C0.-", label="6mW")
NaCsPlot.plot_survival_data(data_100, fmt="C3.-", label="10mW")
grid()
legend()
ylim([0, 1])
title("N=2 Raman")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_6_10")

# figure()
# NaCsPlot.plot_survival_data(data_60, fmt="C0.-", label="6mW")
# NaCsPlot.plot_survival_data(data_73, fmt="C1.-", label="7.3mW")
# NaCsPlot.plot_survival_data(data_86, fmt="C2.-", label="8.6mW")
# NaCsPlot.plot_survival_data(data_100, fmt="C3.-", label="10mW")
# NaCsPlot.plot_survival_data(data_120, fmt="C4.-", label="12mW")
# grid()
# legend()
# ylim([0, 1])
# title("N=2 Raman")
# xlabel("Detuning (kHz)")
# ylabel("Survival")
# NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
