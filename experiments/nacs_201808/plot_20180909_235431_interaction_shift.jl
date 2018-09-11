#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180909_235431.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((-1, 2), (-3, 4)))
data_nacs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1, 2,), (4,)))

const spec_a = OrderedDict(
    :spectrum=>-40:2.0:80,
    :p1=>0:10.0:300,
    :n0=>0:10.0:300,
    :m1=>0:10.0:300,
)

const split_cs_a = NaCsData.split_data(data_cs_a, spec_a)
const split_nacs_a = NaCsData.split_data(data_nacs_a, spec_a)

data_cs_spectrum = split_cs_a[:spectrum]
data_cs_n0 = split_cs_a[:n0]
data_nacs_spectrum = split_nacs_a[:spectrum]
data_nacs_p1 = split_nacs_a[:p1]
data_nacs_n0 = split_nacs_a[:n0]
data_nacs_m1 = split_nacs_a[:m1]

const prefix = joinpath(@__DIR__, "imgs", "data_20180909_235431_interaction_shift")

figure()
NaCsPlot.plot_survival_data(data_cs_spectrum, fmt="C0.-", label="Cs only")
NaCsPlot.plot_survival_data(data_nacs_spectrum, fmt="C1.-", label="Na+Cs")
grid()
legend()
ylim([0, 0.65])
title("4+2 -> 3+2")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_spectrum")

figure()
NaCsPlot.plot_survival_data(data_cs_n0, fmt="C0.-", label="Cs only")
NaCsPlot.plot_survival_data(data_nacs_p1, fmt="C1.-", label="+ shift")
NaCsPlot.plot_survival_data(data_nacs_n0, fmt="C2.-", label="0 shift")
NaCsPlot.plot_survival_data(data_nacs_m1, fmt="C3.-", label="- shift")
grid()
legend()
ylim([0, 1])
title("4+2 -> 3+2 Rabi flopping")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_time")

NaCsPlot.maybe_show()
