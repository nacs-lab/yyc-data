#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180812_161422.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_cs_a = NaCsData.select_count(params_a, logicals_a,
                                  NaCsData.select_single((-1, 2), (4,)), 178000)
data_nacs_a = NaCsData.select_count(params_a, logicals_a,
                                    NaCsData.select_single((1, 2), (4,)), 178000)

const spec = OrderedDict(
    :coprop=>(7 .+ linspace(-30, 30, 61)),
    :rx=>(-200.0:100.0),
)

const split_cs_a = NaCsData.split_data(data_cs_a, spec)
const split_nacs_a = NaCsData.split_data(data_nacs_a, spec)

const prefix = joinpath(@__DIR__, "imgs", "data_20180812_161422_interaction_shift3_rx")

data_cs_coprop = split_cs_a[:coprop]
data_cs_rx = split_cs_a[:rx]
data_nacs_coprop = split_nacs_a[:coprop]
data_nacs_rx = split_nacs_a[:rx]

figure()
NaCsPlot.plot_survival_data(data_cs_coprop, fmt="C0.-", label="Cs only")
NaCsPlot.plot_survival_data(data_nacs_coprop, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("Co-prop interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_coprop")

figure()
NaCsPlot.plot_survival_data(data_cs_rx, fmt="C0.-", label="Cs only")
NaCsPlot.plot_survival_data(data_nacs_rx, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("Radial X interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

figure()
NaCsPlot.plot_survival_data(data_cs_coprop, fmt="C0.-", label="Co-prop")
NaCsPlot.plot_survival_data(data_cs_rx, fmt="C1.-", label="Radial X")
grid()
legend()
ylim([0, 0.8])
title("Cs only")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

figure()
NaCsPlot.plot_survival_data(data_nacs_coprop, fmt="C0.-", label="Co-prop")
NaCsPlot.plot_survival_data(data_nacs_rx, fmt="C1.-", label="Radial X")
grid()
legend()
ylim([0, 0.8])
title("Na + Cs")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_nacs")

NaCsPlot.maybe_show()
