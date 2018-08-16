#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180811_183359.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180811_201412.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)

data_na_a = NaCsData.select_count(params_a, logicals_a,
                                  NaCsData.select_single((1, -2), (3,)))
data_nacs_a = NaCsData.select_count(params_a, logicals_a,
                                    NaCsData.select_single((1, 2), (3,)))
data_na_b = NaCsData.select_count(params_b, logicals_b,
                                  NaCsData.select_single((1, -2), (3,)))
data_nacs_b = NaCsData.select_count(params_b, logicals_b,
                                    NaCsData.select_single((1, 2), (3,)))

const spec = OrderedDict(
    :coprop=>-70.0:-11.0,
    :rx=>-150.0:-0.0,
)

const split_na_a = NaCsData.split_data(data_na_a, spec)
const split_nacs_a = NaCsData.split_data(data_nacs_a, spec)
const split_na_b = NaCsData.split_data(data_na_b, spec)
const split_nacs_b = NaCsData.split_data(data_nacs_b, spec)

const prefix = joinpath(@__DIR__, "imgs", "data_20180811_interaction_shift_rx")

data_na_coprop = [split_na_a[:coprop]; split_na_b[:coprop]]
data_na_rx = [split_na_a[:rx]; split_na_b[:rx]]
data_nacs_coprop = [split_nacs_a[:coprop]; split_nacs_b[:coprop]]
data_nacs_rx = [split_nacs_a[:rx]; split_nacs_b[:rx]]

figure()
NaCsPlot.plot_survival_data(data_na_coprop, fmt="C0.-", label="Na only")
NaCsPlot.plot_survival_data(data_nacs_coprop, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("Co-prop interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_coprop")

figure()
NaCsPlot.plot_survival_data(data_na_rx, fmt="C0.-", label="Na only")
NaCsPlot.plot_survival_data(data_nacs_rx, fmt="C1.-", label="Na + Cs")
grid()
legend()
ylim([0, 0.8])
title("Radial X interaction shift")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

figure()
NaCsPlot.plot_survival_data(data_na_coprop, fmt="C0.-", label="Co-prop")
NaCsPlot.plot_survival_data(data_na_rx, fmt="C1.-", label="Radial X")
grid()
legend()
ylim([0, 0.8])
title("Na only")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

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
