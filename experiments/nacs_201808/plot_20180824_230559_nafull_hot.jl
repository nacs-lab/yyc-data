#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180824_230559.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))

const spec_a = OrderedDict(
    :az=>-140:4.0:500,
    :rx=>-500:10.0:1100,
    :ry=>-500:10.0:1100,
)

const split_na_a = NaCsData.split_data(data_na_a, spec_a)

data_na_az = split_na_a[:az]
data_na_rx = split_na_a[:rx]
data_na_ry = split_na_a[:ry]

const prefix = joinpath(@__DIR__, "imgs", "data_20180824_230559_nafull_hot")

figure()
NaCsPlot.plot_survival_data(data_na_az, fmt="C0.-")
grid()
ylim([0, 0.65])
title("Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az")

figure()
NaCsPlot.plot_survival_data(data_na_rx, fmt="C0.-")
grid()
ylim([0, 0.65])
title("Radial X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

figure()
NaCsPlot.plot_survival_data(data_na_ry, fmt="C0.-")
grid()
ylim([0, 0.65])
title("Radial Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ry")

NaCsPlot.maybe_show()
