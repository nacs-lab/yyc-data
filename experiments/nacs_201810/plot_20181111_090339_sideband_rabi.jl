#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181111_090339.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))

const spec_na_a = OrderedDict(
    :az=>linspace(0, 180, 21),
    :rx=>linspace(0, 80, 21),
    :ry=>linspace(0, 80, 21)
)

const spec_cs_a = OrderedDict(
    :az=>linspace(0, 800, 21),
    :rx=>linspace(0, 400, 21),
    :ry=>linspace(0, 400, 21)
)

const split_na_a = NaCsData.split_data(data_na_a, spec_na_a)
const split_cs_a = NaCsData.split_data(data_cs_a, spec_cs_a)

data_na_az = split_na_a[:az]
data_na_rx = split_na_a[:rx]
data_na_ry = split_na_a[:ry]
data_cs_az = split_cs_a[:az]
data_cs_rx = split_cs_a[:rx]
data_cs_ry = split_cs_a[:ry]

const prefix = joinpath(@__DIR__, "imgs", "data_20181111_090339_sideband_rabi")

figure()
NaCsPlot.plot_survival_data(data_na_rx, fmt="C0.-", label="X")
NaCsPlot.plot_survival_data(data_na_ry, fmt="C1.-", label="Y")
NaCsPlot.plot_survival_data(data_na_az, fmt="C2.-", label="Z")
grid()
legend()
ylim([0, 1])
title("Na heating Rabi flopping")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

figure()
NaCsPlot.plot_survival_data(data_cs_rx, fmt="C0.-", label="X")
NaCsPlot.plot_survival_data(data_cs_ry, fmt="C1.-", label="Y")
NaCsPlot.plot_survival_data(data_cs_az, fmt="C2.-", label="Z")
grid()
legend()
ylim([0, 1])
title("Cs heating Rabi flopping")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
