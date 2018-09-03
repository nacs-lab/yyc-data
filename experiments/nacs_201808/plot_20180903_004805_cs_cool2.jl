#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180903_004805.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))

const spec_cs_a = OrderedDict(
    :az=>(-21:2.0:9, 21:2.0:51),
    :rx=>(-156:6.0:-66, 45:10.0:195),
    :ry=>(-158:6.0:-68, 50:8.0:170),
    :az_rabi=>linspace(0, 800, 21),
    :rx_rabi=>linspace(0, 400, 21),
    :ry_rabi=>linspace(0, 400, 21),
)

const split_cs_a = NaCsData.split_data(data_cs_a, spec_cs_a)

data_cs_az = split_cs_a[:az]
data_cs_rx = split_cs_a[:rx]
data_cs_ry = split_cs_a[:ry]
data_cs_az_rabi = split_cs_a[:az_rabi]
data_cs_rx_rabi = split_cs_a[:rx_rabi]
data_cs_ry_rabi = split_cs_a[:ry_rabi]

const prefix = joinpath(@__DIR__, "imgs", "data_20180903_004805_cs_cool2")

figure()
NaCsPlot.plot_survival_data(data_cs_az[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_az[2], fmt="C0.-")
grid()
ylim([0, 1])
title("Cs Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az")

figure()
NaCsPlot.plot_survival_data(data_cs_rx[1], fmt="C0.-", label="X")
NaCsPlot.plot_survival_data(data_cs_rx[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_ry[1], fmt="C1.-", label="Y")
NaCsPlot.plot_survival_data(data_cs_ry[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Cs Radial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_r")

figure()
NaCsPlot.plot_survival_data(data_cs_az_rabi, fmt="C0.-")
grid()
ylim([0, 1])
title("Cs Axial Heating")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az_rabi")

figure()
NaCsPlot.plot_survival_data(data_cs_rx_rabi, fmt="C0.-", label="X")
NaCsPlot.plot_survival_data(data_cs_ry_rabi, fmt="C1.-", label="Y")
grid()
legend()
ylim([0, 1])
title("Cs Radial Heating")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_r_rabi")

NaCsPlot.maybe_show()
