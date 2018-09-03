#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180902_175939.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))

const spec_na_a = OrderedDict(
    :az=>(-139:4.0:-79, 34:4.0:94),
    :rx=>(-376:6.0:-286, 393:12.0:573),
    :ry=>(-385:6.0:-295, 402:12.0:582)
)

const spec_cs_a = OrderedDict(
    :az=>(-21:2.0:9, 21:2.0:51),
    :rx=>(-156:6.0:-66, 45:10.0:195),
    :ry=>(-158:6.0:-68, 50:8.0:170)
)

const split_na_a = NaCsData.split_data(data_na_a, spec_na_a)
const split_cs_a = NaCsData.split_data(data_cs_a, spec_cs_a)

data_na_az = split_na_a[:az]
data_na_rx = split_na_a[:rx]
data_na_ry = split_na_a[:ry]
data_cs_az = split_cs_a[:az]
data_cs_rx = split_cs_a[:rx]
data_cs_ry = split_cs_a[:ry]

const prefix = joinpath(@__DIR__, "imgs", "data_20180902_175939_cs_cool")

figure()
NaCsPlot.plot_survival_data(data_na_az[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_az[2], fmt="C0.-")
grid()
ylim([0, 1])
title("Na Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
NaCsPlot.plot_survival_data(data_cs_az[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_az[2], fmt="C0.-")
grid()
ylim([0, 1])
title("Cs Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
NaCsPlot.plot_survival_data(data_na_rx[1], fmt="C0.-", label="X")
NaCsPlot.plot_survival_data(data_na_rx[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_ry[1], fmt="C1.-", label="Y")
NaCsPlot.plot_survival_data(data_na_ry[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Na Radial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_r")

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
NaCsPlot.maybe_save("$(prefix)_cs_r")

NaCsPlot.maybe_show()
