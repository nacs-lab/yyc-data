#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181112_095349.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))

const spec_na_a = OrderedDict(
    :az=>(-126 .+ (-30:4.0:30), 76 .+ (-30:4.0:30)),
    :rx=>(-479 .+ (-45:6.0:45), 614 .+ (-90:12.0:90)),
    :ry=>(-490 .+ (-45:6.0:45), 631 .+ (-90:12.0:90))
)

const spec_cs_a = OrderedDict(
    :az=>(-20 .+ (-15:2.0:15), 26 .+ (-15:2.0:15)),
    :rx=>(-135 .+ (-45:6.0:45), 109 .+ (-75:10.0:75)),
    :ry=>(-129 .+ (-45:6.0:45), 112 .+ (-60:8.0:60))
)

const split_na_a = NaCsData.split_data(data_na_a, spec_na_a)
const split_cs_a = NaCsData.split_data(data_cs_a, spec_cs_a)

data_na_az = split_na_a[:az]
data_na_rx = split_na_a[:rx]
data_na_ry = split_na_a[:ry]
data_cs_az = split_cs_a[:az]
data_cs_rx = split_cs_a[:rx]
data_cs_ry = split_cs_a[:ry]

const prefix = joinpath(@__DIR__, "imgs", "data_20181112_095349_sideband_3d")

figure()
NaCsPlot.plot_survival_data(data_na_az[1], fmt="C1.-")
NaCsPlot.plot_survival_data(data_na_az[2], fmt="C1.-")
grid()
ylim([0, 1])
title("Na Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
NaCsPlot.plot_survival_data(data_cs_az[1], fmt="C1.-")
NaCsPlot.plot_survival_data(data_cs_az[2], fmt="C1.-")
grid()
ylim([0, 1])
title("Cs Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
NaCsPlot.plot_survival_data(data_na_rx[1], fmt="C1.-")
NaCsPlot.plot_survival_data(data_na_rx[2], fmt="C1.-")
grid()
ylim([0, 1])
title("Na X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_rx")

figure()
NaCsPlot.plot_survival_data(data_cs_rx[1], fmt="C1.-")
NaCsPlot.plot_survival_data(data_cs_rx[2], fmt="C1.-")
grid()
ylim([0, 1])
title("Cs X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_rx")

figure()
NaCsPlot.plot_survival_data(data_na_ry[1], fmt="C1.-")
NaCsPlot.plot_survival_data(data_na_ry[2], fmt="C1.-")
grid()
ylim([0, 1])
title("Na Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_ry")

figure()
NaCsPlot.plot_survival_data(data_cs_ry[1], fmt="C1.-")
NaCsPlot.plot_survival_data(data_cs_ry[2], fmt="C1.-")
grid()
ylim([0, 1])
title("Cs Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ry")

NaCsPlot.maybe_show()
