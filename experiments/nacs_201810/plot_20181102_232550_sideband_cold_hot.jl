#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181102_232550.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))

const spec_na_a = OrderedDict(
    :cool=>OrderedDict(
        :az=>(-126 .+ (-30:4.0:30), 76 .+ (-30:4.0:30)),
        :rx=>(-479 .+ (-45:6.0:45), 614 .+ (-90:12.0:90)),
        :ry=>(-490 .+ (-45:6.0:45), 631 .+ (-90:12.0:90))
    ),
    :nocool=>OrderedDict(
        :az=>(-126 .+ (-30:4.0:30), 76 .+ (-30:4.0:30)),
        :rx=>(-479 .+ (-45:6.0:45), 614 .+ (-90:12.0:90)),
        :ry=>(-490 .+ (-45:6.0:45), 631 .+ (-90:12.0:90))
    ),
)

const spec_cs_a = OrderedDict(
    :cool=>OrderedDict(
        :az=>(-10 .+ (-15:2.0:15), 33 .+ (-15:2.0:15)),
        :rx=>(-135 .+ (-45:6.0:45), 109 .+ (-75:10.0:75)),
        :ry=>(-129 .+ (-45:6.0:45), 112 .+ (-60:8.0:60))
    ),
    :nocool=>OrderedDict(
        :az=>(-10 .+ (-15:2.0:15), 33 .+ (-15:2.0:15)),
        :rx=>(-135 .+ (-45:6.0:45), 109 .+ (-75:10.0:75)),
        :ry=>(-129 .+ (-45:6.0:45), 112 .+ (-60:8.0:60))
    ),
)

const split_na_a = NaCsData.split_data(data_na_a, spec_na_a)
const split_cs_a = NaCsData.split_data(data_cs_a, spec_cs_a)

data_na_nocool_az = split_na_a[:nocool][:az]
data_na_nocool_rx = split_na_a[:nocool][:rx]
data_na_nocool_ry = split_na_a[:nocool][:ry]
data_cs_nocool_az = split_cs_a[:nocool][:az]
data_cs_nocool_rx = split_cs_a[:nocool][:rx]
data_cs_nocool_ry = split_cs_a[:nocool][:ry]

data_na_cool_az = split_na_a[:cool][:az]
data_na_cool_rx = split_na_a[:cool][:rx]
data_na_cool_ry = split_na_a[:cool][:ry]
data_cs_cool_az = split_cs_a[:cool][:az]
data_cs_cool_rx = split_cs_a[:cool][:rx]
data_cs_cool_ry = split_cs_a[:cool][:ry]

const prefix = joinpath(@__DIR__, "imgs", "data_20181102_232550_sideband_cold_hot")

figure()
NaCsPlot.plot_survival_data(data_na_nocool_az[1], fmt="C0.-", label="No cool")
NaCsPlot.plot_survival_data(data_na_nocool_az[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_cool_az[1], fmt="C1.-", label="Cool")
NaCsPlot.plot_survival_data(data_na_cool_az[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Na Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
NaCsPlot.plot_survival_data(data_cs_nocool_az[1], fmt="C0.-", label="No cool")
NaCsPlot.plot_survival_data(data_cs_nocool_az[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_cool_az[1], fmt="C1.-", label="Cool")
NaCsPlot.plot_survival_data(data_cs_cool_az[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Cs Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
NaCsPlot.plot_survival_data(data_na_nocool_rx[1], fmt="C0.-", label="No cool")
NaCsPlot.plot_survival_data(data_na_nocool_rx[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_cool_rx[1], fmt="C1.-", label="Cool")
NaCsPlot.plot_survival_data(data_na_cool_rx[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Na X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_rx")

figure()
NaCsPlot.plot_survival_data(data_cs_nocool_rx[1], fmt="C0.-", label="No cool")
NaCsPlot.plot_survival_data(data_cs_nocool_rx[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_cool_rx[1], fmt="C1.-", label="Cool")
NaCsPlot.plot_survival_data(data_cs_cool_rx[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Cs X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_rx")

figure()
NaCsPlot.plot_survival_data(data_na_nocool_ry[1], fmt="C0.-", label="No cool")
NaCsPlot.plot_survival_data(data_na_nocool_ry[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_cool_ry[1], fmt="C1.-", label="Cool")
NaCsPlot.plot_survival_data(data_na_cool_ry[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Na Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_ry")

figure()
NaCsPlot.plot_survival_data(data_cs_nocool_ry[1], fmt="C0.-", label="No cool")
NaCsPlot.plot_survival_data(data_cs_nocool_ry[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_cool_ry[1], fmt="C1.-", label="Cool")
NaCsPlot.plot_survival_data(data_cs_cool_ry[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Cs Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ry")

NaCsPlot.maybe_show()
