#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180914_015357.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))

const spec_na_a = OrderedDict(
    :nomerge=>OrderedDict(
        :az=>(-146:4.0:-86, 46:4.0:106),
        :rx=>(-493:6.0:-403, 518:12.0:698),
        :ry=>(-497:6.0:-407, 524:12.0:704)
    ),
    :merge=>OrderedDict(
        :az=>(-146:4.0:-86, 46:4.0:106),
        :rx=>(-493:6.0:-403, 518:12.0:698),
        :ry=>(-497:6.0:-407, 524:12.0:704)
    ),
)

const spec_cs_a = OrderedDict(
    :nomerge=>OrderedDict(
        :az=>(-13:2.0:17, 29:2.0:59),
        :rx=>(-160:6.0:-70, 62:10.0:212),
        :ry=>(-158:6.0:-68, 75:8.0:195)
    ),
    :merge=>OrderedDict(
        :az=>(-13:2.0:17, 29:2.0:59),
        :rx=>(-160:6.0:-70, 62:10.0:212),
        :ry=>(-158:6.0:-68, 75:8.0:195)
    ),
)

const split_na_a = NaCsData.split_data(data_na_a, spec_na_a)
const split_cs_a = NaCsData.split_data(data_cs_a, spec_cs_a)

data_na_nomerge_az = split_na_a[:nomerge][:az]
data_na_nomerge_rx = split_na_a[:nomerge][:rx]
data_na_nomerge_ry = split_na_a[:nomerge][:ry]
data_cs_nomerge_az = split_cs_a[:nomerge][:az]
data_cs_nomerge_rx = split_cs_a[:nomerge][:rx]
data_cs_nomerge_ry = split_cs_a[:nomerge][:ry]

data_na_merge_az = split_na_a[:merge][:az]
data_na_merge_rx = split_na_a[:merge][:rx]
data_na_merge_ry = split_na_a[:merge][:ry]
data_cs_merge_az = split_cs_a[:merge][:az]
data_cs_merge_rx = split_cs_a[:merge][:rx]
data_cs_merge_ry = split_cs_a[:merge][:ry]

const prefix = joinpath(@__DIR__, "imgs", "data_20180914_015357_merge_cmp")

figure()
NaCsPlot.plot_survival_data(data_na_nomerge_az[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_na_nomerge_az[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_merge_az[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_na_merge_az[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Na Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
NaCsPlot.plot_survival_data(data_cs_nomerge_az[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_cs_nomerge_az[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_merge_az[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_cs_merge_az[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Cs Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
NaCsPlot.plot_survival_data(data_na_nomerge_rx[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_na_nomerge_rx[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_merge_rx[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_na_merge_rx[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Na X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_rx")

figure()
NaCsPlot.plot_survival_data(data_cs_nomerge_rx[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_cs_nomerge_rx[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_merge_rx[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_cs_merge_rx[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Cs X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_rx")

figure()
NaCsPlot.plot_survival_data(data_na_nomerge_ry[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_na_nomerge_ry[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_merge_ry[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_na_merge_ry[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Na Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_ry")

figure()
NaCsPlot.plot_survival_data(data_cs_nomerge_ry[1], fmt="C0.-", label="Before")
NaCsPlot.plot_survival_data(data_cs_nomerge_ry[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_merge_ry[1], fmt="C1.-", label="After")
NaCsPlot.plot_survival_data(data_cs_merge_ry[2], fmt="C1.-")
grid()
legend()
ylim([0, 1])
title("Cs Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ry")

NaCsPlot.maybe_show()
