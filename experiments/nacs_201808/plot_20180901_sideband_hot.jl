#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180901_225946.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180902_075701.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
const iname_c = joinpath(@__DIR__, "data", "data_20180902_092224.mat")
const params_c, logicals_c = NaCsData.load_striped_mat(iname_c)
const iname_d = joinpath(@__DIR__, "data", "data_20180902_143525.mat")
const params_d, logicals_d = NaCsData.load_striped_mat(iname_d)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))
data_na_b = NaCsData.select_count(params_b, logicals_b, NaCsData.select_single((1,), (3,)))
data_cs_b = NaCsData.select_count(params_b, logicals_b, NaCsData.select_single((2,), (4,)))
data_na_c = NaCsData.select_count(params_c, logicals_c, NaCsData.select_single((1,), (3,)))
data_cs_c = NaCsData.select_count(params_c, logicals_c, NaCsData.select_single((2,), (4,)))
data_na_d = NaCsData.select_count(params_d, logicals_d, NaCsData.select_single((1,), (3,)))
data_cs_d = NaCsData.select_count(params_d, logicals_d, NaCsData.select_single((2,), (4,)))

const spec_na_a = OrderedDict(
    :az=>-140:4.0:500,
    :rx=>(-466:6.0:-376, 489:12.0:669),
    :ry=>(-475:6.0:-385, 498:12.0:678)
)

const spec_cs_a = OrderedDict(
    :az=>-80:1.0:80,
    :rx=>(-198:6.0:-108, 25:10.0:175),
    :ry=>(-195:6.0:-104, 50:8.0:170)
)

const spec_na_c = OrderedDict(
    :rx=>(-466:6.0:-376, 489:12.0:669),
    :ry=>(-475:6.0:-385, 498:12.0:678)
)

const spec_cs_c = OrderedDict(
    :rx=>(-198:6.0:-108, 25:10.0:175),
    :ry=>(-195:6.0:-104, 50:8.0:170)
)

const spec_na_d = OrderedDict(
    :rx=>(-376:6.0:-286, 393:12.0:573),
    :ry=>(-385:6.0:-295, 402:12.0:582)
)

const spec_cs_d = OrderedDict(
    :rx=>(-156:6.0:-66, 45:10.0:195),
    :ry=>(-158:6.0:-68, 50:8.0:170)
)

const split_na_a = NaCsData.split_data(data_na_a, spec_na_a)
const split_cs_a = NaCsData.split_data(data_cs_a, spec_cs_a)
const split_na_b = NaCsData.split_data(data_na_b, spec_na_a)
const split_cs_b = NaCsData.split_data(data_cs_b, spec_cs_a)
const split_na_c = NaCsData.split_data(data_na_c, spec_na_c)
const split_cs_c = NaCsData.split_data(data_cs_c, spec_cs_c)
const split_na_d = NaCsData.split_data(data_na_d, spec_na_d)
const split_cs_d = NaCsData.split_data(data_cs_d, spec_cs_d)

data_na_az = [split_na_a[:az]; split_na_b[:az]]
data_na_rx = ([split_na_a[:rx][1]; split_na_b[:rx][1]; split_na_c[:rx][1]; split_na_d[:rx][1]],
              [split_na_a[:rx][2]; split_na_b[:rx][2]; split_na_c[:rx][2]; split_na_d[:rx][2]])
data_na_ry = ([split_na_a[:ry][1]; split_na_b[:ry][1]; split_na_c[:ry][1]; split_na_d[:ry][1]],
              [split_na_a[:ry][2]; split_na_b[:ry][2]; split_na_c[:ry][2]; split_na_d[:ry][2]])
data_cs_az = [split_cs_a[:az]; split_cs_b[:az]]
data_cs_rx = ([split_cs_a[:rx][1]; split_cs_b[:rx][1]; split_cs_c[:rx][1]; split_cs_d[:rx][1]],
              [split_cs_a[:rx][2]; split_cs_b[:rx][2]; split_cs_c[:rx][2]; split_cs_d[:rx][2]])
data_cs_ry = ([split_cs_a[:ry][1]; split_cs_b[:ry][1]; split_cs_c[:ry][1]; split_cs_d[:ry][1]],
              [split_cs_a[:ry][2]; split_cs_b[:ry][2]; split_cs_c[:ry][2]; split_cs_d[:ry][2]])

const prefix = joinpath(@__DIR__, "imgs", "data_20180901_sideband_hot")

figure()
NaCsPlot.plot_survival_data(data_na_az, fmt="C0.-")
grid()
ylim([0, 1])
title("Na Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
NaCsPlot.plot_survival_data(data_cs_az, fmt="C0.-")
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
