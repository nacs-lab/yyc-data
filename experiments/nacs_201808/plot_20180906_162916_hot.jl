#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180906_162916.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))

const spec_na_a = OrderedDict(
    :az=>-140:4.0:500,
    :rx=>-500:10.:1100,
    :ry=>-500:10.:1100,
)

const spec_cs_a = OrderedDict(
    :az=>-80:1.0:80,
    :rx=>-160:3.0:320,
    :ry=>-160:3.0:320,
)

const split_na_a = NaCsData.split_data(data_na_a, spec_na_a)
const split_cs_a = NaCsData.split_data(data_cs_a, spec_cs_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20180906_162916_hot")

figure()
NaCsPlot.plot_survival_data(split_na_a[:az], fmt="C0.-")
grid()
ylim([0, 0.6])
title("Na Hot Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
NaCsPlot.plot_survival_data(split_na_a[:rx], fmt="C0.-")
grid()
ylim([0, 0.6])
title("Na Hot Radial X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_rx")

figure()
NaCsPlot.plot_survival_data(split_na_a[:ry], fmt="C0.-")
grid()
ylim([0, 0.6])
title("Na Hot Radial Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_ry")

figure()
NaCsPlot.plot_survival_data(split_cs_a[:az], fmt="C0.-")
grid()
ylim([0, 0.6])
title("Cs Hot Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
NaCsPlot.plot_survival_data(split_cs_a[:rx], fmt="C0.-")
grid()
ylim([0, 0.65])
title("Cs Hot Radial X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_rx")

figure()
NaCsPlot.plot_survival_data(split_cs_a[:ry], fmt="C0.-")
grid()
ylim([0, 0.65])
title("Cs Hot Radial Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ry")

NaCsPlot.maybe_show()
