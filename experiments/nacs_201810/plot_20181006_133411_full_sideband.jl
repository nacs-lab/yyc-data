#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20181006_133411.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))
data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))

const spec_na_a = -140:4.0:500
const spec_cs_a = -80:1.0:80

data_na_az = NaCsData.split_data(data_na_a, spec_na_a)
data_cs_az = NaCsData.split_data(data_cs_a, spec_cs_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20181006_133411_full_sideband")

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

NaCsPlot.maybe_show()
