#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const iname_a = joinpath(@__DIR__, "data", "data_20180912_002538.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

data_cs_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((2,), (4,)))
data_na_a = NaCsData.select_count(params_a, logicals_a, NaCsData.select_single((1,), (3,)))

const spec_a = OrderedDict(
    :x=>1.7:0.1:4.1,
    :y=>1.7:0.1:4.1,
    :z=>1.7:0.1:4.1,
    :total=>1.7:0.1:4.1,
    :align=>2.3:0.1:3.5,
)

const split_cs_a = NaCsData.split_data(data_cs_a, spec_a)
const split_na_a = NaCsData.split_data(data_na_a, spec_a)

data_cs_align = split_cs_a[:align]
data_cs_total = split_cs_a[:total]
data_cs_x = split_cs_a[:x]
data_cs_y = split_cs_a[:y]
data_cs_z = split_cs_a[:z]
data_na_total = split_na_a[:total]
data_na_x = split_na_a[:x]
data_na_y = split_na_a[:y]
data_na_z = split_na_a[:z]

const prefix = joinpath(@__DIR__, "imgs", "data_20180912_002538_merge_temp")

figure()
NaCsPlot.plot_survival_data(data_cs_total, fmt="C0.-", label="Cs")
NaCsPlot.plot_survival_data(data_na_total, fmt="C1.-", label="Na")
NaCsPlot.plot_survival_data(data_cs_align, fmt="C2.-", label="Cs Align")
grid()
legend()
ylim([0, 1])
title("Total Survival")
xlabel("DMerge (\$\\mu m\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_survival")

figure()
NaCsPlot.plot_survival_data(data_cs_x, fmt="C0.-", label="X")
NaCsPlot.plot_survival_data(data_cs_y, fmt="C1.-", label="Y")
NaCsPlot.plot_survival_data(data_cs_z, fmt="C2.-", label="Z")
grid()
legend()
ylim([0, 0.15])
title("Cs Cooling sideband")
xlabel("DMerge (\$\\mu m\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

figure()
NaCsPlot.plot_survival_data(data_na_x, fmt="C0.-", label="X")
NaCsPlot.plot_survival_data(data_na_y, fmt="C1.-", label="Y")
NaCsPlot.plot_survival_data(data_na_z, fmt="C2.-", label="Z")
grid()
legend()
ylim([0, 0.2])
title("Na Cooling sideband")
xlabel("DMerge (\$\\mu m\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

NaCsPlot.maybe_show()
