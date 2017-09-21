#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20170920_125743.csv")

const data_a = NaCsData.load_count_csv(iname_a)

const spec_a = OrderedDict(
    :base=>[0, 1],
    :rx=>[0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14],
)

const split_a = NaCsData.split_data(NaCsData.map_params((i, v)->i, data_a), spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20170920_125743_merge_r")

figure()
NaCsPlot.plot_survival_data([split_a[:rx];], fmt="C0.-")
grid()
ylim([0, ylim()[2]])
title("Total X sideband")
xlabel("Pre-merge depth (V)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

NaCsPlot.maybe_show()
