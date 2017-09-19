#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20170919_110453.csv")

const data_a = NaCsData.load_count_csv(iname_a)

const spec_a = OrderedDict(
    :base=>[0],
    :lower05=>[0, 10, 20, 40, 80],
    :lower10=>[0, 10, 20, 40, 80],
    :lower15=>[0, 10, 20, 40, 80],
)

const split_a = NaCsData.split_data(NaCsData.map_params((i, v)->i, data_a), spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20170919_110453_lower_hold")

@show split_a[:base]

figure()
NaCsPlot.plot_survival_data(split_a[:lower05], fmt=".-", label="Lower=0.05")
NaCsPlot.plot_survival_data(split_a[:lower10], fmt=".-", label="Lower=0.10")
NaCsPlot.plot_survival_data(split_a[:lower15], fmt=".-", label="Lower=0.15")
grid()
ylim([0, 0.25])
title("Axis Z (axial) sideband")
xlabel("Hold time (ms)")
ylabel("Survival")
legend()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
