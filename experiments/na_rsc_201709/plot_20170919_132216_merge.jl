#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20170919_132216.csv")

const data_a = NaCsData.load_count_csv(iname_a)

const spec_a = OrderedDict(
    :base_az=>[0],
    :base_survival=>[0],
    :survival=>[0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14],
    :az=>[0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14],
)

const split_a = NaCsData.split_data(NaCsData.map_params((i, v)->i, data_a), spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20170919_132216_merge")

@show split_a[:base_az]
@show split_a[:base_survival]

base_az = NaCsData.get_values(split_a[:base_az])[2][2]
base_survival = NaCsData.get_values(split_a[:base_survival])[2][2]

figure()
axhline(base_survival, color="C1")
NaCsPlot.plot_survival_data(split_a[:survival], fmt="C0.-")
grid()
ylim([0, 1])
title("Total survival")
xlabel("Pre-merge depth (V)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_survival")

figure()
axhline(base_az, color="C1")
NaCsPlot.plot_survival_data(split_a[:az], fmt="C0.-")
grid()
ylim([0, ylim()[2]])
title("Total axial sideband")
xlabel("Pre-merge depth (V)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az")

NaCsPlot.maybe_show()
