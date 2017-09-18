#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20170918_125622.csv")

const data_a = NaCsData.load_count_csv(iname_a)

const spec_a = OrderedDict(
    :with_cs=>linspace(18.66, 18.76, 21),
    :without_cs=>linspace(18.66, 18.76, 21)
)

const split_a = NaCsData.split_data(data_a, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20170918_125622")

figure()
NaCsPlot.plot_survival_data(split_a[:with_cs], fmt="C3.-", label="With 976")
NaCsPlot.plot_survival_data(split_a[:without_cs], fmt="C0.-", label="Without 976")
grid()
ylim([0, 1])
title("Axis Z (axial) sideband")
xlabel("Detuning from 0 field")
ylabel("Survival")
legend()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
