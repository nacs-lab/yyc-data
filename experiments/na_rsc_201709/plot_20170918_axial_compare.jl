#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20170918_125622.csv")
const iname_b = joinpath(@__DIR__, "data", "data_20170918_145905.csv")

const data_a = NaCsData.load_count_csv(iname_a)
const data_b = NaCsData.load_count_csv(iname_b)

const spec_a = OrderedDict(
    :with_cs=>linspace(18.66, 18.76, 21),
    :without_cs=>linspace(18.66, 18.76, 21)
)
const spec_b = OrderedDict(
    :radial_x=>(linspace(18, 18.4, 21), linspace(19, 19.4, 21), linspace(19.5, 19.8, 16)),
    :radial_y=>(linspace(18, 18.4, 21), linspace(19, 19.4, 21), linspace(19.5, 19.8, 16)),
    :axial_z=>linspace(18.5, 18.85, 71),
    :survival=>[0]
)

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)

const prefix = joinpath(@__DIR__, "imgs", "data_20170918_axial_compare")

figure()
NaCsPlot.plot_survival_data(split_a[:with_cs], fmt="C3.-", label="With 976")
NaCsPlot.plot_survival_data(split_a[:without_cs], fmt="C0.-", label="Without 976")
NaCsPlot.plot_survival_data(split_b[:axial_z], fmt="C2.-", label="Old sequence")
grid()
ylim([0, 1])
title("Axis Z (axial) sideband")
xlabel("Detuning from 0 field")
ylabel("Survival")
legend()
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
