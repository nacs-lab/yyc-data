#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20170916_155350.csv")

const data_a = NaCsData.load_count_csv(iname_a)

const spec_a = OrderedDict(
    :radial_x=>(linspace(18, 18.4, 21), linspace(19, 19.4, 21), linspace(19.5, 19.8, 16)),
    :radial_y=>(linspace(18, 18.4, 21), linspace(19, 19.4, 21), linspace(19.5, 19.8, 16)),
    :axial_z=>linspace(18.5, 19.1, 121),
    :survival=>[0]
)

const split_a = NaCsData.split_data(data_a, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20170916_155350")

to_sideband(f) = (i, v)->(v - f) * 1000

data_rx = NaCsData.map_params(to_sideband(18.705), split_a[:radial_x])
data_ry = NaCsData.map_params(to_sideband(18.713), split_a[:radial_y])
data_az = NaCsData.map_params(to_sideband(18.619), split_a[:axial_z])

figure()
NaCsPlot.plot_survival_data(data_rx[1], fmt="C0o-")
NaCsPlot.plot_survival_data(data_rx[2], fmt="C0o-")
NaCsPlot.plot_survival_data(data_rx[3], fmt="C0o-")
grid()
ylim([0, 1])
title("Axis X (radial)")
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

figure()
NaCsPlot.plot_survival_data(data_ry[1], fmt="C0o-")
NaCsPlot.plot_survival_data(data_ry[2], fmt="C0o-")
NaCsPlot.plot_survival_data(data_ry[3], fmt="C0o-")
grid()
ylim([0, 1])
title("Axis Y (radial)")
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ry")

figure()
NaCsPlot.plot_survival_data(data_az, fmt="C0o-")
grid()
ylim([0, 1])
title("Axis Z (axial)")
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az")

NaCsPlot.maybe_show()
