#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20170913_160122.csv")
const iname_b = joinpath(@__DIR__, "data", "data_20170913_183044.csv")

const data_a = NaCsData.load_count_csv(iname_a)
const data_b = NaCsData.load_count_csv(iname_b)

const spec_a = OrderedDict(
    :radial_y=>(linspace(18, 18.4, 20), linspace(19, 19.4, 20), linspace(19.5, 19.8, 15)),
    :radial_x=>(linspace(18, 18.4, 20), linspace(19, 19.4, 20), linspace(19.5, 19.8, 15)),
)
const spec_b = OrderedDict(
    :radial_y=>(linspace(18, 18.4, 21), linspace(19, 19.4, 21), linspace(19.5, 19.8, 16)),
    :radial_x=>(linspace(18, 18.4, 21), linspace(19, 19.4, 21), linspace(19.5, 19.8, 16)),
)

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)

const prefix = joinpath(@__DIR__, "imgs", "data_20170913_radial_cool")

to_sideband(f) = (i, v)->(v - f) * 1000

data_nocool_rx = NaCsData.map_params(to_sideband(18.705), split_a[:radial_x])
data_nocool_ry = NaCsData.map_params(to_sideband(18.713), split_a[:radial_y])
data_cool_rx = NaCsData.map_params(to_sideband(18.705), split_b[:radial_x])
data_cool_ry = NaCsData.map_params(to_sideband(18.713), split_b[:radial_y])

figure()
# Without cooling
NaCsPlot.plot_survival_data(data_nocool_rx[1], fmt="C1o-", label="Initial")
NaCsPlot.plot_survival_data(data_nocool_rx[2], fmt="C1o-")
NaCsPlot.plot_survival_data(data_nocool_rx[3], fmt="C1o-")
# With cooling
NaCsPlot.plot_survival_data(data_cool_rx[1], fmt="C0o-", label="Cooled")
NaCsPlot.plot_survival_data(data_cool_rx[2], fmt="C0o-")
NaCsPlot.plot_survival_data(data_cool_rx[3], fmt="C0o-")
grid()
ylim([0, 1])
title("Axis X (radial)")
legend()
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

figure()
# Without cooling
NaCsPlot.plot_survival_data(data_nocool_ry[1], fmt="C1o-", label="Initial")
NaCsPlot.plot_survival_data(data_nocool_ry[2], fmt="C1o-")
NaCsPlot.plot_survival_data(data_nocool_ry[3], fmt="C1o-")
# With cooling
NaCsPlot.plot_survival_data(data_cool_ry[1], fmt="C0o-", label="Cooled")
NaCsPlot.plot_survival_data(data_cool_ry[2], fmt="C0o-")
NaCsPlot.plot_survival_data(data_cool_ry[3], fmt="C0o-")
grid()
ylim([0, 1])
title("Axis Y (radial)")
legend()
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ry")

NaCsPlot.maybe_show()
