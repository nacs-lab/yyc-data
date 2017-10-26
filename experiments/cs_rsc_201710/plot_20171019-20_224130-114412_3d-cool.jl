#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20171019_224130.csv")
const iname_b = joinpath(@__DIR__, "data", "data_20171020_114412.csv")

const data_a = NaCsData.load_count_csv(iname_a)
const data_b = NaCsData.load_count_csv(iname_b)

const spec = OrderedDict(
    :radial_x=>(-125:2:-85, 152:2:192),
    :radial_y=>(-130:2:-90, 146:2:186),
    :axial_z=>(5:1.5:35, 55:1.5:85),
    :survival=>[1],
)

const split_a = NaCsData.split_data([data_a; data_b], spec)

const prefix = joinpath(@__DIR__, "imgs", "data_20171019-20_224130-114412_3d-cool")
const sorted_prefix = joinpath(@__DIR__, "sorted", "data_20171019-20_224130-114412_3d-cool")

data_rxp1 = split_a[:radial_x][1]
data_rxm1 = split_a[:radial_x][2]
data_ryp1 = split_a[:radial_y][1]
data_rym1 = split_a[:radial_y][2]
data_azp1 = split_a[:axial_z][1]
data_azm1 = split_a[:axial_z][2]

NaCsData.dump_raw(split_a[:survival])

NaCsData.dump_raw("$(sorted_prefix)_rx.csv", [data_rxp1; data_rxm1])
NaCsData.dump_raw("$(sorted_prefix)_ry.csv", [data_ryp1; data_rym1])
NaCsData.dump_raw("$(sorted_prefix)_az.csv", [data_azp1; data_azm1])

figure()
NaCsPlot.plot_survival_data(data_rxp1, fmt="C0o-")
NaCsPlot.plot_survival_data(data_rxm1, fmt="C0o-")
grid()
ylim([0, 1])
title("Axis X (radial)")
xlabel("Frequency (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

figure()
NaCsPlot.plot_survival_data(data_ryp1, fmt="C0o-")
NaCsPlot.plot_survival_data(data_rym1, fmt="C0o-")
grid()
ylim([0, 1])
title("Axis Y (radial)")
xlabel("Frequency (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ry")

figure()
NaCsPlot.plot_survival_data(data_azp1, fmt="C0o-")
NaCsPlot.plot_survival_data(data_azm1, fmt="C0o-")
grid()
ylim([0, 0.5])
title("Axis Z (axial)")
xlabel("Frequency (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az")

NaCsPlot.maybe_show()
