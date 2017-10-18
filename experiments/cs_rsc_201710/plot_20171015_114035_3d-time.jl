#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20171015_114035.csv")

const data_a = NaCsData.load_count_csv(iname_a)

const spec_a = OrderedDict(
    :axial_z=>0:25:1000,
    :radial_x=>0:5:200,
    :radial_y=>0:5:200,
)

const split_a = NaCsData.split_data(data_a, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20171015_114035_3d-time")

data_rx = split_a[:radial_x]
data_ry = split_a[:radial_y]
data_az = split_a[:axial_z]

figure()
NaCsPlot.plot_survival_data(data_rx, fmt="C0o-")
grid()
ylim([0, 0.5])
title("Axis X (radial)")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

figure()
NaCsPlot.plot_survival_data(data_ry, fmt="C0o-")
grid()
ylim([0, 0.5])
title("Axis Y (radial)")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ry")

figure()
NaCsPlot.plot_survival_data(data_az, fmt="C0o-")
grid()
ylim([0, 0.3])
title("Axis Z (axial)")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az")

NaCsPlot.maybe_show()
