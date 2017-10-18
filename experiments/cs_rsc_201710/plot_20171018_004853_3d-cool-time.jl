#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20171018_004853.csv")

const data_a = NaCsData.load_count_csv(iname_a)

const spec_a = OrderedDict(
    :axial_zp1=>0:400:8000,
    :radial_xp1=>0:40:800,
    :radial_yp1=>0:40:800,
    :axial_z0=>0:20:400,
    :radial_x0=>0:5:100,
    :radial_y0=>0:5:100,
)

const split_a = NaCsData.split_data(data_a, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20171018_004853_3d-cool-time")

data_rxp1 = split_a[:radial_xp1]
data_ryp1 = split_a[:radial_yp1]
data_azp1 = split_a[:axial_zp1]
data_rx0 = split_a[:radial_x0]
data_ry0 = split_a[:radial_y0]
data_az0 = split_a[:axial_z0]

figure()
NaCsPlot.plot_survival_data(data_rxp1, fmt="C0o-")
grid()
ylim([0, 1])
title("Axis X Heating (radial)")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rxp1")

figure()
NaCsPlot.plot_survival_data(data_ryp1, fmt="C0o-")
grid()
ylim([0, 1])
title("Axis Y Heating (radial)")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ryp1")

figure()
NaCsPlot.plot_survival_data(data_azp1, fmt="C0o-")
grid()
ylim([0, 0.6])
title("Axis Z Heating (axial)")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_azp1")

figure()
NaCsPlot.plot_survival_data(data_rx0, fmt="C0o-")
grid()
ylim([0, 1])
title("Axis X Carrier (radial)")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx0")

figure()
NaCsPlot.plot_survival_data(data_ry0, fmt="C0o-")
grid()
ylim([0, 1])
title("Axis Y Carrier (radial)")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ry0")

figure()
NaCsPlot.plot_survival_data(data_az0, fmt="C0o-")
grid()
ylim([0, 1])
title("Axis Z Carrier (axial)")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az0")

NaCsPlot.maybe_show()
