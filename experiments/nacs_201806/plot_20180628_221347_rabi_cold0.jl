#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180628_221347.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)

const spec_a = OrderedDict(
    :axial_z=>linspace(0.1, 300, 21),
    :radial_x=>linspace(0.1, 120, 21),
    :radial_y=>linspace(0.1, 100, 21),
)

const split_a = NaCsData.split_data(data_a, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20180628_221347_rabi")

data_rx = split_a[:radial_x]
data_ry = split_a[:radial_y]
data_az = split_a[:axial_z]

figure()
NaCsPlot.plot_survival_data(data_rx, fmt="C0o-") # 40
grid()
ylim([0, 1])
title("Axis X (radial)")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

figure()
NaCsPlot.plot_survival_data(data_ry, fmt="C0o-") # 35
grid()
ylim([0, 1])
title("Axis Y (radial)")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ry")

figure()
NaCsPlot.plot_survival_data(data_az, fmt="C0o-") # 81
grid()
ylim([0, 1])
title("Axis Z (axial)")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az")

NaCsPlot.maybe_show()
