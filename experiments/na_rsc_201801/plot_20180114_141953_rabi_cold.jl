#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180114_141953.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)

const spec = OrderedDict(
    :x=>linspace(0, 200, 41),
    :y=>linspace(4, 160, 40),
    :z=>linspace(6, 240, 40),
)

const split_a = NaCsData.split_data(data_a, spec)

const data_x = split_a[:x]
const data_y = [data_x[1]; split_a[:y]]
const data_z = [data_x[1]; split_a[:z]]

const prefix = joinpath(@__DIR__, "imgs", "data_20180114_141953_rabi_cold")

figure()
NaCsPlot.plot_survival_data(data_x, fmt="o-")
grid()
ylim([0, 1])
title("X")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x")

figure()
NaCsPlot.plot_survival_data(data_y, fmt="o-")
grid()
ylim([0, 1])
title("Y")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y")

figure()
NaCsPlot.plot_survival_data(data_z, fmt="o-")
grid()
ylim([0, 1])
title("Z")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_z")

NaCsPlot.maybe_show()
