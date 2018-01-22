#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180120_135037.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180120_212214.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)
data_b = NaCsData.select_count(params_b, logicals_b, selector)

# 70, 46, 90, 135

const spec_a = OrderedDict(
    :x=>linspace(0, 100, 26),
    :y=>linspace(3, 75, 25),
    :z=>linspace(5, 125, 25),
)

const spec_b = OrderedDict(
    :x=>linspace(0, 100, 26),
    :y=>linspace(3, 75, 25),
    :z=>linspace(10, 250, 25),
)

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)

const data_x = split_a[:x]
const data_y = [data_x[1]; split_a[:y]]
const data_z = [data_x[1]; split_a[:z]]
const data_x2 = split_b[:x]
const data_z2 = [data_x2[1]; split_b[:z]]

const prefix = joinpath(@__DIR__, "imgs", "data_20180120_rabi_cold")

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
NaCsPlot.plot_survival_data(data_z, fmt="o-", label="Full power")
NaCsPlot.plot_survival_data(data_z2, fmt="o-", label="2/3 power")
grid()
ylim([0, 1])
title("Z")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_z")

NaCsPlot.maybe_show()
