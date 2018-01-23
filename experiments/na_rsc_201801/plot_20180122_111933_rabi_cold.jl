#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180122_111933.mat")
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
    :xp1=>linspace(0, 270, 31),
    :yp1=>linspace(6, 180, 30),
    :zp1=>linspace(17, 510, 30),
    :x0=>linspace(3.5, 140, 40),
    :y0=>linspace(2.5, 100, 40),
    :z0=>linspace(7, 280, 40),
)

const split_a = NaCsData.split_data(data_a, spec)

const data_xp1 = split_a[:xp1]
const data_yp1 = [data_xp1[1]; split_a[:yp1]]
const data_zp1 = [data_xp1[1]; split_a[:zp1]]
const data_x0 = [data_xp1[1]; split_a[:x0]]
const data_y0 = [data_xp1[1]; split_a[:y0]]
const data_z0 = [data_xp1[1]; split_a[:z0]]

const prefix = joinpath(@__DIR__, "imgs", "data_20180122_111933_rabi_cold")

figure()
NaCsPlot.plot_survival_data(data_xp1, fmt="o-", label="+1 order")
NaCsPlot.plot_survival_data(data_x0, fmt="o-", label="0 order")
grid()
ylim([0, 1])
title("X")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x")

figure()
NaCsPlot.plot_survival_data(data_yp1, fmt="o-", label="+1 order")
NaCsPlot.plot_survival_data(data_y0, fmt="o-", label="0 order")
grid()
ylim([0, 1])
title("Y")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y")

figure()
NaCsPlot.plot_survival_data(data_zp1, fmt="o-", label="+1 order")
NaCsPlot.plot_survival_data(data_z0, fmt="o-", label="0 order")
grid()
ylim([0, 1])
title("Z")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_z")

NaCsPlot.maybe_show()
