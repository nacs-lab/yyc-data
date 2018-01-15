#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180113_201830.mat")
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
    :x=>(linspace(18.04, 18.26, 12), linspace(18.98, 19.24, 14),
         linspace(19.48, 19.74, 14)),
    :y=>(linspace(18.04, 18.26, 12), linspace(18.98, 19.24, 14),
         linspace(19.48, 19.74, 14)),
    :z=>linspace(18.46, 18.80, 69),
)

const split_a = NaCsData.split_data(data_a, spec)

const data_x = split_a[:x]
const data_y = split_a[:y]
const data_z = split_a[:z]

const prefix = joinpath(@__DIR__, "imgs", "data_20180113_201830_spectrum_cold")

figure()
NaCsPlot.plot_survival_data(data_x[1], fmt="C0o-")
NaCsPlot.plot_survival_data(data_x[2], fmt="C0o-")
NaCsPlot.plot_survival_data(data_x[3], fmt="C0o-")
grid()
ylim([0, 1])
title("X spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x")

figure()
NaCsPlot.plot_survival_data(data_y[1], fmt="C0o-")
NaCsPlot.plot_survival_data(data_y[2], fmt="C0o-")
NaCsPlot.plot_survival_data(data_y[3], fmt="C0o-")
grid()
ylim([0, 1])
title("Y spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y")

figure()
NaCsPlot.plot_survival_data(data_z, fmt="o-")
grid()
ylim([0, 1])
title("Z spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_z")

NaCsPlot.maybe_show()
