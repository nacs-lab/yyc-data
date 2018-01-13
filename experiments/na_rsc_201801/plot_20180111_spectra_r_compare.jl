#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180111_164819.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180112_153738.mat")
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

const spec = OrderedDict(
    :x=>linspace(18, 20, 101),
    :y=>linspace(18, 20, 101),
)

const split_a = NaCsData.split_data(data_a, spec)
const split_b = NaCsData.split_data(data_b, spec)

const data_x_hot = split_a[:x]
const data_y_hot = split_a[:y]
const data_x_cold = split_b[:x]
const data_y_cold = split_b[:y]

const prefix = joinpath(@__DIR__, "imgs", "data_20180111_spectra_r_compare")

figure()
NaCsPlot.plot_survival_data(data_x_hot, label="Hot")
NaCsPlot.plot_survival_data(data_x_cold, label="Cold")
grid()
legend()
ylim([0, 0.8])
title("X spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x")

figure()
NaCsPlot.plot_survival_data(data_y_hot, label="Hot")
NaCsPlot.plot_survival_data(data_y_cold, label="Cold")
grid()
legend()
ylim([0, 0.8])
title("Y spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y")

NaCsPlot.maybe_show()
