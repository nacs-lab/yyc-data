#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180110_221739.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180111_053724.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
const iname_c = joinpath(@__DIR__, "data", "data_20180111_164819.mat")
const params_c, logicals_c = NaCsData.load_striped_mat(iname_c)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)
data_b = NaCsData.select_count(params_b, logicals_b, selector)
data_c = NaCsData.select_count(params_c, logicals_c, selector)

const spec = OrderedDict(
    :x=>linspace(18, 20, 101),
    :y=>linspace(18, 20, 101),
)

const split_a = NaCsData.split_data(data_a, spec)
const split_b = NaCsData.split_data(data_b, spec)
const split_c = NaCsData.split_data(data_c, spec)

const data_x_old = split_a[:x]
const data_y_old = split_a[:y]
const data_x_sqr = split_b[:x]
const data_y_sqr = split_b[:y]
const data_x_new = split_c[:x]
const data_y_new = split_c[:y]

const prefix = joinpath(@__DIR__, "imgs", "data_20180110_spectra_r_hot")

figure()
NaCsPlot.plot_survival_data(data_x_sqr, label="Square")
NaCsPlot.plot_survival_data(data_x_old, label="Old")
NaCsPlot.plot_survival_data(data_x_new, label="New")
grid()
legend()
ylim([0, 0.8])
title("X spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x")

figure()
NaCsPlot.plot_survival_data(data_y_sqr, label="Square")
NaCsPlot.plot_survival_data(data_y_old, label="Old")
NaCsPlot.plot_survival_data(data_y_new, label="New")
grid()
legend()
ylim([0, 0.8])
title("Y spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y")

NaCsPlot.maybe_show()
