#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180122_214615.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180122_111933.mat")
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

const spec_a = OrderedDict(
    :xp1=>linspace(0, 135, 16),
    :yp1=>linspace(6, 90, 15),
    :zp1=>linspace(17, 255, 15),
    :x0=>linspace(3.5, 52.5, 15),
    :y0=>linspace(2.5, 37.5, 15),
    :z0=>linspace(7, 105, 15),
    :xf=>linspace(18, 20, 101),
    :yf=>linspace(18, 20, 101),
    :zf=>linspace(18.5, 19.1, 121)
)

const spec_b = OrderedDict(
    :xp1=>linspace(0, 270, 31),
    :yp1=>linspace(6, 180, 30),
    :zp1=>linspace(17, 510, 30),
    :x0=>linspace(3.5, 140, 40),
    :y0=>linspace(2.5, 100, 40),
    :z0=>linspace(7, 280, 40),
)

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)

const data_hot_xp1 = split_a[:xp1]
const data_hot_yp1 = [data_hot_xp1[1]; split_a[:yp1]]
const data_hot_zp1 = [data_hot_xp1[1]; split_a[:zp1]]
const data_hot_x0 = [data_hot_xp1[1]; split_a[:x0]]
const data_hot_y0 = [data_hot_xp1[1]; split_a[:y0]]
const data_hot_z0 = [data_hot_xp1[1]; split_a[:z0]]

const data_cold_xp1 = split_b[:xp1]
const data_cold_yp1 = [data_cold_xp1[1]; split_b[:yp1]]
const data_cold_zp1 = [data_cold_xp1[1]; split_b[:zp1]]
const data_cold_x0 = [data_cold_xp1[1]; split_b[:x0]]
const data_cold_y0 = [data_cold_xp1[1]; split_b[:y0]]
const data_cold_z0 = [data_cold_xp1[1]; split_b[:z0]]

const prefix = joinpath(@__DIR__, "imgs", "fit_20180122_rabi")

#### Fitting

# TODO

figure()
NaCsPlot.plot_survival_data(data_cold_xp1, fmt="C0o-", label="Cold")
NaCsPlot.plot_survival_data(data_hot_xp1, fmt="C1o-", label="Hot")
grid()
ylim([0, 1])
title("X")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_xp1")

figure()
NaCsPlot.plot_survival_data(data_cold_yp1, fmt="C0o-", label="Cold")
NaCsPlot.plot_survival_data(data_hot_yp1, fmt="C1o-", label="Hot")
grid()
ylim([0, 1])
title("Y")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_yp1")

figure()
NaCsPlot.plot_survival_data(data_cold_zp1, fmt="C0o-", label="Cold")
NaCsPlot.plot_survival_data(data_hot_zp1, fmt="C1o-", label="Hot")
grid()
ylim([0, 1])
title("Z")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_zp1")

figure()
NaCsPlot.plot_survival_data(data_cold_x0, fmt="C0o-", label="Cold")
NaCsPlot.plot_survival_data(data_hot_x0, fmt="C1o-", label="Hot")
grid()
ylim([0, 1])
title("X")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_x0")

figure()
NaCsPlot.plot_survival_data(data_cold_y0, fmt="C0o-", label="Cold")
NaCsPlot.plot_survival_data(data_hot_y0, fmt="C1o-", label="Hot")
grid()
ylim([0, 1])
title("Y")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_y0")

figure()
NaCsPlot.plot_survival_data(data_cold_z0, fmt="C0o-", label="Cold")
NaCsPlot.plot_survival_data(data_hot_z0, fmt="C1o-", label="Hot")
grid()
ylim([0, 1])
title("Z")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_z0")

NaCsPlot.maybe_show()
