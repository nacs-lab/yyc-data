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
const iname_b = joinpath(@__DIR__, "data", "data_20180123_155129.mat")
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
    :xp1=>linspace(0, 135, 16),
    :yp1=>linspace(6, 90, 15),
    :zp1=>linspace(17, 255, 15),
    :x0=>linspace(3.5, 52.5, 15),
    :y0=>linspace(2.5, 37.5, 15),
    :z0=>linspace(7, 105, 15),
    :xf=>linspace(18, 20, 101),
    :yf=>linspace(18, 20, 101),
    :zf=>linspace(18.4, 18.55, 31),
)

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)

const data_xp1 = split_a[:xp1]
const data_yp1 = [data_xp1[1]; split_a[:yp1]]
const data_zp1 = [data_xp1[1]; split_a[:zp1]]
const data_x0 = [data_xp1[1]; split_a[:x0]]
const data_y0 = [data_xp1[1]; split_a[:y0]]
const data_z0 = [data_xp1[1]; split_a[:z0]]
const data_xf = split_a[:xf]
const data_yf = split_a[:yf]
const data_zf = [split_a[:zf]; split_b[:zf]]

const prefix = joinpath(@__DIR__, "imgs", "data_20180122_214615_hot")

figure()
NaCsPlot.plot_survival_data(data_xp1, fmt="o-", label="+1 order")
NaCsPlot.plot_survival_data(data_x0, fmt="o-", label="0 order")
grid()
ylim([0, 1])
title("X")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_x")

figure()
NaCsPlot.plot_survival_data(data_yp1, fmt="o-", label="+1 order")
NaCsPlot.plot_survival_data(data_y0, fmt="o-", label="0 order")
grid()
ylim([0, 1])
title("Y")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_y")

figure()
NaCsPlot.plot_survival_data(data_zp1, fmt="o-", label="+1 order")
NaCsPlot.plot_survival_data(data_z0, fmt="o-", label="0 order")
grid()
ylim([0, 1])
title("Z")
legend()
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rabi_z")

figure()
NaCsPlot.plot_survival_data(data_xf, fmt="C0o-")
grid()
ylim([0, 1])
title("X spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_spectrum_x")

figure()
NaCsPlot.plot_survival_data(data_yf, fmt="C0o-")
grid()
ylim([0, 1])
title("Y spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_spectrum_y")

figure()
NaCsPlot.plot_survival_data(data_zf, fmt="C0o-")
grid()
ylim([0, 1])
title("Z spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_spectrum_z")

NaCsPlot.maybe_show()
