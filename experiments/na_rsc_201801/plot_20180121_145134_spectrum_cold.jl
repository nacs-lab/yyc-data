#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180121_145134.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)

# 18.163, 19.130
# 18.143, 19.130, 19.615
# 18.515, 18.606, 18.685, 18.773, 18.858, 18.940

const spec = OrderedDict(
    :x=>(linspace(18.11, 18.20, 10), linspace(19.05, 19.20, 16),
         linspace(19.54, 19.69, 16)),
    :y=>(linspace(18.10, 18.20, 11), linspace(19.04, 19.20, 17),
         linspace(19.52, 19.71, 20)),
    :z=>(linspace(18.49, 18.545, 12), linspace(18.57, 18.635, 14),
         linspace(18.66, 18.71, 11), linspace(18.75, 18.80, 11),
         linspace(18.83, 18.89, 13), linspace(18.91, 18.97, 13),
         linspace(18.99, 19.06, 15)),
    :x0=>linspace(18.50, 18.70, 11),
    :y0=>linspace(18.50, 18.70, 11),
)

const split_a = NaCsData.split_data(data_a, spec)

const data_x = split_a[:x]
const data_y = split_a[:y]
const data_z = split_a[:z]
const data_x0 = split_a[:x0]
const data_y0 = split_a[:y0]

const prefix = joinpath(@__DIR__, "imgs", "data_20180121_145134_spectrum_cold")

figure()
NaCsPlot.plot_survival_data(data_x[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_x[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_x[3], fmt="C0.-")
grid()
ylim([0, 1])
title("X spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x")

figure()
NaCsPlot.plot_survival_data(data_y[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_y[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_y[3], fmt="C0.-")
grid()
ylim([0, 1])
title("Y spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y")

figure()
NaCsPlot.plot_survival_data(data_z[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_z[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_z[3], fmt="C0.-")
NaCsPlot.plot_survival_data(data_z[4], fmt="C0.-")
NaCsPlot.plot_survival_data(data_z[5], fmt="C0.-")
NaCsPlot.plot_survival_data(data_z[6], fmt="C0.-")
NaCsPlot.plot_survival_data(data_z[7], fmt="C0.-")
grid()
ylim([0, 1])
title("Z spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_z")

figure()
NaCsPlot.plot_survival_data(data_x0, fmt="C0.-")
grid()
ylim([0, 1])
title("X carrier with square pulse")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x0")

figure()
NaCsPlot.plot_survival_data(data_y0, fmt="C0.-")
grid()
ylim([0, 1])
title("Y carrier with square pulse")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y0")

NaCsPlot.maybe_show()
