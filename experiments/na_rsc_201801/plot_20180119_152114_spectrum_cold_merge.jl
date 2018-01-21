#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180119_152114.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)

# 18.157, 19.140
# 18.140, 19.160
# 18.513, 18.602, 18.684, 18.771

const spec = OrderedDict(
    :x_merge=>(linspace(18.04, 18.26, 18), linspace(18.98, 19.24, 14),
         linspace(19.48, 19.74, 18)),
    :y_merge=>(linspace(18.04, 18.26, 18), linspace(18.98, 19.24, 14),
         linspace(19.48, 19.74, 18)),
    :z_merge=>linspace(18.46, 18.80, 50),
    :x_nomerge=>(linspace(18.04, 18.26, 18), linspace(18.98, 19.24, 14),
         linspace(19.48, 19.74, 18)),
    :y_nomerge=>(linspace(18.04, 18.26, 18), linspace(18.98, 19.24, 14),
         linspace(19.48, 19.74, 18)),
    :z_nomerge=>linspace(18.46, 18.80, 50),
)

const split_a = NaCsData.split_data(data_a, spec)

const data_x_merge = split_a[:x_merge]
const data_y_merge = split_a[:y_merge]
const data_z_merge = split_a[:z_merge]
const data_x_nomerge = split_a[:x_nomerge]
const data_y_nomerge = split_a[:y_nomerge]
const data_z_nomerge = split_a[:z_nomerge]

const prefix = joinpath(@__DIR__, "imgs", "data_20180119_152114_spectrum_cold_merge")

figure()
NaCsPlot.plot_survival_data(data_x_merge[1], fmt="C1o-", label="With merge")
NaCsPlot.plot_survival_data(data_x_merge[2], fmt="C1o-")
NaCsPlot.plot_survival_data(data_x_merge[3], fmt="C1o-")
NaCsPlot.plot_survival_data(data_x_nomerge[1], fmt="C0o-", label="Without merge")
NaCsPlot.plot_survival_data(data_x_nomerge[2], fmt="C0o-")
NaCsPlot.plot_survival_data(data_x_nomerge[3], fmt="C0o-")
grid()
ylim([0, 1])
legend()
title("X spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x")

figure()
NaCsPlot.plot_survival_data(data_y_merge[1], fmt="C1o-", label="With merge")
NaCsPlot.plot_survival_data(data_y_merge[2], fmt="C1o-")
NaCsPlot.plot_survival_data(data_y_merge[3], fmt="C1o-")
NaCsPlot.plot_survival_data(data_y_nomerge[1], fmt="C0o-", label="Without merge")
NaCsPlot.plot_survival_data(data_y_nomerge[2], fmt="C0o-")
NaCsPlot.plot_survival_data(data_y_nomerge[3], fmt="C0o-")
grid()
ylim([0, 1])
legend()
title("Y spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y")

figure()
NaCsPlot.plot_survival_data(data_z_merge, fmt="C1o-", label="With merge")
NaCsPlot.plot_survival_data(data_z_nomerge, fmt="C0o-", label="Without merge")
grid()
ylim([0, 1])
legend()
title("Z spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_z")

NaCsPlot.maybe_show()
