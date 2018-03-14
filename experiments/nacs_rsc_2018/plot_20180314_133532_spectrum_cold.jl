#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180314_133532.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function gen_selector(na::Bool)
    function selector(logicals)
        @assert size(logicals, 2) == 1
        if logicals[2 - na, 1] == 0
            return Int[1, 0, 0]
        end
        return Int[1, 1, logicals[4 - na, 1]]
    end
end

data_na = NaCsData.select_count(params_a, logicals_a, gen_selector(true))
data_cs = NaCsData.select_count(params_a, logicals_a, gen_selector(false))

const spec_na = OrderedDict(
    :z=>(linspace(18.49, 18.545, 16), linspace(18.66, 18.71, 16)),
    :x=>(linspace(18.09, 18.18, 16), linspace(19.00, 19.20, 16)),
    :y=>(linspace(18.05, 18.16, 16), linspace(19.04, 19.20, 16)),
)
const spec_cs = OrderedDict(
    :z=>(linspace(5, 25, 16), linspace(55, 75, 16)),
    :x=>(linspace(-160, -70, 16), linspace(110, 190, 16)),
    :y=>(linspace(-160, -70, 16), linspace(110, 190, 16)),
)

const split_na = NaCsData.split_data(data_na, spec_na)
const split_cs = NaCsData.split_data(data_cs, spec_cs)

const data_na_z = split_na[:z]
const data_na_x = split_na[:x]
const data_na_y = split_na[:y]

const data_cs_z = split_cs[:z]
const data_cs_x = split_cs[:x]
const data_cs_y = split_cs[:y]

const prefix = joinpath(@__DIR__, "imgs", "data_20180314_133532_spectrum_cold")

figure()
NaCsPlot.plot_survival_data(data_na_z[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_z[2], fmt="C0.-")
grid()
ylim([0, 0.9])
title("Na Z spectrum")
xlabel("Detuning (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_z")

figure()
NaCsPlot.plot_survival_data(data_na_x[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_x[2], fmt="C0.-")
grid()
ylim([0, 0.9])
title("Na X spectrum")
xlabel("Detuning (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_x")

figure()
NaCsPlot.plot_survival_data(data_na_y[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_y[2], fmt="C0.-")
grid()
ylim([0, 0.9])
title("Na Y spectrum")
xlabel("Detuning (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_y")

figure()
NaCsPlot.plot_survival_data(data_cs_z[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_z[2], fmt="C0.-")
grid()
ylim([0, 0.9])
title("Cs Z spectrum")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_z")

figure()
NaCsPlot.plot_survival_data(data_cs_x[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_x[2], fmt="C0.-")
grid()
ylim([0, 0.9])
title("Cs X spectrum")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_x")

figure()
NaCsPlot.plot_survival_data(data_cs_y[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_y[2], fmt="C0.-")
grid()
ylim([0, 0.9])
title("Cs Y spectrum")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_y")

NaCsPlot.maybe_show()
