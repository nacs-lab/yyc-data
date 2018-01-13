#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180112_233017.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180113_084031.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)
const iname_c = joinpath(@__DIR__, "data", "data_20180113_114027.mat")
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

const spec_a = OrderedDict(
    :p1_3=>linspace(18.4, 19.1, 141),
    :p2_3=>linspace(18.4, 19.1, 141),
)
const spec_b = OrderedDict(
    :p1_3=>linspace(18.46, 18.80, 69),
    :p2_3=>linspace(18.46, 18.80, 69),
)
const spec_c = OrderedDict(
    :p3_3=>linspace(18.46, 18.80, 69),
)

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)
const split_c = NaCsData.split_data(data_c, spec_c)

const data_hot_1_3 = split_a[:p1_3]
const data_hot_2_3 = split_a[:p2_3]
const data_cold_1_3 = split_b[:p1_3]
const data_cold_2_3 = split_b[:p2_3]
const data_cold_3_3 = split_c[:p3_3]

const prefix = joinpath(@__DIR__, "imgs", "data_20180112_spectra_a_compare")

figure()
NaCsPlot.plot_survival_data(data_hot_1_3, fmt="o-", label="1/3 Power, Hot")
NaCsPlot.plot_survival_data(data_hot_2_3, fmt="o-", label="2/3 Power, Hot")
NaCsPlot.plot_survival_data(data_cold_1_3, fmt="o-", label="1/3 Power, Cold")
NaCsPlot.plot_survival_data(data_cold_2_3, fmt="o-", label="2/3 Power, Cold")
NaCsPlot.plot_survival_data(data_cold_3_3, fmt="o-", label="Full Power, Cold")
grid()
ylim([0, 0.85])
legend()
title("Z spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

figure()
NaCsPlot.plot_survival_data(data_cold_1_3, fmt="o-", label="1/3 Power")
NaCsPlot.plot_survival_data(data_cold_2_3, fmt="o-", label="2/3 Power")
NaCsPlot.plot_survival_data(data_cold_3_3, fmt="o-", label="Full Power")
grid()
ylim([0, 0.85])
legend()
title("Z spectrum (Cold)")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cold")

NaCsPlot.maybe_show()
