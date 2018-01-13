#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180113_084031.mat")
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
    :p1_3=>linspace(18.46, 18.80, 69),
    :p2_3=>linspace(18.46, 18.80, 69),
)

const split_a = NaCsData.split_data(data_a, spec)

const data_1_3 = split_a[:p1_3]
const data_2_3 = split_a[:p2_3]

const prefix = joinpath(@__DIR__, "imgs", "data_20180113_084031_spectrum_a_cold")

figure()
NaCsPlot.plot_survival_data(data_1_3, fmt="o-", label="1/3 Power")
NaCsPlot.plot_survival_data(data_2_3, fmt="o-", label="2/3 Power")
grid()
ylim([0, 0.7])
legend()
title("Z spectrum")
xlabel("Frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
