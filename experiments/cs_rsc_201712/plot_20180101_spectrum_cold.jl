#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180101_155210.mat")
const iname_b = joinpath(@__DIR__, "data", "data_20180101_192255.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)

function gen_selector(site)
    function selector(logicals)
        @assert size(logicals, 2) >= site
        if logicals[1, site] == 0
            return Int[1, 0, 0]
        end
        return Int[1, 1, logicals[2, site]]
    end
end

data_50_a = NaCsData.select_count(params_a, logicals_a, gen_selector(1))
data_62_a = NaCsData.select_count(params_a, logicals_a, gen_selector(2))
data_50_b = NaCsData.select_count(params_b, logicals_b, gen_selector(1))
data_62_b = NaCsData.select_count(params_b, logicals_b, gen_selector(2))

const spec = OrderedDict(
    :z=>linspace(-20, 100, 61),
    :x=>linspace(-180, 220, 51),
    :y=>linspace(-180, 220, 51),
)

const split_50 = NaCsData.split_data([data_50_a; data_50_b], spec)
const split_62 = NaCsData.split_data([data_62_a; data_62_b], spec)

const data_50_z = split_50[:z]
const data_50_x = split_50[:x]
const data_50_y = split_50[:y]

const data_62_z = split_62[:z]
const data_62_x = split_62[:x]
const data_62_y = split_62[:y]

const prefix = joinpath(@__DIR__, "imgs", "data_20180101")

figure()
NaCsPlot.plot_survival_data(data_50_z, fmt="C0o-", label="50MHz")
NaCsPlot.plot_survival_data(data_62_z, fmt="C1o-", label="62MHz")
grid()
ylim([0, 0.7])
legend()
title("Z spectrum")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_z")

figure()
NaCsPlot.plot_survival_data(data_50_x, fmt="C0o-", label="50MHz")
NaCsPlot.plot_survival_data(data_62_x, fmt="C1o-", label="62MHz")
grid()
ylim([0, 0.8])
legend()
title("X spectrum")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x")

figure()
NaCsPlot.plot_survival_data(data_50_y, fmt="C0o-", label="50MHz")
NaCsPlot.plot_survival_data(data_62_y, fmt="C1o-", label="62MHz")
grid()
ylim([0, 0.8])
legend()
title("Y spectrum")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y")

NaCsPlot.maybe_show()
