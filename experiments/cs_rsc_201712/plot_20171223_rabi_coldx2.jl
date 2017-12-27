#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20171223_235301.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20171224_183640.mat")
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
    :z62=>linspace(0, 800, 26),
    :z50=>linspace(32, 800, 25),
    :x62=>linspace(16, 400, 25),
    :x50=>linspace(16, 400, 25),
    :y62=>linspace(16, 400, 25),
    :y50=>linspace(16, 400, 25),
)

const split_50 = NaCsData.split_data([data_50_a; data_50_b], spec)
const split_62 = NaCsData.split_data([data_62_a; data_62_b], spec)

const data_50_z62 = split_50[:z62]
const data_50_z50 = [data_50_z62[1]; split_50[:z50]]
const data_50_x62 = [data_50_z62[1]; split_50[:x62]]
const data_50_x50 = [data_50_z62[1]; split_50[:x50]]
const data_50_y62 = [data_50_z62[1]; split_50[:y62]]
const data_50_y50 = [data_50_z62[1]; split_50[:y50]]

const data_62_z62 = split_62[:z62]
const data_62_z50 = [data_62_z62[1]; split_62[:z50]]
const data_62_x62 = [data_62_z62[1]; split_62[:x62]]
const data_62_x50 = [data_62_z62[1]; split_62[:x50]]
const data_62_y62 = [data_62_z62[1]; split_62[:y62]]
const data_62_y50 = [data_62_z62[1]; split_62[:y50]]

const prefix = joinpath(@__DIR__, "imgs", "data_20171223")

figure()
NaCsPlot.plot_survival_data(data_50_z62, fmt="C0o-", label="50MHz")
NaCsPlot.plot_survival_data(data_62_z62, fmt="C1o-", label="62MHz")
grid()
xlim([0, 810])
ylim([0, 0.5])
legend()
title("Z heating for site 62MHz")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_z62")

figure()
NaCsPlot.plot_survival_data(data_50_z50, fmt="C0o-", label="50MHz")
NaCsPlot.plot_survival_data(data_62_z50, fmt="C1o-", label="62MHz")
grid()
xlim([0, 810])
ylim([0, 0.5])
legend()
title("Z heating for site 50MHz")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_z50")

figure()
NaCsPlot.plot_survival_data(data_50_x62, fmt="C0o-", label="50MHz")
NaCsPlot.plot_survival_data(data_62_x62, fmt="C1o-", label="62MHz")
grid()
xlim([0, 410])
ylim([0, 1])
legend()
title("X heating for site 62MHz")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x62")

figure()
NaCsPlot.plot_survival_data(data_50_x50, fmt="C0o-", label="50MHz")
NaCsPlot.plot_survival_data(data_62_x50, fmt="C1o-", label="62MHz")
grid()
xlim([0, 410])
ylim([0, 1])
legend()
title("X heating for site 50MHz")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_x50")

figure()
NaCsPlot.plot_survival_data(data_50_y62, fmt="C0o-", label="50MHz")
NaCsPlot.plot_survival_data(data_62_y62, fmt="C1o-", label="62MHz")
grid()
xlim([0, 410])
ylim([0, 1])
legend()
title("Y heating for site 62MHz")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y62")

figure()
NaCsPlot.plot_survival_data(data_50_y50, fmt="C0o-", label="50MHz")
NaCsPlot.plot_survival_data(data_62_y50, fmt="C1o-", label="62MHz")
grid()
xlim([0, 410])
ylim([0, 1])
legend()
title("Y heating for site 50MHz")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_y50")

NaCsPlot.maybe_show()
