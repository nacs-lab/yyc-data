#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180711_080904.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0 || logicals[2, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1] != 0 && logicals[4, 1] != 0]
end

params_a .= params_a .* 1000

split_a2 = NaCsData.select_count(params_a, logicals_a, selector, 2000)
split_a4 = NaCsData.select_count(params_a, logicals_a, selector, 4000)
split_a6 = NaCsData.select_count(params_a, logicals_a, selector, 6000)
split_a8 = NaCsData.select_count(params_a, logicals_a, selector, 8000)
split_a10 = NaCsData.select_count(params_a, logicals_a, selector)

const prefix = joinpath(@__DIR__, "imgs", "data_20180711_080904_pa_raman_times")

figure()
NaCsPlot.plot_survival_data(split_a2, fmt="o-", label="2k")
NaCsPlot.plot_survival_data(split_a4, fmt="o-", label="4k")
NaCsPlot.plot_survival_data(split_a6, fmt="o-", label="6k")
NaCsPlot.plot_survival_data(split_a8, fmt="o-", label="8k")
NaCsPlot.plot_survival_data(split_a10, fmt="o-", label="full")
grid()
ylim([0.5, 1])
legend()
xlabel("Raman time (ms)")
ylabel("Two body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
