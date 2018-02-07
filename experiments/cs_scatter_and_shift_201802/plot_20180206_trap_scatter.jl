#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
import NaCsCalc.Format: Unc, Sci

const iname_a = joinpath(@__DIR__, "data", "data_20180206_150235.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)
const iname_b = joinpath(@__DIR__, "data", "data_20180206_164731.mat")
const params_b, logicals_b = NaCsData.load_striped_mat(iname_b)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[2, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[4, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)
data_b = NaCsData.select_count(params_b, logicals_b, selector)

const spec_a = OrderedDict(
    :f3=>[0.001, 50, 100, 200, 400, 800],
    :f44=>[0.001, 50, 100, 200, 400, 800],
    :all=>[0.001, 50, 100, 200, 400, 800],
)

const spec_b = OrderedDict(
    :f3=>[1600.0, 2400.0],
    :f44=>[1600.0, 2400.0],
    :all=>[1600.0, 2400.0],
)

const split_a = NaCsData.split_data(data_a, spec_a)
const split_b = NaCsData.split_data(data_b, spec_b)

const data_f3 = [split_a[:f3]; split_b[:f3]]
const data_f44 = [split_a[:f44]; split_b[:f44]]
const data_all = [split_a[:all]; split_b[:all]]

const prefix = joinpath(@__DIR__, "imgs", "data_20180206_trap_scatter")

function get_survival_data(data)
    params, ratios, uncs = NaCsData.get_values(data)
    perm = sortperm(params)
    params = params[perm]
    return params, ratios[perm, 2], uncs[perm, 2]
end

function get_survival_ratios(base, sub)
    params_b, ratios_b, uncs_b = get_survival_data(base)
    params_s, ratios_s, uncs_s = get_survival_data(sub)
    params_b, ratios_s ./ ratios_b, sqrt.((uncs_s ./ ratios_b).^2 .+
                                          (uncs_b ./ ratios_b.^2 .* ratios_s).^2)
end

function plot_survival_ratios(base, sub; kws...)
    params, ratios, uncs = get_survival_ratios(base, sub)
    errorbar(params, ratios, uncs; kws...)
end

figure()
NaCsPlot.plot_survival_data(data_f3, fmt="o-", label="F=3")
NaCsPlot.plot_survival_data(data_f44, fmt="o-", label="4, 4")
NaCsPlot.plot_survival_data(data_all, fmt="o-", label="Total")
grid()
xlim([0, 2500])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_raw")

figure()
plot_survival_ratios(data_all, data_f3, fmt="o-", label="F=3")
plot_survival_ratios(data_all, data_f44, fmt="o-", label="4, 4")
grid()
xlim([0, 2500])
ylim([0, 1])
legend()
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ratio")

NaCsPlot.maybe_show()
