#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180706_000629.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0 || logicals[2, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1] != 0 && logicals[4, 1] != 0]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)

const spec_a = OrderedDict(
    :with_cool=>[0, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200],
    :without_cool=>[0, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200],
)

const split_a = NaCsData.split_data(data_a, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20180706_000629_pa_times")

data_with_cool = split_a[:with_cool]
data_without_cool = split_a[:without_cool]

figure()
NaCsPlot.plot_survival_data(data_with_cool, fmt="C0.-", label="With cooling")
NaCsPlot.plot_survival_data(data_without_cool, fmt="C1.-", label="Without cooling")
grid()
ylim([0, 1])
xlim([-5, 210])
legend()
xlabel("PA time (ms)")
ylabel("Two body survival")
NaCsPlot.maybe_save("$(prefix)")

figure()
NaCsPlot.plot_survival_data(data_with_cool, fmt="C0.-", label="With cooling")
NaCsPlot.plot_survival_data(data_without_cool, fmt="C1.-", label="Without cooling")
gca()[:set_xscale]("log", nonposx="clip")
grid()
ylim([0, 1])
xticks([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
legend()
xlabel("PA time (ms)")
ylabel("Two body survival")
NaCsPlot.maybe_save("$(prefix)_log")

NaCsPlot.maybe_show()
