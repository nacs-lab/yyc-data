#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180711_142041.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0 || logicals[2, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1] != 0 && logicals[4, 1] != 0]
end

@show length(params_a)

data_a5 = NaCsData.select_count(params_a, logicals_a, selector, 5000)
data_a10 = NaCsData.select_count(params_a, logicals_a, selector, 10000)
data_a = NaCsData.select_count(params_a, logicals_a, selector)

const spec = OrderedDict(
    :p16_5=>93 .+ 25 .* linspace(-1, 1, 26),
    :p12=>93 .+ 25 .* linspace(-1, 1, 26),
)

const split5 = NaCsData.split_data(data_a5, spec)
const split10 = NaCsData.split_data(data_a10, spec)
const split = NaCsData.split_data(data_a, spec)

const prefix = joinpath(@__DIR__, "imgs", "data_20180711_pa_raman_shift")

data5_p16_5 = split5[:p16_5]
data5_p12 = split5[:p12]
data10_p16_5 = split10[:p16_5]
data10_p12 = split10[:p12]
data_p16_5 = split[:p16_5]
data_p12 = split[:p12]


figure()
NaCsPlot.plot_survival_data(data5_p16_5, fmt="C0.-", label="5k")
NaCsPlot.plot_survival_data(data10_p16_5, fmt="C1.-", label="10k")
NaCsPlot.plot_survival_data(data_p16_5, fmt="C2.-", label="full")
grid()
ylim([0, 1])
legend()
xlabel("Raman frequency (298XXX) kHz")
ylabel("Two body survival")
NaCsPlot.maybe_save("$(prefix)")

figure()
NaCsPlot.plot_survival_data(data_p16_5, fmt="C0.-", label="16.5mW")
NaCsPlot.plot_survival_data(data_p12, fmt="C1.-", label="12mW")
grid()
ylim([0, 1])
legend()
xlabel("Raman frequency (298XXX) kHz")
ylabel("Two body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
