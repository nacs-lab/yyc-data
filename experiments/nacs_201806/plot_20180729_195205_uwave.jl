#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180729_195205.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0 || logicals[2, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1] != 0 && logicals[4, 1] != 0]
end

params_a .= params_a ./ 1e3 .+ 298e3

split_a = NaCsData.select_count(params_a, logicals_a, selector)

const prefix = joinpath(@__DIR__, "imgs", "data_20180729_195205_uwave")

figure()
NaCsPlot.plot_survival_data(split_a, fmt="o-")
grid()
ylim([0.63, 0.82])
xlabel("Detuning from -298 MHz (kHz)")
ylabel("Two body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
