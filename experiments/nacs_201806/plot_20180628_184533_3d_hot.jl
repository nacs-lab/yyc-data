#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180628_184533.mat")
const params_a, logicals_a = NaCsData.load_striped_mat(iname_a)

function selector(logicals)
    @assert size(logicals, 2) == 1
    if logicals[1, 1] == 0
        return Int[1, 0, 0]
    end
    return Int[1, 1, logicals[3, 1]]
end

data_a = NaCsData.select_count(params_a, logicals_a, selector)

const spec_a = OrderedDict(
    :axial_z=>(-96 .+ 2 .* 30 .* linspace(-1, 1, 31), 75 .+ 2 .* 30 .* linspace(-1, 1, 31)),
    :radial_x=>(-481 .+ 2 .* 45 .* linspace(-1, 1, 31), 489 .+ 2 .* 90 .* linspace(-1, 1, 31)),
    :radial_y=>(-499 .+ 2 .* 45 .* linspace(-1, 1, 31),  528 .+ 2 .* 90 .* linspace(-1, 1, 31)),
)

const split_a = NaCsData.split_data(data_a, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20180628_184533")

data_rx = split_a[:radial_x]
data_ry = split_a[:radial_y]
data_az = split_a[:axial_z]

figure()
NaCsPlot.plot_survival_data(data_rx[1], fmt="C0o-") # -421
NaCsPlot.plot_survival_data(data_rx[2], fmt="C0o-") # 579
grid()
ylim([0, 1])
title("Axis X (radial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

figure()
NaCsPlot.plot_survival_data(data_ry[1], fmt="C0o-") # -430
NaCsPlot.plot_survival_data(data_ry[2], fmt="C0o-") # 588
grid()
ylim([0, 1])
title("Axis Y (radial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ry")

figure()
NaCsPlot.plot_survival_data(data_az[1], fmt="C0o-") # -109
NaCsPlot.plot_survival_data(data_az[2], fmt="C0o-") # 64
grid()
ylim([0, 1])
title("Axis Z (axial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az")

NaCsPlot.maybe_show()
