#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20180628_225959.mat")
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
    :axial_z=>(-109 .+ 30 .* linspace(-1, 1, 16), 64 .+ 30 .* linspace(-1, 1, 16)),
    :radial_x=>(-421 .+ 45 .* linspace(-1, 1, 16), 579 .+ 90 .* linspace(-1, 1, 16)),
    :radial_y=>(-430 .+ 45 .* linspace(-1, 1, 16),  588 .+ 90 .* linspace(-1, 1, 16)),
    :axial_z2=>(-109 .+ 30 .* linspace(-1, 1, 16), 64 .+ 30 .* linspace(-1, 1, 16)),
    :radial_x2=>(-421 .+ 45 .* linspace(-1, 1, 16), 579 .+ 90 .* linspace(-1, 1, 16)),
    :radial_y2=>(-430 .+ 45 .* linspace(-1, 1, 16),  588 .+ 90 .* linspace(-1, 1, 16)),
)

const split_a = NaCsData.split_data(data_a, spec_a)

const prefix = joinpath(@__DIR__, "imgs", "data_20180628_225959_cold1")

data_rx = split_a[:radial_x]
data_ry = split_a[:radial_y]
data_az = split_a[:axial_z]

data_rx2 = split_a[:radial_x2]
data_ry2 = split_a[:radial_y2]
data_az2 = split_a[:axial_z2]

figure()
NaCsPlot.plot_survival_data(data_rx[1], fmt="C0.-", label="1")
NaCsPlot.plot_survival_data(data_rx[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_rx2[1], fmt="C1.-", label="2")
NaCsPlot.plot_survival_data(data_rx2[2], fmt="C1.-")
grid()
ylim([0, 1])
legend()
title("Axis X (radial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

figure()
NaCsPlot.plot_survival_data(data_ry[1], fmt="C0.-", label="1")
NaCsPlot.plot_survival_data(data_ry[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_ry2[1], fmt="C1.-", label="2")
NaCsPlot.plot_survival_data(data_ry2[2], fmt="C1.-")
grid()
ylim([0, 1])
legend()
title("Axis Y (radial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ry")

figure()
NaCsPlot.plot_survival_data(data_az[1], fmt="C0.-", label="1")
NaCsPlot.plot_survival_data(data_az[2], fmt="C0.-")
NaCsPlot.plot_survival_data(data_az2[1], fmt="C1.-", label="2")
NaCsPlot.plot_survival_data(data_az2[2], fmt="C1.-")
grid()
ylim([0, 1])
legend()
title("Axis Z (axial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az")

function freq_arrow(x, y, len, color)
    xl = xlim()
    arrow(x, y, 0, -len, color=color,
          width=(xl[2] - xl[1]) / 80, length_includes_head=true, head_length=len / 3)
end

figure()
NaCsPlot.plot_survival_data(data_rx[2], fmt="C0.-", label="1")
NaCsPlot.plot_survival_data(data_rx2[2], fmt="C1.-", label="2")
freq_arrow(579, 0.05, 0.006, "C0")
freq_arrow(591, 0.05, 0.006, "C1")
grid()
ylim([0, 0.06])
legend()
title("Axis X (radial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx_-1")

figure()
NaCsPlot.plot_survival_data(data_ry[2], fmt="C0.-", label="1")
NaCsPlot.plot_survival_data(data_ry2[2], fmt="C1.-", label="2")
freq_arrow(588, 0.04, 0.006, "C0")
freq_arrow(605, 0.04, 0.006, "C1")
grid()
ylim([0, 0.06])
legend()
title("Axis Y (radial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ry_-1")

figure()
NaCsPlot.plot_survival_data(data_az[2], fmt="C0.-", label="1")
NaCsPlot.plot_survival_data(data_az2[2], fmt="C1.-", label="2")
freq_arrow(71, 0.04, 0.006, "C0")
freq_arrow(62, 0.04, 0.006, "C1")
grid()
ylim([0, 0.06])
legend()
title("Axis Z (axial)")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az_-1")

NaCsPlot.maybe_show()
