#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures

const iname_a = joinpath(@__DIR__, "data", "data_20170926_193126.csv")

const data_a = NaCsData.load_count_csv(iname_a)

const spec_ap = OrderedDict(
    :radial_x=>(linspace(18, 18.4, 21), linspace(19, 19.4, 21), linspace(19.5, 19.8, 16)),
    :radial_y=>(linspace(18, 18.4, 21), linspace(19, 19.4, 21), linspace(19.5, 19.8, 16)),
    :axial_z=>linspace(18.5, 18.85, 71),
    :survival=>[0]
)
const spec_am = OrderedDict(
    :radial_x=>(linspace(18, 18.4, 21), linspace(19, 19.4, 21), linspace(19.5, 19.8, 16)),
    :radial_y=>(linspace(18, 18.4, 21), linspace(19, 19.4, 21), linspace(19.5, 19.8, 16)),
    :axial_z=>linspace(18.5, 18.85, 71),
    :survival=>[0]
)

const data_ap = data_a[data_a.params .>= 0]
const data_am = NaCsData.map_params((i, v)->-v, data_a[data_a.params .<= 0])

const split_ap = NaCsData.split_data(data_ap, spec_ap)
const split_am = NaCsData.split_data(data_am, spec_am)

const prefix = joinpath(@__DIR__, "imgs", "data_20170926_193126_cool_merge")

to_sideband(f) = (i, v)->(v - f) * 1000

nomerge_rx = NaCsData.map_params(to_sideband(18.705), split_ap[:radial_x])
nomerge_ry = NaCsData.map_params(to_sideband(18.713), split_ap[:radial_y])
nomerge_az = NaCsData.map_params(to_sideband(18.619), split_ap[:axial_z])
merge_rx = NaCsData.map_params(to_sideband(18.705), split_am[:radial_x])
merge_ry = NaCsData.map_params(to_sideband(18.713), split_am[:radial_y])
merge_az = NaCsData.map_params(to_sideband(18.619), split_am[:axial_z])

# nomerge_rx = split_ap[:radial_x]
# nomerge_ry = split_ap[:radial_y]
# nomerge_az = split_ap[:axial_z]
# merge_rx = split_am[:radial_x]
# merge_ry = split_am[:radial_y]
# merge_az = split_am[:axial_z]

@show split_ap[:survival]
@show split_am[:survival]

figure()
NaCsPlot.plot_survival_data(nomerge_rx[1], fmt="C0.-", label="No merge")
NaCsPlot.plot_survival_data(nomerge_rx[2], fmt="C0.-")
NaCsPlot.plot_survival_data(nomerge_rx[3], fmt="C0.-")
NaCsPlot.plot_survival_data(merge_rx[1], fmt="C1.-", label="Merge")
NaCsPlot.plot_survival_data(merge_rx[2], fmt="C1.-")
NaCsPlot.plot_survival_data(merge_rx[3], fmt="C1.-")
axhline(NaCsData.get_values(split_ap[:survival])[2][2], color="C0", ls="--")
axhline(NaCsData.get_values(split_am[:survival])[2][2], color="C1", ls="--")
grid()
ylim([0, 1])
legend()
title("Axis X (radial)")
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_rx")

figure()
NaCsPlot.plot_survival_data(nomerge_ry[1], fmt="C0.-", label="No merge")
NaCsPlot.plot_survival_data(nomerge_ry[2], fmt="C0.-")
NaCsPlot.plot_survival_data(nomerge_ry[3], fmt="C0.-")
NaCsPlot.plot_survival_data(merge_ry[1], fmt="C1.-", label="Merge")
NaCsPlot.plot_survival_data(merge_ry[2], fmt="C1.-")
NaCsPlot.plot_survival_data(merge_ry[3], fmt="C1.-")
axhline(NaCsData.get_values(split_ap[:survival])[2][2], color="C0", ls="--")
axhline(NaCsData.get_values(split_am[:survival])[2][2], color="C1", ls="--")
grid()
ylim([0, 1])
legend()
title("Axis Y (radial)")
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_ry")

figure()
NaCsPlot.plot_survival_data(nomerge_az, fmt="C0.-", label="No merge")
NaCsPlot.plot_survival_data(merge_az, fmt="C1.-", label="Merge")
axhline(NaCsData.get_values(split_ap[:survival])[2][2], color="C0", ls="--")
axhline(NaCsData.get_values(split_am[:survival])[2][2], color="C1", ls="--")
grid()
ylim([0, 1])
legend()
title("Axis Z (axial)")
xlabel("Detuning from carrier (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_az")

NaCsPlot.maybe_show()
