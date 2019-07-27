#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190726_124910.mat",
                "data_20190726_154625.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int), typemax(Int)]
const specs_na = [OrderedDict(
    :az=>(-109 .+ (-30:4.0:30), 68 .+ (-30:4.0:30)),
    :rx=>(-507 .+ (-45:6.0:45), 495 .+ (-90:12.0:90)),
    :ry=>(-496 .+ (-45:6.0:45), 492 .+ (-90:12.0:90))
), OrderedDict(
    :az=>(-109 .+ (-30:4.0:30), 68 .+ (-30:4.0:30)),
    :rx=>(-507 .+ (-45:6.0:45), 495 .+ (-90:12.0:90)),
    :ry=>(-496 .+ (-45:6.0:45), 492 .+ (-90:12.0:90))
)]
const specs_cs = [OrderedDict(
    :az=>(-20 .+ (-15:2.0:15), 23 .+ (-15:2.0:15)),
    :rx=>(-131 .+ (-45:6.0:45), 115 .+ (-75:10.0:75)),
    :ry=>(-134 .+ (-45:6.0:45), 114 .+ (-60:8.0:60))
), OrderedDict(
    :az=>(-20 .+ (-15:2.0:15), 23 .+ (-15:2.0:15)),
    :rx=>(-131 .+ (-45:6.0:45), 115 .+ (-75:10.0:75)),
    :ry=>(-134 .+ (-45:6.0:45), 114 .+ (-60:8.0:60))
)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]
const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_na_az = ([datas_na[1][:az][1]; datas_na[2][:az][1]],
              [datas_na[1][:az][2]; datas_na[2][:az][2]])
data_na_rx = ([datas_na[1][:rx][1]; datas_na[2][:rx][1]],
              [datas_na[1][:rx][2]; datas_na[2][:rx][2]])
data_na_ry = ([datas_na[1][:ry][1]; datas_na[2][:ry][1]],
              [datas_na[1][:ry][2]; datas_na[2][:ry][2]])
data_cs_az = ([datas_cs[1][:az][1]; datas_cs[2][:az][1]],
              [datas_cs[1][:az][2]; datas_cs[2][:az][2]])
data_cs_rx = ([datas_cs[1][:rx][1]; datas_cs[2][:rx][1]],
              [datas_cs[1][:rx][2]; datas_cs[2][:rx][2]])
data_cs_ry = ([datas_cs[1][:ry][1]; datas_cs[2][:ry][1]],
              [datas_cs[1][:ry][2]; datas_cs[2][:ry][2]])

const prefix = joinpath(@__DIR__, "imgs", "data_20190726_124910_temp")

figure(figsize=[11.2, 8.4])
subplot(2, 2, 1)
NaCsPlot.plot_survival_data(data_na_az[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_az[2], fmt="C0.-")
grid()
ylim([0, ylim()[2]])
title("Na Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")

subplot(2, 2, 2)
NaCsPlot.plot_survival_data(data_na_rx[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_rx[2], fmt="C0.-")
grid()
ylim([0, ylim()[2]])
title("Na X")
xlabel("Detuning (kHz)")
ylabel("Survival")

subplot(2, 2, 4)
NaCsPlot.plot_survival_data(data_na_ry[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_ry[2], fmt="C0.-")
grid()
ylim([0, ylim()[2]])
title("Na Y")
xlabel("Detuning (kHz)")
ylabel("Survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_na")

figure(figsize=[11.2, 8.4])
subplot(2, 2, 1)
NaCsPlot.plot_survival_data(data_cs_az[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_az[2], fmt="C0.-")
grid()
ylim([0, ylim()[2]])
title("Cs Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")

subplot(2, 2, 2)
NaCsPlot.plot_survival_data(data_cs_rx[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_rx[2], fmt="C0.-")
grid()
ylim([0, ylim()[2]])
title("Cs X")
xlabel("Detuning (kHz)")
ylabel("Survival")

subplot(2, 2, 4)
NaCsPlot.plot_survival_data(data_cs_ry[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_ry[2], fmt="C0.-")
grid()
ylim([0, ylim()[2]])
title("Cs Y")
xlabel("Detuning (kHz)")
ylabel("Survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
