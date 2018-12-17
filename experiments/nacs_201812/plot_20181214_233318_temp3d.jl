#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181214_233318.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [OrderedDict(
    :az=>(-124 .+ (-30:4.0:30), 72 .+ (-30:4.0:30)),
    :rx=>(-483 .+ (-45:6.0:45), 579 .+ (-90:12.0:90)),
    :ry=>(-492 .+ (-45:6.0:45), 589 .+ (-90:12.0:90))
)]
const specs_cs = [OrderedDict(
    :az=>(-4 .+ (-15:2.0:15), 43 .+ (-15:2.0:15)),
    :rx=>(-124 .+ (-45:6.0:45), 104 .+ (-75:10.0:75)),
    :ry=>(-116 .+ (-45:6.0:45), 124 .+ (-60:8.0:60))
)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]
const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_na_az = datas_na[1][:az]
data_na_rx = datas_na[1][:rx]
data_na_ry = datas_na[1][:ry]
data_cs_az = datas_cs[1][:az]
data_cs_rx = datas_cs[1][:rx]
data_cs_ry = datas_cs[1][:ry]

const prefix = joinpath(@__DIR__, "imgs", "data_20181214_233318_temp")

figure()
NaCsPlot.plot_survival_data(data_na_az[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_az[2], fmt="C0.-")
grid()
ylim([0, 1])
title("Na Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
NaCsPlot.plot_survival_data(data_cs_az[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_az[2], fmt="C0.-")
grid()
ylim([0, 1])
title("Cs Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
NaCsPlot.plot_survival_data(data_na_rx[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_rx[2], fmt="C0.-")
grid()
ylim([0, 1])
title("Na X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_rx")

figure()
NaCsPlot.plot_survival_data(data_cs_rx[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_rx[2], fmt="C0.-")
grid()
ylim([0, 1])
title("Cs X")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_rx")

figure()
NaCsPlot.plot_survival_data(data_na_ry[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_na_ry[2], fmt="C0.-")
grid()
ylim([0, 1])
title("Na Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_ry")

figure()
NaCsPlot.plot_survival_data(data_cs_ry[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_ry[2], fmt="C0.-")
grid()
ylim([0, 1])
title("Cs Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ry")

NaCsPlot.maybe_show()
