#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190209_150517.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_cs = [OrderedDict(
    :az=>(-100:2:20),
    :rx=>(-182 .+ (-45:6.0:45), 50 .+ (-75:10.0:75)),
    :ry=>(-190 .+ (-45:6.0:45), 60 .+ (-60:8.0:60))
)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_cs_az = datas_cs[1][:az]
data_cs_rx = datas_cs[1][:rx]
data_cs_ry = datas_cs[1][:ry]

const prefix = joinpath(@__DIR__, "imgs", "data_20190209_150517_temp")

figure()
NaCsPlot.plot_survival_data(data_cs_az, fmt="C0.-")
grid()
ylim([0, 1])
title("Cs Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_az")

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
NaCsPlot.plot_survival_data(data_cs_ry[1], fmt="C0.-")
NaCsPlot.plot_survival_data(data_cs_ry[2], fmt="C0.-")
grid()
ylim([0, 1])
title("Cs Y")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ry")

NaCsPlot.maybe_show()
