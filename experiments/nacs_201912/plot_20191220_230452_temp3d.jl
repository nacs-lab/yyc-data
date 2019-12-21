#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20191220_230452.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [OrderedDict(
    :az=>(-101 .+ (-30:4.0:30), 80 .+ (-30:4.0:30)),
    :rx=>(-502 .+ (-45:6.0:45), 474 .+ (-90:12.0:90)),
    :ry=>(-487 .+ (-45:6.0:45), 507 .+ (-90:12.0:90))
)]
const specs_cs = [OrderedDict(
    :az=>(5 .+ (-15:2.0:15), 40 .+ (-15:2.0:15)),
    :rx=>(-99 .+ (-45:6.0:45), 110 .+ (-75:10.0:75)),
    :ry=>(-99 .+ (-45:6.0:45), 112 .+ (-60:8.0:60))
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

const prefix = joinpath(@__DIR__, "imgs", "data_20191220_230452_temp")

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
