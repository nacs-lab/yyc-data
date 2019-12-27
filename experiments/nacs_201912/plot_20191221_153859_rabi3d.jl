#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20191221_153859.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [OrderedDict(
    :az=>0:20.0:400,
    :rx=>0:10.0:200,
    :ry=>0:10.0:200
)]
const specs_cs = [OrderedDict(
    :az=>0:60.0:1200,
    :rx=>0:40.0:800,
    :ry=>0:40.0:800
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

const prefix = joinpath(@__DIR__, "imgs", "data_20191221_153859_rabi")

figure(figsize=[11.2, 8.4])
subplot(2, 2, 1)
NaCsPlot.plot_survival_data(data_na_az, fmt="C0.-")
grid()
xlim([0, xlim()[2]])
ylim([0, ylim()[2]])
title("Na Axial")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")

subplot(2, 2, 2)
NaCsPlot.plot_survival_data(data_na_rx, fmt="C0.-")
grid()
xlim([0, xlim()[2]])
ylim([0, ylim()[2]])
title("Na X")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")

subplot(2, 2, 4)
NaCsPlot.plot_survival_data(data_na_ry, fmt="C0.-")
grid()
xlim([0, xlim()[2]])
ylim([0, ylim()[2]])
title("Na Y")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_na")

figure(figsize=[11.2, 8.4])
subplot(2, 2, 1)
NaCsPlot.plot_survival_data(data_cs_az, fmt="C0.-")
grid()
xlim([0, xlim()[2]])
ylim([0, ylim()[2]])
title("Cs Axial")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")

subplot(2, 2, 2)
NaCsPlot.plot_survival_data(data_cs_rx, fmt="C0.-")
grid()
xlim([0, xlim()[2]])
ylim([0, ylim()[2]])
title("Cs X")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")

subplot(2, 2, 4)
NaCsPlot.plot_survival_data(data_cs_ry, fmt="C0.-")
grid()
xlim([0, xlim()[2]])
ylim([0, ylim()[2]])
title("Cs Y")
xlabel("Time (\$\\mu\$s)")
ylabel("Survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
