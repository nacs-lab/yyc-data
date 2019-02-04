#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190126_223602.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [OrderedDict(
    :az=>0:10.0:400,
    :rx=>0:5.0:200,
    :ry=>0:5.0:200
)]
const specs_cs = [OrderedDict(
    :az=>0:30.0:1200,
    :rx=>0:20.0:800,
    :ry=>0:20.0:800
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

const prefix = joinpath(@__DIR__, "imgs", "data_20190126_223602_rabi")

figure()
NaCsPlot.plot_survival_data(data_na_az, fmt="C0.-")
grid()
ylim([0, 1])
title("Na Axial")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_az")

figure()
NaCsPlot.plot_survival_data(data_cs_az, fmt="C0.-")
grid()
ylim([0, 1])
title("Cs Axial")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_az")

figure()
NaCsPlot.plot_survival_data(data_na_rx, fmt="C0.-")
grid()
ylim([0, 1])
title("Na X")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_rx")

figure()
NaCsPlot.plot_survival_data(data_cs_rx, fmt="C0.-")
grid()
ylim([0, 1])
title("Cs X")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_rx")

figure()
NaCsPlot.plot_survival_data(data_na_ry, fmt="C0.-")
grid()
ylim([0, 1])
title("Na Y")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na_ry")

figure()
NaCsPlot.plot_survival_data(data_cs_ry, fmt="C0.-")
grid()
ylim([0, 1])
title("Cs Y")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs_ry")

NaCsPlot.maybe_show()
