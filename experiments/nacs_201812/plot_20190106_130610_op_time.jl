#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit
using NaCsCalc.Atomic: all_scatter_D

const inames = ["data_20190106_130610.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [([0, 0.05, 0.1, 0.2, 0.5, 1, 2, 4] .* 0.1,
                   [0, 0.05, 0.1, 0.2, 0.5, 1, 2, 4] .* 0.1,
                   [0, 0.05, 0.1, 0.2, 0.5, 1, 2, 4] .* 0.1)]
const specs_cs = [([0, 0.05, 0.1, 0.2, 0.5, 1, 2, 4],
                   [0, 0.05, 0.1, 0.2, 0.5, 1, 2, 4],
                   [0, 0.05, 0.1, 0.2, 0.5, 1, 2, 4])]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_na = datas_na[1]
data_cs = datas_cs[1]

const prefix = joinpath(@__DIR__, "imgs", "data_20190106_130610_op_time")

figure()
NaCsPlot.plot_survival_data(data_na[1], fmt="C0.-", label="Total")
NaCsPlot.plot_survival_data(data_na[2], fmt="C1.-", label="F=1")
NaCsPlot.plot_survival_data(data_na[3], fmt="C2.-", label="2, 2")
legend()
grid()
ylim([0, 1])
title("OP survival (Na)")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

figure()
NaCsPlot.plot_survival_data(data_cs[1], fmt="C0.-", label="Total")
NaCsPlot.plot_survival_data(data_cs[2], fmt="C1.-", label="F=3")
NaCsPlot.plot_survival_data(data_cs[3], fmt="C2.-", label="4, 4")
legend()
grid()
ylim([0, 1])
title("OP survival (Cs)")
xlabel("Time (\$ms\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
