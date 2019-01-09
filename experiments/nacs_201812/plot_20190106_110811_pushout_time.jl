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

const inames = ["data_20190106_110811.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [[0, 0.5, 1, 2, 5, 10, 20]]
const specs_cs = [[0, 0.5, 1, 2, 5, 10, 20]]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_na = datas_na[1]
data_cs = datas_cs[1]

const prefix = joinpath(@__DIR__, "imgs", "data_20190106_110811_pushout_time")

figure()
NaCsPlot.plot_survival_data(data_na, fmt="C0.-", label="Na")
NaCsPlot.plot_survival_data(data_cs, fmt="C1.-", label="Cs")
legend()
grid()
ylim([0, 1])
title("Pushout")
xlabel("Time (\$\\mu s\$)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
