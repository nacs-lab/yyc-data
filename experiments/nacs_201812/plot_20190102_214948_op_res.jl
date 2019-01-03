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

const inames = ["data_20190102_214948.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_cs = [-205.0916 .+ (0:2:40)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_cs = datas_cs[1]

const prefix = joinpath(@__DIR__, "imgs", "data_20190102_214948_op_res")

figure()
NaCsPlot.plot_survival_data(data_cs, fmt="C0.-")
grid()
ylim([0, 0.6])
title("Cs F4 OP")
xlabel("Detuning (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
