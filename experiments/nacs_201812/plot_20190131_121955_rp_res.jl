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

const inames = ["data_20190131_121955.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_cs = [(80:2.0:100,
                   80:2.0:100,
                   80:2.0:100)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_cs = datas_cs[1]

const prefix = joinpath(@__DIR__, "imgs", "data_20190131_121955_rp_res")

figure()
NaCsPlot.plot_survival_data(data_cs[1], fmt="C1.-", label="10mW")
NaCsPlot.plot_survival_data(data_cs[2], fmt="C2.-", label="20mW")
NaCsPlot.plot_survival_data(data_cs[3], fmt="C3.-", label="30mW")
legend()
grid()
ylim([0, 1])
title("Cs RP resonance")
xlabel("CS RPDet (double pass) (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
