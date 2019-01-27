#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20190122_065658.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs_na = [OrderedDict(
    :cold=>-140:6.0:220,
    :hot=>-140:6.0:220,
)]
const specs_cs = [OrderedDict(
    :cold=>-40:2.0:80,
    :hot=>-40:2.0:80,
)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]
const datas_na = select_datas(datas, NaCsData.select_single((1,), (3,)), maxcnts, specs_na)
const datas_cs = select_datas(datas, NaCsData.select_single((2,), (4,)), maxcnts, specs_cs)

data_na_cold = datas_na[1][:cold]
data_na_hot = datas_na[1][:hot]
data_cs_cold = datas_cs[1][:cold]
data_cs_hot = datas_cs[1][:hot]

const prefix = joinpath(@__DIR__, "imgs", "data_20190122_065658_axial")

figure()
NaCsPlot.plot_survival_data(data_na_cold, fmt="C0.-", label="Cold")
NaCsPlot.plot_survival_data(data_na_hot, fmt="C1.-", label="Hot")
legend()
grid()
ylim([0, 1])
title("Na Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_na")

figure()
NaCsPlot.plot_survival_data(data_cs_cold, fmt="C0.-", label="Cold")
NaCsPlot.plot_survival_data(data_cs_hot, fmt="C1.-", label="Hot")
legend()
grid()
ylim([0, 1])
title("Cs Axial")
xlabel("Detuning (kHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)_cs")

NaCsPlot.maybe_show()
