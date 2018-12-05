#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20181205_095818.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [OrderedDict(:r50=>311.3 .+ (-400:50:400) .* 1e-3,
                           :r200=>311.3 .+ (-400:50:400) .* 1e-3)]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2,), (3, 4,)), maxcnts, specs)

const prefix = joinpath(@__DIR__, "imgs", "data_20181205_095818_raman_ratios")

figure()
NaCsPlot.plot_survival_data(datas_nacs[1][:r50], fmt="C0.-", label="0.02")
NaCsPlot.plot_survival_data(datas_nacs[1][:r200], fmt="C1.-", label="0.005")
grid()
legend()
xlabel("Raman frequency (MHz)")
ylabel("Survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
