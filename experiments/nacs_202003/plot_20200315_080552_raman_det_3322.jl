#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200315_080552.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [770.449 .+ (-10:2:10) .* 1e-3]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

const data = datas_nacs[1]

const prefix = joinpath(@__DIR__, "imgs", "data_20200315_080552_raman_det_3322")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
fit = fit_survival(model_lorentzian, data, [0.6, 0.2, 770.449, 0.05])

figure()
NaCsPlot.plot_survival_data(data, fmt="C0.")
plot(fit.plotx, fit.ploty, "C0-")
title("288560 GHz, 9 mW, 2 ms")
text(770.438, 0.47, "\$f=$(fit.uncs[3])\$ MHz", color="C0", fontsize="small")
text(770.438, 0.495, "\$\\Gamma=$(fit.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
