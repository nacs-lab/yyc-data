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

const inames = ["data_20200313_175808.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [[0, 0.8, 1.0, 2.0, 3.0, 5.0, 10.0] .* 0.2]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

const data_t = datas_nacs[1]

const prefix = joinpath(@__DIR__, "imgs", "data_20200313_175808_raman_time_3322")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
function model_expoff(x, p)
    p[1] .+ p[2] .* exp.(.-x ./ p[3])
end
fit_t = fit_survival(model_expoff, data_t, [0.45, 0.25, 0.5])

figure()
NaCsPlot.plot_survival_data(data_t, fmt="C0.")
plot(fit_t.plotx, fit_t.ploty, "C0-")
title("288668.35 GHz, 12 mW, 771.35 MHz")
xlim([0, 2.2])
text(0.7, 0.56, "\$\\tau=$(fit_t.uncs[3])\$ ms", color="C0")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
