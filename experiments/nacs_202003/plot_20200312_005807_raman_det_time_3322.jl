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

const inames = ["data_20200312_005807.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [(772.1 .+ linspace(-0.1, 0.175, 45),
                [0, 0.8, 1.0, 2.0, 3.0, 5.0, 10.0] .* 0.1)]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

const data_f = datas_nacs[1][1]
const data_t = datas_nacs[1][2]

const prefix = joinpath(@__DIR__, "imgs", "data_20200312_005807_raman_3322")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
fit_f = fit_survival(model_lorentzian, data_f, [0.6, 0.2, 772.1, 0.1])

figure()
NaCsPlot.plot_survival_data(data_f, fmt="C0.")
plot(fit_f.plotx, fit_f.ploty, "C0-")
title("288668.35 GHz, 21 mW, 0.15 ms")
text(771.98, 0.45, "\$f=$(fit_f.uncs[3])\$ MHz", color="C0", fontsize="small")
text(771.98, 0.47, "\$\\Gamma=$(fit_f.uncs[4])\$ MHz", color="C0", fontsize="small")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_det")

figure()
NaCsPlot.plot_survival_data(data_t, fmt="C0.-")
title("288668.35 GHz, 21 mW, 772.1 MHz")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_time")

NaCsPlot.maybe_show()
