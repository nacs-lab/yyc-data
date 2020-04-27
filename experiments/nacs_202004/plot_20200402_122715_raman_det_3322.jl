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

const inames = ["data_20200402_122715.mat",
                "data_20200402_154116.mat",
                "data_20200402_173658.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [614.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.14 ms
               377.5 .+ [-10; -2.5:0.5:2.5; 10], # 6 mW, 0.48 ms
               292.0 .+ [-5; -0.9:0.18:0.9; 5], # 3 mW, 1.3 ms
               ]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

const prefix = joinpath(@__DIR__, "imgs", "data_20200402_122715_raman_det_3322")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
fit15 = fit_survival(model_lorentzian, datas_nacs[1], [0.4, 0.3, 614.0, 10])
fit6 = fit_survival(model_lorentzian, datas_nacs[2], [0.3, 0.25, 377.5, 5])
fit3 = fit_survival(model_lorentzian, datas_nacs[3], [0.4, 0.3, 292.0, 2])

figure(figsize=[12.6, 11.2])

subplot(2, 2, 1)
NaCsPlot.plot_survival_data(datas_nacs[1], fmt="C0.")
plot(fit15.plotx, fit15.ploty, "C0")
text(586.3, 0.305, "\$f=$(770 + fit15.uncs[3] / 1000)\\ MHz\$")
title("288565 GHz, 15 mW, 0.13 ms")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")

subplot(2, 2, 2)
NaCsPlot.plot_survival_data(datas_nacs[2], fmt="C0.")
plot(fit6.plotx, fit6.ploty, "C0")
text(368, 0.305, "\$f=$(770 + fit6.uncs[3] / 1000)\\ MHz\$")
title("288565 GHz, 6 mW, 0.47 ms")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")

subplot(2, 2, 3)
NaCsPlot.plot_survival_data(datas_nacs[3], fmt="C0.")
plot(fit3.plotx, fit3.ploty, "C0")
text(287.5, 0.305, "\$f=$(770 + fit3.uncs[3] / 1000)\\ MHz\$")
title("288565 GHz, 3 mW, 1.29 ms")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
