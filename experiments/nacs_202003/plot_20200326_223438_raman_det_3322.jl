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

const inames = ["data_20200326_101325.mat",
                "data_20200326_223438.mat",
                "data_20200327_004910.mat",
                "data_20200330_033004.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [(593.0 .+ [-20; -6:2:6; 20], # 15 mW, 0.09 ms
                522.7 .+ [-15; -4.5:1.5:4.5; 15], # 12 mW, 0.12 ms
                447.5 .+ [-10; -3:0.6:3; 5], # 9 mW, 0.16 ms
                369.5 .+ [-5; -2:0.5:2; 5], # 6 mW, 0.27 ms
                287.5 .+ [-5; -0.6:0.15:0.6; 5], # 3 mW, 0.9 ms
                ),
               (593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.1 ms
                593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.14 ms
                ),
               ([0], # 15 mW, 0 ms
                593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.1 ms
                593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.14 ms
                ),
               (593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.32 ms
                369.5 .+ [-10; -2.5:0.5:2.5; 10], # 6 mW, 0.96 ms
                )
               ]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function get_ratio_val(data)
    params, ratios, uncs = NaCsData.get_values(data)
    return ratios[2], uncs[2]
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const data_nacs_00 = datas_nacs[3][1]
const data_nacs_05 = datas_nacs[1][1]
const data_nacs_06 = [datas_nacs[2][1]; datas_nacs[3][2]]
const data_nacs_10 = [datas_nacs[2][2]; datas_nacs[3][3]]
const data_nacs_28 = datas_nacs[4][1]

const prefix = joinpath(@__DIR__, "imgs", "data_20200326_223438_raman_det_3322")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
fit_05 = fit_survival(model_lorentzian, data_nacs_05, [0.35, 0.15, 594, 5])
fit_06 = fit_survival(model_lorentzian, data_nacs_06, [0.35, 0.15, 594, 5])
fit_10 = fit_survival(model_lorentzian, data_nacs_10, [0.35, 0.25, 594, 5])
fit_28 = fit_survival(model_lorentzian, data_nacs_28, [0.35, 0.25, 594, 5])

# @show fit_05.uncs
# @show fit_06.uncs
# @show fit_10.uncs
# @show fit_28.uncs

figure()
ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
errorbar([560, 625], [ratio_00, ratio_00], [uncs_00, uncs_00], fmt="C0.-", label="0ms")
NaCsPlot.plot_survival_data(data_nacs_05, fmt="C1.", label="0.05 ms")
plot(fit_05.plotx, fit_05.ploty, "C1")
NaCsPlot.plot_survival_data(data_nacs_06, fmt="C2.", label="0.06 ms")
plot(fit_06.plotx, fit_06.ploty, "C2")
NaCsPlot.plot_survival_data(data_nacs_10, fmt="C3.", label="0.10 ms")
plot(fit_10.plotx, fit_10.ploty, "C3")
NaCsPlot.plot_survival_data(data_nacs_28, fmt="C4.", label="0.28 ms")
plot(fit_28.plotx, fit_28.ploty, "C4")
text(557.5, 0.055, ("\$\\Gamma_{0.28}=\\!$(fit_28.uncs[4]) kHz\$\n" *
                    "\$f_{0.28}=\\!$(770 + fit_28.uncs[3] / 1000) MHz\$"),
     color="C4", fontsize="x-small", linespacing=0.8)
text(557.5, 0.090, ("\$\\Gamma_{0.10}=\\!$(fit_10.uncs[4]) kHz\$\n" *
                    "\$f_{0.10}=\\!$(770 + fit_10.uncs[3] / 1000) MHz\$"),
     color="C3", fontsize="x-small", linespacing=0.8)
text(557.5, 0.125, ("\$\\Gamma_{0.06}=\\!$(fit_06.uncs[4]) kHz\$\n" *
                    "\$f_{0.06}=\\!$(770 + fit_06.uncs[3] / 1000) MHz\$"),
     color="C2", fontsize="x-small", linespacing=0.8)
text(557.5, 0.160, ("\$\\Gamma_{0.05}=\\!$(fit_05.uncs[4]) kHz\$\n" *
                    "\$f_{0.05}=\\!$(770 + fit_05.uncs[3] / 1000) MHz\$"),
     color="C1", fontsize="x-small", linespacing=0.8)
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")
title("288560 GHz, 15 mW")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
