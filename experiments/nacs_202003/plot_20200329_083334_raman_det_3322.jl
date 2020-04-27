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
                "data_20200329_083334.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [(593.0 .+ [-20; -6:2:6; 20], # 15 mW, 0.09 ms
                522.7 .+ [-15; -4.5:1.5:4.5; 15], # 12 mW, 0.12 ms
                447.5 .+ [-10; -3:0.6:3; 5], # 9 mW, 0.16 ms
                369.5 .+ [-5; -2:0.5:2; 5], # 6 mW, 0.27 ms
                287.5 .+ [-5; -0.6:0.15:0.6; 5], # 3 mW, 0.9 ms
                ),
               ([0], # 3 mW, 0 ms
                287.2 .+ [-5; -0.9:0.18:0.9; 5], # 3 mW, 0.75 ms
                287.2 .+ [-5; -0.9:0.18:0.9; 5], # 3 mW, 1.3 ms
                287.2 .+ [-5; -0.9:0.18:0.9; 5], # 3 mW, 2.0 ms
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
const data_nacs_00 = datas_nacs[2][1]
const data_nacs_71 = datas_nacs[2][2]
const data_nacs_86 = datas_nacs[1][5]
const data_nacs_126 = datas_nacs[2][3]
const data_nacs_196 = datas_nacs[2][4]

const prefix = joinpath(@__DIR__, "imgs", "data_20200329_083334_raman_det_3322")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
fit_71 = fit_survival(model_lorentzian, data_nacs_71, [0.35, 0.15, 287.2, 1])
fit_86 = fit_survival(model_lorentzian, data_nacs_86, [0.35, 0.15, 287.2, 1])
fit_126 = fit_survival(model_lorentzian, data_nacs_126, [0.35, 0.25, 287.2, 1])
fit_196 = fit_survival(model_lorentzian, data_nacs_196, [0.35, 0.25, 287.2, 1])

# @show fit_71.uncs
# @show fit_86.uncs
# @show fit_126.uncs
# @show fit_196.uncs

figure()
ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
errorbar([281.6, 293.0], [ratio_00, ratio_00], [uncs_00, uncs_00], fmt="C0.-", label="0ms")
NaCsPlot.plot_survival_data(data_nacs_71, fmt="C1.", label="0.71 ms")
plot(fit_71.plotx, fit_71.ploty, "C1")
NaCsPlot.plot_survival_data(data_nacs_86, fmt="C2.", label="0.86 ms")
plot(fit_86.plotx, fit_86.ploty, "C2")
NaCsPlot.plot_survival_data(data_nacs_126, fmt="C3.", label="1.26 ms")
plot(fit_126.plotx, fit_126.ploty, "C3")
NaCsPlot.plot_survival_data(data_nacs_196, fmt="C4.", label="1.96 ms")
plot(fit_196.plotx, fit_196.ploty, "C4")
text(281.131, 0.0348, ("\$\\Gamma_{1.96}=\\!$(fit_196.uncs[4]) kHz\$\n" *
                       "\$f_{1.96}=\\!$(770 + fit_196.uncs[3] / 1000) MHz\$"),
     color="C4", fontsize="x-small", linespacing=0.8)
text(281.131, 0.0668, ("\$\\Gamma_{1.26}=\\!$(fit_126.uncs[4]) kHz\$\n" *
                       "\$f_{1.26}=\\!$(770 + fit_126.uncs[3] / 1000) MHz\$"),
     color="C3", fontsize="x-small", linespacing=0.8)
text(281.131, 0.0988, ("\$\\Gamma_{0.86}=\\!$(fit_86.uncs[4]) kHz\$\n" *
                       "\$f_{0.86}=\\!$(770 + fit_86.uncs[3] / 1000) MHz\$"),
     color="C2", fontsize="x-small", linespacing=0.8)
text(281.131, 0.1308, ("\$\\Gamma_{0.71}=\\!$(fit_71.uncs[4]) kHz\$\n" *
                       "\$f_{0.71}=\\!$(770 + fit_71.uncs[3] / 1000) MHz\$"),
     color="C1", fontsize="x-small", linespacing=0.8)
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")
title("288560 GHz, 3 mW")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
