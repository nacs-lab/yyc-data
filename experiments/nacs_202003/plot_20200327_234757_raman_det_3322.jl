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
                "data_20200327_234757.mat",
                "data_20200328_015629.mat",
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
               (369.5 .+ [-10; -2.5:0.5:2.5; 10], # 6 mW, 0.48 ms
                ),
               ([0], # 6 mW, 0 ms
                369.5 .+ [-10; -2.5:0.5:2.5; 10], # 6 mW, 0.48 ms
                369.5 .+ [-10; -2.5:0.5:2.5; 10], # 6 mW, 0.24 ms
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
const data_nacs_20 = datas_nacs[3][3]
const data_nacs_23 = datas_nacs[1][4]
const data_nacs_44 = [datas_nacs[2][1]; datas_nacs[3][2]]
const data_nacs_92 = datas_nacs[4][2]

const prefix = joinpath(@__DIR__, "imgs", "data_20200327_234757_raman_det_3322")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
fit_20 = fit_survival(model_lorentzian, data_nacs_20, [0.35, 0.15, 369.5, 5])
fit_23 = fit_survival(model_lorentzian, data_nacs_23, [0.35, 0.15, 369.5, 5])
fit_44 = fit_survival(model_lorentzian, data_nacs_44, [0.35, 0.25, 369.5, 5])
fit_92 = fit_survival(model_lorentzian, data_nacs_92, [0.35, 0.25, 369.5, 5])

# @show fit_20.uncs
# @show fit_23.uncs
# @show fit_44.uncs
# @show fit_92.uncs

figure()
ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
errorbar([358.6, 380.3], [ratio_00, ratio_00], [uncs_00, uncs_00], fmt="C0.-", label="0ms")
NaCsPlot.plot_survival_data(data_nacs_20, fmt="C1.", label="0.20 ms")
plot(fit_20.plotx, fit_20.ploty, "C1")
NaCsPlot.plot_survival_data(data_nacs_23, fmt="C2.", label="0.23 ms")
plot(fit_23.plotx, fit_23.ploty, "C2")
NaCsPlot.plot_survival_data(data_nacs_44, fmt="C3.", label="0.44 ms")
plot(fit_44.plotx, fit_44.ploty, "C3")
NaCsPlot.plot_survival_data(data_nacs_92, fmt="C4.", label="0.92 ms")
plot(fit_92.plotx, fit_92.ploty, "C4")
text(357.6, 0.017, ("\$\\Gamma_{0.92}=\\!$(fit_92.uncs[4]) kHz\$\n" *
                    "\$f_{0.92}=\\!$(770 + fit_92.uncs[3] / 1000) MHz\$"),
     color="C4", fontsize="x-small", linespacing=0.8)
text(357.6, 0.052, ("\$\\Gamma_{0.44}=\\!$(fit_44.uncs[4]) kHz\$\n" *
                    "\$f_{0.44}=\\!$(770 + fit_44.uncs[3] / 1000) MHz\$"),
     color="C3", fontsize="x-small", linespacing=0.8)
text(357.6, 0.087, ("\$\\Gamma_{0.23}=\\!$(fit_23.uncs[4]) kHz\$\n" *
                    "\$f_{0.23}=\\!$(770 + fit_23.uncs[3] / 1000) MHz\$"),
     color="C2", fontsize="x-small", linespacing=0.8)
text(357.6, 0.122, ("\$\\Gamma_{0.20}=\\!$(fit_20.uncs[4]) kHz\$\n" *
                    "\$f_{0.20}=\\!$(770 + fit_20.uncs[3] / 1000) MHz\$"),
     color="C1", fontsize="x-small", linespacing=0.8)
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")
title("288560 GHz, 6 mW")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
