#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../../lib"))

include("../../../experiments/nacs_202003/molecular_raman_model.jl")

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const expdir = joinpath(@__DIR__, "../../../experiments")

const inames = ["nacs_202003/data/data_20200326_101325.mat",
                "nacs_202003/data/data_20200326_223438.mat",
                "nacs_202003/data/data_20200327_004910.mat",
                "nacs_202003/data/data_20200330_033004.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(expdir, iname)) for iname in inames]
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
                ),
               ]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

const data_nacs_00 = datas_nacs[3][1]
const data_nacs_05 = datas_nacs[1][1]
const data_nacs_06 = [datas_nacs[2][1]; datas_nacs[3][2]]
const data_nacs_10 = [datas_nacs[2][2]; datas_nacs[3][3]]
const data_nacs_28 = datas_nacs[4][1]

const prefix = joinpath(@__DIR__, "../figures/raman_transfer_spectrum")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
fit_l10 = fit_survival(model_lorentzian, data_nacs_10, [0.35, 0.25, 594, 5])

const data_nacs_10_shift = NaCsData.map_params((i, v)->v - fit_l10.param[3], data_nacs_10)

figure()
NaCsPlot.plot_survival_data(data_nacs_10_shift, fmt="C0.", label="0.10 ms")
plot(fit_l10.plotx .- fit_l10.param[3], fit_l10.ploty, "C0")
legend(fontsize="x-small", loc="lower right")
text(-24, 0.32, "\$f_{\\mathrm{Raman}}=$(770 + fit_l10.uncs[3] / 1000) \\mathrm{MHz}\$", color="C0", fontsize="small")
grid()
xlabel("Detuning from resonance (kHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

figure()
NaCsPlot.plot_survival_data(data_nacs_10_shift, fmt="C0.", label="0.10 ms")
plot(fit_l10.plotx .- fit_l10.param[3], fit_l10.ploty, "C0")
x₋ = fit_l10.param[3] - fit_l10.param[4] / 2
x₊ = fit_l10.param[3] + fit_l10.param[4] / 2
y_pm = (model_lorentzian(x₋, fit_l10.param) + model_lorentzian(x₊, fit_l10.param)) / 2
ax = gca()
ax.annotate("", (-fit_l10.param[4] / 2 - 1, y_pm), xytext=(-fit_l10.param[4] / 2 - 11, y_pm),
            arrowprops=Dict(:color=>"C3"))
ax.annotate("", (fit_l10.param[4] / 2 + 1, y_pm), xytext=(fit_l10.param[4] / 2 + 11, y_pm),
            arrowprops=Dict(:color=>"C3"))
text(-21, 0.16, "\$\\mathbf{\\Gamma_{FWHM}=$(fit_l10.uncs[4]) kHz}\$",
     fontsize=17, color="C3")
text(-24, 0.32, "\$f_{\\mathrm{Raman}}=$(770 + fit_l10.uncs[3] / 1000) \\mathrm{MHz}\$", color="C0", fontsize="small")
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("Detuning from resonance (kHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_text")

NaCsPlot.maybe_show()
