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

const inames = ["data_20200317_110129.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
# const specs = [[0, 0.06, 0.11, 0.16, 0.21, 0.26, 0.31]]
const specs = [[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

const data = datas_nacs[1]

const prefix = joinpath(@__DIR__, "imgs", "data_20200317_110129_raman_time_3322")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
function model_expoff(x, p)
    p[1] .+ p[2] .* exp.(.-x ./ p[3])
end
function model_expsin(x, p)
    p[1] .* cos.(x .* 2π .* p[2]).^2 .* exp.(.-x ./ p[3]) .+ p[4]
end
fit = fit_survival(model_expsin, data, [0.25, 3, 0.25, 0.02])

figure()
NaCsPlot.plot_survival_data(data, fmt="C0.")
plot(fit.plotx, fit.ploty)
title("288560 GHz, 15 mW, 770.5945 MHz")
text(0.01, 0.26, "\$p_0\\cdot \\cos^2(\\Omega\\cdot t)\\cdot\\exp(-t/\\tau)+p_1\$", color="C0", fontsize="small")
text(0.1, 0.19, "\$\\Omega=2\\pi\\cdot$(fit.uncs[2])\$ kHz", color="C0", fontsize="small")
text(0.1, 0.16, "\$\\tau=$(fit.uncs[3])\$ ms", color="C0", fontsize="small")
xlim([0, 0.26])
grid()
xlabel("Raman time (ms)\${}_{\\mathrm{(with\\ 0.06\\ ms\\ offset)}}\$")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
