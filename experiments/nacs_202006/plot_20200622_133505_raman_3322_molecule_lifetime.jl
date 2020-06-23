#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200622_133505.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [[0, 0.01, 0.02, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0],
               ]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const data_nacs = datas_nacs[1]

const prefix = joinpath(@__DIR__, "imgs", "data_20200622_133505_raman_3322_molecule_lifetime")

function model_exp(x, p)
    if length(p) > 2
        p[1] .* exp.(.- x .* p[2]) .+ p[3]
    else
        p[1] .* exp.(.- x .* p[2])
    end
end

fit = fit_survival(model_exp, data_nacs, [0.13, 5, 0.02])

@show fit.uncs

figure()
NaCsPlot.plot_survival_data(data_nacs, fmt="C0.")
plot(fit.plotx, fit.ploty, "C0")
text(0.35, 0.08, "\$\\Gamma_{molecule}=2\\pi\\times$(fit.uncs[2] / 2π)\$ kHz")
xlim([0, 2.2])
title("288625 GHz, 6 mW")
grid()
xlabel("Time (ms)")
ylabel("Molecule survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
