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

const inames = ["data_20200401_060025.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int)]
const specs = [([0],
                [0],
                [0, 0.08, 0.16, 0.24, 0.3],
                )]
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

const prefix = joinpath(@__DIR__, "imgs", "data_20200401_060025_molecule_lifetime_3322")

function model_expoff(x, p)
    p[3] .+ p[1] .* exp.(.- x .* p[2])
end
function model_exp(x, p)
    p[1] .* exp.(.- x .* p[2])
end
fit1 = fit_survival(model_expoff, datas_nacs[1][3], [0.1, 6, 0.01])

# @show fit1.uncs

figure()
NaCsPlot.plot_survival_data(datas_nacs[1][3], fmt="C0s")
plot(fit1.plotx, fit1.ploty, "C0-")
xlim([0, 0.32])
ylim([0.01, 0.11])
text(0.1, 0.06, "\$\\gamma=2\\pi\\times$(fit1.uncs[2] / 2π)\$ kHz", color="C0")
title("306560 GHz, 15 mW")
grid()
xlabel("Molecule Hold Time (ms)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
