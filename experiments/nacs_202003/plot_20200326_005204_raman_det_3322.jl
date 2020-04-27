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

const inames = ["data_20200326_005204.mat",
                "data_20200326_101325.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [(593.0 .+ [-20; -6:2:6; 20], # 15 mW, 0.09 ms
                522.7 .+ [-15; -4.5:1.5:4.5; 15], # 12 mW, 0.12 ms
                447.5 .+ [-10; -3:0.6:3; 5], # 9 mW, 0.16 ms
                369.5 .+ [-5; -2:0.5:2; 5], # 6 mW, 0.27 ms
                287.5 .+ [-5; -0.6:0.15:0.6; 5], # 3 mW, 1 ms
                ),
               (593.0 .+ [-20; -6:2:6; 20], # 15 mW, 0.09 ms
                522.7 .+ [-15; -4.5:1.5:4.5; 15], # 12 mW, 0.12 ms
                447.5 .+ [-10; -3:0.6:3; 5], # 9 mW, 0.16 ms
                369.5 .+ [-5; -2:0.5:2; 5], # 6 mW, 0.27 ms
                287.5 .+ [-5; -0.6:0.15:0.6; 5], # 3 mW, 0.9 ms
                )]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

const prefix = joinpath(@__DIR__, "imgs", "data_20200326_005204_raman_det_3322")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
# fit = fit_survival(model_lorentzian, data, [0.6, 0.2, 770.374, 0.02])

powers = [15, 12, 9, 6, 3]
times_old = [0.09, 0.12, 0.16, 0.27, 1]
times_new = [0.09, 0.12, 0.16, 0.27, 0.9]

figure(figsize=[12.6, 16.8])

for i in 1:5
    subplot(3, 2, i)
    time_old = times_old[i]
    time_new = times_new[i]
    if time_old == time_new
        title_suffix = ", $time_old ms"
        label_old_suffix = ""
        label_new_suffix = ""
    else
        title_suffix = ""
        label_old_suffix = ", $time_old ms"
        label_new_suffix = ", $time_new ms"
    end
    NaCsPlot.plot_survival_data(datas_nacs[1][i], fmt="C0.-", label="Old" * label_old_suffix)
    NaCsPlot.plot_survival_data(datas_nacs[2][i], fmt="C1.-", label="New" * label_new_suffix)
    title("288560 GHz, $(powers[i]) mW" * title_suffix)
    legend()
    grid()
    xlabel("2-Photon Detuning (770XXX kHz)")
    ylabel("Two-body survival")
end
tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
