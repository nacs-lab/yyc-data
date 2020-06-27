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

const inames = ["data_20200622_231831.mat",
                "data_20200623_100410.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [[0, 0.01, 0.02, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0] .* 0.1,
               [0, 0.10, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0] .* 0.4]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const data_nacs1 = datas_nacs[1]
const data_nacs2 = datas_nacs[2]
const data_fit = [NaCsData.map_params((i, v) -> (v, 1), data_nacs1);
                  NaCsData.map_params((i, v) -> (v, 2), data_nacs2)]

const prefix = joinpath(@__DIR__, "imgs", "data_20200622_231831_raman_3322_molecule_lifetime")

function model_exp(x, p)
    if length(p) > 2
        p[1] .* exp.(.- x .* p[2]) .+ p[3]
    else
        p[1] .* exp.(.- x .* p[2])
    end
end

function get_model_param(p, idx)
    return (p[2 * idx], p[1], p[2 * idx + 1])
end

function model(xs, p)
    function wrapper(x)
        t, idx = x
        return model_exp(t, get_model_param(p, idx))
    end
    return wrapper.(xs)
end

fit = fit_survival(model, data_fit, [10, 0.13, 0.02, 0.13, 0.02], plotx=false)
@show fit.uncs
const param_1 = get_model_param(fit.param, 1)
const param_2 = get_model_param(fit.param, 2)
const uncs_1 = get_model_param(fit.uncs, 1)

const plot_t = linspace(0, 0.85, 1000)

figure()
NaCsPlot.plot_survival_data(data_nacs1, fmt="C0.")
plot(plot_t, model_exp.(plot_t, (param_1,)), "C0")
NaCsPlot.plot_survival_data(data_nacs2, fmt="C1.")
plot(plot_t, model_exp.(plot_t, (param_2,)), "C1")
text(0.15, 0.08, "\$\\Gamma_{molecule}=2\\pi\\times$(uncs_1[2] / 2π)\$ kHz")
xlim([0, 0.9])
title("288625 GHz, 15 mW")
grid()
xlabel("Time (ms)")
ylabel("Molecule survival")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
