#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

include("../nacs_202003/molecular_raman_model.jl")

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200618_173603.mat",
                "data_20200618_201020.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [8050,
                 typemax(Int)]
const specs = [([0.04, 0.07, 0.085, 0.10, 0.13, 0.145, 0.16, 0.175, 0.19, 0.22, 0.25,
                 0.28] * 6 .- 0.01, # 3 mW, 770.3472 MHz
                [0.0], # 3 mW, 0 ms
                347.2 .+ (-2:0.4:2), # 3 mW, 1.02 ms
                347.2 .+ (-2:0.4:2), # 3 mW, 0.51 ms
                ),
               ([0.04, 0.07, 0.085, 0.10, 0.13, 0.145, 0.16, 0.175, 0.19, 0.22, 0.25,
                 0.28] * 6 .- 0.01, # 3 mW, 770.3472 MHz
                [0.0], # 3 mW, 0 ms
                347.2 .+ (-2:0.4:2), # 3 mW, 1.02 ms
                347.2 .+ (-2:0.4:2), # 3 mW, 0.51 ms
                ),
               ]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function get_ratio_val(data)
    params, ratios, uncs = NaCsData.get_values(data)
    return ratios[2], uncs[2]
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const data_nacs_t = [datas_nacs[1][1]; datas_nacs[2][1]] # Survival 1
const data_nacs_00 = [datas_nacs[1][2]; datas_nacs[2][2]] # Survival 1
const data_nacs_101 = [datas_nacs[1][3]; datas_nacs[2][3]] # Survival 1
const data_nacs_50 = [datas_nacs[1][4]; datas_nacs[2][4]] # Survival 1

const data_fit = [NaCsData.map_params((i, v) -> (347.2, v, 1), data_nacs_00);
                  NaCsData.map_params((i, v) -> (v, 0.50, 1), data_nacs_50);
                  NaCsData.map_params((i, v) -> (v, 1.01, 1), data_nacs_101);
                  NaCsData.map_params((i, v) -> (347.2, v, 1), data_nacs_t)]

const prefix = joinpath(@__DIR__, "imgs", "fit_20200618_173603_raman_3322")

function get_model_param(p, idx)
    f0, Ω, Γ1, Γ2, p0r = p
    p1 = p[5 + idx]
    p0 = p1 * p0r
    return (p0, p1, f0, Ω, Γ1, Γ2)
end

function model(xs, p)
    function wrapper(x)
        f, t, idx = x
        return model_2d(t, f, get_model_param(p, idx))
    end
    return wrapper.(xs)
end
fit = fit_survival(model, data_fit, [347.2, 2π * 1.0, 0, 2π * 0.3, 0.2, 0.30],
                   plotx=false, lower=zeros(6))
@show fit.uncs
const param_1 = get_model_param(fit.param, 1)
const uncs_1 = get_model_param(fit.uncs, 1)

const plot_freq_lo = 347.2 - 2.1
const plot_freq_hi = 347.2 + 2.1
const plot_freq = linspace(plot_freq_lo, plot_freq_hi, 1000)

figure(figsize=[12.6, 5.6])

subplot(1, 2, 1)
ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
v00 = model_2d(0, 0, param_1)
errorbar([param_1[3]], [ratio_00], [uncs_00], fmt="C0.", label="0.00 ms")
plot([plot_freq_lo, plot_freq_hi], [v00, v00], "C0-")
NaCsPlot.plot_survival_data(data_nacs_50, fmt="C1.", label="0.50 ms")
plot(plot_freq, model_2d.(0.50, plot_freq, (param_1,)), "C1")
NaCsPlot.plot_survival_data(data_nacs_101, fmt="C2.", label="1.01 ms")
plot(plot_freq, model_2d.(1.01, plot_freq, (param_1,)), "C2")
legend(fontsize="x-small", loc="lower right")
xlim([plot_freq_lo, plot_freq_hi])
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")
title("288625 GHz, 3 mW")

subplot(1, 2, 2)
const plot_time = linspace(0, 1.70, 1000)
NaCsPlot.plot_survival_data([data_nacs_00; data_nacs_t], fmt="C0.", label="770.3472 MHz")
plot(plot_time, model_2d.(plot_time, 347.2, (param_1,)), "C0")
title("Time")
text(0.35, 0.20, ("\$f_{res}=$(uncs_1[3] / 1000 + 770)\$ MHz\n" *
                  "\$\\mathbf{\\Omega_{Raman}=2\\pi\\times$(uncs_1[4] / 2π) kHz}\$\n" *
                  "\$\\Gamma_{atom}=2\\pi\\times$(uncs_1[5] / 2π * 1000)\$ Hz\n" *
                  "\$\\mathbf{\\Gamma_{molecule}=2\\pi\\times$(uncs_1[6] / 2π) kHz}\$\n"),
     color="C0", fontsize="small")
xlim([0, 1.75])
legend(fontsize="x-small", loc="upper right")
grid()
xlabel("Raman time (ms)\${}_{\\mathrm{(with\\ 0.01\\ ms\\ offset)}}\$")
ylabel("Two-body survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()