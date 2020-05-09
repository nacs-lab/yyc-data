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

const inames = ["data_20200505_023054.mat",
                "data_20200505_204707.mat",
                "data_20200507_172024.mat",
                "data_20200508_205921.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [(547.00 .+ [-45; -7.5:1.5:7.5; 45], #  9 mW, 0.135 ms
                437.60 .+ [-12; -2.0:0.4:2.0; 12], #  6 mW, 0.23 ms
                321.25 .+ [-8; -1.25:0.25:1.25; 8], # 3 mW, 0.5 ms
                ),
               (765.00 .+ [-120; -20:4.0:20; 120], # 15 mW, 0.07 ms
                654.00 .+ [-90; -15:3.0:15; 90], #   12 mW, 0.10 ms
                ),
               ([[0.12, 0.23, 0.33, 0.44, 0.55, 0.73, 0.9] .- 0.01;], # 6 mW, 770.43678 MHz
                [0.0], # 9 mW, 0 ms
                436.78 .+ (-2:0.4:2), # 9 mW, 0.16 ms
                436.78 .+ (-2:0.4:2), # 9 mW, 0.34 ms
                436.78 .+ (-2:0.4:2), # 9 mW, 0.56 ms
                ),
               [0; [0.18, 0.23, 0.33, 0.44, 0.55, 0.73, 0.9] .- 0.01], # 6 mW, 770.437066 MHz
               ]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function get_ratio_val(data)
    params, ratios, uncs = NaCsData.get_values(data)
    return ratios[2], uncs[2]
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const data_nacs_t = datas_nacs[3][1] # Survival 1
const data_nacs_t2 = datas_nacs[4] # Survival 3
const data_nacs_00 = datas_nacs[3][2] # Survival 1
const data_nacs_22 = datas_nacs[1][2] # Survival 2
const data_nacs_15 = datas_nacs[3][3] # Survival 1
const data_nacs_33 = datas_nacs[3][4] # Survival 1
const data_nacs_55 = datas_nacs[3][5] # Survival 1

const data_fit = [NaCsData.map_params((i, v) -> (v, 0.0, 1), data_nacs_00);
                  NaCsData.map_params((i, v) -> (v, 0.15, 1), data_nacs_15);
                  NaCsData.map_params((i, v) -> (v, 0.22, 2), data_nacs_22);
                  NaCsData.map_params((i, v) -> (v, 0.33, 1), data_nacs_33);
                  NaCsData.map_params((i, v) -> (v, 0.55, 1), data_nacs_55);
                  NaCsData.map_params((i, v) -> (436.78, v, 1), data_nacs_t);
                  NaCsData.map_params((i, v) -> (437.066, v, 3), data_nacs_t2)]

const prefix = joinpath(@__DIR__, "imgs", "fit_20200507_172024_raman_3322")

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
fit = fit_survival(model, data_fit, [436.78, 2π * 1.5, 0, 2π / 0.2, 0.1, 0.3, 0.3, 0.3],
                   plotx=false, lower=[0.0, 0.0, 0.0, 0, 0, 0, 0, 0])
@show fit.uncs
const param_1 = get_model_param(fit.param, 1)
const uncs_1 = get_model_param(fit.uncs, 1)
const param_2 = get_model_param(fit.param, 2)
const param_3 = get_model_param(fit.param, 3)

const plot_freq_lo = 436.78 - 13
const plot_freq_hi = 436.78 + 13
const plot_freq = linspace(plot_freq_lo, plot_freq_hi, 1000)

figure(figsize=[12.6, 5.6])

subplot(1, 2, 1)
ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
v00 = model_2d(0, 0, param_1)
errorbar([param_1[3]], [ratio_00], [uncs_00], fmt="C0.", label="0.00 ms")
plot([plot_freq_lo, plot_freq_hi], [v00, v00], "C0-")
NaCsPlot.plot_survival_data(data_nacs_15, fmt="C1.", label="0.15 ms")
plot(plot_freq, model_2d.(0.15, plot_freq, (param_1,)), "C1")
NaCsPlot.plot_survival_data(data_nacs_22, fmt="C2.", label="0.22 ms")
plot(plot_freq, model_2d.(0.22, plot_freq, (param_2,)), "C2")
NaCsPlot.plot_survival_data(data_nacs_33, fmt="C3.", label="0.33 ms")
plot(plot_freq, model_2d.(0.33, plot_freq, (param_1,)), "C3")
NaCsPlot.plot_survival_data(data_nacs_55, fmt="C4.", label="0.55 ms")
plot(plot_freq, model_2d.(0.55, plot_freq, (param_1,)), "C4")
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")
title("288605 GHz, 6 mW")

subplot(1, 2, 2)
const plot_time = linspace(0, 0.95, 1000)
NaCsPlot.plot_survival_data([data_nacs_00; data_nacs_t], fmt="C0.", label="770.43678 MHz")
plot(plot_time, model_2d.(plot_time, 436.78, (param_1,)), "C0")
NaCsPlot.plot_survival_data(data_nacs_t2, fmt="C1.", label="770.437066 MHz")
plot(plot_time, model_2d.(plot_time, 437.066, (param_3,)), "C1")
title("Time")
legend(fontsize="x-small", loc="upper right")
text(0.20, 0.14, ("\$f_{res}=$(uncs_1[3] / 1000 + 770)\$ MHz\n" *
                  "\$\\mathbf{\\Omega_{Raman}=2\\pi\\times$(uncs_1[4] / 2π) kHz}\$\n" *
                  "\$\\Gamma_{atom}=2\\pi\\times$(uncs_1[5] / 2π * 1000)\$ Hz\n" *
                  "\$\\mathbf{\\Gamma_{molecule}=2\\pi\\times$(uncs_1[6] / 2π) kHz}\$\n"),
     color="C0", fontsize="small")
xlim([0, 0.96])
grid()
xlabel("Raman time (ms)\${}_{\\mathrm{(with\\ 0.01\\ ms\\ offset)}}\$")
ylabel("Two-body survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

# figure()
# param2 = copy(fit.param)
# # p0r, p11, p12, p13, f0, Ω, Γ1, Γ2
# param2[5] = 436.78
# param2[6] *= (6 / 15)^(1.29)
# param2[7] *= (6 / 15)^(2.58)
# param2[8] *= (6 / 15)^(1)
# const img2_freq = param2[5] .+ linspace(-2.5, 2.5, 201)
# const img2_time = linspace(0, 1, 201)
# const mol2_2d = [model_2d(t, f, get_model_param(param2, 1), molecule=false)
#                 for f in img2_freq, t in img2_time]
# imshow(mol2_2d, aspect="auto", interpolation="none", origin="lower",
#        extent=[img2_time[1] - step(img2_time) / 2, img2_time[end] + step(img2_time) / 2,
#                img2_freq[1] - step(img2_freq) / 2, img2_freq[end] + step(img2_freq) / 2])
# xlabel("Raman time (ms)\${}_{\\mathrm{(with\\ 0.01\\ ms\\ offset)}}\$")
# ylabel("2-Photon Detuning (770XXX kHz)")
# title("Molecule Population")
# grid()
# colorbar()
# NaCsPlot.maybe_save("$(prefix)_atm6")

NaCsPlot.maybe_show()
