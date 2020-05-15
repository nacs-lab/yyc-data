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
                "data_20200511_233303.mat",
                "data_20200512_094924.mat"]
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
               ([[0.21, 0.41, 0.61, 0.81, 1.01, 1.21, 1.51, 1.91] .- 0.01;], # 3 mW, 770.32140 MHz
                [0.0], # 3 mW, 0 ms
                321.40 .+ (-1.5:0.3:1.5), # 3 mW, 1.25 ms
                321.40 .+ (-1.5:0.3:1.5), # 3 mW, 0.62 ms
                ),
               ([[0.21, 0.41, 0.61, 0.81, 1.01, 1.21, 1.51, 1.91] .- 0.01;], # 3 mW, 770.32140 MHz
                [0.0], # 3 mW, 0 ms
                321.40 .+ (-1.5:0.3:1.5), # 3 mW, 1.25 ms
                321.40 .+ (-1.5:0.3:1.5), # 3 mW, 0.62 ms
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
const data_nacs_t = [datas_nacs[3][1]; datas_nacs[4][1]] # Survival 1
const data_nacs_00 = [datas_nacs[3][2]; datas_nacs[4][2]] # Survival 1
const data_nacs_49 = datas_nacs[1][3] # Survival 2
const data_nacs_61 = [datas_nacs[3][4]; datas_nacs[4][4]] # Survival 1
const data_nacs_124 = [datas_nacs[3][3]; datas_nacs[4][3]] # Survival 1

const data_fit = [NaCsData.map_params((i, v) -> (v, 0.0, 1), data_nacs_00);
                  NaCsData.map_params((i, v) -> (v, 0.49, 2), data_nacs_49);
                  NaCsData.map_params((i, v) -> (v, 0.61, 1), data_nacs_61);
                  NaCsData.map_params((i, v) -> (v, 1.24, 1), data_nacs_124);
                  NaCsData.map_params((i, v) -> (321.40, v, 1), data_nacs_t)]

const prefix = joinpath(@__DIR__, "imgs", "fit_20200511_233303_raman_3322")

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
fit = fit_survival(model, data_fit, [321.40, 2π * 1.5, 0, 2π / 0.2, 0.1, 0.3, 0.3],
                   plotx=false, lower=[0.0, 0.0, 0.0, 0, 0, 0, 0])
@show fit.uncs
const param_1 = get_model_param(fit.param, 1)
const uncs_1 = get_model_param(fit.uncs, 1)
const param_2 = get_model_param(fit.param, 2)

const plot_freq_lo = 321.40 - 8
const plot_freq_hi = 321.40 + 8
const plot_freq = linspace(plot_freq_lo, plot_freq_hi, 1000)

figure(figsize=[12.6, 5.6])

subplot(1, 2, 1)
ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
v00 = model_2d(0, 0, param_1)
errorbar([param_1[3]], [ratio_00], [uncs_00], fmt="C0.", label="0.00 ms")
plot([plot_freq_lo, plot_freq_hi], [v00, v00], "C0-")
NaCsPlot.plot_survival_data(data_nacs_49, fmt="C1.", label="0.49 ms")
plot(plot_freq, model_2d.(0.49, plot_freq, (param_2,)), "C1")
NaCsPlot.plot_survival_data(data_nacs_61, fmt="C2.", label="0.61 ms")
plot(plot_freq, model_2d.(0.61, plot_freq, (param_1,)), "C2")
NaCsPlot.plot_survival_data(data_nacs_124, fmt="C3.", label="1.24 ms")
plot(plot_freq, model_2d.(1.24, plot_freq, (param_1,)), "C3")
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")
title("288605 GHz, 3 mW")

subplot(1, 2, 2)
const plot_time = linspace(0, 1.95, 1000)
NaCsPlot.plot_survival_data([data_nacs_00; data_nacs_t], fmt="C0.", label="770.32140 MHz")
plot(plot_time, model_2d.(plot_time, 321.40, (param_1,)), "C0")
title("Time")
legend(fontsize="x-small", loc="upper right")
text(0.50, 0.14, ("\$f_{res}=$(uncs_1[3] / 1000 + 770)\$ MHz\n" *
                  "\$\\mathbf{\\Omega_{Raman}=2\\pi\\times$(uncs_1[4] / 2π) kHz}\$\n" *
                  "\$\\Gamma_{atom}=2\\pi\\times$(uncs_1[5] / 2π * 1000)\$ Hz\n" *
                  "\$\\mathbf{\\Gamma_{molecule}=2\\pi\\times$(uncs_1[6] / 2π) kHz}\$\n"),
     color="C0", fontsize="small")
xlim([0, 2.0])
grid()
xlabel("Raman time (ms)\${}_{\\mathrm{(with\\ 0.01\\ ms\\ offset)}}\$")
ylabel("Two-body survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

# figure()
# param2 = copy(fit.param)
# # p0r, p11, p12, p13, f0, Ω, Γ1, Γ2
# param2[5] = 321.40
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
