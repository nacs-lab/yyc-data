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

const inames = ["data_20200421_223616.mat",
                "data_20200422_093717.mat",
                "data_20200428_115959.mat",
                "data_20200428_200837.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [(484.00 .+ [-30; -4.0:0.8:4.0; 30], # 15 mW, 0.20 ms
                433.00 .+ [-20; -2.5:0.5:2.5; 20], # 12 mW, 0.25 ms
                379.75 .+ [-10; -1.5:0.3:1.5; 10], #  9 mW, 0.40 ms
                324.5 .+ [-8; -1.25:0.25:1.25; 8], #  6 mW, 0.80 ms
                265.00 .+ [-5; -1.0:0.2:1.0; 5], # 3 mW, 2.0 ms
                ),
               (484.00 .+ [-30; -4.0:0.8:4.0; 30], # 15 mW, 0.20 ms
                433.00 .+ [-20; -2.5:0.5:2.5; 20], # 12 mW, 0.25 ms
                379.75 .+ [-10; -1.5:0.3:1.5; 10], #  9 mW, 0.40 ms
                324.5 .+ [-8; -1.25:0.25:1.25; 8], #  6 mW, 0.80 ms
                265.00 .+ [-5; -1.0:0.2:1.0; 5], # 3 mW, 2.0 ms
                ),
               ([0; [0.06, 0.11, 0.16, 0.31, 0.46, 0.61] .* 2 .- 0.01], # 9 mW, 770.379740 MHz
                [0.0], # 9 mW, 0 ms
                379.740 .+ [-10; -1.5:0.3:1.5; 10], # 9 mW, 0.21 ms
                379.740 .+ [-10; -1.5:0.3:1.5; 10], # 9 mW, 0.61 ms
                379.740 .+ [-10; -1.5:0.3:1.5; 10], # 9 mW, 0.96 ms
                ),
               ([0; [0.06, 0.11, 0.16, 0.31, 0.46, 0.61] .* 2 .- 0.01], # 9 mW, 770.379740 MHz
                [0.0], # 9 mW, 0 ms
                379.740 .+ [-10; -1.5:0.3:1.5; 10], # 9 mW, 0.21 ms
                379.740 .+ [-10; -1.5:0.3:1.5; 10], # 9 mW, 0.61 ms
                379.740 .+ [-10; -1.5:0.3:1.5; 10], # 9 mW, 0.96 ms
                )]

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
const data_nacs_39 = [datas_nacs[1][3]; datas_nacs[2][3]] # Survival 2
const data_nacs_20 = [datas_nacs[3][3]; datas_nacs[4][3]] # Survival 1
const data_nacs_60 = [datas_nacs[3][4]; datas_nacs[4][4]] # Survival 1
const data_nacs_95 = [datas_nacs[3][5]; datas_nacs[4][5]] # Survival 1

const data_fit = [NaCsData.map_params((i, v) -> (v, 0.0, 1), data_nacs_00);
                  NaCsData.map_params((i, v) -> (v, 0.20, 1), data_nacs_20);
                  NaCsData.map_params((i, v) -> (v, 0.39, 2), data_nacs_39);
                  NaCsData.map_params((i, v) -> (v, 0.60, 1), data_nacs_60);
                  NaCsData.map_params((i, v) -> (v, 0.95, 1), data_nacs_95);
                  NaCsData.map_params((i, v) -> (379.740, v, 1), data_nacs_t)]

const prefix = joinpath(@__DIR__, "imgs", "fit_20200428_115959_raman_3322")

function get_model_param(p, idx)
    p0r, p11, p12, f0, Ω, Γ1, Γ2 = p
    p1 = (p11, p12)[idx]
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
fit = fit_survival(model, data_fit, [0.1, 0.3, 0.3, 379.740, 2π * 1.5, 0, 2π / 0.2], plotx=false)
@show fit.uncs

const plot_freq_lo = 379.740 - 11
const plot_freq_hi = 379.740 + 11
const plot_freq = linspace(plot_freq_lo, plot_freq_hi, 1000)

figure()
ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
v00 = model_2d(0, 0, get_model_param(fit.param, 1))
errorbar([plot_freq_lo + 1, plot_freq_hi - 1], [ratio_00, ratio_00], [uncs_00, uncs_00], fmt="C0.", label="0.00 ms")
plot([plot_freq_lo, plot_freq_hi], [v00, v00], "C0-")
NaCsPlot.plot_survival_data(data_nacs_20, fmt="C1.", label="0.20 ms")
plot(plot_freq, model_2d.(0.20, plot_freq, (get_model_param(fit.param, 1),)), "C1")
NaCsPlot.plot_survival_data(data_nacs_39, fmt="C2.", label="0.39 ms")
plot(plot_freq, model_2d.(0.39, plot_freq, (get_model_param(fit.param, 2),)), "C2")
NaCsPlot.plot_survival_data(data_nacs_60, fmt="C3.", label="0.60 ms")
plot(plot_freq, model_2d.(0.60, plot_freq, (get_model_param(fit.param, 1),)), "C3")
NaCsPlot.plot_survival_data(data_nacs_95, fmt="C4.", label="0.95 ms")
plot(plot_freq, model_2d.(0.95, plot_freq, (get_model_param(fit.param, 1),)), "C4")
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")
title("288503 GHz, 9 mW")
NaCsPlot.maybe_save("$(prefix)_f")

const plot_time = linspace(0, 1.3, 1000)

figure()
NaCsPlot.plot_survival_data(data_nacs_t, fmt="C0.")
plot(plot_time, model_2d.(plot_time, 379.740, (get_model_param(fit.param, 1),)), "C0")
title("288503 GHz, 9 mW, 770.379740 MHz")
text(0.15, 0.19, ("\$f_{res}=$(fit.uncs[4] / 1000 + 770)\$ MHz\n" *
                  "\$\\mathbf{\\Omega_{Raman}=2\\pi\\times$(fit.uncs[5] / 2π) kHz}\$\n" *
                  "\$\\Gamma_{atom}=2\\pi\\times$(fit.uncs[6] / 2π)\$ kHz\n" *
                  "\$\\mathbf{\\Gamma_{molecule}=2\\pi\\times$(fit.uncs[7] / 2π) kHz}\$\n"),
     color="C0", fontsize="small")
xlim([0, 1.35])
grid()
xlabel("Raman time (ms)\${}_{\\mathrm{(with\\ 0.01\\ ms\\ offset)}}\$")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_t")

# figure()
# const img_freq = fit.param[4] .+ linspace(-10, 10, 201)
# const img_time = linspace(0, 0.6, 201)
# const mol_2d = [model_2d(t, f, get_model_param(fit.param, 1), molecule=true)
#                 for f in img_freq, t in img_time]
# imshow(mol_2d, aspect="auto", interpolation="none", origin="lower",
#        extent=[img_time[1] - step(img_time) / 2, img_time[end] + step(img_time) / 2,
#                img_freq[1] - step(img_freq) / 2, img_freq[end] + step(img_freq) / 2])
# xlabel("Raman time (ms)\${}_{\\mathrm{(with\\ 0.01\\ ms\\ offset)}}\$")
# ylabel("2-Photon Detuning (770XXX kHz)")
# title("Molecule Population")
# grid()
# colorbar()
# NaCsPlot.maybe_save("$(prefix)_mol")

NaCsPlot.maybe_show()
