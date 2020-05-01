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
                "data_20200424_112943.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
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
               ([0; [0.26, 0.51, 0.81, 1.01, 1.41] .* 1.5 .- 0.01], # 6 mW, 770.324501 MHz
                [0.0], # 6 mW, 0 ms
                324.501 .+ [-8; -1.25:0.25:1.25; 8], # 6 mW, 0.5 ms
                324.501 .+ [-8; -1.25:0.25:1.25; 8], # 6 mW, 1.0 ms
                324.501 .+ [-8; -1.25:0.25:1.25; 8], # 6 mW, 1.8 ms
                )]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function get_ratio_val(data)
    params, ratios, uncs = NaCsData.get_values(data)
    return ratios[2], uncs[2]
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const data_nacs_t = datas_nacs[3][1] # Survival 1
const data_nacs_00 = datas_nacs[3][2] # Survival 1
const data_nacs_79 = [datas_nacs[1][4]; datas_nacs[2][4]] # Survival 2
const data_nacs_49 = datas_nacs[3][3] # Survival 1
const data_nacs_99 = datas_nacs[3][4] # Survival 1
const data_nacs_179 = datas_nacs[3][5] # Survival 1

const data_fit = [NaCsData.map_params((i, v) -> (v, 0.0, 1), data_nacs_00);
                  NaCsData.map_params((i, v) -> (v, 0.49, 1), data_nacs_49);
                  NaCsData.map_params((i, v) -> (v, 0.79, 2), data_nacs_79);
                  NaCsData.map_params((i, v) -> (v, 0.99, 1), data_nacs_99);
                  NaCsData.map_params((i, v) -> (v, 1.79, 1), data_nacs_179);
                  NaCsData.map_params((i, v) -> (324.501, v, 1), data_nacs_t)]

const prefix = joinpath(@__DIR__, "imgs", "fit_20200424_112943_raman_3322")

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
fit = fit_survival(model, data_fit, [0.1, 0.3, 0.3, 324.501, 2π * 1.5, 0, 2π / 0.2], plotx=false)
@show fit.uncs

const plot_freq_lo = 324.501 - 10
const plot_freq_hi = 324.501 + 10
const plot_freq = linspace(plot_freq_lo, plot_freq_hi, 1000)

figure()
ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
v00 = model_2d(0, 0, get_model_param(fit.param, 1))
errorbar([plot_freq_lo + 2, plot_freq_hi - 2], [ratio_00, ratio_00], [uncs_00, uncs_00], fmt="C0.", label="0.00 ms")
plot([plot_freq_lo, plot_freq_hi], [v00, v00], "C0-")
NaCsPlot.plot_survival_data(data_nacs_49, fmt="C1.", label="0.49 ms")
plot(plot_freq, model_2d.(0.49, plot_freq, (get_model_param(fit.param, 1),)), "C1")
NaCsPlot.plot_survival_data(data_nacs_79, fmt="C2.", label="0.79 ms")
plot(plot_freq, model_2d.(0.79, plot_freq, (get_model_param(fit.param, 2),)), "C2")
NaCsPlot.plot_survival_data(data_nacs_99, fmt="C3.", label="0.99 ms")
plot(plot_freq, model_2d.(0.99, plot_freq, (get_model_param(fit.param, 1),)), "C3")
NaCsPlot.plot_survival_data(data_nacs_179, fmt="C4.", label="1.79 ms")
plot(plot_freq, model_2d.(1.79, plot_freq, (get_model_param(fit.param, 1),)), "C4")
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")
title("288503 GHz, 6 mW")
NaCsPlot.maybe_save("$(prefix)_f")

const plot_time = linspace(0, 2.5, 1000)

figure()
NaCsPlot.plot_survival_data(data_nacs_t, fmt="C0.")
plot(plot_time, model_2d.(plot_time, 324.501, (get_model_param(fit.param, 1),)), "C0")
title("288503 GHz, 6 mW, 770.324501 MHz")
text(0.4, 0.19, ("\$f_{res}=$(fit.uncs[4] / 1000 + 770)\$ MHz\n" *
                  "\$\\mathbf{\\Omega_{Raman}=2\\pi\\times$(fit.uncs[5] / 2π) kHz}\$\n" *
                  "\$\\Gamma_{atom}=2\\pi\\times$(fit.uncs[6] / 2π)\$ kHz\n" *
                  "\$\\mathbf{\\Gamma_{molecule}=2\\pi\\times$(fit.uncs[7] / 2π) kHz}\$\n"),
     color="C0", fontsize="small")
xlim([0, 2.6])
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
