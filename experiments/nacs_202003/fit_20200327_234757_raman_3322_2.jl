#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

include("molecular_raman_model.jl")

import NaCsCalc.Format: Unc, Sci
using NaCsCalc
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsData.Fitting: fit_data, fit_survival
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200326_101325.mat",
                "data_20200327_234757.mat",
                "data_20200328_015629.mat",
                "data_20200330_033004.mat",
                "data_20200328_200418.mat",
                "data_20200330_214529.mat",
                "data_20200331_020659.mat",
                "../../nacs_202004/data/data_20200401_182808.mat",
                "../../nacs_202004/data/data_20200401_220309.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [(593.0 .+ [-20; -6:2:6; 20], # 15 mW, 0.09 ms
                522.7 .+ [-15; -4.5:1.5:4.5; 15], # 12 mW, 0.12 ms
                447.5 .+ [-10; -3:0.6:3; 5], # 9 mW, 0.16 ms
                369.5 .+ [-5; -2:0.5:2; 5], # 6 mW, 0.27 ms
                287.5 .+ [-5; -0.6:0.15:0.6; 5], # 3 mW, 0.9 ms
                ),
               (369.5 .+ [-10; -2.5:0.5:2.5; 10], # 6 mW, 0.48 ms
                ),
               ([0], # 6 mW, 0 ms
                369.5 .+ [-10; -2.5:0.5:2.5; 10], # 6 mW, 0.48 ms
                369.5 .+ [-10; -2.5:0.5:2.5; 10], # 6 mW, 0.24 ms
                ),
               (593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.32 ms
                369.5 .+ [-10; -2.5:0.5:2.5; 10], # 6 mW, 0.96 ms
                ),
               [0; [0.08, 0.16, 0.24, 0.32, 0.4, 0.48] .* 3 .- 0.04], # 6 mW, 770.369079 MHz
               [0, 20, 40, 80, 120] .* 1.0,
               ([0.0, 2, 4, 8, 12, 20],
                [0.0, 2, 4, 8, 12, 20],
                [0.0, 20, 40, 80, 160],
                [0.0, 20, 40, 80, 160],
                [0.0, 20, 40, 80, 160] .* 2,
                [0.0, 20, 40, 80, 160] .* 2),
               [0, 0.08, 0.16, 0.24, 0.3] .* 2,
               [0, 0.08, 0.16, 0.24, 0.3] .* 2
               ]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function get_ratio_val(data)
    params, ratios, uncs = NaCsData.get_values(data)
    return ratios[2], uncs[2]
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const datas_cs = select_datas(datas, NaCsData.select_single((-1, 2), (-3, 4,)), maxcnts, specs)
const datas_na = select_datas(datas, NaCsData.select_single((1, -2), (3, -4,)), maxcnts, specs)

const data_nacs_00 = datas_nacs[3][1]
const data_nacs_20 = datas_nacs[3][3]
const data_nacs_23 = datas_nacs[1][4]
const data_nacs_44 = [datas_nacs[2][1]; datas_nacs[3][2]]
const data_nacs_92 = datas_nacs[4][2]
const data_nacs_t = datas_nacs[5]

const data_pa_na1 = datas_na[6]
const data_pa_m1 = datas_nacs[6]

const data_pa_a2 = datas_nacs[7][3]
const data_pa_m2 = datas_nacs[7][4]
const data_pa_cs2 = datas_cs[7][3]
const data_pa_na2 = [datas_na[7][3]; datas_na[7][4]]

const data_nacs_m1 = datas_nacs[8]
const data_nacs_m2 = datas_nacs[9]

const data_fit = [NaCsData.map_params((i, v) -> (1, v, 0.0, 1), data_nacs_00);
                  NaCsData.map_params((i, v) -> (1, v, 0.20, 1), data_nacs_20);
                  NaCsData.map_params((i, v) -> (1, v, 0.23, 2), data_nacs_23);
                  NaCsData.map_params((i, v) -> (1, v, 0.44, 1), data_nacs_44);
                  NaCsData.map_params((i, v) -> (1, v, 0.92, 1), data_nacs_92);
                  NaCsData.map_params((i, v) -> (1, 369.079, v, 1), data_nacs_t);
                  NaCsData.map_params((i, v) -> (2, v, 1, 3), data_pa_na1);
                  NaCsData.map_params((i, v) -> (2, v, 4, 4), data_pa_m1);
                  NaCsData.map_params((i, v) -> (2, v, 1, 5), data_pa_na2);
                  NaCsData.map_params((i, v) -> (2, v, 2, 6), data_pa_cs2);
                  NaCsData.map_params((i, v) -> (2, v, 3, 7), data_pa_a2);
                  NaCsData.map_params((i, v) -> (2, v, 4, 8), data_pa_m2);
                  NaCsData.map_params((i, v) -> (3, v, 5, 9, 2), data_nacs_m1);
                  NaCsData.map_params((i, v) -> (3, v, 5, 10, 3), data_nacs_m2);
                  ]

const prefix = joinpath(@__DIR__, "imgs", "fit_20200327_234757_raman_3322_2")

function model_exp(x, p)
    p[1] .* exp.(.- x .* p[2])
end

function model_expoff(x, p)
    p[3] .+ p[1] .* exp.(.- x .* p[2])
end

function get_model_param(p, idx)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0r = p
    p1 = p[9 + idx]
    p0 = p1 * p0r
    return (p0, p1, f0, Ω, Γ1 + Γ_na + Γ_cs, Γ2)
end

function gen_exp_param(p, typ, idx)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs = p
    p1 = p[9 + idx]
    if typ == 1
        return (p1, Γ_na)
    elseif typ == 2
        return (p1, Γ_cs)
    elseif typ == 3
        return (p1, Γ_na + Γ_cs)
    elseif typ == 4
        return (p1, Γ1 + Γ_na + Γ_cs)
    elseif typ == 5
        return (p1, Γ2)
    end
end

function gen_exp_param(p, typ, idx, idx0)
    f0, Ω, Γ1, Γ2, Γ_na, Γ_cs = p
    p0 = p[6 + idx0]
    p1 = p[9 + idx]
    if typ == 1
        return (p1, Γ_na, p0)
    elseif typ == 2
        return (p1, Γ_cs, p0)
    elseif typ == 3
        return (p1, Γ_na + Γ_cs, p0)
    elseif typ == 4
        return (p1, Γ1 + Γ_na + Γ_cs, p0)
    elseif typ == 5
        return (p1, Γ2, p0)
    end
end

function model(xs, p)
    # p: f0, Ω, Γ1, Γ2, Γ_na, Γ_cs, p0s..., p1s...
    function wrapper(x)
        typ = x[1]
        if typ == 1
            _, f, t, idx = x
            return model_2d(t, f, get_model_param(p, idx))
        elseif typ == 2
            t = x[2]
            _, t, etyp, idx = x
            return model_exp(t, gen_exp_param(p, etyp, idx))
        else
            _, t, etyp, idx, idx0 = x
            return model_expoff(t, gen_exp_param(p, etyp, idx, idx0))
        end
    end
    return wrapper.(xs)
end
fit = fit_survival(model, data_fit, [369.079, 2π * 0.7, 0, 2π / 0.4, 0, 0, 0.1, 0.01, 0.01,
                                     0.3, 0.3, 0.8, 0.3, 0.8, 0.3, 0.3, 0.3, 0.1, 0.1],
                   plotx=false)
@show fit.uncs
const param_1 = get_model_param(fit.param, 1)
const uncs_1 = get_model_param(fit.uncs, 1)
const param_2 = get_model_param(fit.param, 2)

const plot_freq_lo = 369.079 - 10
const plot_freq_hi = 369.079 + 10
const plot_freq = linspace(plot_freq_lo, plot_freq_hi, 1000)

figure(figsize=[12.6, 11.2])

subplot(2, 2, 1)
ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
v00 = model_2d(0, 0, param_1)
errorbar([param_1[3]], [ratio_00], [uncs_00], fmt="C0.", label="0.00 ms")
plot([plot_freq_lo, plot_freq_hi], [v00, v00], "C0-")
NaCsPlot.plot_survival_data(data_nacs_20, fmt="C1.", label="0.20 ms")
plot(plot_freq, model_2d.(0.20, plot_freq, (param_1,)), "C1")
NaCsPlot.plot_survival_data(data_nacs_23, fmt="C2.", label="0.23 ms")
plot(plot_freq, model_2d.(0.23, plot_freq, (param_2,)), "C2")
NaCsPlot.plot_survival_data(data_nacs_44, fmt="C3.", label="0.44 ms")
plot(plot_freq, model_2d.(0.44, plot_freq, (param_1,)), "C3")
NaCsPlot.plot_survival_data(data_nacs_92, fmt="C4.", label="0.92 ms")
plot(plot_freq, model_2d.(0.92, plot_freq, (param_1,)), "C4")
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")
title("288560 GHz, 6 mW")

subplot(2, 2, 2)
const plot_time = linspace(0, 1.44, 1000)
NaCsPlot.plot_survival_data(data_nacs_t, fmt="C0.")
plot(plot_time, model_2d.(plot_time, 369.079, (param_1,)), "C0")
title("770.369079 MHz")
text(0.26, 0.12, ("\$f_{res}=$(uncs_1[3] / 1000 + 770)\$ MHz\n" *
                  "\$\\Omega_{Raman}=2\\pi\\times$(uncs_1[4] / 2π)\$ kHz\n" *
                  "\$\\Gamma_{atom}=2\\pi\\times$(uncs_1[5] / 2π * 1000)\$ Hz\n" *
                  "\$\\Gamma_{molecule}=2\\pi\\times$(uncs_1[6] / 2π)\$ kHz\n"),
     color="C0", fontsize="small")
xlim([0, 1.46])
grid()
xlabel("Raman time (ms)\${}_{\\mathrm{(with\\ 0.04\\ ms\\ offset)}}\$")
ylabel("Two-body survival")

subplot(2, 2, 3)
const plot_pa_time = linspace(0, 160, 1000)
NaCsPlot.plot_survival_data(data_pa_na1, fmt="C0s", label="Na 1")
plot(plot_pa_time, model_exp(plot_pa_time, gen_exp_param(fit.param, 1, 3)), "C0-")
NaCsPlot.plot_survival_data(data_pa_na2, fmt="C1s", label="Na 2")
plot(plot_pa_time, model_exp(plot_pa_time, gen_exp_param(fit.param, 1, 5)), "C1-")
NaCsPlot.plot_survival_data(data_pa_cs2, fmt="C2s", label="Cs")
plot(plot_pa_time, model_exp(plot_pa_time, gen_exp_param(fit.param, 2, 6)), "C2-")
NaCsPlot.plot_survival_data(data_pa_a2, fmt="C3s", label="Hot")
plot(plot_pa_time, model_exp(plot_pa_time, gen_exp_param(fit.param, 3, 7)), "C3-")
NaCsPlot.plot_survival_data(data_pa_m1, fmt="C4s", label="Cold 1")
plot(plot_pa_time, model_exp(plot_pa_time, gen_exp_param(fit.param, 4, 4)), "C4-")
NaCsPlot.plot_survival_data(data_pa_m2, fmt="C5s", label="Cold 2")
plot(plot_pa_time, model_exp(plot_pa_time, gen_exp_param(fit.param, 4, 8)), "C5-")
xlim([0, 170])
ylim([0, 0.9])
legend(fontsize="small", ncol=3, borderpad=0.2, labelspacing=0.2,
       handletextpad=0.3, columnspacing=0.2, borderaxespad=0.4, loc="center left")
title("Atomic loss/PA")
grid()
xlabel("Time (ms)")
ylabel("Survival")

subplot(2, 2, 4)
const plot_m_time = linspace(0, 0.62, 1000)
NaCsPlot.plot_survival_data(data_nacs_m1, fmt="C0s")
plot(plot_m_time, model_expoff(plot_m_time, gen_exp_param(fit.param, 5, 9, 2)), "C0-")
NaCsPlot.plot_survival_data(data_nacs_m2, fmt="C1s")
plot(plot_m_time, model_expoff(plot_m_time, gen_exp_param(fit.param, 5, 10, 3)), "C1-")
xlim([0, 0.62])
title("Molecular lifetime")
grid()
xlabel("Molecule Hold Time (ms)")
ylabel("Two-body survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
