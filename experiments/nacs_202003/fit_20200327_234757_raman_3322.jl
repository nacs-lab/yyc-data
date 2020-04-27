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
                "data_20200328_200418.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
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
               ]

select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

function get_ratio_val(data)
    params, ratios, uncs = NaCsData.get_values(data)
    return ratios[2], uncs[2]
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const data_nacs_00 = datas_nacs[3][1]
const data_nacs_20 = datas_nacs[3][3]
const data_nacs_23 = datas_nacs[1][4]
const data_nacs_44 = [datas_nacs[2][1]; datas_nacs[3][2]]
const data_nacs_92 = datas_nacs[4][2]
const data_nacs_t = datas_nacs[5]

function gen_data_all(datas, freqs, times)
    local data_all
    freq_all = Float64[]
    time_all = Float64[]
    for i in 1:length(datas)
        data = datas[i]
        freq = freqs[i]
        time = times[i]
        nd = size(data, 1)
        if isnan(freq)
            append!(freq_all, data.params)
            append!(time_all, linspace(time, time, nd))
        else
            append!(freq_all, linspace(freq, freq, nd))
            append!(time_all, data.params)
        end
        if @isdefined(data_all)
            data_all = [data_all; NaCsData.map_params((i, v)->i + size(data_all, 1), data)]
        else
            data_all = NaCsData.map_params((i, v)->i, data)
        end
    end
    return data_all, freq_all, time_all
end

const data_all, freq_all, time_all =
    gen_data_all([data_nacs_00, data_nacs_20, data_nacs_23,
                  data_nacs_44, data_nacs_92, data_nacs_t],
                 [NaN, NaN, NaN, NaN, NaN, 369.079],
                 [0, 0.20, 0.23, 0.44, 0.92, NaN])

const prefix = joinpath(@__DIR__, "imgs", "fit_20200327_234757_raman_3322")

function model(i, p)
    function wrapper(i)
        t = time_all[i]
        f = freq_all[i]
        return model_2d(t, f, p)
    end
    return wrapper.(i)
end
fit = fit_survival(model, data_all, [0.03, 0.3, 369.079, 2π * 0.7, 0, 2π / 0.4], plotx=false)
@show fit.uncs

const plot_freq = linspace(358.6, 380.3, 1000)


figure()
ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
v00 = model_2d(0, 0, fit.param)
errorbar([358.6, 380.3], [ratio_00, ratio_00], [uncs_00, uncs_00], fmt="C0.", label="0.00 ms")
plot([358.6, 380.3], [v00, v00], "C0-")
NaCsPlot.plot_survival_data(data_nacs_20, fmt="C1.", label="0.20 ms")
plot(plot_freq, model_2d.(0.20, plot_freq, (fit.param,)), "C1")
NaCsPlot.plot_survival_data(data_nacs_23, fmt="C2.", label="0.23 ms")
plot(plot_freq, model_2d.(0.23, plot_freq, (fit.param,)), "C2")
NaCsPlot.plot_survival_data(data_nacs_44, fmt="C3.", label="0.44 ms")
plot(plot_freq, model_2d.(0.44, plot_freq, (fit.param,)), "C3")
NaCsPlot.plot_survival_data(data_nacs_92, fmt="C4.", label="0.92 ms")
plot(plot_freq, model_2d.(0.92, plot_freq, (fit.param,)), "C4")
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")
title("288560 GHz, 6 mW")
NaCsPlot.maybe_save("$(prefix)_f")

const plot_time = linspace(0, 1.44, 1000)

figure()
NaCsPlot.plot_survival_data(data_nacs_t, fmt="C0.")
plot(plot_time, model_2d.(plot_time, 369.079, (fit.param,)), "C0")
title("288560 GHz, 6 mW, 770.369079 MHz")
text(0.26, 0.12, ("\$p_0=$(fit.uncs[1])\$\n" *
                  "\$p_1=$(fit.uncs[2])\$\n" *
                  "\$f_{res}=$(fit.uncs[3] / 1000 + 770)\$ MHz\n" *
                  "\$\\mathbf{\\Omega_{Raman}=2\\pi\\times$(fit.uncs[4] / 2π) kHz}\$\n" *
                  "\$\\Gamma_{atom}=2\\pi\\times$(fit.uncs[5] / 2π)\$ kHz\n" *
                  "\$\\mathbf{\\Gamma_{molecule}=2\\pi\\times$(fit.uncs[6] / 2π) kHz}\$\n"),
     color="C0", fontsize="small")
xlim([0, 1.46])
grid()
xlabel("Raman time (ms)\${}_{\\mathrm{(with\\ 0.04\\ ms\\ offset)}}\$")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_t")

NaCsPlot.maybe_show()
