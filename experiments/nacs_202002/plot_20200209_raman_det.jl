#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200209_004034.mat",
                "data_20200209_084201.mat",
                "data_20200211_002123.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [(303.36 .+ [-2000; -300:50:300; 2000] .* 1e-3, # 8 mW, 0.28 ms
                300.66 .+ [-500; -80:20:80; 500] .* 1e-3, # 4 mW, 0.46 ms
                299.25 .+ [-300; -40:10:40; 300] .* 1e-3, # 2 mW, 0.82 ms
                298.505 .+ [-200; -50; -25:5:25; 50; 200] .* 1e-3), # 1 mW, 2 ms
               (303.36 .+ [-2000; -300:50:300; 2000] .* 1e-3, # 8 mW, 0.23 ms
                300.66 .+ [-500; -80:20:80; 500] .* 1e-3, # 4 mW, 0.4 ms
                299.25 .+ [-300; -20:5:20; 300] .* 1e-3, # 2 mW, 0.82 ms
                298.505 .+ [-200; -20; -15:2.5:15; 20; 200] .* 1e-3, # 1 mW, 1.6 ms
                298.132 .+ [-100; -10; -7.5:1.5:7.5; 10; 100] .* 1e-3), # 0.5 mW, 4.6 ms
               (297.939 .+ [-50; -6:1.5:6; 50] .* 1e-3, # 0.25 mW, 50 ms
                297.939 .+ [-50; -6:1.5:6; 50] .* 1e-3, # 0.25 mW, 30 ms
                297.842 .+ [-30; -9:1.5:9; 30] .* 1e-3, # 0.125 mW, 100 ms
                297.842 .+ [-30; -9:1.5:9; 30] .* 1e-3)] # 0.125 mW, 150 ms
select_datas(datas, selector, maxcnts, specs) =
    [NaCsData.split_data(NaCsData.select_count(data..., selector, maxcnt), spec)
     for (data, maxcnt, spec) in zip(datas, maxcnts, specs)]

fit_data(model, x, y, p0; kws...) =
    fit_data(model, x, y, nothing, p0; kws...)

function fit_data(model, params, ratios, uncs, p0;
                  plotx=nothing, plot_lo=nothing, plot_hi=nothing, plot_scale=1.1)
    use_unc = uncs !== nothing
    if plotx === nothing
        lo = minimum(params)
        hi = maximum(params)
        span = hi - lo
        mid = (hi + lo) / 2
        if plot_lo === nothing
            plot_lo = mid - span * plot_scale / 2
            if plot_lo * lo <= 0
                plot_lo = 0
            end
        end
        if plot_hi === nothing
            plot_hi = mid + span * plot_scale / 2
            if plot_hi * hi <= 0
                plot_hi = 0
            end
        end
        plotx = linspace(plot_lo, plot_hi, 10000)
    end
    if use_unc
        fit = curve_fit(model, params, ratios, uncs.^-(2/3), p0)
    else
        fit = curve_fit(model, params, ratios, p0)
    end
    param = fit.param
    unc = estimate_errors(fit)
    return (param=param, unc=unc,
            uncs=Unc.(param, unc, Sci),
            plotx=plotx, ploty=model.(plotx, (fit.param,)))
end

function fit_survival(model, data, p0; use_unc=true, kws...)
    if use_unc
        params, ratios, uncs = NaCsData.get_values(data)
        return fit_data(model, params, ratios[:, 2], uncs[:, 2], p0; kws...)
    else
        params, ratios, uncs = NaCsData.get_values(data, 0.0)
        return fit_data(model, params, ratios[:, 2], p0; kws...)
    end
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)

const data_8_28 = datas_nacs[1][1]
const data_8_23 = datas_nacs[2][1]
const data_4_46 = datas_nacs[1][2]
const data_4_40 = datas_nacs[2][2]
const data_2_82 = [datas_nacs[1][3]; datas_nacs[2][3]]
const data_1_20 = datas_nacs[1][4]
const data_1_16 = datas_nacs[2][4]
const data_05_46 = datas_nacs[2][5]
const data_025_5 = datas_nacs[3][1]
const data_025_3 = datas_nacs[3][2]
const data_0125_15 = datas_nacs[3][4]
const data_0125_10 = datas_nacs[3][3]

const prefix = joinpath(@__DIR__, "imgs", "data_20200209_raman_det")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end

const fit_8_28 = fit_survival(model_lorentzian, data_8_28, [0.7, 0.3, 303.36, 0.2])
const fit_8_23 = fit_survival(model_lorentzian, data_8_23, [0.7, 0.3, 303.36, 0.2])
const fit_4_46 = fit_survival(model_lorentzian, data_4_46, [0.7, 0.3, 300.66, 0.1])
const fit_4_40 = fit_survival(model_lorentzian, data_4_40, [0.7, 0.3, 300.66, 0.1])
const fit_2_82 = fit_survival(model_lorentzian, data_2_82, [0.7, 0.3, 299.25, 0.05])
const fit_1_20 = fit_survival(model_lorentzian, data_1_20, [0.7, 0.3, 298.505, 0.02])
const fit_1_16 = fit_survival(model_lorentzian, data_1_16, [0.7, 0.3, 298.505, 0.02])
const fit_05_46 = fit_survival(model_lorentzian, data_05_46, [0.7, 0.3, 298.132, 0.01])
const fit_025_5 = fit_survival(model_lorentzian, data_025_5, [0.7, 0.3, 297.939, 0.005])
const fit_025_3 = fit_survival(model_lorentzian, data_025_3, [0.7, 0.3, 297.939, 0.005])
const fit_0125_15 = fit_survival(model_lorentzian, data_0125_15, [0.7, 0.3, 297.842, 0.002])
const fit_0125_10 = fit_survival(model_lorentzian, data_0125_10, [0.7, 0.3, 297.842, 0.002])

figure(figsize=[12.6, 22.4])

subplot(4, 2, 1)
NaCsPlot.plot_survival_data(data_8_28, fmt="C0.", label="0.28 ms")
plot(fit_8_28.plotx, fit_8_28.ploty, "C0-")
NaCsPlot.plot_survival_data(data_8_23, fmt="C1.", label="0.23 ms")
plot(fit_8_23.plotx, fit_8_23.ploty, "C1-")
title("306491 GHz, 8 mW")
# ylim([0.44, 0.66])
text(303.0, 0.38, "\$f_{0.28}=$(fit_8_28.uncs[3])\$ MHz", color="C0", fontsize="small")
text(303.5, 0.42, "\$\\Gamma_{0.28}=$(fit_8_28.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
text(301, 0.485, "\$f_{0.23}=$(fit_8_23.uncs[3])\$ MHz", color="C1", fontsize="small")
text(301, 0.525, "\$\\Gamma_{0.23}=$(fit_8_23.uncs[4] * 1000)\$ kHz", color="C1", fontsize="small")
legend(fontsize="small")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")

subplot(4, 2, 2)
NaCsPlot.plot_survival_data(data_4_46, fmt="C0.", label="0.46 ms")
plot(fit_4_46.plotx, fit_4_46.ploty, "C0-")
NaCsPlot.plot_survival_data(data_4_40, fmt="C1.", label="0.4 ms")
plot(fit_4_40.plotx, fit_4_40.ploty, "C1-")
title("306491 GHz, 4 mW")
# ylim([0.44, 0.66])
text(300.53, 0.405, "\$f_{0.46}=$(fit_4_46.uncs[3])\$ MHz", color="C0", fontsize="small")
text(300.7, 0.435, "\$\\Gamma_{0.46}=$(fit_4_46.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
text(300.08, 0.5, "\$f_{0.4}=$(fit_4_40.uncs[3])\$ MHz", color="C1", fontsize="small")
text(300.08, 0.53, "\$\\Gamma_{0.4}=$(fit_4_40.uncs[4] * 1000)\$ kHz", color="C1", fontsize="small")
legend(fontsize="small")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")

subplot(4, 2, 3)
NaCsPlot.plot_survival_data(data_2_82, fmt="C0.")
plot(fit_2_82.plotx, fit_2_82.ploty, "C0-")
title("306491 GHz, 2 mW, 0.82 ms")
# ylim([0.44, 0.66])
text(298.908, 0.47, "\$f=$(fit_2_82.uncs[3])\$ MHz", color="C0", fontsize="small")
text(298.908, 0.50, "\$\\Gamma=$(fit_2_82.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")

subplot(4, 2, 4)
NaCsPlot.plot_survival_data(data_1_20, fmt="C0.", label="2 ms")
plot(fit_1_20.plotx, fit_1_20.ploty, "C0-")
NaCsPlot.plot_survival_data(data_1_16, fmt="C1.", label="1.6 ms")
plot(fit_1_16.plotx, fit_1_16.ploty, "C1-")
title("306491 GHz, 1 mW")
# ylim([0.44, 0.66])
text(298.45, 0.468, "\$f_{2}=$(fit_1_20.uncs[3])\$ MHz", color="C0", fontsize="small")
text(298.515, 0.498, "\$\\Gamma_{2}=$(fit_1_20.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
text(298.27, 0.55, "\$f_{1.6}=$(fit_1_16.uncs[3])\$ MHz", color="C1", fontsize="small")
text(298.27, 0.58, "\$\\Gamma_{1.6}=$(fit_1_16.uncs[4] * 1000)\$ kHz", color="C1", fontsize="small")
legend(fontsize="small")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")

subplot(4, 2, 5)
NaCsPlot.plot_survival_data(data_05_46, fmt="C0.")
plot(fit_05_46.plotx, fit_05_46.ploty, "C0-")
title("306491 GHz, 0.5 mW, 4.6 ms")
# ylim([0.44, 0.66])
text(298.015, 0.53, "\$f=$(fit_05_46.uncs[3])\$ MHz", color="C0", fontsize="small")
text(298.015, 0.56, "\$\\Gamma=$(fit_05_46.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")

subplot(4, 2, 6)
NaCsPlot.plot_survival_data(data_025_5, fmt="C0.", label="50 ms")
plot(fit_025_5.plotx, fit_025_5.ploty, "C0-")
NaCsPlot.plot_survival_data(data_025_3, fmt="C1.", label="30 ms")
plot(fit_025_3.plotx, fit_025_3.ploty, "C1-")
title("306491 GHz, 0.25 mW")
# ylim([0.44, 0.66])
text(297.927, 0.439, "\$f_{50}=$(fit_025_5.uncs[3])\$ MHz", color="C0", fontsize="small")
text(297.945, 0.452, "\$\\Gamma_{50}=$(fit_025_5.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
text(297.88, 0.465, "\$f_{30}=$(fit_025_3.uncs[3])\$ MHz", color="C1", fontsize="small")
text(297.88, 0.478, "\$\\Gamma_{30}=$(fit_025_3.uncs[4] * 1000)\$ kHz", color="C1", fontsize="small")
legend(fontsize="small", loc="center right")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")

subplot(4, 2, 7)
NaCsPlot.plot_survival_data(data_0125_15, fmt="C0.", label="150 ms")
plot(fit_0125_15.plotx, fit_0125_15.ploty, "C0-")
NaCsPlot.plot_survival_data(data_0125_10, fmt="C1.", label="100 ms")
plot(fit_0125_10.plotx, fit_0125_10.ploty, "C1-")
title("306491 GHz, 0.125 mW")
# ylim([0.44, 0.66])
text(297.836, 0.44, "\$f_{150}=$(fit_0125_15.uncs[3])\$ MHz", color="C0", fontsize="small")
text(297.845, 0.455, "\$\\Gamma_{150}=$(fit_0125_15.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
text(297.806, 0.470, "\$f_{100}=$(fit_0125_10.uncs[3])\$ MHz", color="C1", fontsize="small")
text(297.806, 0.485, "\$\\Gamma_{100}=$(fit_0125_10.uncs[4] * 1000)\$ kHz", color="C1", fontsize="small")
legend(fontsize="small", loc="center right")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
