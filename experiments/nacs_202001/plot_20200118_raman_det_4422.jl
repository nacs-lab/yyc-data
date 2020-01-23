#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200118_213737.mat",
                "data_20200119_082135.mat",
                "data_20200119_232425.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [(298.227 .+ [-60; -10:2:10; 60] .* 1e-3,
                299.535 .+ [-100; -30:6:30; 100] .* 1e-3,
                300 .+ [-150; -40:8:40; 150] .* 1e-3),
               (298.227 .+ [-60; -10:2:10; 60] .* 1e-3,
                299.535 .+ [-100; -30:6:30; 100] .* 1e-3,
                300 .+ [-230; -120:16:40; 150] .* 1e-3,
                300 .+ [-230; -120:16:40; 150] .* 1e-3),
               (298.227 .+ [-60; -10:2:10; 60] .* 1e-3,
                299.535 .+ [-100; -30:6:30; 100] .* 1e-3,
                300 .+ [-420; -120:16:40; 340] .* 1e-3)]
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

const data_nacs_4 = [datas_nacs[1][1]; datas_nacs[2][1]; datas_nacs[3][1]]
const data_nacs_16 = [datas_nacs[1][2]; datas_nacs[2][2]; datas_nacs[3][2]]
const data_nacs_24 = [datas_nacs[1][3]; datas_nacs[2][3]; datas_nacs[2][4]; datas_nacs[3][3]]

const prefix = joinpath(@__DIR__, "imgs", "data_20200118_raman_det_4422")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
fit_4 = fit_survival(model_lorentzian, data_nacs_4, [0.7, 0.2, 298.227, 0.05])
fit_16 = fit_survival(model_lorentzian, data_nacs_16, [0.7, 0.2, 299.535, 0.05])
fit_24 = fit_survival(model_lorentzian, data_nacs_24, [0.7, 0.2, 299.960, 0.05])

# @show fit1.uncs

figure(figsize=[13.2, 10.0])
subplot(2, 2, 1)
NaCsPlot.plot_survival_data(data_nacs_4, fmt="C0.")
plot(fit_4.plotx, fit_4.ploty, "C0-")
title("4 mW, 1.5 ms")
# ylim([0.18, 0.42])
text(298.16, 0.530, "\$f=$(fit_4.uncs[3])\$ MHz", color="C0", fontsize="small")
text(298.16, 0.555, "\$\\Gamma=$(fit_4.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")

subplot(2, 2, 2)
NaCsPlot.plot_survival_data(data_nacs_16, fmt="C0.")
plot(fit_16.plotx, fit_16.ploty, "C0-")
title("16 mW, 0.5 ms, ratio 0.05")
# ylim([0.18, 0.42])
text(299.43, 0.550, "\$f=$(fit_16.uncs[3])\$ MHz", color="C0", fontsize="small")
text(299.43, 0.575, "\$\\Gamma=$(fit_16.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")

subplot(2, 2, 4)
NaCsPlot.plot_survival_data(data_nacs_24, fmt="C0.")
plot(fit_24.plotx, fit_24.ploty, "C0-")
title("24 mW, 1 ms, ratio 0.02")
# ylim([0.18, 0.42])
text(299.53, 0.550, "\$f=$(fit_24.uncs[3])\$ MHz", color="C0", fontsize="small")
text(299.53, 0.575, "\$\\Gamma=$(fit_24.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")

tight_layout(pad=0.6)
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
