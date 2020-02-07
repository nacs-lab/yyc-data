#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200204_232730.mat",
                "data_20200205_145529.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int)]
const specs = [(367.9375 .+ [-20; -6; -5; -4:0.5:4; 5; 6; 20] .* 1e-3,
                367.8625 .+ [-10; -5; -2.5; -2:0.2:2; 2.5; 5; 10] .* 1e-3),
               (367.9375 .+ [-20; -7:1:5; 20] .* 1e-3,
                367.8613 .+ [-10; -5; -2.5:0.5:2.5; 5; 10] .* 1e-3)]
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

const data_1 = [datas_nacs[1][1]; datas_nacs[2][1]]
const data_05_80 = datas_nacs[1][2]
const data_05_40 = datas_nacs[2][2]

const prefix = joinpath(@__DIR__, "imgs", "data_20200204_232730_raman_det_3311")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
fit_1 = fit_survival(model_lorentzian, data_1, [0.7, 0.2, 367.9375, 0.01])
fit_05_80 = fit_survival(model_lorentzian, data_05_80, [0.7, 0.2, 367.8613, 0.01])
fit_05_40 = fit_survival(model_lorentzian, data_05_40, [0.7, 0.2, 367.8613, 0.01])

# @show fit_1.uncs
# @show fit_05_80.uncs
# @show fit_05_40.uncs

figure()
NaCsPlot.plot_survival_data(data_1, fmt="C0.")
plot(fit_1.plotx, fit_1.ploty, "C0-")
title("306490.7 GHz, 1 mW, 20 ms")
# ylim([0.44, 0.66])
text(367.914, 0.51, "\$f=$(fit_1.uncs[3])\$ MHz", color="C0", fontsize="small")
text(367.914, 0.525, "\$\\Gamma=$(fit_1.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
# legend(fontsize="small")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_1")

figure()
NaCsPlot.plot_survival_data(data_05_80, fmt="C0.", label="80 ms")
plot(fit_05_80.plotx, fit_05_80.ploty, "C0-")
NaCsPlot.plot_survival_data(data_05_40, fmt="C1.", label="40 ms")
plot(fit_05_40.plotx, fit_05_40.ploty, "C1-")
title("306490.7 GHz, 0.5 mW")
ylim([0.43, 0.69])
text(367.859, 0.432, "\$f_{80}=$(fit_05_80.uncs[3])\$ MHz", color="C0", fontsize="small")
text(367.863, 0.451, "\$\\Gamma_{80}=$(fit_05_80.uncs[4] * 1000)\$ kHz", color="C0", fontsize="small")
text(367.859, 0.65, "\$f_{40}=$(fit_05_40.uncs[3])\$ MHz", color="C1", fontsize="small")
text(367.859, 0.67, "\$\\Gamma_{40}=$(fit_05_40.uncs[4] * 1000)\$ kHz", color="C1", fontsize="small")
legend(fontsize="small")
grid()
xlabel("2-Photon Detuning (MHz)")
ylabel("Two-body survival")
NaCsPlot.maybe_save("$(prefix)_05")

NaCsPlot.maybe_show()
