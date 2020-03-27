#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const inames = ["data_20200326_101325.mat",
                "data_20200326_223438.mat",
                "data_20200327_004910.mat"]
const datas = [NaCsData.load_striped_mat(joinpath(@__DIR__, "data", iname)) for iname in inames]
const maxcnts = [typemax(Int),
                 typemax(Int),
                 typemax(Int)]
const specs = [(593.0 .+ [-20; -6:2:6; 20], # 15 mW, 0.09 ms
                522.7 .+ [-15; -4.5:1.5:4.5; 15], # 12 mW, 0.12 ms
                447.5 .+ [-10; -3:0.6:3; 5], # 9 mW, 0.16 ms
                369.5 .+ [-5; -2:0.5:2; 5], # 6 mW, 0.27 ms
                287.5 .+ [-5; -0.6:0.15:0.6; 5], # 3 mW, 0.9 ms
                ),
               (593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.1 ms
                593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.14 ms
                ),
               ([0], # 15 mW, 0 ms
                593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.1 ms
                593.0 .+ [-30; -7.5:1.5:7.5; 30], # 15 mW, 0.14 ms
                )
               ]

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

function get_ratio_val(data)
    params, ratios, uncs = NaCsData.get_values(data)
    return ratios[2], uncs[2]
end

const datas_nacs = select_datas(datas, NaCsData.select_single((1, 2), (3, 4,)), maxcnts, specs)
const data_nacs_00 = datas_nacs[3][1]
const data_nacs_08 = datas_nacs[1][1]
const data_nacs_09 = [datas_nacs[2][1]; datas_nacs[3][2]]
const data_nacs_13 = [datas_nacs[2][2]; datas_nacs[3][3]]

const prefix = joinpath(@__DIR__, "imgs", "data_20200326_005204_raman_det_3322")

function model_lorentzian(x, p)
    p[1] .- p[2] ./ (1 .+ ((x .- p[3]) ./ (p[4] / 2)).^2)
end
function model_gaussian(x, p)
    p[1] .- p[2] ./ exp.(((x .- p[3]) ./ p[4]).^2)
end
fit_08 = fit_survival(model_lorentzian, data_nacs_08, [0.35, 0.15, 594, 5])
fit_09 = fit_survival(model_lorentzian, data_nacs_09, [0.35, 0.15, 594, 5])
fit_13 = fit_survival(model_lorentzian, data_nacs_13, [0.35, 0.25, 594, 5])

# @show fit_08.uncs
# @show fit_09.uncs
# @show fit_13.uncs

figure()
ratio_00, uncs_00 = get_ratio_val(data_nacs_00)
errorbar([560, 625], [ratio_00, ratio_00], [uncs_00, uncs_00], fmt="C0.-", label="0ms")
NaCsPlot.plot_survival_data(data_nacs_08, fmt="C1.", label="0.08 ms")
plot(fit_08.plotx, fit_08.ploty, "C1")
NaCsPlot.plot_survival_data(data_nacs_09, fmt="C2.", label="0.09 ms")
plot(fit_09.plotx, fit_09.ploty, "C2")
NaCsPlot.plot_survival_data(data_nacs_13, fmt="C3.", label="0.13 ms")
plot(fit_13.plotx, fit_13.ploty, "C3")
text(557.5, 0.055, ("\$\\Gamma_{0.13}=\\!$(fit_13.uncs[4]) kHz\$\n" *
                    "\$f_{0.13}=\\!$(770 + fit_13.uncs[3] / 1000) MHz\$"),
     color="C3", fontsize="x-small")
text(557.5, 0.105, ("\$\\Gamma_{0.09}=\\!$(fit_09.uncs[4]) kHz\$\n" *
                    "\$f_{0.09}=\\!$(770 + fit_09.uncs[3] / 1000) MHz\$"),
     color="C2", fontsize="x-small")
text(557.5, 0.155, ("\$\\Gamma_{0.08}=\\!$(fit_08.uncs[4]) kHz\$\n" *
                    "\$f_{0.08}=\\!$(770 + fit_08.uncs[3] / 1000) MHz\$"),
     color="C1", fontsize="x-small")
legend(fontsize="x-small", loc="lower right")
grid()
xlabel("2-Photon Detuning (770XXX kHz)")
ylabel("Two-body survival")
title("288560 GHz, 15 mW")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
