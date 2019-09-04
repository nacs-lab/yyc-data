#!/usr/bin/julia

push!(LOAD_PATH, joinpath(@__DIR__, "../../lib"))

import NaCsCalc.Format: Unc, Sci
using NaCsCalc.Utils: interactive
using NaCsData
using NaCsPlot
using PyPlot
using DataStructures
using LsqFit

const pwrs = [15, 10, 5]
const freqs = [916.8, 840, 670]
const freqs_s = [6.0, 20, 15]

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

model_lin0(x, p) = x .* p[1]
model_lin1(x, p) = x .* p[1] .+ p[2]

fit_freq = fit_data(model_lin1, pwrs, freqs, freqs_s, [6.0, 630], plot_lo=0)

const prefix = joinpath(@__DIR__, "imgs", "data_20190903_pa_res")

figure()
errorbar(pwrs, freqs, freqs_s, fmt="C0.")
plot(fit_freq.plotx, fit_freq.ploty, "C0")
text(1, 850, "\$f = f_0 + b \\cdot P\$", color="C0", fontsize="large")
text(5, 550, ("\$f_0 = $(fit_freq.uncs[2]) \\mathrm{MHz}\$\n" *
              "\$b = $(fit_freq.uncs[1]) \\mathrm{MHz\\cdot mW^{-1}}\$"), color="C0")
xlim([0, 16])
grid()
xlabel("Tweezer Power (mW)")
ylabel("306496XXX MHz")
title("PA resonance")
NaCsPlot.maybe_save("$(prefix)")

NaCsPlot.maybe_show()
